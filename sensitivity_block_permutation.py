"""
Pipeline: Block-Shuffle Permutation Null — Upgraded to 1000 Permutations
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Addresses Kill Shot #2: The original permutation shuffles time bin labels
    between conditions (label-shuffle null). This destroys temporal autocorrelation
    structure shared by both conditions (e.g., slow fluctuations in brain state)
    and may produce an overly conservative null (artificially low null WJ),
    inflating p-values.

    Implements TWO null models for comparison:
    1. Label-shuffle (current): concatenate A+B bins, random split at 1000 perms
    2. Block-shuffle within condition A: shuffle time blocks of ~10 bins within
       condition A, compute WJ(shuffled_A, B). Preserves temporal autocorrelation
       distribution while destroying spatial correlation structure.

    The block-shuffle is the more principled null for the question:
    "Is the correlation structure in condition A distinguishable from a
    temporally scrambled version of itself?" This tests whether the correlation
    architecture is meaningful (time-locked to stimulus) rather than an artifact
    of temporal autocorrelation.

    Also upgrades from 200 to 1000 permutations per comparison to improve
    p-value resolution below p=0.005.

Dependencies: numpy, scipy, pandas, matplotlib, h5py
Input:
    G:/My Drive/inner_architecture_research/neuropixels_wj/data/sub-699733573_ses-715093703.nwb
Output:
    results/block_permutation_comparison.csv
    figures/block_permutation_null_comparison.png
"""
import os
import sys
import time
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata, mannwhitneyu

warnings.filterwarnings('ignore')

sys.path.insert(0, r'G:\My Drive\inner_architecture_research')
from wj_utils import (weighted_jaccard, implementation_divergence,
                      fast_spearman_matrix)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py

# ============================================================================
# CONFIGURATION
# ============================================================================
RANDOM_SEED = 42
FORCE_RECOMPUTE = True
N_PERMUTATIONS = 1000   # Upgraded from 200; provides p<0.001 resolution
BLOCK_SIZE = 10         # Consecutive bins per block for block-shuffle
BIN_SIZE_SEC = 0.1
MIN_FIRING_RATE = 0.5

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

NWB_FILE = os.path.join(DATA_DIR, 'sub-699733573_ses-715093703.nwb')

# Test these representative condition pairs
PAIRS_TO_TEST = [
    ('natural_scenes_presentations_9.0', 'natural_scenes_presentations_10.0'),   # same-stim
    ('natural_scenes_presentations_9.0', 'drifting_gratings_presentations_2.0'), # diff-stim
    ('spontaneous_presentations_spontaneous', 'natural_scenes_presentations_9.0'), # diff-stim
]

np.random.seed(RANDOM_SEED)
START = time.time()
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def log(msg):
    print(f"[{(time.time()-START)/60:5.1f}m] {msg}", flush=True)


def load_spike_data(nwb_path):
    """Load unit spike times and stimulus intervals from NWB."""
    f = h5py.File(nwb_path, 'r')
    units = f['units']
    spike_times_idx = units['spike_times_index'][:]
    spike_times_data = units['spike_times'][:]
    n_units = len(spike_times_idx)
    unit_spikes = []
    prev = 0
    for i in range(n_units):
        end = spike_times_idx[i]
        unit_spikes.append(spike_times_data[prev:end].copy())
        prev = end
    firing_rates = units['firing_rate'][:] if 'firing_rate' in units else None

    stim_info = {}
    if 'intervals' in f:
        for ik in f['intervals'].keys():
            if ik == 'units':
                continue
            try:
                interval = f['intervals'][ik]
                starts = interval['start_time'][:]
                stops = interval['stop_time'][:]
                name_keys = [k for k in interval.keys()
                             if 'stimulus' in k.lower() or 'name' in k.lower()]
                if name_keys:
                    names = interval[name_keys[0]][:]
                    if len(names) > 0 and isinstance(names[0], bytes):
                        names = [n.decode() for n in names]
                    for stype in sorted(set(names)):
                        mask = np.array(names) == stype
                        stim_info[f"{ik}_{stype}"] = {
                            'starts': starts[mask],
                            'stops': stops[mask],
                        }
            except Exception:
                pass
    f.close()
    return unit_spikes, firing_rates, stim_info


def build_spike_matrix(unit_spikes_list, starts, stops, bin_size=BIN_SIZE_SEC):
    n_units = len(unit_spikes_list)
    n_bins_per_win = [max(0, int((s1 - s0) / bin_size)) for s0, s1 in zip(starts, stops)]
    total_bins = sum(n_bins_per_win)
    if total_bins < 10:
        return None
    counts = np.zeros((n_units, total_bins), dtype=np.float32)
    offset = 0
    for t0, t1, n_b in zip(starts, stops, n_bins_per_win):
        if n_b == 0:
            continue
        for i, spikes in enumerate(unit_spikes_list):
            mask = (spikes >= t0) & (spikes < t1)
            if mask.any():
                rel = spikes[mask] - t0
                bidx = np.clip((rel / bin_size).astype(int), 0, n_b - 1) + offset
                np.add.at(counts[i], bidx, 1)
        offset += n_b
    return counts


def block_shuffle(counts, block_size=BLOCK_SIZE, rng=None):
    """
    Block-shuffle columns of counts matrix.
    Divides columns into blocks of `block_size` and shuffles block order.
    Preserves within-block temporal autocorrelation while destroying
    global temporal structure.
    """
    if rng is None:
        rng = np.random.default_rng(RANDOM_SEED)
    n_bins = counts.shape[1]
    n_blocks = int(np.ceil(n_bins / block_size))
    block_indices = np.array_split(np.arange(n_bins), n_blocks)
    shuffled_order = rng.permutation(len(block_indices))
    new_cols = np.concatenate([block_indices[k] for k in shuffled_order])
    return counts[:, new_cols]


# ============================================================================
# PERMUTATION FUNCTIONS
# ============================================================================
def label_shuffle_null(counts_a, counts_b, n_perm=N_PERMUTATIONS, seed=RANDOM_SEED):
    """Current null: concatenate A+B bins, randomly split, compute WJ."""
    rng = np.random.default_rng(seed)
    all_counts = np.hstack([counts_a, counts_b])
    n_a = counts_a.shape[1]
    null_wj = np.zeros(n_perm)
    for p in range(n_perm):
        idx = rng.permutation(all_counts.shape[1])
        pa = fast_spearman_matrix(all_counts[:, idx[:n_a]])
        pb = fast_spearman_matrix(all_counts[:, idx[n_a:]])
        null_wj[p] = weighted_jaccard(pa, pb)
        if (p + 1) % 200 == 0:
            log(f"    label-shuffle {p+1}/{n_perm}")
    return null_wj


def block_shuffle_null(counts_a, counts_b, n_perm=N_PERMUTATIONS,
                       block_size=BLOCK_SIZE, seed=RANDOM_SEED):
    """
    New null: block-shuffle within condition A, compare to original B.
    Tests whether observed WJ is significantly different from what you'd
    get with a temporally scrambled version of condition A's activity.
    """
    rng = np.random.default_rng(seed)
    corr_b = fast_spearman_matrix(counts_b)
    obs_corr_a = fast_spearman_matrix(counts_a)
    obs_wj = weighted_jaccard(obs_corr_a, corr_b)
    null_wj = np.zeros(n_perm)
    for p in range(n_perm):
        shuffled_a = block_shuffle(counts_a, block_size=block_size, rng=rng)
        perm_corr_a = fast_spearman_matrix(shuffled_a)
        null_wj[p] = weighted_jaccard(perm_corr_a, corr_b)
        if (p + 1) % 200 == 0:
            log(f"    block-shuffle {p+1}/{n_perm}")
    return obs_wj, null_wj


# ============================================================================
# MAIN
# ============================================================================
def main():
    log("=" * 70)
    log("BLOCK-SHUFFLE PERMUTATION COMPARISON")
    log(f"N_PERMUTATIONS = {N_PERMUTATIONS}, BLOCK_SIZE = {BLOCK_SIZE} bins")
    log("=" * 70)

    unit_spikes, firing_rates, stim_info = load_spike_data(NWB_FILE)
    n_units = len(unit_spikes)

    all_spikes = np.concatenate([s for s in unit_spikes if len(s) > 0])
    t_total = all_spikes.max() - all_spikes.min()
    if firing_rates is None:
        firing_rates = np.array([len(s) / t_total for s in unit_spikes])
    active_mask = firing_rates >= MIN_FIRING_RATE
    active_spikes = [unit_spikes[i] for i in range(n_units) if active_mask[i]]
    log(f"Active units: {len(active_spikes)}/{n_units}")

    rows = []
    figures_data = []

    for cond_a, cond_b in PAIRS_TO_TEST:
        if cond_a not in stim_info or cond_b not in stim_info:
            log(f"SKIP: {cond_a} or {cond_b} not found")
            continue

        stim_a = cond_a.split('_presentations_')[0] if '_presentations_' in cond_a else cond_a
        stim_b = cond_b.split('_presentations_')[0] if '_presentations_' in cond_b else cond_b
        comp_type = 'same-stimulus' if stim_a == stim_b else 'diff-stimulus'
        log(f"\n{'='*60}")
        log(f"  {comp_type}: {cond_a} vs {cond_b}")

        counts_a = build_spike_matrix(active_spikes,
                                      stim_info[cond_a]['starts'],
                                      stim_info[cond_a]['stops'])
        counts_b = build_spike_matrix(active_spikes,
                                      stim_info[cond_b]['starts'],
                                      stim_info[cond_b]['stops'])
        if counts_a is None or counts_b is None:
            log("  SKIP: insufficient data")
            continue

        # Remove constant units
        var_mask = (np.std(counts_a, axis=1) > 0) & (np.std(counts_b, axis=1) > 0)
        counts_a = counts_a[var_mask]
        counts_b = counts_b[var_mask]
        log(f"  Variable units: {var_mask.sum()} | bins A={counts_a.shape[1]}, B={counts_b.shape[1]}")

        # Observed WJ
        corr_a = fast_spearman_matrix(counts_a)
        corr_b = fast_spearman_matrix(counts_b)
        obs_wj = weighted_jaccard(corr_a, corr_b)
        log(f"  Observed WJ: {obs_wj:.4f}")

        # Null model 1: label-shuffle (current, 1000 perms)
        log(f"  Running label-shuffle null ({N_PERMUTATIONS} perms)...")
        null_label = label_shuffle_null(counts_a, counts_b, N_PERMUTATIONS)
        p_label = float(np.mean(null_label <= obs_wj))
        log(f"  Label-shuffle: null mean={null_label.mean():.4f}, "
            f"std={null_label.std():.4f}, p={p_label:.4f}")

        # Null model 2: block-shuffle (new, 1000 perms)
        log(f"  Running block-shuffle null ({N_PERMUTATIONS} perms, block={BLOCK_SIZE})...")
        obs_wj_block, null_block = block_shuffle_null(counts_a, counts_b, N_PERMUTATIONS)
        p_block = float(np.mean(null_block <= obs_wj))
        log(f"  Block-shuffle: null mean={null_block.mean():.4f}, "
            f"std={null_block.std():.4f}, p={p_block:.4f}")

        # Compare null distributions
        mw_stat, mw_p = mannwhitneyu(null_label, null_block, alternative='two-sided')
        log(f"  Null distributions differ? MW p={mw_p:.4e}")

        rows.append({
            'condition_a': cond_a,
            'condition_b': cond_b,
            'comparison_type': comp_type,
            'n_units': int(var_mask.sum()),
            'n_bins_a': counts_a.shape[1],
            'n_bins_b': counts_b.shape[1],
            'obs_wj': round(obs_wj, 4),
            'label_null_mean': round(float(null_label.mean()), 4),
            'label_null_std': round(float(null_label.std()), 6),
            'label_p': round(p_label, 4),
            'block_null_mean': round(float(null_block.mean()), 4),
            'block_null_std': round(float(null_block.std()), 6),
            'block_p': round(p_block, 4),
            'null_distributions_differ_p': round(float(mw_p), 6),
        })
        figures_data.append({
            'label': f"{comp_type}\n{cond_a.split('_')[-1]} vs {cond_b.split('_')[-1]}",
            'obs_wj': obs_wj,
            'null_label': null_label,
            'null_block': null_block,
        })

    if not rows:
        log("No comparisons completed.")
        return

    df = pd.DataFrame(rows)
    out_path = os.path.join(RESULTS_DIR, 'block_permutation_comparison.csv')
    df.to_csv(out_path, index=False)
    log(f"\nSaved: {out_path}")

    # Figure: null distribution comparisons
    n_plots = len(figures_data)
    if n_plots > 0:
        fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5))
        if n_plots == 1:
            axes = [axes]
        for ax, fdata in zip(axes, figures_data):
            bins = np.linspace(min(fdata['null_label'].min(), fdata['null_block'].min(),
                                   fdata['obs_wj'] - 0.05),
                               max(fdata['null_label'].max(), fdata['null_block'].max(),
                                   fdata['obs_wj'] + 0.05), 40)
            ax.hist(fdata['null_label'], bins=bins, alpha=0.6, color='#1f77b4',
                    label='Label-shuffle null', density=True)
            ax.hist(fdata['null_block'], bins=bins, alpha=0.6, color='#ff7f0e',
                    label='Block-shuffle null', density=True)
            ax.axvline(fdata['obs_wj'], color='red', linewidth=2.5,
                       linestyle='-', label=f"Observed WJ = {fdata['obs_wj']:.3f}")
            ax.set_xlabel('WJ (null distribution)', fontsize=11)
            ax.set_ylabel('Density', fontsize=11)
            ax.set_title(fdata['label'], fontsize=10)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        plt.suptitle(f'Block-Shuffle vs Label-Shuffle Permutation Null\n'
                     f'({N_PERMUTATIONS} permutations each, block_size={BLOCK_SIZE} bins)',
                     fontsize=13, fontweight='bold', y=1.02)
        plt.tight_layout()
        fig_path = os.path.join(FIGURES_DIR, 'block_permutation_null_comparison.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        log(f"Saved: {fig_path}")

    log("\n" + "=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log(df[['comparison_type', 'obs_wj', 'label_null_mean', 'label_p',
            'block_null_mean', 'block_p']].to_string(index=False))
    log("\nInterpretation:")
    log("  If label_p ≈ block_p → null model choice doesn't matter")
    log("  If block_p > label_p → block-shuffle is more conservative (expected)")
    log("  If both p < 0.05 → finding is robust to null model choice")
    log(f"\nTotal time: {(time.time()-START)/60:.1f}m")


if __name__ == '__main__':
    main()
