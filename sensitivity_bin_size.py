"""
Pipeline: Bin Size Sensitivity Analysis (50ms, 100ms, 200ms, 500ms)
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Addresses Kill Shot #5: The 100ms bin size was not empirically validated
    for Neuropixels data. Different bin sizes change the effective temporal
    averaging of spike counts, altering the correlation structure and potentially
    the WJ values.

    Tests bin sizes: 50ms, 100ms (current), 200ms, 500ms on representative
    condition pairs. Reports WJ, sign_inversion_pct, Pearson WJ, and n_bins
    at each bin size.

    Expected: WJ should be quantitatively similar across bin sizes if the
    finding is robust. Large changes would indicate bin-size dependence.

Dependencies: numpy, scipy, pandas, matplotlib, h5py
Input:
    G:/My Drive/inner_architecture_research/neuropixels_wj/data/sub-699733573_ses-715093703.nwb
Output:
    results/bin_size_sensitivity.csv
    figures/bin_size_sensitivity.png
"""
import os
import sys
import time
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

sys.path.insert(0, r'G:\My Drive\inner_architecture_research')
from wj_utils import (weighted_jaccard, implementation_divergence,
                      fast_spearman_matrix, fast_pearson_matrix)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py

# ============================================================================
# CONFIGURATION
# ============================================================================
RANDOM_SEED = 42
FORCE_RECOMPUTE = True
MIN_FIRING_RATE = 0.5

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

NWB_FILE = os.path.join(DATA_DIR, 'sub-699733573_ses-715093703.nwb')

# Bin sizes to test in seconds
BIN_SIZES = [0.05, 0.10, 0.20, 0.50]
BIN_SIZE_LABELS = ['50ms', '100ms (current)', '200ms', '500ms']

# Condition pairs to test
PAIRS_TO_TEST = [
    # Same-stimulus pairs
    ('natural_scenes_presentations_9.0', 'natural_scenes_presentations_10.0', 'same-stim'),
    ('natural_scenes_presentations_9.0', 'natural_scenes_presentations_13.0', 'same-stim'),
    # Different-stimulus pairs
    ('natural_scenes_presentations_9.0', 'drifting_gratings_presentations_2.0', 'diff-stim'),
    ('spontaneous_presentations_spontaneous', 'natural_scenes_presentations_9.0', 'diff-stim'),
]

np.random.seed(RANDOM_SEED)
START = time.time()
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def log(msg):
    print(f"[{(time.time()-START)/60:5.1f}m] {msg}", flush=True)


def load_spike_data(nwb_path):
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


def build_spike_matrix(unit_spikes_list, starts, stops, bin_size):
    n_units = len(unit_spikes_list)
    n_bins_per_win = [max(0, int((s1 - s0) / bin_size)) for s0, s1 in zip(starts, stops)]
    total_bins = sum(n_bins_per_win)
    if total_bins < 5:
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


# ============================================================================
# MAIN
# ============================================================================
def main():
    log("=" * 70)
    log("BIN SIZE SENSITIVITY ANALYSIS")
    log(f"Bin sizes: {BIN_SIZES}")
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

    for cond_a, cond_b, comp_type in PAIRS_TO_TEST:
        if cond_a not in stim_info or cond_b not in stim_info:
            log(f"SKIP: {cond_a} not found")
            continue

        log(f"\n{'='*60}")
        log(f"  {comp_type}: {cond_a} vs {cond_b}")

        for bin_size, bs_label in zip(BIN_SIZES, BIN_SIZE_LABELS):
            log(f"  Bin size: {bs_label}...")
            counts_a = build_spike_matrix(active_spikes,
                                          stim_info[cond_a]['starts'],
                                          stim_info[cond_a]['stops'],
                                          bin_size=bin_size)
            counts_b = build_spike_matrix(active_spikes,
                                          stim_info[cond_b]['starts'],
                                          stim_info[cond_b]['stops'],
                                          bin_size=bin_size)
            if counts_a is None or counts_b is None:
                log(f"    SKIP: insufficient data at {bs_label}")
                continue

            # Remove constant units
            var_mask = (np.std(counts_a, axis=1) > 0) & (np.std(counts_b, axis=1) > 0)
            counts_a = counts_a[var_mask]
            counts_b = counts_b[var_mask]
            n_var = int(var_mask.sum())

            if n_var < 50:
                log(f"    SKIP: only {n_var} variable units at {bs_label}")
                continue

            corr_a = fast_spearman_matrix(counts_a)
            corr_b = fast_spearman_matrix(counts_b)
            div = implementation_divergence(corr_a, corr_b)

            # Pearson sensitivity
            corr_a_p = fast_pearson_matrix(counts_a)
            corr_b_p = fast_pearson_matrix(counts_b)
            wj_pearson = weighted_jaccard(corr_a_p, corr_b_p)

            log(f"    n_units={n_var}, bins_A={counts_a.shape[1]}, bins_B={counts_b.shape[1]}")
            log(f"    WJ_u={div['wj_unsigned']:.4f}, WJ_s={div['wj_signed']:.4f}, "
                f"sign_inv={div['sign_inversion_pct']:.1f}%, WJ_pearson={wj_pearson:.4f}")

            rows.append({
                'condition_a': cond_a,
                'condition_b': cond_b,
                'comparison_type': comp_type,
                'bin_size_s': bin_size,
                'bin_size_label': bs_label,
                'n_units': n_var,
                'n_bins_a': counts_a.shape[1],
                'n_bins_b': counts_b.shape[1],
                'wj_unsigned': round(div['wj_unsigned'], 4),
                'wj_signed': round(div['wj_signed'], 4),
                'gap': round(div['gap'], 4),
                'sign_inversion_pct': round(div['sign_inversion_pct'], 1),
                'wj_pearson': round(wj_pearson, 4),
            })

    if not rows:
        log("No comparisons completed.")
        return

    df = pd.DataFrame(rows)
    out_path = os.path.join(RESULTS_DIR, 'bin_size_sensitivity.csv')
    df.to_csv(out_path, index=False)
    log(f"\nSaved: {out_path}")

    # Figure: WJ vs bin size for each comparison
    pairs_done = df[['condition_a', 'condition_b', 'comparison_type']].drop_duplicates()
    n_pairs = len(pairs_done)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = plt.cm.tab10(np.linspace(0, 1, n_pairs))

    # Panel A: WJ_unsigned vs bin size
    ax = axes[0]
    for (_, row), color in zip(pairs_done.iterrows(), colors):
        sub = df[(df['condition_a'] == row['condition_a']) &
                 (df['condition_b'] == row['condition_b'])]
        lbl = f"{row['comparison_type']}: {row['condition_a'].split('_')[-1]} vs {row['condition_b'].split('_')[-1]}"
        ax.plot(sub['bin_size_s'] * 1000, sub['wj_unsigned'], 'o-',
                color=color, label=lbl[:45], linewidth=2, markersize=7)
    ax.set_xlabel('Bin Size (ms)', fontsize=12)
    ax.set_ylabel('WJ (unsigned)', fontsize=12)
    ax.set_xscale('log')
    ax.set_xticks([50, 100, 200, 500])
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_title('WJ (Unsigned Spearman) vs Bin Size', fontsize=12)
    ax.legend(fontsize=7, loc='best')
    ax.grid(True, alpha=0.3)

    # Panel B: sign_inversion_pct vs bin size
    ax = axes[1]
    for (_, row), color in zip(pairs_done.iterrows(), colors):
        sub = df[(df['condition_a'] == row['condition_a']) &
                 (df['condition_b'] == row['condition_b'])]
        lbl = f"{row['comparison_type']}: {row['condition_a'].split('_')[-1]} vs {row['condition_b'].split('_')[-1]}"
        ax.plot(sub['bin_size_s'] * 1000, sub['sign_inversion_pct'], 's-',
                color=color, label=lbl[:45], linewidth=2, markersize=7)
    ax.set_xlabel('Bin Size (ms)', fontsize=12)
    ax.set_ylabel('Sign-Inversion % Metric', fontsize=12)
    ax.set_xscale('log')
    ax.set_xticks([50, 100, 200, 500])
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_title('Sign-Inversion Metric vs Bin Size', fontsize=12)
    ax.legend(fontsize=7, loc='best')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Bin Size Sensitivity Analysis\n(50ms, 100ms, 200ms, 500ms)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, 'bin_size_sensitivity.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    log(f"Saved: {fig_path}")

    log("\n" + "=" * 70)
    log("SUMMARY TABLE")
    log("=" * 70)
    log(df[['comparison_type', 'bin_size_label', 'n_units', 'n_bins_a',
            'wj_unsigned', 'sign_inversion_pct']].to_string(index=False))
    log(f"\nTotal time: {(time.time()-START)/60:.1f}m")


if __name__ == '__main__':
    main()
