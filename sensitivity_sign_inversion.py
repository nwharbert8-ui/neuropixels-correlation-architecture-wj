"""
Pipeline: Sign-Inversion Sensitivity Analysis — Stratified by |r| Magnitude
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Addresses reviewer concern that the 92% sign-inversion metric is dominated
    by small-|r| pairs where correlation sign is unreliable due to sampling noise.

    Computes pair-level correlation matrices for representative condition pairs
    (same-stimulus and different-stimulus), then stratifies all pairs by
    minimum(|r_a|, |r_b|) threshold. Reports:
      1. Actual pair-level sign-flip rate at each threshold
      2. WJ-derived sign_inversion_pct metric at each threshold
      3. Chance baseline (50%) and binomial test p-value
      4. Fraction of total pairs surviving each threshold

    If 92% sign-inversion holds at |r|_min >= 0.10 and 0.20, the finding is
    robust to the sampling-noise concern.

Dependencies: numpy, scipy, pandas, matplotlib, h5py
Input:
    G:/My Drive/inner_architecture_research/neuropixels_wj/data/sub-699733573_ses-715093703.nwb
Output:
    results/sign_inversion_sensitivity.csv
    figures/sign_inversion_sensitivity.png
"""
import os
import sys
import time
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from scipy.stats import rankdata

warnings.filterwarnings('ignore')

sys.path.insert(0, r'G:\My Drive\inner_architecture_research')
from wj_utils import (weighted_jaccard, signed_weighted_jaccard,
                      implementation_divergence, fast_spearman_matrix)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py

# ============================================================================
# CONFIGURATION
# ============================================================================
RANDOM_SEED = 42
FORCE_RECOMPUTE = True
BIN_SIZE_SEC = 0.1   # Match main pipeline
MIN_FIRING_RATE = 0.5  # Hz

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

NWB_FILE = os.path.join(DATA_DIR, 'sub-699733573_ses-715093703.nwb')

# |r| magnitude thresholds to test (minimum of |r_a| and |r_b| for a pair)
THRESHOLDS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]

# Representative condition pairs to analyze (same-stim and diff-stim)
PAIRS_TO_ANALYZE = [
    # Same-stimulus: should show sign preservation
    ('natural_scenes_presentations_9.0', 'natural_scenes_presentations_10.0'),
    ('natural_scenes_presentations_9.0', 'natural_scenes_presentations_13.0'),
    # Different-stimulus: the "reorganization" case
    ('natural_scenes_presentations_9.0', 'drifting_gratings_presentations_2.0'),
    ('spontaneous_presentations_spontaneous', 'drifting_gratings_presentations_2.0'),
]

np.random.seed(RANDOM_SEED)
START = time.time()

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def log(msg):
    print(f"[{(time.time()-START)/60:5.1f}m] {msg}", flush=True)


# ============================================================================
# FAST VECTORIZED SPIKE COUNT BUILDER
# ============================================================================
def build_spike_matrix(unit_spikes_list, starts, stops, bin_size=BIN_SIZE_SEC):
    """
    Build spike count matrix from concatenated stimulus windows.
    Returns (n_units, n_bins) float32 array.
    """
    n_units = len(unit_spikes_list)

    # Total bins across all windows
    n_bins_per_window = [max(0, int((s1 - s0) / bin_size)) for s0, s1 in zip(starts, stops)]
    total_bins = sum(n_bins_per_window)
    if total_bins < 10:
        return None

    counts = np.zeros((n_units, total_bins), dtype=np.float32)
    bin_offset = 0

    for win_idx, (t0, t1, n_b) in enumerate(zip(starts, stops, n_bins_per_window)):
        if n_b == 0:
            continue
        for i, spikes in enumerate(unit_spikes_list):
            mask = (spikes >= t0) & (spikes < t1)
            if mask.any():
                rel_times = spikes[mask] - t0
                bin_idxs = (rel_times / bin_size).astype(int)
                bin_idxs = np.clip(bin_idxs, 0, n_b - 1) + bin_offset
                np.add.at(counts[i], bin_idxs, 1)
        bin_offset += n_b

    return counts


def spearman_corr_matrix(counts):
    """Compute Spearman correlation matrix (units x units)."""
    return fast_spearman_matrix(counts)


# ============================================================================
# LOAD NWB SESSION
# ============================================================================
def load_session_data(nwb_path):
    """Load spike times, firing rates, and stimulus intervals."""
    log(f"Loading: {os.path.basename(nwb_path)}")
    f = h5py.File(nwb_path, 'r')
    units = f['units']

    # Spike times
    spike_times_idx = units['spike_times_index'][:]
    spike_times_data = units['spike_times'][:]
    n_units = len(spike_times_idx)

    unit_spikes = []
    prev = 0
    for i in range(n_units):
        end = spike_times_idx[i]
        unit_spikes.append(spike_times_data[prev:end].copy())
        prev = end

    # Firing rates
    firing_rates = units['firing_rate'][:] if 'firing_rate' in units else None

    # Stimulus intervals
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
                        key = f"{ik}_{stype}"
                        stim_info[key] = {
                            'starts': starts[mask],
                            'stops': stops[mask],
                        }
                else:
                    stim_info[ik] = {'starts': starts, 'stops': stops}
            except Exception as e:
                log(f"  Warning: {ik}: {e}")

    f.close()
    log(f"  {n_units} units, {len(stim_info)} stimulus conditions")
    return unit_spikes, firing_rates, stim_info


# ============================================================================
# SIGN INVERSION ANALYSIS
# ============================================================================
def analyze_sign_inversions(r_a, r_b, label, thresholds=THRESHOLDS):
    """
    Stratify pairs by min(|r_a|, |r_b|) and compute sign-flip rate at each threshold.

    Returns DataFrame with columns:
        threshold, n_pairs, pct_total, sign_flip_rate, sign_inv_pct, binom_p
    """
    r_a_flat = r_a.flatten()
    r_b_flat = r_b.flatten()
    total_pairs = len(r_a_flat)

    # Only upper triangle pairs
    n = int(np.sqrt(total_pairs))
    idx_i, idx_j = np.triu_indices(n, k=1)
    r_a_pairs = r_a[idx_i, idx_j]
    r_b_pairs = r_b[idx_i, idx_j]

    # Sign flip: true when signs differ
    sign_flip = np.sign(r_a_pairs) != np.sign(r_b_pairs)

    # Magnitude: use min of |r_a|, |r_b| (conservative — both must be above threshold)
    min_abs_r = np.minimum(np.abs(r_a_pairs), np.abs(r_b_pairs))

    rows = []
    for thresh in thresholds:
        mask = min_abs_r >= thresh
        n_pairs = mask.sum()
        if n_pairs < 10:
            continue

        flip_rate = float(sign_flip[mask].mean())
        n_flip = int(sign_flip[mask].sum())

        # Binomial test against 50% chance baseline
        binom_result = binomtest(n_flip, n_pairs, p=0.5, alternative='two-sided')
        binom_p = binom_result.pvalue

        # WJ sign_inversion_pct restricted to this subset
        ra_sub = r_a_pairs[mask]
        rb_sub = r_b_pairs[mask]
        # Build fake matrices for WJ computation (diagonal = 1 for placeholder)
        wj_u = _wj_from_pairs(ra_sub, rb_sub, unsigned=True)
        wj_s = _wj_from_pairs(ra_sub, rb_sub, unsigned=False)
        gap = wj_s - wj_u
        reorg = 1 - wj_u
        sign_inv_pct = (gap / reorg * 100) if reorg > 1e-9 else 0.0

        rows.append({
            'comparison': label,
            'threshold': thresh,
            'n_pairs': n_pairs,
            'pct_total': 100 * n_pairs / len(r_a_pairs),
            'sign_flip_rate_pct': round(flip_rate * 100, 2),
            'sign_inv_pct_metric': round(sign_inv_pct, 2),
            'wj_unsigned': round(wj_u, 4),
            'wj_signed': round(wj_s, 4),
            'binom_p': round(binom_p, 6),
            'n_sign_flips': n_flip,
        })

    return pd.DataFrame(rows)


def _wj_from_pairs(ra, rb, unsigned=True):
    """Compute WJ from 1D arrays of pair-level r values."""
    if unsigned:
        a = np.abs(ra)
        b = np.abs(rb)
    else:
        a = ra + 1.0  # shift [-1,1] -> [0,2]
        b = rb + 1.0
    mins = np.minimum(a, b)
    maxs = np.maximum(a, b)
    denom = maxs.sum()
    return float(mins.sum() / denom) if denom > 1e-12 else 0.0


# ============================================================================
# MAIN
# ============================================================================
def main():
    log("=" * 70)
    log("SIGN INVERSION SENSITIVITY ANALYSIS")
    log("=" * 70)

    unit_spikes, firing_rates, stim_info = load_session_data(NWB_FILE)
    n_units = len(unit_spikes)

    # Filter by firing rate (match main pipeline)
    all_spikes = np.concatenate([s for s in unit_spikes if len(s) > 0])
    t_total = all_spikes.max() - all_spikes.min()
    if firing_rates is None:
        firing_rates = np.array([len(s) / t_total for s in unit_spikes])
    active_mask = firing_rates >= MIN_FIRING_RATE
    active_spikes = [unit_spikes[i] for i in range(n_units) if active_mask[i]]
    n_active = len(active_spikes)
    log(f"Active units (>={MIN_FIRING_RATE}Hz): {n_active}/{n_units}")

    available = sorted(stim_info.keys())
    log(f"Available conditions ({len(available)}):")
    for k in available[:20]:
        log(f"  {k}")

    all_dfs = []

    for cond_a, cond_b in PAIRS_TO_ANALYZE:
        if cond_a not in stim_info:
            log(f"  SKIP: {cond_a} not found")
            continue
        if cond_b not in stim_info:
            log(f"  SKIP: {cond_b} not found")
            continue

        log(f"\n--- {cond_a} vs {cond_b} ---")
        info_a = stim_info[cond_a]
        info_b = stim_info[cond_b]

        log(f"  Building spike matrix A ({len(info_a['starts'])} windows)...")
        counts_a = build_spike_matrix(active_spikes, info_a['starts'], info_a['stops'])
        log(f"  Building spike matrix B ({len(info_b['starts'])} windows)...")
        counts_b = build_spike_matrix(active_spikes, info_b['starts'], info_b['stops'])

        if counts_a is None or counts_b is None:
            log("  SKIP: insufficient data")
            continue

        # Remove constant units (zero variance across bins)
        std_a = np.std(counts_a, axis=1)
        std_b = np.std(counts_b, axis=1)
        var_mask = (std_a > 0) & (std_b > 0)
        counts_a = counts_a[var_mask]
        counts_b = counts_b[var_mask]
        n_var = var_mask.sum()
        log(f"  Variable units: {n_var}/{n_active} | bins A={counts_a.shape[1]}, B={counts_b.shape[1]}")

        if n_var < 50:
            log("  SKIP: fewer than 50 variable units")
            continue

        log("  Computing Spearman correlation matrices...")
        corr_a = spearman_corr_matrix(counts_a)
        corr_b = spearman_corr_matrix(counts_b)

        # Global WJ for reference
        div = implementation_divergence(corr_a, corr_b)
        log(f"  Global: WJ_u={div['wj_unsigned']:.4f}, WJ_s={div['wj_signed']:.4f}, "
            f"sign_inv_pct={div['sign_inversion_pct']:.1f}%")

        # Pair-level analysis
        n = corr_a.shape[0]
        idx_i, idx_j = np.triu_indices(n, k=1)
        ra_pairs = corr_a[idx_i, idx_j]
        rb_pairs = corr_b[idx_i, idx_j]
        n_pairs_total = len(ra_pairs)

        # Sign flip overall
        sign_flip_all = np.sign(ra_pairs) != np.sign(rb_pairs)
        log(f"  Overall sign flip rate: {sign_flip_all.mean()*100:.1f}% of {n_pairs_total:,} pairs")

        label = f"{cond_a.split('_')[-2]}{cond_a.split('_')[-1]}_vs_{cond_b.split('_')[-2]}{cond_b.split('_')[-1]}"
        # Determine if same or different stimulus class
        stim_class_a = cond_a.split('_presentations_')[0] if '_presentations_' in cond_a else cond_a
        stim_class_b = cond_b.split('_presentations_')[0] if '_presentations_' in cond_b else cond_b
        comparison_type = 'same-stimulus' if stim_class_a == stim_class_b else 'diff-stimulus'
        label = f"{comparison_type}: {cond_a.split('_')[-1]} vs {cond_b.split('_')[-1]}"

        log("\n  Stratification by min(|r_a|, |r_b|) threshold:")
        log(f"  {'Threshold':>10} {'N pairs':>10} {'% total':>8} {'Sign flip%':>12} {'Chance':>8} {'Binom-p':>12}")

        df = analyze_sign_inversions(corr_a, corr_b, label)
        for _, row in df.iterrows():
            log(f"  {row['threshold']:>10.2f} {row['n_pairs']:>10,} {row['pct_total']:>8.1f}% "
                f"{row['sign_flip_rate_pct']:>10.1f}% {'50.0%':>8} {row['binom_p']:>12.2e}")

        all_dfs.append(df)

    if not all_dfs:
        log("ERROR: No comparisons completed. Check condition names vs available list.")
        return

    out_df = pd.concat(all_dfs, ignore_index=True)
    out_path = os.path.join(RESULTS_DIR, 'sign_inversion_sensitivity.csv')
    out_df.to_csv(out_path, index=False)
    log(f"\nSaved: {out_path}")

    # Figure
    comparisons = out_df['comparison'].unique()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(comparisons)))

    # Panel A: Pair-level sign flip rate vs threshold
    ax = axes[0]
    for (comp, color) in zip(comparisons, colors):
        sub = out_df[out_df['comparison'] == comp]
        ax.plot(sub['threshold'], sub['sign_flip_rate_pct'], 'o-',
                label=comp[:40], color=color, linewidth=2, markersize=7)
    ax.axhline(50, color='red', linestyle='--', linewidth=1.5, label='Chance (50%)', zorder=5)
    ax.set_xlabel('Min |r| Threshold', fontsize=12)
    ax.set_ylabel('Sign Flip Rate (%)', fontsize=12)
    ax.set_title('Pair-Level Sign Flip Rate\nby Correlation Magnitude Threshold', fontsize=12)
    ax.legend(fontsize=7, loc='upper right')
    ax.set_ylim(0, 110)
    ax.grid(True, alpha=0.3)

    # Panel B: WJ sign_inversion_pct metric vs threshold
    ax = axes[1]
    for (comp, color) in zip(comparisons, colors):
        sub = out_df[out_df['comparison'] == comp]
        ax.plot(sub['threshold'], sub['sign_inv_pct_metric'], 's-',
                label=comp[:40], color=color, linewidth=2, markersize=7)
    ax.axhline(0, color='gray', linestyle='-', linewidth=0.8, alpha=0.5)
    ax.axhline(50, color='orange', linestyle='--', linewidth=1.5,
               label='50% reference', zorder=5)
    ax.set_xlabel('Min |r| Threshold', fontsize=12)
    ax.set_ylabel('WJ Sign-Inversion % Metric', fontsize=12)
    ax.set_title('WJ-Derived Sign-Inversion Metric\nby Correlation Magnitude Threshold', fontsize=12)
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Sign-Inversion Sensitivity — Stratified by |r| Magnitude\n'
                 '(Robustness check: does 92% metric hold at high |r|?)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, 'sign_inversion_sensitivity.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    log(f"Saved: {fig_path}")

    # Print summary
    log("\n" + "=" * 70)
    log("SUMMARY: Sign Inversion Sensitivity")
    log("=" * 70)
    log(out_df[['comparison', 'threshold', 'n_pairs', 'sign_flip_rate_pct',
                'sign_inv_pct_metric', 'binom_p']].to_string(index=False))
    log("\nInterpretation guide:")
    log("  sign_flip_rate_pct > 50% AND binom_p < 0.05 -> sign inversion is real at that threshold")
    log("  sign_flip_rate_pct ~= 50% -> sign pattern is consistent with random noise")
    log("  sign_inv_pct_metric ~= 92% at high thresholds -> metric holds for strong correlations")
    log(f"\nTotal time: {(time.time()-START)/60:.1f}m")


if __name__ == '__main__':
    main()
