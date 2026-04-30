"""
Pipeline: Spike Sorting Quality Filter — Allen QC Criteria
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Addresses Kill Shot #1: The main pipeline filters units only by firing rate
    (>= 0.5 Hz), excluding the Allen-standard spike sorting quality criteria.
    Multi-unit noise and poorly isolated units could introduce spurious
    correlations that inflate WJ differences between conditions.

    Applies Allen Brain Observatory QC thresholds:
      - quality == 'good' (vs 'noise')
      - isi_violations < 0.5  (contamination from refractory period violations)
      - presence_ratio >= 0.9 (unit present throughout most of recording)
      - amplitude_cutoff < 0.1 (minimal spike clipping)

    Compares n_units before/after filtering and re-runs WJ analysis on
    representative condition pairs to verify findings are robust to unit quality.

Dependencies: numpy, scipy, pandas, matplotlib, h5py
Input:
    G:/My Drive/inner_architecture_research/neuropixels_wj/data/sub-699733573_ses-715093703.nwb
Output:
    results/quality_filter_effect.csv
    figures/quality_filter_effect.png
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
MIN_FIRING_RATE = 0.5

# Allen QC thresholds
QC_QUALITY = 'good'
QC_ISI_MAX = 0.5
QC_PRESENCE_MIN = 0.9
QC_AMPLITUDE_CUTOFF_MAX = 0.1

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

NWB_FILE = os.path.join(DATA_DIR, 'sub-699733573_ses-715093703.nwb')

PAIRS_TO_TEST = [
    ('natural_scenes_presentations_9.0', 'natural_scenes_presentations_10.0', 'same-stim'),
    ('natural_scenes_presentations_9.0', 'drifting_gratings_presentations_2.0', 'diff-stim'),
    ('spontaneous_presentations_spontaneous', 'natural_scenes_presentations_9.0', 'diff-stim'),
]

np.random.seed(RANDOM_SEED)
START = time.time()
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def log(msg):
    print(f"[{(time.time()-START)/60:5.1f}m] {msg}", flush=True)


def load_all_unit_data(nwb_path):
    """
    Load spike times and all quality metrics from NWB.
    Returns: unit_spikes list, quality_df DataFrame with QC columns.
    """
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

    # Build quality DataFrame
    qc_cols = {}
    for col in ['firing_rate', 'isi_violations', 'presence_ratio',
                 'amplitude_cutoff', 'quality', 'snr', 'd_prime']:
        if col in units:
            try:
                vals = units[col][:]
                if vals.dtype.kind in ('S', 'O'):
                    vals = np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
                qc_cols[col] = vals
            except Exception as e:
                log(f"  Could not read {col}: {e}")

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

    quality_df = pd.DataFrame(qc_cols)
    quality_df.index.name = 'unit_idx'
    return unit_spikes, quality_df, stim_info


def build_spike_matrix(unit_spikes_list, starts, stops, bin_size=0.1):
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


def compute_wj_for_units(unit_indices, unit_spikes_all, stim_info, cond_a, cond_b):
    """Compute WJ for a given subset of unit indices."""
    spikes_subset = [unit_spikes_all[i] for i in unit_indices]
    counts_a = build_spike_matrix(spikes_subset, stim_info[cond_a]['starts'],
                                   stim_info[cond_a]['stops'])
    counts_b = build_spike_matrix(spikes_subset, stim_info[cond_b]['starts'],
                                   stim_info[cond_b]['stops'])
    if counts_a is None or counts_b is None:
        return None
    var_mask = (np.std(counts_a, axis=1) > 0) & (np.std(counts_b, axis=1) > 0)
    counts_a = counts_a[var_mask]
    counts_b = counts_b[var_mask]
    if var_mask.sum() < 30:
        return None
    corr_a = fast_spearman_matrix(counts_a)
    corr_b = fast_spearman_matrix(counts_b)
    div = implementation_divergence(corr_a, corr_b)
    return div, int(var_mask.sum())


# ============================================================================
# MAIN
# ============================================================================
def main():
    log("=" * 70)
    log("SPIKE SORTING QUALITY FILTER EFFECT")
    log(f"Allen QC: quality={QC_QUALITY!r}, ISI<{QC_ISI_MAX}, "
        f"presence>={QC_PRESENCE_MIN}, amplitude_cutoff<{QC_AMPLITUDE_CUTOFF_MAX}")
    log("=" * 70)

    unit_spikes, qdf, stim_info = load_all_unit_data(NWB_FILE)
    n_total = len(unit_spikes)
    log(f"Total units: {n_total}")
    log(f"QC columns available: {list(qdf.columns)}")

    # ---- Define filter masks ----
    all_idx = np.arange(n_total)

    # Firing rate filter (current pipeline)
    all_spikes_cat = np.concatenate([s for s in unit_spikes if len(s) > 0])
    t_total = all_spikes_cat.max() - all_spikes_cat.min()
    if 'firing_rate' in qdf.columns:
        fr = qdf['firing_rate'].values
    else:
        fr = np.array([len(s) / t_total for s in unit_spikes])
    fr_mask = fr >= MIN_FIRING_RATE

    # Allen QC mask
    qc_mask = np.ones(n_total, dtype=bool)
    if 'quality' in qdf.columns:
        qc_mask &= (qdf['quality'].values == QC_QUALITY)
        n_good = qc_mask.sum()
        log(f"  quality=='good': {n_good}/{n_total} units")
    if 'isi_violations' in qdf.columns:
        isi_pass = qdf['isi_violations'].values < QC_ISI_MAX
        log(f"  isi_violations<{QC_ISI_MAX}: {isi_pass.sum()}/{n_total} units")
        qc_mask &= isi_pass
    if 'presence_ratio' in qdf.columns:
        pr_pass = qdf['presence_ratio'].values >= QC_PRESENCE_MIN
        log(f"  presence_ratio>={QC_PRESENCE_MIN}: {pr_pass.sum()}/{n_total} units")
        qc_mask &= pr_pass
    if 'amplitude_cutoff' in qdf.columns:
        ac_pass = qdf['amplitude_cutoff'].values < QC_AMPLITUDE_CUTOFF_MAX
        log(f"  amplitude_cutoff<{QC_AMPLITUDE_CUTOFF_MAX}: {ac_pass.sum()}/{n_total} units")
        qc_mask &= ac_pass

    # Combined mask: FR + Allen QC
    combined_mask = fr_mask & qc_mask

    log(f"\nFilter comparison:")
    log(f"  Original (FR >= {MIN_FIRING_RATE} Hz):          {fr_mask.sum()}/{n_total}")
    log(f"  Allen QC only:                     {qc_mask.sum()}/{n_total}")
    log(f"  Combined (FR + Allen QC):           {combined_mask.sum()}/{n_total}")

    # Unit indices for each filter
    idx_fr = np.where(fr_mask)[0]
    idx_combined = np.where(combined_mask)[0]

    # QC column distributions
    log("\nQC metric distributions for FR-passing units:")
    for col in ['isi_violations', 'presence_ratio', 'amplitude_cutoff']:
        if col in qdf.columns:
            vals = qdf[col].values[fr_mask]
            log(f"  {col}: mean={vals.mean():.3f}, median={np.median(vals):.3f}, "
                f"std={vals.std():.3f}, range=[{vals.min():.3f}, {vals.max():.3f}]")

    # ---- Compare WJ results with and without QC filter ----
    rows = []
    for cond_a, cond_b, comp_type in PAIRS_TO_TEST:
        if cond_a not in stim_info or cond_b not in stim_info:
            log(f"\nSKIP: {cond_a} or {cond_b} not in stim_info")
            continue

        log(f"\n--- {comp_type}: {cond_a} vs {cond_b} ---")

        for filter_label, unit_idx in [('FR_only', idx_fr), ('FR+AllenQC', idx_combined)]:
            log(f"  Computing WJ with {filter_label} ({len(unit_idx)} units)...")
            result = compute_wj_for_units(unit_idx, unit_spikes, stim_info, cond_a, cond_b)
            if result is None:
                log(f"  SKIP: insufficient units or data for {filter_label}")
                continue
            div, n_var = result
            log(f"    n_variable_units={n_var}, WJ_u={div['wj_unsigned']:.4f}, "
                f"WJ_s={div['wj_signed']:.4f}, sign_inv={div['sign_inversion_pct']:.1f}%")
            rows.append({
                'condition_a': cond_a,
                'condition_b': cond_b,
                'comparison_type': comp_type,
                'filter': filter_label,
                'n_units_filter': len(unit_idx),
                'n_units_variable': n_var,
                'wj_unsigned': round(div['wj_unsigned'], 4),
                'wj_signed': round(div['wj_signed'], 4),
                'gap': round(div['gap'], 4),
                'sign_inversion_pct': round(div['sign_inversion_pct'], 1),
            })

    if not rows:
        log("No comparisons completed.")
        return

    df = pd.DataFrame(rows)
    out_path = os.path.join(RESULTS_DIR, 'quality_filter_effect.csv')
    df.to_csv(out_path, index=False)
    log(f"\nSaved: {out_path}")

    # Figure: before/after QC for each comparison
    pairs_done = df[['condition_a', 'condition_b', 'comparison_type']].drop_duplicates()
    n_pairs = len(pairs_done)

    fig, axes = plt.subplots(1, 2, figsize=(12, max(5, 3 * n_pairs)))
    bar_w = 0.35
    x = np.arange(n_pairs)
    labels_short = []

    wj_fr = []
    wj_qc = []
    si_fr = []
    si_qc = []

    for _, row in pairs_done.iterrows():
        sub = df[(df['condition_a'] == row['condition_a']) &
                 (df['condition_b'] == row['condition_b'])]
        fr_row = sub[sub['filter'] == 'FR_only']
        qc_row = sub[sub['filter'] == 'FR+AllenQC']
        lbl = f"{row['comparison_type']}\n{row['condition_a'].split('_')[-1]} vs {row['condition_b'].split('_')[-1]}"
        labels_short.append(lbl)
        wj_fr.append(float(fr_row['wj_unsigned'].iloc[0]) if len(fr_row) else np.nan)
        wj_qc.append(float(qc_row['wj_unsigned'].iloc[0]) if len(qc_row) else np.nan)
        si_fr.append(float(fr_row['sign_inversion_pct'].iloc[0]) if len(fr_row) else np.nan)
        si_qc.append(float(qc_row['sign_inversion_pct'].iloc[0]) if len(qc_row) else np.nan)

    ax = axes[0]
    ax.bar(x - bar_w/2, wj_fr, bar_w, label='FR only (current)', color='#1f77b4', alpha=0.8)
    ax.bar(x + bar_w/2, wj_qc, bar_w, label='FR + Allen QC', color='#2ca02c', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_short, fontsize=8)
    ax.set_ylabel('WJ (unsigned)', fontsize=12)
    ax.set_title('WJ Unsigned\nBefore vs After Quality Filter', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    ax = axes[1]
    ax.bar(x - bar_w/2, si_fr, bar_w, label='FR only (current)', color='#1f77b4', alpha=0.8)
    ax.bar(x + bar_w/2, si_qc, bar_w, label='FR + Allen QC', color='#2ca02c', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_short, fontsize=8)
    ax.set_ylabel('Sign-Inversion % Metric', fontsize=12)
    ax.set_title('Sign-Inversion Metric\nBefore vs After Quality Filter', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    plt.suptitle(f'Quality Filter Effect on WJ Analysis\n'
                 f'(Allen QC: quality=good, ISI<{QC_ISI_MAX}, '
                 f'presence>={QC_PRESENCE_MIN}, amplitude_cutoff<{QC_AMPLITUDE_CUTOFF_MAX})',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, 'quality_filter_effect.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    log(f"Saved: {fig_path}")

    log("\n" + "=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log(df[['comparison_type', 'filter', 'n_units_variable',
            'wj_unsigned', 'sign_inversion_pct']].to_string(index=False))
    log(f"\nTotal time: {(time.time()-START)/60:.1f}m")


if __name__ == '__main__':
    main()
