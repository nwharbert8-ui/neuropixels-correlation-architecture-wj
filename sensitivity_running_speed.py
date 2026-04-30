"""
Pipeline: Running Speed Confound Analysis — Per-Stimulus Condition
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Addresses Kill Shot #4: Running speed modulates V1 correlation structure
    (Niell & Stryker 2010, McGinley 2015, Vinck 2015). If animals run at
    systematically different speeds during different stimulus conditions, then
    condition differences in WJ could reflect locomotion state rather than
    stimulus-driven reorganization.

    Extracts running speed from Allen NWB files during each stimulus window.
    Reports mean ± SD per condition and tests for significant differences
    across conditions (Kruskal-Wallis). Computes Spearman correlation between
    mean running speed and WJ values across comparisons.

Dependencies: numpy, scipy, pandas, matplotlib, h5py
Input:
    G:/My Drive/inner_architecture_research/neuropixels_wj/data/*.nwb
    G:/My Drive/inner_architecture_research/neuropixels_wj/results/ses-715093703_checkpoint.json
Output:
    results/running_speed_by_condition.csv
    results/running_speed_wj_correlation.csv
    figures/running_speed_analysis.png
"""
import os
import sys
import time
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr

warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import h5py

# ============================================================================
# CONFIGURATION
# ============================================================================
RANDOM_SEED = 42
FORCE_RECOMPUTE = True

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

np.random.seed(RANDOM_SEED)
START = time.time()
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def log(msg):
    print(f"[{(time.time()-START)/60:5.1f}m] {msg}", flush=True)


# ============================================================================
# RUNNING SPEED EXTRACTION
# ============================================================================
def extract_running_speed(nwb_path):
    """
    Extract running speed time series from NWB file.
    Path: processing/running/running_speed/data + timestamps
    Returns (timestamps, speed_cm_per_s) arrays.
    """
    with h5py.File(nwb_path, 'r') as f:
        speed = f['processing/running/running_speed/data'][:]
        times = f['processing/running/running_speed/timestamps'][:]
    return times, speed


def extract_stimulus_conditions(nwb_path):
    """Extract stimulus condition intervals from NWB file."""
    stim_info = {}
    with h5py.File(nwb_path, 'r') as f:
        if 'intervals' not in f:
            return stim_info
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
                            'n_presentations': mask.sum(),
                        }
            except Exception as e:
                log(f"  Warning {ik}: {e}")
    return stim_info


def mean_speed_in_windows(speed_times, speed_vals, starts, stops):
    """Compute mean running speed during each stimulus window."""
    window_means = []
    for t0, t1 in zip(starts, stops):
        mask = (speed_times >= t0) & (speed_times <= t1)
        if mask.sum() >= 1:
            window_means.append(float(np.mean(speed_vals[mask])))
    return np.array(window_means)


# ============================================================================
# MAIN
# ============================================================================
def main():
    log("=" * 70)
    log("RUNNING SPEED CONFOUND ANALYSIS")
    log("=" * 70)

    nwb_files = sorted([f for f in os.listdir(DATA_DIR) if f.endswith('.nwb')])
    log(f"Found {len(nwb_files)} NWB files")

    all_condition_rows = []
    wj_corr_rows = []

    for nwb_file in nwb_files:
        nwb_path = os.path.join(DATA_DIR, nwb_file)
        session_id = nwb_file.replace('.nwb', '').split('_')[-1]
        log(f"\n--- Session {session_id} ---")

        try:
            speed_times, speed_vals = extract_running_speed(nwb_path)
        except Exception as e:
            log(f"  ERROR extracting running speed: {e}")
            continue

        log(f"  Running speed: {len(speed_vals)} samples, "
            f"range {speed_vals.min():.1f} - {speed_vals.max():.1f} cm/s, "
            f"mean={speed_vals.mean():.1f}")

        stim_info = extract_stimulus_conditions(nwb_path)
        log(f"  Stimulus conditions: {len(stim_info)}")

        condition_data = {}
        for cond_key, info in stim_info.items():
            means = mean_speed_in_windows(speed_times, speed_vals,
                                          info['starts'], info['stops'])
            if len(means) < 2:
                continue
            condition_data[cond_key] = means
            stim_class = cond_key.split('_presentations_')[0] if '_presentations_' in cond_key else cond_key

            all_condition_rows.append({
                'session': session_id,
                'condition': cond_key,
                'stim_class': stim_class,
                'n_presentations': len(means),
                'mean_speed_cm_s': round(float(np.mean(means)), 3),
                'std_speed_cm_s': round(float(np.std(means)), 3),
                'median_speed_cm_s': round(float(np.median(means)), 3),
                'pct_running': round(float(np.mean(means > 5.0)) * 100, 1),
            })

        # Kruskal-Wallis test across conditions for this session
        if len(condition_data) >= 2:
            groups = list(condition_data.values())
            try:
                stat, p = kruskal(*groups)
                log(f"  Kruskal-Wallis across {len(groups)} conditions: "
                    f"H={stat:.2f}, p={p:.4e}")
                if p < 0.05:
                    log("  WARNING: Running speed differs significantly across conditions")
                else:
                    log("  Running speed does not significantly differ across conditions")
            except Exception as e:
                log(f"  KW test error: {e}")

        # Load checkpoint WJ values to correlate with running speed
        ckpt_file = os.path.join(RESULTS_DIR, f'{session_id}_checkpoint.json')
        if os.path.exists(ckpt_file):
            log("  Loading WJ checkpoint for correlation analysis...")
            with open(ckpt_file) as f_ckpt:
                ckpt = json.load(f_ckpt)

            for r in ckpt.get('results', []):
                cond_a = r['condition_a']
                cond_b = r['condition_b']
                if cond_a in condition_data and cond_b in condition_data:
                    speed_a = float(np.mean(condition_data[cond_a]))
                    speed_b = float(np.mean(condition_data[cond_b]))
                    delta_speed = abs(speed_a - speed_b)
                    stim_a = cond_a.split('_presentations_')[0] if '_presentations_' in cond_a else cond_a
                    stim_b = cond_b.split('_presentations_')[0] if '_presentations_' in cond_b else cond_b
                    wj_corr_rows.append({
                        'session': session_id,
                        'condition_a': cond_a,
                        'condition_b': cond_b,
                        'comparison_type': 'same' if stim_a == stim_b else 'diff',
                        'wj_unsigned': r['wj_unsigned'],
                        'sign_inversion_pct': r.get('sign_inversion_pct', r.get('sign_inv_pct', np.nan)),
                        'mean_speed_a': round(speed_a, 3),
                        'mean_speed_b': round(speed_b, 3),
                        'delta_speed': round(delta_speed, 3),
                    })

    # Save condition-level summary
    df_cond = pd.DataFrame(all_condition_rows)
    if not df_cond.empty:
        out_cond = os.path.join(RESULTS_DIR, 'running_speed_by_condition.csv')
        df_cond.to_csv(out_cond, index=False)
        log(f"\nSaved: {out_cond}")

        log("\nRunning speed by stimulus class (aggregated across sessions):")
        log(df_cond.groupby('stim_class')[['mean_speed_cm_s', 'pct_running']].mean().to_string())

    # WJ vs running speed correlation
    if wj_corr_rows:
        df_wj = pd.DataFrame(wj_corr_rows)
        out_wj = os.path.join(RESULTS_DIR, 'running_speed_wj_correlation.csv')
        df_wj.to_csv(out_wj, index=False)
        log(f"\nSaved: {out_wj}")

        # Correlation between delta_speed and WJ
        if len(df_wj) >= 10:
            rho, p = spearmanr(df_wj['delta_speed'], df_wj['wj_unsigned'])
            log(f"\nSpearman: delta_speed vs WJ_unsigned: rho={rho:.3f}, p={p:.4e}")
            rho2, p2 = spearmanr(df_wj['delta_speed'], df_wj['sign_inversion_pct'])
            log(f"Spearman: delta_speed vs sign_inv_pct: rho={rho2:.3f}, p={p2:.4e}")

            if abs(rho) < 0.3 and p > 0.05:
                log("RESULT: Running speed delta does NOT significantly predict WJ — "
                    "locomotion is not a confound")
            else:
                log("RESULT: Running speed delta correlates with WJ — "
                    "consider adding running speed as covariate")

    # Figure
    if not df_cond.empty:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Panel A: running speed distribution per stim class
        ax = axes[0]
        stim_classes = sorted(df_cond['stim_class'].unique())
        means_by_class = [df_cond[df_cond['stim_class'] == sc]['mean_speed_cm_s'].values
                          for sc in stim_classes]
        ax.boxplot(means_by_class, labels=[sc[:20] for sc in stim_classes],
                   patch_artist=True, notch=False)
        ax.set_ylabel('Mean Running Speed (cm/s)', fontsize=12)
        ax.set_title('Running Speed per Stimulus Class\n(presentations mean across sessions)', fontsize=11)
        ax.tick_params(axis='x', rotation=30)
        ax.grid(True, alpha=0.3, axis='y')

        # Panel B: delta_speed vs WJ
        if wj_corr_rows:
            ax = axes[1]
            df_wj_plot = pd.DataFrame(wj_corr_rows)
            colors = {'same': '#1f77b4', 'diff': '#d62728'}
            for ctype, grp in df_wj_plot.groupby('comparison_type'):
                ax.scatter(grp['delta_speed'], grp['wj_unsigned'],
                           c=colors.get(ctype, 'gray'), label=f'{ctype}-stimulus',
                           alpha=0.5, s=20)
            ax.set_xlabel('|ΔRunning Speed| (cm/s)', fontsize=12)
            ax.set_ylabel('WJ (unsigned)', fontsize=12)
            ax.set_title('Running Speed Delta vs WJ\n(confound assessment)', fontsize=11)
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.3)
        else:
            axes[1].text(0.5, 0.5, 'No WJ checkpoint data', ha='center', va='center',
                         transform=axes[1].transAxes)

        plt.suptitle('Running Speed Analysis\n(Niell & Stryker 2010 / McGinley 2015 confound check)',
                     fontsize=13, fontweight='bold', y=1.02)
        plt.tight_layout()
        fig_path = os.path.join(FIGURES_DIR, 'running_speed_analysis.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        log(f"Saved: {fig_path}")

    log(f"\nTotal time: {(time.time()-START)/60:.1f}m")
    log("COMPLETE")


if __name__ == '__main__':
    main()
