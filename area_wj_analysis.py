"""
Pipeline: Neuropixels Per-Area WJ Aggregate Analysis
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Aggregates per-area WJ values from all session checkpoints.
    Computes same-stim vs diff-stim WJ by brain area.
    Outputs summary CSV and publication-grade figures.

Dependencies: numpy, scipy, pandas, matplotlib, seaborn
Input: checkpoint JSONs in results/ (with area_wj populated)
Output: area_wj_summary.csv, area_wj_figure*.png in results/
"""

import os
import json
import glob
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIG
# ============================================================================
RANDOM_SEED = 42
BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
MIN_AREA_OBSERVATIONS = 10  # min number of comparison pairs per area for inclusion

np.random.seed(RANDOM_SEED)
START = time.time()


def log(msg):
    print(f"[{(time.time()-START)/60:6.1f}m] {msg}", flush=True)


def stim_type(condition_name):
    """Extract stimulus type by stripping trailing _<value> suffix."""
    parts = condition_name.rsplit('_', 1)
    return parts[0] if len(parts) == 2 else condition_name


def load_all_checkpoints(results_dir):
    """Load all session checkpoints and return flat list of result records."""
    checkpoints = sorted(glob.glob(os.path.join(results_dir, '*_checkpoint.json')))
    log(f"Loading {len(checkpoints)} checkpoints...")

    records = []
    sessions_with_area = 0
    sessions_skipped = 0

    for ckpt_path in checkpoints:
        session_id = os.path.basename(ckpt_path).replace('_checkpoint.json', '')
        with open(ckpt_path) as f:
            ckpt = json.load(f)

        results = ckpt.get('results', [])
        if not results:
            continue

        # Check area_wj populated
        sample = results[0].get('area_wj', {})
        if not sample or list(sample.keys()) == ['unknown']:
            sessions_skipped += 1
            continue

        sessions_with_area += 1
        for r in results:
            area_wj = r.get('area_wj', {})
            if not area_wj:
                continue
            cond_a = r.get('condition_a', '')
            cond_b = r.get('condition_b', '')
            # Derive same_stim from condition names if not stored
            same_stim = r.get('same_stim')
            if same_stim is None:
                same_stim = (stim_type(cond_a) == stim_type(cond_b))
            base = {
                'session_id': session_id,
                'condition_a': cond_a,
                'condition_b': cond_b,
                'same_stim': bool(same_stim),
                'wj_global': r.get('wj_unsigned', np.nan),
                'wj_signed_global': r.get('wj_signed', np.nan),
                'sign_inv_pct': r.get('sign_inversion_pct', np.nan),
                'p_value': r.get('perm_p', np.nan),
            }
            for area, area_data in area_wj.items():
                rec = dict(base)
                rec['area'] = area
                # area_wj values are dicts with wj_unsigned, wj_signed, etc.
                if isinstance(area_data, dict):
                    rec['wj_area'] = area_data.get('wj_unsigned', np.nan)
                    rec['wj_area_signed'] = area_data.get('wj_signed', np.nan)
                    rec['sign_inv_pct_area'] = area_data.get('sign_inversion_pct', np.nan)
                else:
                    rec['wj_area'] = float(area_data)
                    rec['wj_area_signed'] = np.nan
                    rec['sign_inv_pct_area'] = np.nan
                records.append(rec)

    log(f"Sessions with area data: {sessions_with_area} | Skipped (no area): {sessions_skipped}")
    log(f"Total records: {len(records)}")
    return pd.DataFrame(records)


def compute_area_stats(df):
    """Per-area same-stim vs diff-stim WJ statistics."""
    rows = []
    areas = df['area'].unique()

    for area in sorted(areas):
        if area == 'unknown' or area == '':
            continue
        sub = df[df['area'] == area]

        same = sub[sub['same_stim'] == True]['wj_area'].dropna()
        diff = sub[sub['same_stim'] == False]['wj_area'].dropna()

        if len(same) < MIN_AREA_OBSERVATIONS or len(diff) < MIN_AREA_OBSERVATIONS:
            continue

        t, p = stats.ttest_ind(same, diff)
        d = (same.mean() - diff.mean()) / np.sqrt(
            ((len(same)-1)*same.std()**2 + (len(diff)-1)*diff.std()**2) /
            (len(same) + len(diff) - 2)
        )

        rows.append({
            'area': area,
            'n_same': len(same),
            'n_diff': len(diff),
            'wj_same_mean': same.mean(),
            'wj_same_std': same.std(),
            'wj_diff_mean': diff.mean(),
            'wj_diff_std': diff.std(),
            'delta_wj': same.mean() - diff.mean(),
            't_stat': t,
            'p_raw': p,
            'cohens_d': d,
        })

    result = pd.DataFrame(rows)
    if len(result) == 0:
        return result

    _, p_fdr, _, _ = multipletests(result['p_raw'], method='fdr_bh')
    result['p_fdr'] = p_fdr
    result = result.sort_values('wj_same_mean', ascending=False)
    return result


def figure_area_wj_comparison(area_stats, out_path):
    """Bar chart: same-stim vs diff-stim WJ by area."""
    if len(area_stats) == 0:
        log("  No area stats — skipping figure")
        return

    sig = area_stats[area_stats['p_fdr'] < 0.05]
    n_areas = len(area_stats)

    fig, ax = plt.subplots(figsize=(max(10, n_areas * 0.6), 6))
    x = np.arange(n_areas)
    w = 0.35

    palette = sns.color_palette('colorblind')
    bars_same = ax.bar(x - w/2, area_stats['wj_same_mean'], w,
                       yerr=area_stats['wj_same_std'] / np.sqrt(area_stats['n_same']),
                       label='Same stimulus', color=palette[0], capsize=3)
    bars_diff = ax.bar(x + w/2, area_stats['wj_diff_mean'], w,
                       yerr=area_stats['wj_diff_std'] / np.sqrt(area_stats['n_diff']),
                       label='Different stimulus', color=palette[1], capsize=3)

    # Mark FDR-significant areas
    for i, row in area_stats.reset_index(drop=True).iterrows():
        if row['p_fdr'] < 0.05:
            y_max = max(row['wj_same_mean'], row['wj_diff_mean']) + 0.02
            ax.text(i, y_max, '*', ha='center', va='bottom', fontsize=14,
                    fontweight='bold', color='black')

    ax.set_xticks(x)
    ax.set_xticklabels(area_stats['area'], rotation=45, ha='right', fontsize=11)
    ax.set_ylabel('WJ Implementation Divergence', fontsize=12)
    ax.set_xlabel('Brain Area', fontsize=12)
    ax.set_title('Per-Area WJ: Same vs. Different Stimulus Conditions\n'
                 '(* = FDR q < 0.05; error bars = SEM)', fontsize=13)
    ax.legend(fontsize=11)
    ax.set_ylim(bottom=0)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {os.path.basename(out_path)}")


def figure_area_violin(df, area_stats, out_path):
    """Violin plot of WJ distribution by area, split same/diff."""
    if len(area_stats) == 0:
        return

    # Only areas in area_stats (pass MIN_AREA_OBSERVATIONS)
    valid_areas = list(area_stats['area'])
    sub = df[df['area'].isin(valid_areas)].copy()
    sub['condition'] = sub['same_stim'].map({True: 'Same', False: 'Different'})

    # Sort by same-stim WJ mean
    order = list(area_stats['area'])

    fig, ax = plt.subplots(figsize=(max(12, len(valid_areas) * 0.7), 6))
    palette = {'Same': sns.color_palette('colorblind')[0],
               'Different': sns.color_palette('colorblind')[1]}
    sns.violinplot(data=sub, x='area', y='wj_area', hue='condition',
                   order=order, palette=palette, split=False, inner='quartile',
                   ax=ax, cut=0)
    ax.set_xlabel('Brain Area', fontsize=12)
    ax.set_ylabel('WJ Implementation Divergence', fontsize=12)
    ax.set_title('WJ Distribution by Brain Area and Stimulus Condition', fontsize=13)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    ax.legend(title='Condition', fontsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {os.path.basename(out_path)}")


def figure_delta_wj_ranked(area_stats, out_path):
    """Ranked delta WJ (same - diff) across areas."""
    if len(area_stats) == 0:
        return

    ranked = area_stats.sort_values('delta_wj', ascending=True).reset_index(drop=True)
    colors = ['#d73027' if row['p_fdr'] < 0.05 else '#636363'
              for _, row in ranked.iterrows()]

    fig, ax = plt.subplots(figsize=(8, max(5, len(ranked) * 0.35)))
    bars = ax.barh(ranked['area'], ranked['delta_wj'], color=colors)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel('ΔWJ (Same − Different stimulus)', fontsize=12)
    ax.set_title('WJ Reorganization by Brain Area\n(red = FDR q < 0.05)', fontsize=13)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {os.path.basename(out_path)}")


def main():
    log("=" * 70)
    log("AREA WJ AGGREGATE ANALYSIS")
    log("=" * 70)

    # Load all checkpoint data
    df = load_all_checkpoints(RESULTS_DIR)
    if df.empty:
        log("ERROR: No records loaded — run area_decomposition.py first")
        return

    # Overall area coverage
    area_counts = df.groupby('area')['session_id'].nunique().sort_values(ascending=False)
    log(f"\nArea coverage (sessions with data):")
    for area, n in area_counts.items():
        if area and area != 'unknown':
            n_pairs = len(df[df['area'] == area])
            log(f"  {area:12s}: {n:2d} sessions, {n_pairs:5d} comparison pairs")

    # Compute per-area stats
    log("\nComputing per-area same vs. diff-stim statistics...")
    area_stats = compute_area_stats(df)

    if area_stats.empty:
        log("No areas passed minimum observation threshold")
        return

    # Save CSV
    csv_path = os.path.join(RESULTS_DIR, 'area_wj_summary.csv')
    area_stats.to_csv(csv_path, index=False, float_format='%.6f')
    log(f"\nSaved: area_wj_summary.csv ({len(area_stats)} areas)")

    # Print summary table
    log("\nArea WJ Summary (sorted by same-stim WJ mean):")
    log(f"{'Area':12s} {'n_same':>7} {'n_diff':>7} {'WJ_same':>8} {'WJ_diff':>8} {'Delta':>8} {'Cohen_d':>8} {'p_FDR':>10}")
    log("-" * 75)
    for _, row in area_stats.iterrows():
        sig = '*' if row['p_fdr'] < 0.05 else ' '
        log(f"{row['area']:12s} {int(row['n_same']):>7} {int(row['n_diff']):>7} "
            f"{row['wj_same_mean']:>8.4f} {row['wj_diff_mean']:>8.4f} "
            f"{row['delta_wj']:>8.4f} {row['cohens_d']:>8.3f} {row['p_fdr']:>10.4f}{sig}")

    # Figures
    log("\nGenerating figures...")
    figure_area_wj_comparison(area_stats,
        os.path.join(RESULTS_DIR, 'area_wj_bar_comparison.png'))
    figure_area_violin(df, area_stats,
        os.path.join(RESULTS_DIR, 'area_wj_violin.png'))
    figure_delta_wj_ranked(area_stats,
        os.path.join(RESULTS_DIR, 'area_wj_delta_ranked.png'))

    elapsed = (time.time() - START) / 60
    log(f"\n{'='*70}")
    log(f"ANALYSIS COMPLETE — {elapsed:.1f} minutes")
    log(f"Outputs in: {RESULTS_DIR}")


if __name__ == '__main__':
    main()
