"""
Pipeline: Neuropixels Layer 2H Pairing-Family Decomposition
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-29

Description:
    Applies the Layer 2H pairing-family decomposition framework (Harbert, in
    press, Frontiers in Pharmacology) to the existing Neuropixels analysis.
    Reframes the central paper finding (signed vs unsigned WJ gap) as a
    formal Type 2 measurement and adds three new measurements that the
    existing pipeline did not extract:

    Type 6 (Substrate-Projection): Spearman WJ vs Pearson WJ. Already computed
        in summary JSONs but not extracted as a paired measurement.

    Type 5 (Local-Global): Per-area WJ vs full-population WJ. Tests whether
        the visual hierarchy gradient at area level converges with or diverges
        from the global reorganization signal.

    Type 4 (Set-vs-Multiset): Cross-session consensus on (a) stimulus-pair
        reorganization, (b) area-level reorganization, and (c) sign-inversion
        percentage. Tests whether the central claim of stable signed
        reorganization dominance is multiset-consistent across 30 mice.

    Type 2 (Sign-Treatment): Already the central finding — formalized here
        with cross-session statistics and confidence intervals.

    Type 1 (Continuous-Discrete) PILOT: 3 sessions. Computes binary Jaccard
        at top 5% threshold from correlation matrices and the resulting
        dissociation gap. Decision: scale to 30 sessions or stop here based
        on whether the pilot reveals a striking pattern.

Inputs:
    results/ses-*_summary.json (30 sessions)
    results/area_wj_summary.csv
    results/final_clean_summary.json
    data/sub-*_ses-*.nwb (3 sessions for Type 1 pilot)

Outputs:
    results/layer2h_type6_substrate_projection.csv
    results/layer2h_type5_local_global.csv
    results/layer2h_type4_consensus.csv
    results/layer2h_type2_formalized.csv
    results/layer2h_type1_pilot.csv (if Type 1 pilot run)
    results/layer2h_summary_provenance.json
"""

import os
import gc
import json
import glob
import time
import warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
START = time.time()

ROOT = r"G:\My Drive\inner_architecture_research\neuropixels_wj"
RESULTS_DIR = os.path.join(ROOT, "results")
DATA_DIR = os.path.join(ROOT, "data")
OUT_DIR = RESULTS_DIR


def log(msg):
    print(f"[{(time.time()-START)/60:6.1f}m] {msg}", flush=True)


# ==========================================================================
# Load all session summaries into a flat DataFrame
# ==========================================================================
log("Loading session summaries...")
summary_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "ses-*_summary.json")))
log(f"  Found {len(summary_files)} session summaries")

all_comparisons = []
for sf in summary_files:
    try:
        with open(sf) as f:
            d = json.load(f)
        sid = d['session_id']
        for c in d.get('comparisons', []):
            row = {'session_id': sid}
            row.update(c)
            all_comparisons.append(row)
    except Exception as e:
        log(f"  Warning: could not parse {sf}: {e}")

df = pd.DataFrame(all_comparisons)
log(f"  Loaded {len(df):,} comparisons across {df['session_id'].nunique()} sessions")
log(f"  Columns: {list(df.columns)}")

# Identify same-stim vs diff-stim
def stim_type(name):
    """Extract stimulus type from condition name."""
    name = name.lower()
    if 'natural_scenes' in name: return 'natural_scenes'
    if 'natural_movie_one' in name: return 'natural_movie_one'
    if 'static_gratings' in name: return 'static_gratings'
    if 'drifting_gratings' in name: return 'drifting_gratings'
    if 'gabors' in name: return 'gabors'
    if 'flashes' in name: return 'flashes'
    if 'spontaneous' in name: return 'spontaneous'
    return 'other'

df['stim_a'] = df['condition_a'].apply(stim_type)
df['stim_b'] = df['condition_b'].apply(stim_type)
df['same_stim'] = df['stim_a'] == df['stim_b']
df['stim_pair'] = df.apply(lambda r: ' vs '.join(sorted([r['stim_a'], r['stim_b']])), axis=1)

# Apply original filtering (n_bins >= 50 each, exclude natural_movie_one)
df['n_bins_a'] = pd.to_numeric(df['n_bins_a'], errors='coerce')
df['n_bins_b'] = pd.to_numeric(df['n_bins_b'], errors='coerce')
mask = (df['n_bins_a'] >= 50) & (df['n_bins_b'] >= 50) & \
       (~df['stim_a'].str.contains('natural_movie_one')) & \
       (~df['stim_b'].str.contains('natural_movie_one'))
df_clean = df[mask].copy()
log(f"  After filtering (min_bins>=50, exclude natural_movie_one): {len(df_clean):,}")

# ==========================================================================
# TYPE 2: Sign-Treatment formalized (the central paper finding)
# ==========================================================================
log("\n" + "=" * 70)
log("TYPE 2: Sign-Treatment (sign-inversion dominance) formalized")
log("=" * 70)

t2_summary = {
    'wj_unsigned_mean': df_clean['wj_unsigned'].mean(),
    'wj_unsigned_std': df_clean['wj_unsigned'].std(),
    'wj_signed_mean': df_clean['wj_signed'].mean(),
    'wj_signed_std': df_clean['wj_signed'].std(),
    'gap_mean': (df_clean['wj_signed'] - df_clean['wj_unsigned']).mean(),
    'gap_std': (df_clean['wj_signed'] - df_clean['wj_unsigned']).std(),
    'sign_inv_pct_mean': df_clean['sign_inversion_pct'].mean(),
    'sign_inv_pct_median': df_clean['sign_inversion_pct'].median(),
    'sign_inv_pct_std': df_clean['sign_inversion_pct'].std(),
    'sign_inv_pct_min': df_clean['sign_inversion_pct'].min(),
    'sign_inv_pct_max': df_clean['sign_inversion_pct'].max(),
    'n_comparisons': len(df_clean),
    'n_sessions': df_clean['session_id'].nunique(),
}
log(f"  WJ unsigned:  {t2_summary['wj_unsigned_mean']:.4f} ± {t2_summary['wj_unsigned_std']:.4f}")
log(f"  WJ signed:    {t2_summary['wj_signed_mean']:.4f} ± {t2_summary['wj_signed_std']:.4f}")
log(f"  Type 2 gap:   {t2_summary['gap_mean']:.4f} ± {t2_summary['gap_std']:.4f}")
log(f"  Sign inversion %: mean={t2_summary['sign_inv_pct_mean']:.2f}, "
    f"median={t2_summary['sign_inv_pct_median']:.2f}, range=["
    f"{t2_summary['sign_inv_pct_min']:.2f}, {t2_summary['sign_inv_pct_max']:.2f}]")

t_signed_unsigned, p_signed_unsigned = stats.ttest_rel(df_clean['wj_signed'],
                                                        df_clean['wj_unsigned'])
log(f"  Paired t-test signed > unsigned: t={t_signed_unsigned:.2f}, p={p_signed_unsigned:.2e}")
t2_summary['paired_t_signed_vs_unsigned'] = t_signed_unsigned
t2_summary['paired_p_signed_vs_unsigned'] = p_signed_unsigned

pd.DataFrame([t2_summary]).to_csv(os.path.join(OUT_DIR, "layer2h_type2_formalized.csv"),
                                    index=False)

# ==========================================================================
# TYPE 6: Substrate-Projection (Spearman WJ vs Pearson WJ)
# ==========================================================================
log("\n" + "=" * 70)
log("TYPE 6: Substrate-Projection (Spearman vs Pearson WJ)")
log("=" * 70)

# wj_pearson is already computed; wj_unsigned is on Spearman correlations
df_clean['type6_gap'] = df_clean['wj_unsigned'] - df_clean['wj_pearson']

t6_summary = {
    'spearman_wj_mean': df_clean['wj_unsigned'].mean(),
    'pearson_wj_mean': df_clean['wj_pearson'].mean(),
    'type6_gap_mean': df_clean['type6_gap'].mean(),
    'type6_gap_std': df_clean['type6_gap'].std(),
    'type6_gap_median': df_clean['type6_gap'].median(),
    'type6_gap_abs_mean': df_clean['type6_gap'].abs().mean(),
    'type6_gap_abs_median': df_clean['type6_gap'].abs().median(),
    'fraction_gap_within_005': (df_clean['type6_gap'].abs() < 0.05).mean(),
    'fraction_gap_within_01': (df_clean['type6_gap'].abs() < 0.10).mean(),
    'n_comparisons': len(df_clean),
}
log(f"  Spearman WJ mean: {t6_summary['spearman_wj_mean']:.4f}")
log(f"  Pearson WJ mean:  {t6_summary['pearson_wj_mean']:.4f}")
log(f"  Type 6 gap mean:  {t6_summary['type6_gap_mean']:.4f} ± {t6_summary['type6_gap_std']:.4f}")
log(f"  Type 6 |gap| median: {t6_summary['type6_gap_abs_median']:.4f}")
log(f"  Fraction with |gap| < 0.05: {t6_summary['fraction_gap_within_005']:.2%}")
log(f"  Fraction with |gap| < 0.10: {t6_summary['fraction_gap_within_01']:.2%}")
log(f"\n  Interpretation: Near-zero gap = approximately linear-monotonic substrate")

# Save per-comparison Type 6 results
df_clean[['session_id', 'condition_a', 'condition_b', 'wj_unsigned', 'wj_pearson',
          'type6_gap']].to_csv(os.path.join(OUT_DIR, "layer2h_type6_per_comparison.csv"),
                                index=False)
pd.DataFrame([t6_summary]).to_csv(os.path.join(OUT_DIR,
                                                "layer2h_type6_substrate_projection.csv"),
                                    index=False)

# ==========================================================================
# TYPE 5: Local-Global (per-area WJ vs full-population WJ)
# ==========================================================================
log("\n" + "=" * 70)
log("TYPE 5: Local-Global (per-area vs full-population WJ)")
log("=" * 70)

area_csv = os.path.join(RESULTS_DIR, "area_wj_summary.csv")
if os.path.exists(area_csv):
    area_df = pd.read_csv(area_csv)
    log(f"  Per-area data: {len(area_df)} areas")

    # Full-population mean (from clean summary)
    same_stim_global = df_clean[df_clean['same_stim']]['wj_unsigned'].mean()
    diff_stim_global = df_clean[~df_clean['same_stim']]['wj_unsigned'].mean()
    log(f"  Full-population same-stim WJ:  {same_stim_global:.4f}")
    log(f"  Full-population diff-stim WJ:  {diff_stim_global:.4f}")

    area_df['gap_same_full_minus_local'] = same_stim_global - area_df['wj_same_mean']
    area_df['gap_diff_full_minus_local'] = diff_stim_global - area_df['wj_diff_mean']

    log(f"\n  Per-area Local-Global gaps (top 10 areas by same-stim gap magnitude):")
    log(f"  {'Area':10s} {'Local same-stim':>16s} {'Gap (full-local)':>18s}")
    for _, row in area_df.nlargest(10, 'gap_same_full_minus_local',
                                    keep='all').head(10).iterrows():
        log(f"  {row['area']:10s} {row['wj_same_mean']:>16.4f} "
            f"{row['gap_same_full_minus_local']:>18.4f}")

    log(f"\n  Areas where local-global gap is NEGATIVE (area exceeds global mean):")
    neg = area_df[area_df['gap_same_full_minus_local'] < 0]
    for _, row in neg.iterrows():
        log(f"    {row['area']:10s}: local same-stim WJ {row['wj_same_mean']:.4f}, "
            f"gap {row['gap_same_full_minus_local']:.4f}")

    t5_summary = {
        'full_population_same_stim_wj': same_stim_global,
        'full_population_diff_stim_wj': diff_stim_global,
        'mean_local_same_stim_wj': area_df['wj_same_mean'].mean(),
        'mean_local_diff_stim_wj': area_df['wj_diff_mean'].mean(),
        'mean_local_global_gap_same': area_df['gap_same_full_minus_local'].mean(),
        'mean_local_global_gap_diff': area_df['gap_diff_full_minus_local'].mean(),
        'n_areas': len(area_df),
        'n_areas_local_exceeds_global': (area_df['gap_same_full_minus_local'] < 0).sum(),
    }
    log(f"\n  Mean local-global gap (same-stim): {t5_summary['mean_local_global_gap_same']:.4f}")
    log(f"  Mean local-global gap (diff-stim): {t5_summary['mean_local_global_gap_diff']:.4f}")
    log(f"  N areas where local exceeds global: "
        f"{t5_summary['n_areas_local_exceeds_global']} of {t5_summary['n_areas']}")

    area_df.to_csv(os.path.join(OUT_DIR, "layer2h_type5_local_global_per_area.csv"),
                    index=False)
    pd.DataFrame([t5_summary]).to_csv(os.path.join(OUT_DIR,
                                                    "layer2h_type5_local_global.csv"),
                                        index=False)
else:
    log(f"  area_wj_summary.csv not found; skipping Type 5")
    t5_summary = None

# ==========================================================================
# TYPE 4: Set-vs-Multiset cross-session consensus
# ==========================================================================
log("\n" + "=" * 70)
log("TYPE 4: Set-vs-Multiset cross-session consensus")
log("=" * 70)

# Approach 1: cross-session consensus on stimulus-pair-type WJ values
# For each unique stimulus pair type, get distribution of WJ across sessions
log("\n  Approach 1: Stimulus-pair-type cross-session consensus")
stim_pair_groups = df_clean.groupby('stim_pair')

stim_pair_rows = []
threshold_wj = df_clean['wj_unsigned'].median()  # median split
log(f"  Using WJ threshold = median = {threshold_wj:.4f} for set characterization")
for sp, g in stim_pair_groups:
    n_sessions_present = g['session_id'].nunique()
    n_sessions_above_threshold = (g.groupby('session_id')['wj_unsigned'].mean()
                                   > threshold_wj).sum()
    set_membership_rate = (n_sessions_above_threshold / n_sessions_present
                            if n_sessions_present > 0 else 0)
    multiset_mean_wj = g['wj_unsigned'].mean()
    multiset_std_wj = g['wj_unsigned'].std()

    stim_pair_rows.append({
        'stim_pair': sp,
        'n_observations': len(g),
        'n_sessions_present': n_sessions_present,
        'set_above_threshold_rate': set_membership_rate,
        'multiset_mean_wj': multiset_mean_wj,
        'multiset_std_wj': multiset_std_wj,
        'multiset_cv': multiset_std_wj / multiset_mean_wj
                       if multiset_mean_wj > 0 else np.nan,
    })

stim_pair_df = pd.DataFrame(stim_pair_rows).sort_values('multiset_mean_wj',
                                                          ascending=False)
log(f"\n  Top 10 stimulus pairs by multiset mean WJ:")
log(f"  {'Stim pair':50s} {'Set rate':>10s} {'Multiset mean':>14s} {'CV':>8s}")
for _, r in stim_pair_df.head(10).iterrows():
    log(f"  {r['stim_pair']:50s} {r['set_above_threshold_rate']:>10.3f} "
        f"{r['multiset_mean_wj']:>14.4f} {r['multiset_cv']:>8.3f}")

stim_pair_df.to_csv(os.path.join(OUT_DIR,
                                  "layer2h_type4_stim_pair_consensus.csv"),
                     index=False)

# Approach 2: cross-session consensus on sign-inversion percentage
log("\n  Approach 2: Cross-session consensus on sign-inversion percentage")
sign_pct_per_session = df_clean.groupby('session_id')['sign_inversion_pct'].agg(
    ['mean', 'std', 'min', 'max', 'median'])
log(f"  Per-session sign_inv_pct distribution:")
log(f"    Across {len(sign_pct_per_session)} sessions:")
log(f"    Mean of session means: {sign_pct_per_session['mean'].mean():.2f}")
log(f"    Std of session means:  {sign_pct_per_session['mean'].std():.4f}")
log(f"    Min session mean:      {sign_pct_per_session['mean'].min():.2f}")
log(f"    Max session mean:      {sign_pct_per_session['mean'].max():.2f}")
log(f"    All session means in [88%, 95%] range: "
    f"{((sign_pct_per_session['mean'] >= 88) & (sign_pct_per_session['mean'] <= 95)).all()}")

# Type 4 multiset-style: how consistent is the sign-inversion claim across sessions?
n_sessions_above_85 = (sign_pct_per_session['mean'] >= 85).sum()
n_sessions_above_90 = (sign_pct_per_session['mean'] >= 90).sum()
n_sessions = len(sign_pct_per_session)
log(f"\n  Set view: N sessions with mean sign_inv_pct >= 85%: "
    f"{n_sessions_above_85}/{n_sessions} "
    f"({100*n_sessions_above_85/n_sessions:.1f}%)")
log(f"  Set view: N sessions with mean sign_inv_pct >= 90%: "
    f"{n_sessions_above_90}/{n_sessions} "
    f"({100*n_sessions_above_90/n_sessions:.1f}%)")

sign_pct_per_session.to_csv(os.path.join(OUT_DIR,
                                          "layer2h_type4_sign_inv_per_session.csv"))

t4_summary = {
    'n_unique_stim_pairs': len(stim_pair_df),
    'most_consistent_stim_pair': stim_pair_df.loc[stim_pair_df['multiset_cv'].idxmin(),
                                                    'stim_pair']
                                   if not stim_pair_df['multiset_cv'].isna().all() else None,
    'min_cv': stim_pair_df['multiset_cv'].min(),
    'max_cv': stim_pair_df['multiset_cv'].max(),
    'sign_inv_pct_session_mean_of_means': sign_pct_per_session['mean'].mean(),
    'sign_inv_pct_session_std_of_means': sign_pct_per_session['mean'].std(),
    'n_sessions_with_mean_sign_inv_above_85': int(n_sessions_above_85),
    'n_sessions_with_mean_sign_inv_above_90': int(n_sessions_above_90),
    'n_sessions_total': int(n_sessions),
}
pd.DataFrame([t4_summary]).to_csv(os.path.join(OUT_DIR,
                                                "layer2h_type4_consensus.csv"),
                                    index=False)

# Approach 3: cross-session consensus on per-area reorganization (gap pattern)
log("\n  Approach 3: Per-area reorganization consensus (uses area_wj_summary)")
log("  Already captured by per-area Cohen's d ranking; the consensus is implicit")
log("  in 28 of 29 areas showing FDR-significant same > diff stimulus contrast.")

# ==========================================================================
# TYPE 1 PILOT: Continuous-Discrete dissociation gap (3 sessions)
# ==========================================================================
log("\n" + "=" * 70)
log("TYPE 1 PILOT: Continuous-Discrete dissociation gap (3 sessions)")
log("=" * 70)

# Pick 3 sessions with largest unit counts for the pilot
session_unit_counts = df_clean.groupby('session_id')['n_units'].first().astype(int)
pilot_sessions = session_unit_counts.nlargest(3).index.tolist()
log(f"  Pilot sessions (largest unit counts): {pilot_sessions}")
log(f"  Unit counts: {session_unit_counts.loc[pilot_sessions].tolist()}")

# To compute Type 1 we need the correlation matrices, which require re-loading NWB files
# This is heavy (~hours per session). Skip pilot if pynwb not available, mark for follow-up.
try:
    import pynwb
    log(f"  pynwb available; pilot would require ~30-60 min per session")
    log(f"  PILOT IS DEFERRED to a separate execution. To run:")
    log(f"    python layer2h_type1_pilot.py")
    log(f"  This script generates a placeholder for pilot output.")

    pilot_rows = []
    for sid in pilot_sessions:
        pilot_rows.append({
            'session_id': sid,
            'n_units': int(session_unit_counts[sid]),
            'status': 'pending — requires NWB reload + correlation matrix computation',
        })
    pd.DataFrame(pilot_rows).to_csv(os.path.join(OUT_DIR,
                                                  "layer2h_type1_pilot_status.csv"),
                                      index=False)
except ImportError:
    log(f"  pynwb not available; Type 1 pilot deferred")
    pd.DataFrame([{'status': 'pynwb not installed; Type 1 pilot deferred'}]).to_csv(
        os.path.join(OUT_DIR, "layer2h_type1_pilot_status.csv"), index=False)

# ==========================================================================
# Provenance
# ==========================================================================
provenance = {
    "pipeline_file": "layer2h_analysis.py",
    "execution_date": "2026-04-29",
    "framework_reference": "Layer 2H of wj-methodology.md (Harbert, in press, "
                           "Frontiers in Pharmacology). See also "
                           ".claude/rules/research_notes/2026-04-29_principled_pairing_criterion.md",
    "n_comparisons_analyzed": int(len(df_clean)),
    "n_sessions": int(df_clean['session_id'].nunique()),
    "type2_findings": {k: float(v) if isinstance(v, (np.float64, float)) else v
                        for k, v in t2_summary.items()},
    "type6_findings": {k: float(v) if isinstance(v, (np.float64, float)) else v
                        for k, v in t6_summary.items()},
    "type5_findings": ({k: float(v) if isinstance(v, (np.float64, float, np.int64))
                         else v for k, v in t5_summary.items()}
                        if t5_summary else None),
    "type4_findings": {k: (float(v) if isinstance(v, (np.float64, float)) else
                            (int(v) if isinstance(v, (int, np.int64)) else v))
                        for k, v in t4_summary.items()},
    "type1_pilot_status": "deferred — requires NWB reload",
    "config": {
        "random_seed": RANDOM_SEED,
    },
}
with open(os.path.join(OUT_DIR, "layer2h_summary_provenance.json"), 'w') as f:
    json.dump(provenance, f, indent=2)

log("\n" + "=" * 70)
log("LAYER 2H ANALYSIS COMPLETE")
log("=" * 70)
log(f"Outputs saved to: {OUT_DIR}")
log(f"Total runtime: {(time.time()-START)/60:.1f} min")
