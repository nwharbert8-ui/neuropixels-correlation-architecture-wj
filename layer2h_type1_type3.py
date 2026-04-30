"""
Pipeline: Layer 2H Types 1 + 3 — Full 30-Session Scale-Up
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-29

Description:
    Full-scale Layer 2H Type 1 (Continuous-Discrete) and Type 3 (Asymmetric-
    Symmetric Tversky) decomposition on the Neuropixels Visual Coding dataset.
    Runs on all NWB sessions present in data/. For each pairwise stimulus
    condition comparison, computes:
        - Continuous WJ unsigned on full pairwise Spearman correlations
        - Binary Jaccard at top 5% threshold (high-coupled pair overlap)
        - Type 1 dissociation gap = continuous WJ - binary Jaccard
        - Tversky T(A,B; alpha=1, beta=0) — fraction of A preserved in B
        - Tversky T(A,B; alpha=0, beta=1) — fraction of B preserved in A
        - Type 3 asymmetric gap = T_A_weighted - T_B_weighted

    Type 1 localizes whether sign-driven reorganization concentrates at top-
    coupled pairs vs spreads broadly. Type 3 captures directional asymmetry
    in architectural overlap between stimulus categories.

Inputs:
    data/sub-*_ses-*.nwb (all sessions present)

Outputs:
    results/layer2h_type1_type3_per_comparison.csv
    results/layer2h_type1_type3_summary.csv
    results/layer2h_type1_type3_provenance.json
"""

import os
import sys
import time
import json
import gc
import warnings
import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata
warnings.filterwarnings('ignore')

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
START = time.time()

BIN_SIZE_SEC = 0.1
MIN_FIRING_RATE = 0.5
TOP_PCT = 5
MIN_BINS = 50

ROOT = r"G:\My Drive\inner_architecture_research\neuropixels_wj"
DATA_DIR = os.path.join(ROOT, "data")
RESULTS_DIR = os.path.join(ROOT, "results")

# Run on all sessions present in data folder (replaces 3-session pilot)
PILOT_SESSIONS = None  # None = auto-detect from NWB files


def log(msg):
    print(f"[{(time.time()-START)/60:6.1f}m] {msg}", flush=True)


# ==========================================================================
# Utilities
# ==========================================================================
def load_nwb_session(nwb_path):
    """Load spike times and stimulus intervals from NWB file."""
    log(f"  Loading {os.path.basename(nwb_path)} "
        f"({os.path.getsize(nwb_path)/(1024**3):.1f} GB)")
    f = h5py.File(nwb_path, 'r')

    units = f['units']
    spike_times_idx = units['spike_times_index'][:]
    spike_times_data = units['spike_times'][:]

    n_units = len(spike_times_idx)
    log(f"    Total units: {n_units}")

    unit_spikes = []
    prev_idx = 0
    for i in range(n_units):
        end_idx = spike_times_idx[i]
        unit_spikes.append(spike_times_data[prev_idx:end_idx])
        prev_idx = end_idx

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
                stim_name_keys = [k for k in interval.keys() if
                                  'stimulus' in k.lower() or 'name' in k.lower() or
                                  'type' in k.lower()]
                if stim_name_keys:
                    stim_names = interval[stim_name_keys[0]][:]
                    if len(stim_names) > 0 and isinstance(stim_names[0], bytes):
                        stim_names = [s.decode() for s in stim_names]
                    stim_info[ik] = {
                        'starts': starts, 'stops': stops, 'names': stim_names,
                        'unique_stim': sorted(set(stim_names)),
                    }
                else:
                    stim_info[ik] = {'starts': starts, 'stops': stops}
            except Exception:
                pass

    f.close()
    return unit_spikes, stim_info


def get_stimulus_windows(stim_info):
    """Identify distinct stimulus conditions matching existing pipeline naming."""
    conditions = {}
    for key, info in stim_info.items():
        if 'names' in info:
            for stim_type in info.get('unique_stim', []):
                mask = np.array(info['names']) == stim_type
                starts = info['starts'][mask]
                stops = info['stops'][mask]
                total_dur = np.sum(stops - starts)
                if total_dur > 10:
                    conditions[f"{key}_{stim_type}"] = {
                        'starts': starts, 'stops': stops,
                        'total_duration': total_dur,
                    }
    return conditions


def compute_spike_counts(unit_spikes, starts, stops, bin_size=BIN_SIZE_SEC):
    """Concatenated binned spike counts across multiple windows."""
    all_bins = []
    for s, e in zip(starts, stops):
        n_bins = int((e - s) / bin_size)
        if n_bins < 1:
            continue
        bins_in_window = np.zeros((len(unit_spikes), n_bins), dtype=np.float32)
        for i, spikes in enumerate(unit_spikes):
            mask = (spikes >= s) & (spikes < e)
            if mask.any():
                rel = (spikes[mask] - s) / bin_size
                rel_int = np.clip(rel.astype(int), 0, n_bins - 1)
                np.add.at(bins_in_window[i], rel_int, 1)
        all_bins.append(bins_in_window)
    if not all_bins:
        return None
    return np.concatenate(all_bins, axis=1)


def fast_spearman(X):
    """Spearman correlation matrix on N units x B bins. Returns NxN matrix."""
    n_units = X.shape[0]
    # Rank each row
    Xr = np.zeros_like(X)
    for i in range(n_units):
        Xr[i] = rankdata(X[i])
    # Center and normalize
    Xr = Xr - Xr.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(Xr, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Xn = Xr / norms
    C = Xn @ Xn.T
    return C


def weighted_jaccard(a, b):
    """Standard unsigned weighted Jaccard on |r| values."""
    a_abs = np.abs(a)
    b_abs = np.abs(b)
    s_min = np.minimum(a_abs, b_abs).sum()
    s_max = np.maximum(a_abs, b_abs).sum()
    return s_min / s_max if s_max > 0 else np.nan


def binary_jaccard_top_pct(a, b, pct=5):
    """Binary Jaccard at top-pct threshold using |r|."""
    a_abs = np.abs(a)
    b_abs = np.abs(b)
    thr_a = np.percentile(a_abs, 100 - pct)
    thr_b = np.percentile(b_abs, 100 - pct)
    in_a = a_abs >= thr_a
    in_b = b_abs >= thr_b
    inter = (in_a & in_b).sum()
    union = (in_a | in_b).sum()
    return inter / union if union > 0 else np.nan


def tversky_asymmetric(a, b, alpha, beta, pct=5):
    """Asymmetric Tversky on top-pct thresholded sets.
    T(A,B; alpha, beta) = |A intersect B| / (|A intersect B| + alpha*|A-B| + beta*|B-A|)
    alpha=1, beta=0: weighted toward A — fraction of A preserved in B
    alpha=0, beta=1: weighted toward B — fraction of B inherited from A
    """
    a_abs = np.abs(a)
    b_abs = np.abs(b)
    thr_a = np.percentile(a_abs, 100 - pct)
    thr_b = np.percentile(b_abs, 100 - pct)
    in_a = a_abs >= thr_a
    in_b = b_abs >= thr_b
    inter = (in_a & in_b).sum()
    a_only = (in_a & ~in_b).sum()
    b_only = (~in_a & in_b).sum()
    denom = inter + alpha * a_only + beta * b_only
    return inter / denom if denom > 0 else np.nan


# ==========================================================================
# Main
# ==========================================================================
def main():
    log("=" * 70)
    log(f"LAYER 2H TYPES 1 + 3: FULL 30-SESSION SCALE-UP")
    log("=" * 70)

    nwb_files = []
    if PILOT_SESSIONS is None:
        # Auto-detect from NWB files in data dir
        for f in sorted(os.listdir(DATA_DIR)):
            if f.endswith('.nwb'):
                # Extract ses-NNNNN from filename
                parts = f.replace('.nwb', '').split('_')
                sid = next((p for p in parts if p.startswith('ses-')), None)
                if sid:
                    nwb_files.append((sid, os.path.join(DATA_DIR, f)))
    else:
        for sid in PILOT_SESSIONS:
            matches = [f for f in os.listdir(DATA_DIR)
                       if sid in f and f.endswith('.nwb')]
            if not matches:
                log(f"  WARNING: no NWB file for {sid}")
                continue
            nwb_files.append((sid, os.path.join(DATA_DIR, matches[0])))

    log(f"\nSessions to process: {len(nwb_files)}")
    for sid, _ in nwb_files:
        log(f"  {sid}")

    all_results = []

    for sid, nwb_path in nwb_files:
        log(f"\n{'='*70}")
        log(f"SESSION: {sid}")
        log(f"{'='*70}")

        try:
            unit_spikes, stim_info = load_nwb_session(nwb_path)

            # Filter units by firing rate
            all_t = np.concatenate([s for s in unit_spikes if len(s) > 0])
            t_total = all_t.max() - all_t.min()
            firing_rates = np.array([len(s) / t_total if t_total > 0 else 0
                                      for s in unit_spikes])
            active_mask = firing_rates >= MIN_FIRING_RATE
            active_spikes = [unit_spikes[i] for i in range(len(unit_spikes))
                             if active_mask[i]]
            log(f"    Active units: {active_mask.sum()}/{len(unit_spikes)}")

            conditions = get_stimulus_windows(stim_info)
            filtered = {k: v for k, v in conditions.items()
                        if 'natural_movie_one' not in k.lower()}
            log(f"    Conditions: {len(filtered)}")

            # FIX: Compute spike count matrices first, then identify
            # session-common active units (non-zero variance in ALL conditions)
            # before computing correlation matrices. This ensures all
            # correlation vectors use the same unit set and are directly
            # comparable.
            spike_count_matrices = {}
            valid_conditions = []
            for cond_name, info in filtered.items():
                counts = compute_spike_counts(active_spikes, info['starts'],
                                                info['stops'])
                if counts is None or counts.shape[1] < MIN_BINS:
                    log(f"    SKIP {cond_name}: insufficient bins "
                        f"({counts.shape[1] if counts is not None else 0})")
                    continue
                spike_count_matrices[cond_name] = counts
                valid_conditions.append(cond_name)

            if len(valid_conditions) < 2:
                log(f"    SKIP session: fewer than 2 valid conditions")
                continue

            # Find session-common active units: non-zero variance in EVERY condition
            n_total_units = spike_count_matrices[valid_conditions[0]].shape[0]
            common_active = np.ones(n_total_units, dtype=bool)
            for cond_name in valid_conditions:
                stds = spike_count_matrices[cond_name].std(axis=1)
                common_active &= (stds > 0)

            n_common = int(common_active.sum())
            log(f"    Session-common active units (non-zero variance in ALL "
                f"{len(valid_conditions)} conditions): {n_common}")

            if n_common < 100:
                log(f"    SKIP session: too few common units ({n_common})")
                continue

            # Compute correlation matrices on the common unit set
            corr_matrices = {}
            for cond_name in valid_conditions:
                t0 = time.time()
                counts = spike_count_matrices[cond_name][common_active]
                C = fast_spearman(counts)
                iu = np.triu_indices(C.shape[0], k=1)
                r_vec = C[iu]
                corr_matrices[cond_name] = (r_vec, counts.shape[1])
                log(f"    {cond_name}: {counts.shape[1]} bins, "
                    f"{n_common} units, |r|>0.3 fraction "
                    f"{(np.abs(r_vec) > 0.3).mean():.3f}, "
                    f"computed in {time.time()-t0:.1f}s")
                del C
                gc.collect()

            del spike_count_matrices
            gc.collect()

            log(f"\n  Computing pairwise Type 1 gaps...")
            cond_names = list(corr_matrices.keys())
            n_pairs = 0
            for i, ca in enumerate(cond_names):
                ra, nbins_a = corr_matrices[ca]
                for cb in cond_names[i+1:]:
                    rb, nbins_b = corr_matrices[cb]

                    wj_unsigned = weighted_jaccard(ra, rb)
                    bj_top5 = binary_jaccard_top_pct(ra, rb, pct=TOP_PCT)
                    type1_gap = wj_unsigned - bj_top5

                    # Type 3 asymmetric Tversky
                    t_a_weighted = tversky_asymmetric(ra, rb, alpha=1.0,
                                                       beta=0.0, pct=TOP_PCT)
                    t_b_weighted = tversky_asymmetric(ra, rb, alpha=0.0,
                                                       beta=1.0, pct=TOP_PCT)
                    type3_gap = t_a_weighted - t_b_weighted

                    all_results.append({
                        'session_id': sid,
                        'condition_a': ca,
                        'condition_b': cb,
                        'n_units': n_common,
                        'n_bins_a': int(nbins_a),
                        'n_bins_b': int(nbins_b),
                        'wj_unsigned': float(wj_unsigned),
                        'binary_jaccard_top5': float(bj_top5),
                        'type1_gap': float(type1_gap),
                        'tversky_a_weighted': float(t_a_weighted),
                        'tversky_b_weighted': float(t_b_weighted),
                        'type3_asymmetry_gap': float(type3_gap),
                    })
                    n_pairs += 1

            log(f"  Computed Type 1 for {n_pairs} comparisons in this session")
        except Exception as e:
            log(f"  ERROR processing {sid}: {e}")
            import traceback
            log(traceback.format_exc())

    # Save
    df = pd.DataFrame(all_results)
    df.to_csv(os.path.join(RESULTS_DIR,
                            "layer2h_type1_type3_per_comparison.csv"),
                index=False)

    # Summary
    log(f"\n{'='*70}")
    log(f"PILOT SUMMARY (Layer 2H Types 1 + 3 on {df['session_id'].nunique() if len(df) > 0 else 0} sessions)")
    log(f"{'='*70}")
    log(f"  Total comparisons: {len(df)}")
    if len(df) > 0:
        log(f"\n  TYPE 1 (Continuous-Discrete):")
        log(f"    WJ unsigned:         {df['wj_unsigned'].mean():.4f} ± "
            f"{df['wj_unsigned'].std():.4f}")
        log(f"    Binary Jaccard top5: {df['binary_jaccard_top5'].mean():.4f} ± "
            f"{df['binary_jaccard_top5'].std():.4f}")
        log(f"    Type 1 gap:          {df['type1_gap'].mean():.4f} ± "
            f"{df['type1_gap'].std():.4f}")
        log(f"    Type 1 gap range:    [{df['type1_gap'].min():.4f}, "
            f"{df['type1_gap'].max():.4f}]")
        log(f"\n  TYPE 3 (Asymmetric-Symmetric Tversky, top 5%):")
        log(f"    T_A_weighted (frac of A in B):     "
            f"{df['tversky_a_weighted'].mean():.4f} ± "
            f"{df['tversky_a_weighted'].std():.4f}")
        log(f"    T_B_weighted (frac of B in A):     "
            f"{df['tversky_b_weighted'].mean():.4f} ± "
            f"{df['tversky_b_weighted'].std():.4f}")
        log(f"    Type 3 asymmetry gap:              "
            f"{df['type3_asymmetry_gap'].mean():.4f} ± "
            f"{df['type3_asymmetry_gap'].std():.4f}")
        log(f"    Type 3 |gap| mean:                 "
            f"{df['type3_asymmetry_gap'].abs().mean():.4f}")
        log(f"    Type 3 gap range:                  "
            f"[{df['type3_asymmetry_gap'].min():.4f}, "
            f"{df['type3_asymmetry_gap'].max():.4f}]")

    summary = {
        'n_sessions_completed': df['session_id'].nunique() if len(df) > 0 else 0,
        'n_comparisons': len(df),
        'wj_unsigned_mean': float(df['wj_unsigned'].mean()) if len(df) > 0 else None,
        'binary_jaccard_top5_mean': (float(df['binary_jaccard_top5'].mean())
                                       if len(df) > 0 else None),
        'type1_gap_mean': float(df['type1_gap'].mean()) if len(df) > 0 else None,
        'type1_gap_std': float(df['type1_gap'].std()) if len(df) > 0 else None,
        'type1_gap_min': float(df['type1_gap'].min()) if len(df) > 0 else None,
        'type1_gap_max': float(df['type1_gap'].max()) if len(df) > 0 else None,
        'tversky_a_weighted_mean': (float(df['tversky_a_weighted'].mean())
                                      if len(df) > 0 else None),
        'tversky_b_weighted_mean': (float(df['tversky_b_weighted'].mean())
                                      if len(df) > 0 else None),
        'type3_asymmetry_gap_mean': (float(df['type3_asymmetry_gap'].mean())
                                       if len(df) > 0 else None),
        'type3_asymmetry_gap_std': (float(df['type3_asymmetry_gap'].std())
                                      if len(df) > 0 else None),
        'type3_asymmetry_gap_abs_mean': (float(df['type3_asymmetry_gap'].abs().mean())
                                           if len(df) > 0 else None),
        'type3_asymmetry_gap_min': (float(df['type3_asymmetry_gap'].min())
                                      if len(df) > 0 else None),
        'type3_asymmetry_gap_max': (float(df['type3_asymmetry_gap'].max())
                                      if len(df) > 0 else None),
    }
    pd.DataFrame([summary]).to_csv(
        os.path.join(RESULTS_DIR, "layer2h_type1_type3_summary.csv"),
        index=False)

    provenance = {
        'pipeline_file': 'layer2h_type1_pilot.py',
        'execution_date': '2026-04-29',
        'sessions_processed': sorted(df['session_id'].unique().tolist())
                                if len(df) > 0 else [],
        'config': {
            'random_seed': RANDOM_SEED,
            'bin_size_sec': BIN_SIZE_SEC,
            'min_firing_rate': MIN_FIRING_RATE,
            'top_pct': TOP_PCT,
            'min_bins': MIN_BINS,
        },
        'summary': summary,
        'runtime_minutes': (time.time() - START) / 60,
    }
    with open(os.path.join(RESULTS_DIR,
                            "layer2h_type1_type3_provenance.json"), 'w') as f:
        json.dump(provenance, f, indent=2)

    log(f"\nFiles saved:")
    log(f"  layer2h_type1_type3_per_comparison.csv")
    log(f"  layer2h_type1_type3_summary.csv")
    log(f"  layer2h_type1_type3_provenance.json")
    log(f"\nTotal runtime: {(time.time()-START)/60:.1f} min")


if __name__ == "__main__":
    main()
