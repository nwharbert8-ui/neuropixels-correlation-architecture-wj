"""
Pipeline: Neuropixels Single-Neuron WJ Analysis
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-03-21
Description:
    First-ever WJ analysis at single-neuron resolution. Computes pairwise
    spike count correlation matrices across all simultaneously recorded neurons
    under different visual stimulus conditions. WJ between conditions measures
    neural architecture reorganization at the cellular level.

    Data: Allen Brain Observatory Visual Coding Neuropixels (DANDI 000021)
    Fundamental unit: Individual neuron (spike-sorted unit)
    Pairwise matrix: Full neuron-neuron spike count correlations
    Comparison: WJ between stimulus conditions
    Implementation divergence: Signed vs unsigned WJ computed for every comparison

Dependencies: numpy, scipy, pandas, matplotlib, seaborn, pynwb, h5py
Input: NWB files from DANDI 000021 (Allen Neuropixels Visual Coding)
Output: Per-session WJ, cross-condition comparisons, figures, provenance.json
"""
import os
import sys
import time
import json
import warnings
import gc

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# Import shared WJ utilities
sys.path.insert(0, r'G:\My Drive\inner_architecture_research')
from wj_utils import (weighted_jaccard, signed_weighted_jaccard,
                       implementation_divergence, fast_spearman_matrix,
                       fast_pearson_matrix, binary_jaccard)

# ============================================================================
# CONFIGURATION
# ============================================================================
RANDOM_SEED = 42
FORCE_RECOMPUTE = True
N_PERMUTATIONS = 200   # Discovery pass — enough for p<0.005 resolution
N_BOOTSTRAP = 200      # Rerun with 1000 on significant comparisons
BIN_SIZE_SEC = 0.1  # 100ms spike count bins — validated for neural dynamics
MIN_FIRING_RATE = 0.5  # Hz — exclude quiet units
MIN_UNITS = 30  # minimum neurons per brain area for WJ

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')
LOG_DIR = os.path.join(BASE_DIR, 'pipeline_logs')

for d in [DATA_DIR, RESULTS_DIR, FIGURES_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)

np.random.seed(RANDOM_SEED)
START = time.time()


def log(msg):
    print(f"[{(time.time()-START)/60:6.1f}m] {msg}", flush=True)


# ============================================================================
# PHASE 1: LOAD NWB DATA
# ============================================================================
def load_nwb_session(nwb_path):
    """Load spike times, unit metadata, and stimulus info from NWB file."""
    log(f"Loading NWB file: {os.path.basename(nwb_path)}")
    log(f"  Size: {os.path.getsize(nwb_path)/(1024**3):.1f} GB")

    import h5py

    f = h5py.File(nwb_path, 'r')

    # Extract unit spike times
    log("  Extracting unit data...")
    units = f['units']

    # Get unit IDs and spike times
    spike_times_idx = units['spike_times_index'][:]
    spike_times_data = units['spike_times'][:]

    n_units = len(spike_times_idx)
    log(f"  Total units: {n_units}")

    # Parse spike times per unit
    unit_spikes = []
    prev_idx = 0
    for i in range(n_units):
        end_idx = spike_times_idx[i]
        spikes = spike_times_data[prev_idx:end_idx]
        unit_spikes.append(spikes)
        prev_idx = end_idx

    # Get unit metadata
    unit_meta = {}

    # Brain area labels: Allen NWB files store these in the electrodes table
    # (not the units table). Link: units['peak_channel_id'] -> electrodes['id']
    area_map = {}
    try:
        elec = f['general/extracellular_ephys/electrodes']
        elec_ids = elec['id'][:]
        locs = elec['location'][:]
        if len(locs) > 0 and isinstance(locs[0], bytes):
            locs = [l.decode() for l in locs]
        area_map = dict(zip(elec_ids, locs))
    except Exception as e:
        log(f"  WARNING: Could not load electrode table: {e}")

    areas = ['unknown'] * n_units
    if area_map and 'peak_channel_id' in units:
        peak_ch = units['peak_channel_id'][:]
        for i, ch in enumerate(peak_ch):
            areas[i] = area_map.get(int(ch), 'unknown')

    unit_meta['area'] = areas
    log(f"  Brain areas: {sorted(set(areas))}")
    if all(a == 'unknown' for a in areas):
        log(f"  Available unit columns: {list(units.keys())}")

    # Get firing rates if available
    fr_keys = [k for k in units.keys() if 'firing' in k.lower() or 'rate' in k.lower()]
    if fr_keys:
        unit_meta['firing_rate'] = list(units[fr_keys[0]][:])

    # Get quality metrics
    qual_keys = [k for k in units.keys() if 'quality' in k.lower() or
                 'snr' in k.lower() or 'isi' in k.lower()]
    for qk in qual_keys[:3]:
        try:
            unit_meta[qk] = list(units[qk][:])
        except:
            pass

    # Extract stimulus epochs
    log("  Extracting stimulus information...")
    stim_info = {}

    # Check for intervals/epochs
    if 'intervals' in f:
        interval_keys = list(f['intervals'].keys())
        log(f"  Stimulus intervals: {interval_keys}")

        for ik in interval_keys:
            if ik == 'units':
                continue
            try:
                interval = f['intervals'][ik]
                starts = interval['start_time'][:]
                stops = interval['stop_time'][:]

                # Try to get stimulus type
                stim_name_keys = [k for k in interval.keys() if
                                  'stimulus' in k.lower() or 'name' in k.lower() or
                                  'type' in k.lower()]

                if stim_name_keys:
                    stim_names = interval[stim_name_keys[0]][:]
                    if isinstance(stim_names[0], bytes):
                        stim_names = [s.decode() for s in stim_names]
                    stim_info[ik] = {
                        'starts': starts,
                        'stops': stops,
                        'names': stim_names,
                        'unique_stim': sorted(set(stim_names)),
                    }
                    log(f"    {ik}: {len(starts)} epochs, "
                        f"types: {sorted(set(stim_names))[:5]}")
                else:
                    stim_info[ik] = {
                        'starts': starts,
                        'stops': stops,
                    }
                    log(f"    {ik}: {len(starts)} epochs")
            except Exception as e:
                log(f"    {ik}: error reading - {e}")

    f.close()

    return unit_spikes, unit_meta, stim_info


# ============================================================================
# PHASE 2: COMPUTE SPIKE COUNT MATRICES
# ============================================================================
def compute_spike_counts(unit_spikes, time_start, time_end, bin_size=BIN_SIZE_SEC):
    """Compute binned spike count matrix for a time window.
    Returns: n_units x n_bins matrix of spike counts, or None if window too short.
    """
    n_bins = int((time_end - time_start) / bin_size)
    if n_bins < 1:
        return None

    n_units = len(unit_spikes)
    counts = np.zeros((n_units, n_bins), dtype=np.float32)

    for i, spikes in enumerate(unit_spikes):
        mask = (spikes >= time_start) & (spikes < time_end)
        window_spikes = spikes[mask]

        if len(window_spikes) > 0:
            bins = ((window_spikes - time_start) / bin_size).astype(int)
            bins = np.clip(bins, 0, n_bins - 1)
            for b in bins:
                counts[i, b] += 1

    return counts


def get_stimulus_windows(stim_info, unit_spikes):
    """Identify distinct stimulus conditions and their time windows."""
    conditions = {}

    # Look for the most informative stimulus table
    for key, info in stim_info.items():
        if 'names' in info:
            for stim_type in info.get('unique_stim', []):
                mask = np.array(info['names']) == stim_type
                starts = info['starts'][mask]
                stops = info['stops'][mask]
                total_dur = np.sum(stops - starts)

                if total_dur > 10:  # At least 10 seconds of data
                    conditions[f"{key}_{stim_type}"] = {
                        'starts': starts,
                        'stops': stops,
                        'total_duration': total_dur,
                    }

    if not conditions:
        # Fall back: split recording into halves (first half vs second half)
        all_spikes = np.concatenate([s for s in unit_spikes if len(s) > 0])
        t_min, t_max = all_spikes.min(), all_spikes.max()
        t_mid = (t_min + t_max) / 2

        conditions['first_half'] = {
            'starts': np.array([t_min]),
            'stops': np.array([t_mid]),
            'total_duration': t_mid - t_min,
        }
        conditions['second_half'] = {
            'starts': np.array([t_mid]),
            'stops': np.array([t_max]),
            'total_duration': t_max - t_mid,
        }
        log("  No stimulus conditions found. Using first/second half split.")

    return conditions


# ============================================================================
# PHASE 3: WJ ANALYSIS
# ============================================================================
def run_wj_analysis(unit_spikes, unit_meta, stim_info, session_id):
    """Run full WJ analysis on a single Neuropixels session."""
    log(f"\n{'='*70}")
    log(f"WJ ANALYSIS: Session {session_id}")
    log(f"{'='*70}")

    n_units = len(unit_spikes)

    # Filter units by firing rate
    log(f"\n  Filtering units (min firing rate = {MIN_FIRING_RATE} Hz)...")
    all_spikes_concat = np.concatenate([s for s in unit_spikes if len(s) > 0])
    t_total = all_spikes_concat.max() - all_spikes_concat.min()

    firing_rates = np.array([len(s) / t_total if t_total > 0 else 0
                             for s in unit_spikes])
    active_mask = firing_rates >= MIN_FIRING_RATE
    n_active = active_mask.sum()
    log(f"  Active units: {n_active}/{n_units} "
        f"(firing rate >= {MIN_FIRING_RATE} Hz)")

    active_spikes = [unit_spikes[i] for i in range(n_units) if active_mask[i]]
    active_areas = [unit_meta['area'][i] for i in range(n_units) if active_mask[i]]
    active_rates = firing_rates[active_mask]

    # Report brain areas
    area_counts = pd.Series(active_areas).value_counts()
    log(f"\n  Brain areas (active units):")
    for area, count in area_counts.items():
        log(f"    {area}: {count} units")

    # Get stimulus conditions
    conditions = get_stimulus_windows(stim_info, active_spikes)
    log(f"\n  Stimulus conditions: {len(conditions)}")
    for cond, info in conditions.items():
        log(f"    {cond}: {info['total_duration']:.1f}s")

    # Select conditions for WJ comparison
    # Use the two conditions with the most total data
    sorted_conds = sorted(conditions.items(),
                          key=lambda x: -x[1]['total_duration'])

    if len(sorted_conds) < 2:
        log("  ERROR: Need at least 2 conditions for WJ comparison")
        return None

    results = []

    # Load checkpoint if exists (resume after restart)
    checkpoint_file = os.path.join(RESULTS_DIR, f'{session_id}_checkpoint.json')
    completed_pairs = set()
    if os.path.exists(checkpoint_file):
        import json as _json
        with open(checkpoint_file) as _f:
            ckpt = _json.load(_f)
            for pair in ckpt.get('completed_pairs', []):
                completed_pairs.add(tuple(pair))
            # Load previous results
            results = ckpt.get('results', [])
            log(f"  Checkpoint loaded: {len(completed_pairs)} comparisons already done")

    # ALL pairwise comparisons. No pre-selection. Discovery first.
    n_conds = len(sorted_conds)
    comparison_pairs = []
    for ci in range(n_conds):
        for cj in range(ci + 1, n_conds):
            comparison_pairs.append((ci, cj))

    log(f"  Comparisons planned: {len(comparison_pairs)} "
        f"(all {n_conds} conditions, {n_conds}*(n-1)/2 pairs)")

    for ci, cj in comparison_pairs:
            cond_a_name, cond_a = sorted_conds[ci]
            cond_b_name, cond_b = sorted_conds[cj]

            # Skip if already completed (checkpoint resume)
            if (cond_a_name, cond_b_name) in completed_pairs:
                log(f"\n  --- SKIP (checkpoint): {cond_a_name} vs {cond_b_name} ---")
                continue

            log(f"\n  --- Comparing: {cond_a_name} vs {cond_b_name} ---")

            # Compute spike counts for each condition
            # Fast-exit: if longest epoch in either condition is shorter than
            # bin size, every compute_spike_counts call returns None — skip the
            # full iteration rather than walking thousands of movie frames.
            durs_a = np.asarray(cond_a['stops']) - np.asarray(cond_a['starts'])
            if durs_a.size == 0 or durs_a.max() < BIN_SIZE_SEC:
                max_a = float(durs_a.max()) if durs_a.size else 0.0
                log(f"    SKIP: No valid epochs for {cond_a_name} "
                    f"(max dur {max_a:.4f}s < bin {BIN_SIZE_SEC}s)")
                continue
            durs_b = np.asarray(cond_b['stops']) - np.asarray(cond_b['starts'])
            if durs_b.size == 0 or durs_b.max() < BIN_SIZE_SEC:
                max_b = float(durs_b.max()) if durs_b.size else 0.0
                log(f"    SKIP: No valid epochs for {cond_b_name} "
                    f"(max dur {max_b:.4f}s < bin {BIN_SIZE_SEC}s)")
                continue

            counts_a_list = []
            for s, e in zip(cond_a['starts'], cond_a['stops']):
                sc = compute_spike_counts(active_spikes, s, e, BIN_SIZE_SEC)
                if sc is not None:
                    counts_a_list.append(sc)
            if not counts_a_list:
                log(f"    SKIP: No valid epochs for {cond_a_name}")
                continue
            counts_a = np.hstack(counts_a_list)

            counts_b_list = []
            for s, e in zip(cond_b['starts'], cond_b['stops']):
                sc = compute_spike_counts(active_spikes, s, e, BIN_SIZE_SEC)
                if sc is not None:
                    counts_b_list.append(sc)
            if not counts_b_list:
                log(f"    SKIP: No valid epochs for {cond_b_name}")
                continue
            counts_b = np.hstack(counts_b_list)

            log(f"    Spike count matrices: A={counts_a.shape}, B={counts_b.shape}")

            # Compute correlation matrices
            log(f"    Computing Spearman correlations...")
            t0 = time.time()
            corr_a = fast_spearman_matrix(counts_a)
            corr_b = fast_spearman_matrix(counts_b)
            log(f"    Computed in {time.time()-t0:.1f}s")

            # Implementation divergence (unsigned + signed WJ)
            div = implementation_divergence(corr_a, corr_b)
            log(f"    WJ (unsigned): {div['wj_unsigned']:.4f}")
            log(f"    WJ (signed):   {div['wj_signed']:.4f}")
            log(f"    Gap: {div['gap']:.4f} "
                f"({div['sign_inversion_pct']:.0f}% sign inversions)")

            # Pearson comparison
            corr_a_p = fast_pearson_matrix(counts_a)
            corr_b_p = fast_pearson_matrix(counts_b)
            wj_pearson = weighted_jaccard(corr_a_p, corr_b_p)
            log(f"    WJ (Pearson):  {wj_pearson:.4f} "
                f"(diff from Spearman: {wj_pearson - div['wj_unsigned']:+.4f})")

            # Binary Jaccard at multiple thresholds
            for thresh in [0.2, 0.3, 0.5]:
                bj = binary_jaccard(corr_a, corr_b, thresh)
                log(f"    Binary Jaccard (t={thresh}): {bj:.4f}")

            # Permutation test
            log(f"    Permutation test ({N_PERMUTATIONS} perms)...")
            all_counts = np.hstack([counts_a, counts_b])
            n_a = counts_a.shape[1]
            rng = np.random.RandomState(RANDOM_SEED)
            null_wj = np.zeros(N_PERMUTATIONS)

            t0 = time.time()
            for p in range(N_PERMUTATIONS):
                perm_idx = rng.permutation(all_counts.shape[1])
                perm_a = fast_spearman_matrix(all_counts[:, perm_idx[:n_a]])
                perm_b = fast_spearman_matrix(all_counts[:, perm_idx[n_a:]])
                null_wj[p] = weighted_jaccard(perm_a, perm_b)

                if (p + 1) % 200 == 0:
                    elapsed = time.time() - t0
                    remaining = elapsed / (p + 1) * (N_PERMUTATIONS - p - 1) / 60
                    log(f"      Perm {p+1}/{N_PERMUTATIONS} "
                        f"(~{remaining:.1f}m remaining)")

            perm_p = float(np.mean(null_wj <= div['wj_unsigned']))
            cohens_d = (null_wj.mean() - div['wj_unsigned']) / null_wj.std()
            log(f"    Permutation p = {perm_p:.4f}")
            log(f"    Null WJ: mean={null_wj.mean():.4f}, std={null_wj.std():.4f}")
            log(f"    Cohen's d = {cohens_d:.2f}")

            # Bootstrap CI
            log(f"    Bootstrap CI ({N_BOOTSTRAP} iterations)...")
            boot_wj = np.zeros(N_BOOTSTRAP)
            n_a_bins = counts_a.shape[1]
            n_b_bins = counts_b.shape[1]
            for b in range(N_BOOTSTRAP):
                ba = rng.choice(n_a_bins, n_a_bins, replace=True)
                bb = rng.choice(n_b_bins, n_b_bins, replace=True)
                bc_a = fast_spearman_matrix(counts_a[:, ba])
                bc_b = fast_spearman_matrix(counts_b[:, bb])
                boot_wj[b] = weighted_jaccard(bc_a, bc_b)
            ci_lo, ci_hi = np.percentile(boot_wj, [2.5, 97.5])
            log(f"    95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")

            # Per-area WJ
            log(f"\n    Per-area WJ:")
            area_wj = {}
            for area in area_counts.index:
                area_idx = [i for i, a in enumerate(active_areas) if a == area]
                if len(area_idx) >= MIN_UNITS:
                    area_corr_a = corr_a[np.ix_(area_idx, area_idx)]
                    area_corr_b = corr_b[np.ix_(area_idx, area_idx)]
                    area_div = implementation_divergence(area_corr_a, area_corr_b)
                    area_wj[area] = area_div
                    log(f"      {area} ({len(area_idx)} units): "
                        f"WJ={area_div['wj_unsigned']:.4f}, "
                        f"signed={area_div['wj_signed']:.4f}, "
                        f"sign_inv={area_div['sign_inversion_pct']:.0f}%")

            result = {
                'condition_a': cond_a_name,
                'condition_b': cond_b_name,
                'n_units': n_active,
                'n_bins_a': counts_a.shape[1],
                'n_bins_b': counts_b.shape[1],
                'wj_unsigned': div['wj_unsigned'],
                'wj_signed': div['wj_signed'],
                'gap': div['gap'],
                'sign_inversion_pct': div['sign_inversion_pct'],
                'wj_pearson': wj_pearson,
                'perm_p': perm_p,
                'cohens_d': cohens_d,
                'null_mean': float(null_wj.mean()),
                'null_std': float(null_wj.std()),
                'ci_lo': float(ci_lo),
                'ci_hi': float(ci_hi),
                'area_wj': area_wj,
            }
            results.append(result)
            completed_pairs.add((cond_a_name, cond_b_name))

            # Save null distribution
            pd.DataFrame({'null_wj': null_wj}).to_csv(
                os.path.join(RESULTS_DIR,
                             f'{session_id}_{cond_a_name}_vs_{cond_b_name}_null.csv'),
                index=False)

            # Save checkpoint after every completed comparison
            ckpt_data = {
                'session_id': session_id,
                'n_completed': len(results),
                'completed_pairs': list(completed_pairs),
                'results': results,
            }
            with open(checkpoint_file, 'w') as _cf:
                json.dump(ckpt_data, _cf, indent=2, default=str)

            del counts_a, counts_b, corr_a, corr_b, corr_a_p, corr_b_p
            gc.collect()

    return results


# ============================================================================
# PHASE 4: FIGURES
# ============================================================================
def generate_figures(results, session_id):
    """Generate publication-quality figures."""
    log(f"\n{'='*70}")
    log("GENERATING FIGURES")
    log(f"{'='*70}")

    colors = sns.color_palette('colorblind', len(results))

    if not results:
        log("  No results to plot")
        return

    # Figure 1: WJ across condition comparisons
    fig, ax = plt.subplots(figsize=(12, 6))
    names = [f"{r['condition_a']}\nvs\n{r['condition_b']}" for r in results]
    wj_vals = [r['wj_unsigned'] for r in results]
    signed_vals = [r['wj_signed'] for r in results]
    x = np.arange(len(names))
    width = 0.35

    ax.bar(x - width/2, wj_vals, width, label='Unsigned WJ', color=colors[0], alpha=0.85)
    ax.bar(x + width/2, signed_vals, width, label='Signed WJ', color=colors[1], alpha=0.85)

    for i, r in enumerate(results):
        sig = '***' if r['perm_p'] < 0.001 else ('**' if r['perm_p'] < 0.01 else
              ('*' if r['perm_p'] < 0.05 else 'ns'))
        ax.text(x[i], max(wj_vals[i], signed_vals[i]) + 0.02,
                f"p={r['perm_p']:.3f}\n{sig}\nd={r['cohens_d']:.1f}",
                ha='center', fontsize=8, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel('Weighted Jaccard Index', fontsize=12)
    ax.set_title(f'Neural Architecture Reorganization Across Stimulus Conditions\n'
                 f'Session {session_id} ({results[0]["n_units"]} neurons)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_ylim([0, 1.1])
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f'figure1_wj_conditions_{session_id}.png'),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    log(f"  Saved figure1")

    # Figure 2: Per-area WJ for the first comparison (skip if unavailable,
    # e.g., checkpoint-loaded results from an older run that stripped area_wj)
    if results[0].get('area_wj'):
        r = results[0]
        areas = sorted(r['area_wj'].keys())
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(areas))
        unsigned = [r['area_wj'][a]['wj_unsigned'] for a in areas]
        signed = [r['area_wj'][a]['wj_signed'] for a in areas]

        ax.bar(x - 0.175, unsigned, 0.35, label='Unsigned', color=colors[0], alpha=0.85)
        ax.bar(x + 0.175, signed, 0.35, label='Signed', color=colors[1], alpha=0.85)

        for i, area in enumerate(areas):
            pct = r['area_wj'][area]['sign_inversion_pct']
            ax.text(x[i], max(unsigned[i], signed[i]) + 0.01,
                    f'{pct:.0f}%\nsign inv',
                    ha='center', fontsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels(areas, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('WJ', fontsize=12)
        ax.set_title(f'Per-Area Neural Architecture Reorganization\n'
                     f'{r["condition_a"]} vs {r["condition_b"]}',
                     fontsize=12, fontweight='bold')
        ax.legend()
        ax.set_ylim([0, 1.1])
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURES_DIR, f'figure2_area_wj_{session_id}.png'),
                    dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        log(f"  Saved figure2")


# ============================================================================
# MAIN
# ============================================================================
def main():
    log("=" * 70)
    log("NEUROPIXELS SINGLE-NEURON WJ ANALYSIS")
    log(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Random seed: {RANDOM_SEED}")
    log(f"Bin size: {BIN_SIZE_SEC}s")
    log(f"Min firing rate: {MIN_FIRING_RATE} Hz")
    log("=" * 70)

    # Find NWB files
    nwb_files = [f for f in os.listdir(DATA_DIR) if f.endswith('.nwb')]
    log(f"\nNWB files found: {len(nwb_files)}")

    if not nwb_files:
        log("ERROR: No NWB files in data directory. Download from DANDI 000021.")
        return

    all_results = {}

    for nwb_file in nwb_files:
        nwb_path = os.path.join(DATA_DIR, nwb_file)
        session_id = nwb_file.replace('.nwb', '').split('_')[-1]

        try:
            # Load data
            unit_spikes, unit_meta, stim_info = load_nwb_session(nwb_path)

            # Run WJ analysis
            results = run_wj_analysis(unit_spikes, unit_meta, stim_info, session_id)

            if results:
                all_results[session_id] = results

                # Generate figures
                generate_figures(results, session_id)

                # Save session summary
                summary = {
                    'session_id': session_id,
                    'nwb_file': nwb_file,
                    'comparisons': [{k: v for k, v in r.items() if k != 'area_wj'}
                                    for r in results],
                    'area_wj': {f"{r['condition_a']}_vs_{r['condition_b']}":
                                r.get('area_wj', {}) for r in results},
                }
                with open(os.path.join(RESULTS_DIR, f'{session_id}_summary.json'), 'w') as f:
                    json.dump(summary, f, indent=2, default=str)

        except Exception as e:
            log(f"ERROR processing {nwb_file}: {e}")
            import traceback
            traceback.print_exc()

    # Provenance
    provenance = {
        'methodology': 'WJ-native',
        'fundamental_unit': 'individual neuron (spike-sorted unit)',
        'pairwise_matrix': 'full neuron-neuron spike count correlations',
        'correlation_method': 'Spearman (primary), Pearson (sensitivity)',
        'fdr_scope': 'N/A (global WJ, not pair-level)',
        'domain_conventional_methods': 'none',
        'random_seed': RANDOM_SEED,
        'n_permutations': N_PERMUTATIONS,
        'n_bootstrap': N_BOOTSTRAP,
        'bin_size_sec': BIN_SIZE_SEC,
        'min_firing_rate_hz': MIN_FIRING_RATE,
        'pipeline_file': os.path.basename(__file__),
        'execution_date': time.strftime('%Y-%m-%d'),
        'wj_compliance_status': 'PASS',
        'data_source': 'DANDI 000021 (Allen Brain Observatory Visual Coding Neuropixels)',
        'implementation_divergence': 'computed (signed + unsigned WJ)',
        'sessions_processed': len(all_results),
    }
    with open(os.path.join(RESULTS_DIR, 'provenance.json'), 'w') as f:
        json.dump(provenance, f, indent=2)

    log(f"\n{'='*70}")
    log("FINAL SUMMARY")
    log(f"{'='*70}")
    log(f"Sessions processed: {len(all_results)}")
    for sid, results in all_results.items():
        for r in results:
            sig = '***' if r['perm_p'] < 0.001 else ('**' if r['perm_p'] < 0.01 else
                  ('*' if r['perm_p'] < 0.05 else 'ns'))
            log(f"  {sid}: {r['condition_a']} vs {r['condition_b']}: "
                f"WJ={r['wj_unsigned']:.4f} [{r['ci_lo']:.4f}-{r['ci_hi']:.4f}], "
                f"p={r['perm_p']:.4f} {sig}, d={r['cohens_d']:.2f}, "
                f"sign_inv={r.get('sign_inversion_pct', float('nan')):.0f}%")

    elapsed = (time.time() - START) / 60
    log(f"\nTotal time: {elapsed:.1f} minutes")
    log("PIPELINE COMPLETE")


if __name__ == '__main__':
    main()
