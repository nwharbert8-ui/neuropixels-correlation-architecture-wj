"""
Pipeline: Neuropixels Per-Area WJ Decomposition (Post-Processing)
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-04-26
Description:
    Retroactively adds per-area WJ decomposition to all existing session
    checkpoints. The main pipeline returned 'unknown' for all brain areas
    because it searched the units table for area labels; Allen NWB files
    store area labels in the electrodes table, linked via peak_channel_id.

    This script:
    1. Loads each session NWB file
    2. Builds electrode_id -> brain_area mapping from the electrodes table
    3. Assigns each unit its brain area via peak_channel_id
    4. For each completed comparison pair, re-bins spike counts and recomputes
       per-area WJ (no permutation testing - global WJ already done)
    5. Updates checkpoint JSONs with area_wj data

Dependencies: numpy, scipy, pandas, h5py
Input: NWB files in data/, checkpoint JSONs in results/
Output: Updated checkpoint JSONs with area_wj entries
"""

import os
import sys
import time
import json
import warnings
import glob

import numpy as np
import pandas as pd
from scipy.stats import rankdata

warnings.filterwarnings('ignore')

sys.path.insert(0, r'G:\My Drive\inner_architecture_research')
from wj_utils import implementation_divergence, fast_spearman_matrix

# ============================================================================
# CONFIG
# ============================================================================
RANDOM_SEED = 42
BIN_SIZE_SEC = 0.1
MIN_FIRING_RATE = 0.5
MIN_UNITS = 30

BASE_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')

np.random.seed(RANDOM_SEED)
START = time.time()


def log(msg):
    print(f"[{(time.time()-START)/60:6.1f}m] {msg}", flush=True)


def load_nwb_with_areas(nwb_path):
    """Load spike times, area labels (via electrode table), and stimulus info."""
    import h5py

    log(f"  Loading: {os.path.basename(nwb_path)}")
    f = h5py.File(nwb_path, 'r')

    units = f['units']
    n_units = len(units['spike_times_index'][:])

    # --- Area labels via electrode table ---
    area_map = {}
    try:
        elec = f['general/extracellular_ephys/electrodes']
        elec_ids = elec['id'][:]
        locs = elec['location'][:]
        if len(locs) > 0 and isinstance(locs[0], bytes):
            locs = [l.decode() for l in locs]
        else:
            locs = list(locs)
        area_map = dict(zip(elec_ids, locs))
    except Exception as e:
        log(f"    WARNING: Could not load electrode table: {e}")

    areas = ['unknown'] * n_units
    if area_map and 'peak_channel_id' in units:
        peak_ch = units['peak_channel_id'][:]
        for i, ch in enumerate(peak_ch):
            areas[i] = area_map.get(int(ch), 'unknown')

    # --- Spike times ---
    spike_times_idx = units['spike_times_index'][:]
    spike_times_data = units['spike_times'][:]
    unit_spikes = []
    prev = 0
    for i in range(n_units):
        end = spike_times_idx[i]
        unit_spikes.append(spike_times_data[prev:end])
        prev = end

    # --- Stimulus intervals ---
    stim_info = {}
    if 'intervals' in f:
        for ik in f['intervals'].keys():
            if ik == 'units':
                continue
            try:
                interval = f['intervals'][ik]
                starts = interval['start_time'][:]
                stops = interval['stop_time'][:]
                name_keys = [k for k in interval.keys() if
                             'stimulus' in k.lower() or 'name' in k.lower()
                             or 'type' in k.lower()]
                if name_keys:
                    snames = interval[name_keys[0]][:]
                    if isinstance(snames[0], bytes):
                        snames = [s.decode() for s in snames]
                    stim_info[ik] = {
                        'starts': starts, 'stops': stops,
                        'names': snames,
                        'unique_stim': sorted(set(snames)),
                    }
            except Exception:
                pass

    f.close()
    return unit_spikes, areas, stim_info


def apply_firing_rate_filter(unit_spikes, areas):
    """Return active spike list + active areas list."""
    if len(unit_spikes) == 0:
        return [], []
    all_sp = np.concatenate([s for s in unit_spikes if len(s) > 0])
    if len(all_sp) == 0:
        return [], []
    t_total = all_sp.max() - all_sp.min()
    if t_total <= 0:
        return [], []
    rates = np.array([len(s) / t_total for s in unit_spikes])
    mask = rates >= MIN_FIRING_RATE
    active_spikes = [unit_spikes[i] for i in range(len(unit_spikes)) if mask[i]]
    active_areas = [areas[i] for i in range(len(areas)) if mask[i]]
    return active_spikes, active_areas


def build_conditions(stim_info):
    """Rebuild conditions dict exactly as the main pipeline does."""
    conditions = {}
    for key, info in stim_info.items():
        if 'names' in info:
            for stim_type in info.get('unique_stim', []):
                mask = np.array(info['names']) == stim_type
                starts = info['starts'][mask]
                stops = info['stops'][mask]
                total_dur = float(np.sum(stops - starts))
                if total_dur > 10:
                    conditions[f"{key}_{stim_type}"] = {
                        'starts': starts, 'stops': stops,
                        'total_duration': total_dur,
                    }
    return conditions


def build_spike_index(unit_spikes):
    """Pre-sort all spikes with unit labels for fast searchsorted-based binning.
    Returns (sorted_times, sorted_unit_ids) — call once per session.
    """
    times_list = []
    uid_list = []
    for i, spikes in enumerate(unit_spikes):
        if len(spikes) > 0:
            times_list.append(spikes)
            uid_list.append(np.full(len(spikes), i, dtype=np.int32))
    if not times_list:
        return np.array([]), np.array([], dtype=np.int32)
    all_times = np.concatenate(times_list)
    all_uids = np.concatenate(uid_list)
    order = np.argsort(all_times)
    return all_times[order], all_uids[order]


def compute_spike_counts(sorted_times, sorted_uids, n_units, starts, stops,
                         bin_size=BIN_SIZE_SEC):
    """Bin spike counts with no Python loops over units.
    Uses searchsorted to slice each epoch, then 2D bincount via flat index.
    """
    all_cols = []
    for t_start, t_stop in zip(starts, stops):
        dur = t_stop - t_start
        if dur < bin_size:
            continue
        n_bins = int(dur / bin_size)
        lo = np.searchsorted(sorted_times, t_start)
        hi = np.searchsorted(sorted_times, t_stop)
        if hi <= lo:
            all_cols.append(np.zeros((n_units, n_bins), dtype=np.float32))
            continue
        epoch_times = sorted_times[lo:hi]
        epoch_uids = sorted_uids[lo:hi]
        tb = np.clip(((epoch_times - t_start) / bin_size).astype(np.int32), 0, n_bins - 1)
        flat = epoch_uids * n_bins + tb
        col = np.bincount(flat, minlength=n_units * n_bins).reshape(n_units, n_bins)
        all_cols.append(col.astype(np.float32))
    if not all_cols:
        return None
    return np.concatenate(all_cols, axis=1)


def compute_area_wj(corr_a, corr_b, active_areas):
    """Compute per-area WJ for areas with >= MIN_UNITS units."""
    area_series = pd.Series(active_areas)
    area_counts = area_series.value_counts()
    area_wj = {}
    for area, count in area_counts.items():
        if area == 'unknown' or count < MIN_UNITS:
            continue
        idx = [i for i, a in enumerate(active_areas) if a == area]
        sub_a = corr_a[np.ix_(idx, idx)]
        sub_b = corr_b[np.ix_(idx, idx)]
        div = implementation_divergence(sub_a, sub_b)
        area_wj[area] = div
    return area_wj


def process_session(ckpt_path, nwb_path):
    """Add area_wj to all results in one session checkpoint."""
    session_id = os.path.basename(ckpt_path).replace('_checkpoint.json', '')
    log(f"\nSession: {session_id}")

    with open(ckpt_path) as f:
        ckpt = json.load(f)

    results = ckpt.get('results', [])
    if not results:
        log(f"  No results — skipping")
        return

    # Check if area_wj is already populated with real areas
    sample_awj = results[0].get('area_wj', {})
    if sample_awj and list(sample_awj.keys()) != ['unknown']:
        log(f"  Area WJ already populated — skipping")
        return

    # Load NWB
    unit_spikes, areas, stim_info = load_nwb_with_areas(nwb_path)
    unique_areas = sorted(set(areas))
    log(f"  Areas: {unique_areas}")

    if all(a == 'unknown' for a in unique_areas):
        log(f"  WARNING: All areas still unknown — no electrode table match")
        return

    # Apply firing rate filter
    active_spikes, active_areas = apply_firing_rate_filter(unit_spikes, areas)
    log(f"  Active units: {len(active_spikes)}")
    area_counts = pd.Series(active_areas).value_counts()
    eligible = {a: c for a, c in area_counts.items()
                if a != 'unknown' and c >= MIN_UNITS}
    log(f"  Eligible areas (>={MIN_UNITS} units): {eligible}")

    if not eligible:
        log(f"  No area meets MIN_UNITS threshold — skipping")
        return

    # Rebuild conditions dict
    conditions = build_conditions(stim_info)

    # Build sorted spike index once for this session
    log(f"  Building spike index...")
    sorted_times, sorted_uids = build_spike_index(active_spikes)
    n_active = len(active_spikes)

    # Pre-compute correlation matrix for each unique condition ONCE
    needed_conds = set()
    for r in results:
        needed_conds.add(r.get('condition_a'))
        needed_conds.add(r.get('condition_b'))
    needed_conds.discard(None)

    log(f"  Pre-computing matrices for {len(needed_conds)} unique conditions...")
    corr_cache = {}
    for cond in sorted(needed_conds):
        if cond not in conditions:
            continue
        info = conditions[cond]
        durs = np.asarray(info['stops']) - np.asarray(info['starts'])
        if durs.size == 0 or durs.max() < BIN_SIZE_SEC:
            continue
        counts = compute_spike_counts(sorted_times, sorted_uids, n_active,
                                      info['starts'], info['stops'])
        if counts is None:
            continue
        corr_cache[cond] = fast_spearman_matrix(counts)
        log(f"    {cond}: matrix {corr_cache[cond].shape}")

    updated = 0
    for r in results:
        cond_a = r.get('condition_a')
        cond_b = r.get('condition_b')
        if cond_a not in corr_cache or cond_b not in corr_cache:
            continue
        area_wj = compute_area_wj(corr_cache[cond_a], corr_cache[cond_b], active_areas)
        r['area_wj'] = area_wj
        updated += 1

    log(f"  Updated {updated}/{len(results)} comparisons with area_wj")

    # Save updated checkpoint
    with open(ckpt_path, 'w') as f:
        json.dump(ckpt, f)
    log(f"  Checkpoint saved")


def main():
    log("=" * 70)
    log("AREA DECOMPOSITION — Retroactive per-area WJ")
    log("=" * 70)

    checkpoints = sorted(glob.glob(os.path.join(RESULTS_DIR, '*_checkpoint.json')))
    log(f"Found {len(checkpoints)} checkpoints")

    # Build session_id -> nwb_path map
    nwb_files = glob.glob(os.path.join(DATA_DIR, '*.nwb'))
    nwb_map = {}
    for nf in nwb_files:
        base = os.path.basename(nf)
        # filename: sub-XXXXX_ses-YYYYYYY.nwb
        # session_id: ses-YYYYYYY
        parts = base.replace('.nwb', '').split('_')
        for p in parts:
            if p.startswith('ses-'):
                nwb_map[p] = nf

    log(f"NWB files indexed: {len(nwb_map)}")

    n_done = 0
    n_skipped = 0
    for ckpt_path in checkpoints:
        session_id = os.path.basename(ckpt_path).replace('_checkpoint.json', '')
        nwb_path = nwb_map.get(session_id)
        if nwb_path is None:
            log(f"\nNo NWB found for {session_id} — skipping")
            n_skipped += 1
            continue
        process_session(ckpt_path, nwb_path)
        n_done += 1

    elapsed = (time.time() - START) / 60
    log(f"\n{'='*70}")
    log(f"AREA DECOMPOSITION COMPLETE")
    log(f"Processed: {n_done} sessions, Skipped: {n_skipped}")
    log(f"Total time: {elapsed:.1f} minutes")


if __name__ == '__main__':
    main()
