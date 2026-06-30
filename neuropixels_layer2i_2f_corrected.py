"""
Pipeline: Neuropixels Layer 2I (direct sign-flip) + Layer 2F (per-unit decomposition), corrected
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-06-30
Description:
    Corrects the sign-inversion headline of the Neuropixels manuscript and adds the
    missing Layer 2F per-unit row decomposition, in a single NWB-reprocessing pass.

    (1) Layer 2I (Hard-Stop #11 compliance): the prior manuscript reported the
        gap-derived sign_inversion_pct = gap/(1-wj_unsigned) as a "sign-flip rate."
        That is prohibited. This pipeline computes the DIRECT pair-level sign-flip
        rate (sign(r_A) != sign(r_B)) across ALL condition comparisons, STRATIFIED
        by min(|r_A|,|r_B|), and tests each stratum against the 50% chance rate
        (binomial). The stratified trajectory is the substrate signature.

    (2) Layer 2F (per-unit row decomposition): for each unit, its coupling-profile
        reorganization (1 - row weighted Jaccard on |r|) and its sign coherence
        (fraction of meaningful-magnitude partners that preserve correlation sign),
        with the four-quadrant classification (stable hub / coherent magnitude /
        split learner / incoherent noise).

    (3) Saves each condition's unit-by-unit Spearman matrix to disk (npz) so future
        analyses never need to reprocess the raw NWB files again.

    Data: Allen Brain Observatory Visual Coding Neuropixels (DANDI 000021).
    Fundamental unit: individual spike-sorted unit. Full pairwise matrix. Spearman.
Dependencies: numpy, pandas, scipy, h5py, wj_utils
Input:  neuropixels_wj/data/*.nwb
Output: neuropixels_wj/results/layer2i_2f/ (per-comparison CSVs, matrices, provenance.json)
"""
import os, sys, time, json, glob, gc, warnings
import numpy as np
import pandas as pd
from scipy.stats import binomtest
warnings.filterwarnings("ignore")

sys.path.insert(0, r"G:\My Drive\inner_architecture_research")
from wj_utils import fast_spearman_matrix

# ============================================================== CONFIG
RANDOM_SEED = 42
FORCE_RECOMPUTE = True
BIN_SIZE_SEC = 0.1            # 100 ms bins (matches main pipeline)
MIN_FIRING_RATE = 0.5        # Hz
STRATA = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]   # min(|r_A|,|r_B|) thresholds
COHERENCE_FLOOR = 0.10       # only count partners with min|r| >= this in sign coherence
REORG_HI = 0.50              # row reorganization split for quadrant classification
COH_HI = 0.90                # sign-coherence high cut
COH_MID = 0.60               # below this (toward 0.5) = split/incoherent
BASE_DIR = r"G:\My Drive\inner_architecture_research\neuropixels_wj"
DATA_DIR = os.path.join(BASE_DIR, "data")
OUT_DIR = os.path.join(BASE_DIR, "results", "layer2i_2f")
MAT_DIR = os.path.join(OUT_DIR, "matrices")
for d in (OUT_DIR, MAT_DIR):
    os.makedirs(d, exist_ok=True)
np.random.seed(RANDOM_SEED)
START = time.time()
def log(m): print(f"[{(time.time()-START)/60:6.1f}m] {m}", flush=True)

# ============================================================== NWB LOADER (from main pipeline)
def load_nwb_session(nwb_path):
    import h5py
    f = h5py.File(nwb_path, "r")
    units = f["units"]
    idx = units["spike_times_index"][:]; data = units["spike_times"][:]
    n_units = len(idx); unit_spikes = []; prev = 0
    for i in range(n_units):
        unit_spikes.append(data[prev:idx[i]]); prev = idx[i]
    stim_info = {}
    if "intervals" in f:
        for ik in f["intervals"].keys():
            if ik == "units":
                continue
            try:
                iv = f["intervals"][ik]; starts = iv["start_time"][:]; stops = iv["stop_time"][:]
                nk = [k for k in iv.keys() if "stimulus" in k.lower() or "name" in k.lower() or "type" in k.lower()]
                if nk:
                    nm = iv[nk[0]][:]
                    if len(nm) and isinstance(nm[0], bytes):
                        nm = [s.decode() for s in nm]
                    stim_info[ik] = {"starts": starts, "stops": stops, "names": nm,
                                     "unique_stim": sorted(set(nm))}
            except Exception as e:
                log(f"    interval {ik}: {e}")
    f.close()
    return unit_spikes

def compute_spike_counts(unit_spikes, t0, t1, bin_size=BIN_SIZE_SEC):
    n_bins = int((t1 - t0) / bin_size)
    if n_bins < 1:
        return None
    counts = np.zeros((len(unit_spikes), n_bins), dtype=np.float32)
    for i, sp in enumerate(unit_spikes):
        w = sp[(sp >= t0) & (sp < t1)]
        if len(w):
            b = np.clip(((w - t0) / bin_size).astype(int), 0, n_bins - 1)
            np.add.at(counts[i], b, 1)
    return counts

def get_conditions(nwb_path):
    import h5py
    f = h5py.File(nwb_path, "r"); stim = {}
    if "intervals" in f:
        for ik in f["intervals"].keys():
            if ik == "units":
                continue
            try:
                iv = f["intervals"][ik]; starts = iv["start_time"][:]; stops = iv["stop_time"][:]
                nk = [k for k in iv.keys() if "stimulus" in k.lower() or "name" in k.lower() or "type" in k.lower()]
                if not nk:
                    continue
                nm = iv[nk[0]][:]
                if len(nm) and isinstance(nm[0], bytes):
                    nm = [s.decode() for s in nm]
                nm = np.array(nm)
                for st in sorted(set(nm.tolist())):
                    m = nm == st; tot = float(np.sum(stops[m] - starts[m]))
                    if tot > 10:
                        stim[f"{ik}_{st}"] = {"starts": starts[m], "stops": stops[m], "dur": tot}
            except Exception:
                pass
    f.close()
    return stim

# ============================================================== WJ HELPERS
def wj_vec(a, b):
    """Unsigned weighted Jaccard on two vectors of correlations."""
    aa, bb = np.abs(a), np.abs(b)
    den = np.maximum(aa, bb).sum()
    return float(np.minimum(aa, bb).sum() / den) if den > 0 else np.nan

def active_filter(unit_spikes):
    cat = np.concatenate([s for s in unit_spikes if len(s)])
    t = cat.max() - cat.min()
    fr = np.array([len(s) / t if t > 0 else 0 for s in unit_spikes])
    keep = fr >= MIN_FIRING_RATE
    return [unit_spikes[i] for i in range(len(unit_spikes)) if keep[i]]

def cond_corr(active, cond):
    parts = []
    for s, e in zip(cond["starts"], cond["stops"]):
        sc = compute_spike_counts(active, s, e)
        if sc is not None:
            parts.append(sc)
    if not parts:
        return None
    counts = np.hstack(parts)
    if counts.shape[1] < 10:
        return None
    return fast_spearman_matrix(counts)

# ============================================================== PER-SESSION
def process_session(nwb_path):
    sid = os.path.basename(nwb_path).split(".")[0]
    log(f"=== {sid} ===")
    unit_spikes = load_nwb_session(nwb_path)
    active = active_filter(unit_spikes)
    n = len(active); log(f"  active units: {n}")
    conds = get_conditions(nwb_path)
    names = sorted(conds, key=lambda k: -conds[k]["dur"])
    if len(names) < 2 or n < 30:
        log("  SKIP (insufficient conditions/units)"); return [], []
    # compute + cache + save each condition's correlation matrix once
    corr = {}
    iu = np.triu_indices(n, k=1)
    for nm in names:
        c = cond_corr(active, conds[nm])
        if c is not None:
            corr[nm] = c
            safe = nm.replace("/", "_").replace(" ", "_")
            np.savez_compressed(os.path.join(MAT_DIR, f"{sid}__{safe}.npz"), corr=c.astype(np.float32))
    usable = [nm for nm in names if nm in corr]
    log(f"  usable conditions: {len(usable)}")
    rows2i, rows2f = [], []
    for i in range(len(usable)):
        for j in range(i + 1, len(usable)):
            a, b = usable[i], usable[j]
            ca, cb = corr[a], corr[b]
            ra, rb = ca[iu], cb[iu]
            same = a.split("_presentations")[0].rsplit("_", 1)[0] == b.split("_presentations")[0].rsplit("_", 1)[0]
            # ---- Layer 2I: direct stratified sign-flip rate ----
            minabs = np.minimum(np.abs(ra), np.abs(rb))
            flip = np.sign(ra) != np.sign(rb)
            for thr in STRATA:
                m = minabs >= thr; npr = int(m.sum())
                if npr < 20:
                    continue
                nfl = int(flip[m].sum()); rate = nfl / npr
                p = binomtest(nfl, npr, 0.5, alternative="two-sided").pvalue
                rows2i.append({"session": sid, "cond_a": a, "cond_b": b, "same_stim": same,
                               "threshold": thr, "n_pairs": npr, "sign_flip_rate": rate,
                               "n_flips": nfl, "binom_p_vs_chance": p,
                               "below_chance": rate < 0.5})
            # ---- Layer 2F: per-unit row decomposition ----
            reorg = np.empty(n); coh = np.full(n, np.nan)
            for u in range(n):
                row_a = np.delete(ca[u], u); row_b = np.delete(cb[u], u)
                reorg[u] = 1.0 - wj_vec(row_a, row_b)
                mm = np.minimum(np.abs(row_a), np.abs(row_b)) >= COHERENCE_FLOOR
                if mm.sum() >= 5:
                    coh[u] = float((np.sign(row_a[mm]) == np.sign(row_b[mm])).mean())
            valid = ~np.isnan(coh)
            if valid.sum() == 0:
                continue
            rv, cv = reorg[valid], coh[valid]
            def frac(mask): return float(mask.mean())
            stable_hub = (rv < REORG_HI) & (cv >= COH_HI)
            coherent_mag = (rv >= REORG_HI) & (cv >= COH_HI)
            split_learner = (rv >= REORG_HI) & (cv < COH_MID)
            incoherent = (rv < REORG_HI) & (cv < COH_MID)
            rows2f.append({"session": sid, "cond_a": a, "cond_b": b, "same_stim": same,
                           "n_units": int(valid.sum()), "mean_row_reorg": float(rv.mean()),
                           "mean_sign_coherence": float(cv.mean()),
                           "frac_stable_hub": frac(stable_hub),
                           "frac_coherent_magnitude": frac(coherent_mag),
                           "frac_split_learner": frac(split_learner),
                           "frac_incoherent": frac(incoherent)})
    del corr; gc.collect()
    return rows2i, rows2f

# ============================================================== MAIN
def main():
    nwbs = sorted(glob.glob(os.path.join(DATA_DIR, "*.nwb")))
    log(f"sessions found: {len(nwbs)}")
    all2i, all2f = [], []
    for p in nwbs:
        try:
            r2i, r2f = process_session(p)
            all2i += r2i; all2f += r2f
            pd.DataFrame(all2i).to_csv(os.path.join(OUT_DIR, "layer2i_signflip_stratified.csv"), index=False)
            pd.DataFrame(all2f).to_csv(os.path.join(OUT_DIR, "layer2f_perunit_quadrants.csv"), index=False)
        except Exception as e:
            log(f"  ERROR {os.path.basename(p)}: {e}")
    d2i = pd.DataFrame(all2i); d2f = pd.DataFrame(all2f)

    # headline summary: direct sign-flip rate by stratum (vs the deprecated gap metric)
    summary = {}
    for thr in STRATA:
        s = d2i[d2i["threshold"] == thr]["sign_flip_rate"]
        if len(s):
            summary[f"thr_{thr}"] = {"mean_sign_flip_rate": round(float(s.mean()), 4),
                                     "median": round(float(s.median()), 4),
                                     "frac_below_chance": round(float((s < 0.5).mean()), 4),
                                     "n_comparisons": int(len(s))}
    quad = {k: round(float(d2f[k].mean()), 4) for k in
            ["frac_stable_hub", "frac_coherent_magnitude", "frac_split_learner", "frac_incoherent"]} if len(d2f) else {}
    prov = {"methodology": "WJ-native; Layer 2I direct stratified sign-flip + Layer 2F per-unit",
            "fundamental_unit": "individual spike-sorted unit", "correlation_method": "Spearman",
            "random_seed": RANDOM_SEED, "n_comparisons_2I": int(d2i[["session","cond_a","cond_b"]].drop_duplicates().shape[0]) if len(d2i) else 0,
            "sign_flip_rate_by_stratum": summary,
            "layer2f_quadrant_fractions_mean": quad,
            "note": ("Direct sign-flip rate replaces the deprecated gap-derived sign_inversion_pct "
                     "(Hard-Stop #11). If rate is at/below 0.5 and falls toward 0 with stratification, "
                     "reorganization is magnitude-driven and sign is preserved.")}
    json.dump(prov, open(os.path.join(OUT_DIR, "provenance.json"), "w"), indent=2)
    log("DONE. Stratified sign-flip summary:")
    print(json.dumps(summary, indent=2))
    log("Layer 2F mean quadrant fractions:")
    print(json.dumps(quad, indent=2))

if __name__ == "__main__":
    main()
