"""
Pipeline: Neuropixels per-unit Layer 2F x waveform shape, threshold-free / training-free
Author: Drake H. Harbert (D.H.H.)
Affiliation: Inner Architecture LLC, Canton, OH
ORCID: 0009-0007-7740-3616
Date: 2026-06-30
Description:
    Runs AFTER neuropixels_layer2i_2f_corrected.py has saved the per-condition
    correlation matrices. Asks whether a neuron's extracellular WAVEFORM SHAPE tracks
    its correlation-architecture reorganization, WITHOUT imposing any category.

    No clustering. No trained model. No imposed cell-type label. No fixed magnitude
    threshold. Everything is continuous and assumption-light:
      - Per unit (continuous): row reorganization (1 - row weighted Jaccard on |r|) and
        sign coherence (fraction of partners preserving sign), the latter computed across
        a SWEEP of magnitude floors (the stratification axis IS the measurement, not a cut).
      - Association: Spearman rank correlation between each raw waveform feature and each
        per-unit reorganization/coherence metric, across all units, with a permutation
        null (shuffle waveform values across units) and a bootstrap CI. Rank-based, no
        training, no parameters fit to the data.

    Categories are findings, not inputs. If spike width tracks sign reorganization, it
    appears as a continuous, orderly relationship, not as two imposed groups.
Dependencies: numpy, pandas, scipy, h5py
Input:  data/*.nwb  +  results/layer2i_2f/matrices/*.npz
Output: results/layer2i_2f/waveform_association/
"""
import os, sys, glob, json, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")

RANDOM_SEED = 42
MIN_FIRING_RATE = 0.5
FLOORS = [0.0, 0.05, 0.10, 0.15, 0.20]     # swept magnitude axis (not a threshold)
N_PERM = 2000
N_BOOT = 1000
WAVE_FEATS = ["waveform_duration", "waveform_halfwidth", "PT_ratio",
              "repolarization_slope", "recovery_slope", "spread",
              "velocity_above", "velocity_below"]
BASE = r"G:\My Drive\inner_architecture_research\neuropixels_wj"
DATA = os.path.join(BASE, "data")
MAT = os.path.join(BASE, "results", "layer2i_2f", "matrices")
OUT = os.path.join(BASE, "results", "layer2i_2f", "waveform_association")
os.makedirs(OUT, exist_ok=True)
rng = np.random.RandomState(RANDOM_SEED)
def log(m): print(m, flush=True)

def load_units(nwb_path):
    import h5py
    f = h5py.File(nwb_path, "r"); u = f["units"]
    idx = u["spike_times_index"][:]; data = u["spike_times"][:]
    n = len(idx); spikes = []; prev = 0
    for i in range(n):
        spikes.append(data[prev:idx[i]]); prev = idx[i]
    cat = np.concatenate([s for s in spikes if len(s)]); tot = cat.max() - cat.min()
    fr = np.array([len(s) / tot if tot > 0 else 0 for s in spikes])
    keep = fr >= MIN_FIRING_RATE
    feats = {k: (u[k][:] if k in u else np.full(n, np.nan)) for k in WAVE_FEATS}
    f.close()
    return keep, pd.DataFrame(feats)

def row_reorg(ca, cb):
    A = np.abs(ca).copy(); B = np.abs(cb).copy()
    np.fill_diagonal(A, 0.0); np.fill_diagonal(B, 0.0)
    den = np.maximum(A, B).sum(1); num = np.minimum(A, B).sum(1)
    return 1.0 - np.where(den > 0, num / den, np.nan)

def row_coh(ca, cb, floor):
    same = np.sign(ca) == np.sign(cb)
    mean = np.minimum(np.abs(ca), np.abs(cb)) >= floor
    np.fill_diagonal(mean, False)
    num = (same & mean).sum(1).astype(float); den = mean.sum(1).astype(float)
    return np.where(den >= 5, num / den, np.nan)

def main():
    nwbs = sorted(glob.glob(os.path.join(DATA, "*.nwb")))
    rows = []
    for nwb in nwbs:
        sid = os.path.basename(nwb).split(".")[0]
        mats = sorted(glob.glob(os.path.join(MAT, f"{sid}__*.npz")))
        if len(mats) < 2:
            continue
        keep, feat_all = load_units(nwb)
        ai = np.where(keep)[0]; n = len(ai)
        C = [np.load(m)["corr"] for m in mats]
        if C[0].shape[0] != n:
            log(f"SKIP {sid}: n={n} != matrix {C[0].shape[0]}"); continue
        feat = feat_all.iloc[ai].reset_index(drop=True)
        acc = {"reorg": np.zeros(n)}
        for fl in FLOORS: acc[f"coh_{fl}"] = np.zeros(n)
        cnt = np.zeros(n); cntc = {fl: np.zeros(n) for fl in FLOORS}
        for i in range(len(C)):
            for j in range(i + 1, len(C)):
                rr = row_reorg(C[i], C[j]); ok = ~np.isnan(rr)
                acc["reorg"][ok] += rr[ok]; cnt[ok] += 1
                for fl in FLOORS:
                    ch = row_coh(C[i], C[j], fl); okc = ~np.isnan(ch)
                    acc[f"coh_{fl}"][okc] += ch[okc]; cntc[fl][okc] += 1
        df = feat.copy(); df["session"] = sid
        df["mean_row_reorg"] = np.where(cnt > 0, acc["reorg"] / np.maximum(cnt, 1), np.nan)
        for fl in FLOORS:
            df[f"sign_coh_floor_{fl}"] = np.where(cntc[fl] > 0, acc[f"coh_{fl}"] / np.maximum(cntc[fl], 1), np.nan)
        df["n_comparisons"] = cnt
        rows.append(df); log(f"{sid}: {int((cnt>0).sum())} units")
    allu = pd.concat(rows, ignore_index=True)
    allu.to_csv(os.path.join(OUT, "perunit_waveform_continuous.csv"), index=False)
    log(f"TOTAL units: {len(allu)}")

    metrics = ["mean_row_reorg"] + [f"sign_coh_floor_{fl}" for fl in FLOORS]
    assoc = []
    for feat in WAVE_FEATS:
        for met in metrics:
            sub = allu[[feat, met]].apply(pd.to_numeric, errors="coerce").dropna()
            if len(sub) < 50:
                continue
            x = sub[feat].values; y = sub[met].values
            rho = spearmanr(x, y).correlation
            # permutation null: shuffle waveform values across units
            cnt_ge = 0
            for _ in range(N_PERM):
                if abs(spearmanr(rng.permutation(x), y).correlation) >= abs(rho):
                    cnt_ge += 1
            p = (cnt_ge + 1) / (N_PERM + 1)
            # bootstrap CI on rho
            boots = []
            idxs = np.arange(len(x))
            for _ in range(N_BOOT):
                bi = rng.choice(idxs, len(idxs), replace=True)
                boots.append(spearmanr(x[bi], y[bi]).correlation)
            lo, hi = np.percentile(boots, [2.5, 97.5])
            assoc.append({"waveform_feature": feat, "metric": met, "n": int(len(sub)),
                          "spearman_rho": round(float(rho), 4), "perm_p": round(float(p), 4),
                          "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4)})
    ad = pd.DataFrame(assoc).sort_values("perm_p")
    ad.to_csv(os.path.join(OUT, "waveform_reorg_association.csv"), index=False)
    prov = {"methodology": "threshold-free, training-free continuous association; no imposed categories",
            "n_units": int(len(allu)), "floors_swept": FLOORS, "n_perm": N_PERM, "n_boot": N_BOOT,
            "random_seed": RANDOM_SEED,
            "note": ("Spearman rho between raw waveform shape and per-unit reorganization/coherence, "
                     "coherence stratified across swept magnitude floors. No clustering, no model fit, "
                     "no cell-type label. Significant continuous association = spike shape tracks "
                     "reorganization; report by measured waveform properties only.")}
    json.dump(prov, open(os.path.join(OUT, "provenance.json"), "w"), indent=2)
    log("Top associations (by perm p):")
    log(ad.head(12).to_string(index=False))

if __name__ == "__main__":
    main()
