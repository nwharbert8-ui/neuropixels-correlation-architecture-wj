"""
Pipeline: Neuropixels finalize battery (per-session inference, scalar baseline,
          reorg-vs-coherence joint, per-neuron effective coupling regimes)
Author: Drake H. Harbert (D.H.H.) | ORCID 0009-0007-7740-3616 | 2026-06-30
Description:
    Runs on the saved correlation matrices + the per-unit waveform CSV. All
    threshold-free / training-free; inference respects session structure.

    A. Per-SESSION waveform association: Spearman(waveform feature, per-unit metric)
       computed within each session, then aggregated across sessions (Wilcoxon of the
       per-session rho vs 0, fraction same sign). Locks the inference that pooling
       55k units inflated.
    B. Reorg-vs-coherence joint: per-session Spearman between a unit's row
       reorganization and its sign coherence (the two-number relationship).
    C. Scalar baseline (framework-family contrast): per condition, the field-standard
       scalar (mean |r|); per condition-pair, the scalar change vs the architecture
       reorganization (1 - unsigned WJ). Shows what the scalar collapses.
    D. Per-neuron effective coupling regimes: participation ratio of a neuron's
       coupling profile across conditions (threshold-free count of distinct coupling
       patterns; ~1 = one stable regime, higher = multiple condition-specific regimes).
Dependencies: numpy, pandas, scipy
Input:  results/layer2i_2f/matrices/*.npz  +  results/layer2i_2f/waveform_association/perunit_waveform_continuous.csv
Output: results/layer2i_2f/finalize/
"""
import os, glob, json, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
warnings.filterwarnings("ignore")
RANDOM_SEED = 42; np.random.seed(RANDOM_SEED)
BASE = r"G:\My Drive\inner_architecture_research\neuropixels_wj"
MAT = os.path.join(BASE, "results", "layer2i_2f", "matrices")
PERUNIT = os.path.join(BASE, "results", "layer2i_2f", "waveform_association", "perunit_waveform_continuous.csv")
OUT = os.path.join(BASE, "results", "layer2i_2f", "finalize"); os.makedirs(OUT, exist_ok=True)
def log(m): print(m, flush=True)

def per_session_assoc(df, feat, metric):
    rhos = []
    for sid, g in df.groupby("session"):
        s = g[[feat, metric]].apply(pd.to_numeric, errors="coerce").dropna()
        if len(s) >= 30:
            r = spearmanr(s[feat], s[metric]).correlation
            if not np.isnan(r): rhos.append(r)
    rhos = np.array(rhos)
    if len(rhos) < 5: return None
    w = wilcoxon(rhos)
    return {"feature": feat, "metric": metric, "n_sessions": len(rhos),
            "mean_rho": round(float(rhos.mean()), 4), "median_rho": round(float(np.median(rhos)), 4),
            "frac_same_sign": round(float((np.sign(rhos) == np.sign(rhos.mean())).mean()), 3),
            "wilcoxon_p": round(float(w.pvalue), 5)}

def upper(m): return m[np.triu_indices(m.shape[0], k=1)]
def wj_unsigned(a, b):
    A, B = np.abs(a), np.abs(b); d = np.maximum(A, B).sum()
    return float(np.minimum(A, B).sum() / d) if d > 0 else np.nan

def main():
    df = pd.read_csv(PERUNIT)
    feats = ["waveform_duration", "waveform_halfwidth", "spread", "PT_ratio",
             "recovery_slope", "repolarization_slope"]
    out = {}

    # ---- A: per-session waveform association ----
    A = []
    for f in feats:
        for m in ["mean_row_reorg", "sign_coh_floor_0.1", "sign_coh_floor_0.2"]:
            r = per_session_assoc(df, f, m)
            if r: A.append(r)
    pd.DataFrame(A).to_csv(os.path.join(OUT, "A_persession_waveform_assoc.csv"), index=False)
    out["A_strongest"] = sorted(A, key=lambda x: -abs(x["mean_rho"]))[:6]

    # ---- B: reorg vs coherence joint ----
    B = [per_session_assoc(df, "mean_row_reorg", m) for m in ["sign_coh_floor_0.1", "sign_coh_floor_0.2"]]
    out["B_reorg_vs_coherence"] = [b for b in B if b]

    # ---- C + D: load matrices per session ----
    sessions = sorted(set(os.path.basename(p).split("__")[0] for p in glob.glob(os.path.join(MAT, "*.npz"))))
    scalar_rows = []; pr_rows = []
    for sid in sessions:
        files = sorted(glob.glob(os.path.join(MAT, f"{sid}__*.npz")))
        if len(files) < 2: continue
        C = [np.load(f)["corr"] for f in files]
        n = C[0].shape[0]
        meanabs = [float(np.abs(upper(c)).mean()) for c in C]
        # C: scalar change vs architecture reorganization, per condition pair
        for i in range(len(C)):
            for j in range(i + 1, len(C)):
                scalar_rows.append({"session": sid,
                                    "scalar_change_mean_absr": abs(meanabs[i] - meanabs[j]),
                                    "architecture_reorg": 1 - wj_unsigned(upper(C[i]), upper(C[j]))})
        # D: per-neuron participation ratio of coupling across conditions
        stack = np.stack([c for c in C], axis=0)   # (n_cond, n, n)
        for u in range(n):
            M = stack[:, u, :]                      # (n_cond, n) coupling profile across conditions
            M = np.delete(M, u, axis=1)             # drop self
            Mc = M - M.mean(0, keepdims=True)
            G = Mc @ Mc.T                            # (n_cond, n_cond)
            ev = np.linalg.eigvalsh(G); ev = ev[ev > 1e-12]
            if len(ev) == 0: continue
            pr = (ev.sum() ** 2) / (np.square(ev).sum())   # participation ratio
            pr_rows.append({"session": sid, "n_conditions": len(C), "participation_ratio": float(pr)})
        log(f"{sid}: scalar+PR done ({n} units, {len(C)} conditions)")
    sc = pd.DataFrame(scalar_rows); pr = pd.DataFrame(pr_rows)
    sc.to_csv(os.path.join(OUT, "C_scalar_vs_architecture.csv"), index=False)
    pr.to_csv(os.path.join(OUT, "D_perneuron_coupling_regimes.csv"), index=False)
    out["C_scalar_baseline"] = {
        "median_scalar_change_mean_absr": round(float(sc["scalar_change_mean_absr"].median()), 4),
        "median_architecture_reorg": round(float(sc["architecture_reorg"].median()), 4),
        "spearman_scalar_vs_architecture": round(float(spearmanr(sc["scalar_change_mean_absr"], sc["architecture_reorg"]).correlation), 4),
        "note": "small scalar change + large architecture reorg + weak coupling = scalar collapses what WJ resolves"}
    out["D_coupling_regimes"] = {
        "n_neurons": int(len(pr)), "median_participation_ratio": round(float(pr["participation_ratio"].median()), 3),
        "mean_participation_ratio": round(float(pr["participation_ratio"].mean()), 3),
        "frac_pr_gt_2": round(float((pr["participation_ratio"] > 2).mean()), 3),
        "max_possible": int(pr["n_conditions"].max()),
        "note": "PR ~1 = single stable coupling regime; >2 = multiple condition-specific regimes (candidate input channels)"}
    json.dump(out, open(os.path.join(OUT, "finalize_summary.json"), "w"), indent=2)
    log(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
