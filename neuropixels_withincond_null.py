"""
Pipeline: Within-condition split-half null for the sign-flip rate (corrects the 50% strawman)
Author: Drake H. Harbert (D.H.H.) | ORCID 0009-0007-7740-3616 | 2026-06-30
Description:
    The headline sign-preservation claim must be tested against a structure-preserving
    null, NOT the 50% chance rate (the strawman corrected in the PTSD/JPR paper and now
    in wj-methodology.md Layer 2I). The correct null is the WITHIN-condition split-half
    sign-flip rate: split each condition's bins into two independent halves, compute a
    correlation matrix for each half, and measure how often a pair reverses sign between
    two samples of the SAME condition (i.e., from sampling noise alone). The observed
    BETWEEN-condition flip rate is then compared to this within-condition null, with
    sessions as the replication unit (Wilcoxon across sessions, per magnitude stratum).
    Observed ~ null => sign preserved beyond sampling; observed > null => real sign
    reorganization. Reuses the loader from neuropixels_layer2i_2f_corrected.py.
Input:  data/*.nwb  +  results/layer2i_2f/layer2i_signflip_stratified.csv (observed)
Output: results/layer2i_2f/withincond_null/
"""
import os, sys, glob, json, gc, warnings
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, r"G:\My Drive\inner_architecture_research\neuropixels_wj")
import neuropixels_layer2i_2f_corrected as base   # __main__-guarded, safe to import

SEED = 42; np.random.seed(SEED)
STRATA = base.STRATA
BASE = r"G:\My Drive\inner_architecture_research\neuropixels_wj"
DATA = os.path.join(BASE, "data")
OBS = os.path.join(BASE, "results", "layer2i_2f", "layer2i_signflip_stratified.csv")
OUT = os.path.join(BASE, "results", "layer2i_2f", "withincond_null"); os.makedirs(OUT, exist_ok=True)
def log(m): print(m, flush=True)

def flip_rates(ca, cb):
    iu = np.triu_indices(ca.shape[0], k=1)
    ra, rb = ca[iu], cb[iu]
    minabs = np.minimum(np.abs(ra), np.abs(rb)); flip = np.sign(ra) != np.sign(rb)
    out = {}
    for thr in STRATA:
        m = minabs >= thr; npr = int(m.sum())
        out[thr] = (flip[m].mean() if npr >= 20 else np.nan)
    return out

def main():
    rng = np.random.RandomState(SEED)
    rows = []
    for nwb in sorted(glob.glob(os.path.join(DATA, "*.nwb"))):
        sid = os.path.basename(nwb).split(".")[0]
        try:
            unit_spikes = base.load_nwb_session(nwb)
            active = base.active_filter(unit_spikes)
            if len(active) < 30:
                log(f"SKIP {sid}"); continue
            conds = base.get_conditions(nwb)
            for nm, cd in conds.items():
                parts = []
                for s, e in zip(cd["starts"], cd["stops"]):
                    sc = base.compute_spike_counts(active, s, e)
                    if sc is not None: parts.append(sc)
                if not parts: continue
                counts = np.hstack(parts)
                nb = counts.shape[1]
                if nb < 40: continue
                perm = rng.permutation(nb); h = nb // 2
                c1 = base.fast_spearman_matrix(counts[:, perm[:h]])
                c2 = base.fast_spearman_matrix(counts[:, perm[h:]])
                fr = flip_rates(c1, c2)
                rows.append({"session": sid, "condition": nm, **{f"null_flip_{t}": fr[t] for t in STRATA}})
            log(f"{sid}: within-condition null done ({len(active)} units, {len(conds)} conditions)")
            del active, unit_spikes; gc.collect()
            pd.DataFrame(rows).to_csv(os.path.join(OUT, "withincond_null_per_condition.csv"), index=False)
        except Exception as e:
            log(f"ERROR {sid}: {e}")
    nulldf = pd.DataFrame(rows)
    # per-session mean within-condition null
    null_sess = nulldf.groupby("session")[[f"null_flip_{t}" for t in STRATA]].mean()
    # observed between-condition flip rate, per session
    obs = pd.read_csv(OBS)
    obs_sess = obs.pivot_table(index="session", columns="threshold", values="sign_flip_rate", aggfunc="mean")
    summary = {}
    for t in STRATA:
        ocol = obs_sess[t] if t in obs_sess.columns else None
        ncol = null_sess[f"null_flip_{t}"]
        merged = pd.concat([ocol.rename("obs"), ncol.rename("null")], axis=1).dropna()
        if len(merged) >= 5:
            w = wilcoxon(merged["obs"], merged["null"])
            summary[f"thr_{t}"] = {
                "n_sessions": int(len(merged)),
                "observed_between_cond_flip": round(float(merged["obs"].mean()), 4),
                "within_cond_null_flip": round(float(merged["null"].mean()), 4),
                "observed_minus_null": round(float((merged["obs"] - merged["null"]).mean()), 4),
                "wilcoxon_p_obs_vs_null": round(float(w.pvalue), 5),
                "observed_exceeds_null": bool((merged["obs"] > merged["null"]).mean() > 0.5)}
    prov = {"null": "within-condition split-half (structure-preserving); 50% binomial retired",
            "replication_unit": "session", "seed": SEED, "by_stratum": summary,
            "note": "observed ~ within-condition null => reorganization preserves sign beyond sampling; "
                    "observed > null => genuine sign reorganization at that magnitude stratum"}
    json.dump(prov, open(os.path.join(OUT, "withincond_null_summary.json"), "w"), indent=2)
    log(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
