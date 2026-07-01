# -*- coding: utf-8 -*-
"""Figures for the Neuropixels Journal of Neurophysiology manuscript.
Author: Drake H. Harbert (D.H.H.) | ORCID 0009-0007-7740-3616 | 2026-06-30
300 DPI, colorblind-safe. Reads results/layer2i_2f/."""
import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
B = r"G:/My Drive/inner_architecture_research/neuropixels_wj/results/layer2i_2f"
OUT = r"G:/My Drive/inner_architecture_research/neuropixels_wj/JNeurophysiol_Submission/Main_Figures"
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#999999"
def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight"); plt.close(fig)

# ---- Fig 1: observed sign-flip rate vs within-condition (sampling) null ----
nj = json.load(open(os.path.join(B, "withincond_null", "withincond_null_summary.json")))["by_stratum"]
thr = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
obs = [nj[f"thr_{t}"]["observed_between_cond_flip"] * 100 for t in thr]
nul = [max(nj[f"thr_{t}"]["within_cond_null_flip"] * 100, 0.005) for t in thr]  # floor for log axis
fig, ax = plt.subplots(figsize=(7.5, 5.5))
ax.plot(thr, obs, "o-", color=ORANGE, lw=2.5, ms=9, mec="black", mew=0.6, label="observed between-condition")
ax.plot(thr, nul, "s--", color=GREY, lw=2, ms=8, mec="black", mew=0.5, label="within-condition null (sampling floor)")
ax.set_yscale("log")
for t, o in zip(thr, obs):
    ax.annotate(f"{o:.2f}%", (t, o), textcoords="offset points", xytext=(4, 8), fontsize=8, color="#8a3b00")
ax.set_xlabel("magnitude stratum  min(|r$_A$|, |r$_B$|)", fontsize=12)
ax.set_ylabel("pairs reversing correlation sign (%, log scale)", fontsize=12)
ax.set_title("Sign reversal is rare but exceeds the sampling floor at every magnitude\n"
             "(observed >> within-condition null; all 30 sessions, Wilcoxon p < 0.001)", fontsize=11)
ax.legend(fontsize=10, loc="upper right"); ax.grid(alpha=0.25, which="both")
save(fig, "Fig1_signflip_stratified")

# ---- Fig 2: per-unit reorganization vs sign coherence (the -0.52 trade-off) ----
pu = pd.read_csv(os.path.join(B, "waveform_association", "perunit_waveform_continuous.csv"))
d = pu[["mean_row_reorg", "sign_coh_floor_0.1"]].dropna()
fig, ax = plt.subplots(figsize=(7, 5.5))
hb = ax.hexbin(d["mean_row_reorg"], d["sign_coh_floor_0.1"], gridsize=45, cmap="viridis", mincnt=1, bins="log")
ax.set_xlabel("per-neuron row reorganization", fontsize=12)
ax.set_ylabel("per-neuron sign coherence (min|r| ≥ 0.1)", fontsize=12)
ax.set_title("Reorganization trades off against sign preservation\n(Spearman rho = -0.52, 30/30 sessions; n = %d neurons)" % len(d), fontsize=11)
cb = fig.colorbar(hb, ax=ax); cb.set_label("neurons (log)", fontsize=10)
save(fig, "Fig2_reorg_vs_coherence")

# ---- Fig 3: scalar baseline vs architecture reorganization ----
sc = pd.read_csv(os.path.join(B, "finalize", "C_scalar_vs_architecture.csv"))
fig, ax = plt.subplots(figsize=(7, 5.5))
parts = ax.violinplot([sc["scalar_change_mean_absr"].values, sc["architecture_reorg"].values],
                      showmedians=True, widths=0.8)
for pc, c in zip(parts["bodies"], [GREY, BLUE]): pc.set_facecolor(c); pc.set_alpha(0.6)
ax.set_xticks([1, 2]); ax.set_xticklabels(["standard scalar\n(Δ mean |r|)", "architecture\nreorganization (1 − WJ)"], fontsize=11)
ax.set_ylabel("change between conditions", fontsize=12)
ax.set_title("A scalar summary collapses the reorganization:\nmedian 0.009 vs 0.58 (~60×) across 2,052 comparisons", fontsize=11)
ax.grid(alpha=0.25, axis="y")
save(fig, "Fig3_scalar_vs_architecture")

# ---- Fig 4: per-neuron coupling regimes vs null ----
pr = pd.read_csv(os.path.join(B, "finalize", "D_perneuron_coupling_regimes.csv"))
nullj = json.load(open(os.path.join(B, "finalize", "D_null_participation_ratio.json")))
fig, ax = plt.subplots(figsize=(7, 5.5))
ax.hist(pr["participation_ratio"], bins=40, color=GREEN, alpha=0.8, label="observed", density=True)
ax.axvline(nullj["observed_median_PR"], color=GREEN, ls="-", lw=2, label=f"observed median {nullj['observed_median_PR']}")
ax.axvline(nullj["null_median_PR"], color=GREY, ls="--", lw=2, label=f"condition-shuffled null median {nullj['null_median_PR']}")
ax.set_xlabel("effective coupling regimes per neuron (participation ratio)", fontsize=12)
ax.set_ylabel("density", fontsize=12)
ax.set_title("Neurons reuse a small, structured coupling repertoire\n(~4 regimes vs null ~11; 100% of neurons below null)", fontsize=11)
ax.legend(fontsize=10)
save(fig, "Fig4_coupling_regimes")

print("Figures written to", OUT)
print(sorted(os.listdir(OUT)))
