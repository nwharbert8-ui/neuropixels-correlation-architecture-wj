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

# ---- Fig 1: direct sign-flip rate vs magnitude stratum ----
prov = json.load(open(os.path.join(B, "provenance.json")))
st = prov["sign_flip_rate_by_stratum"]
thr = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
rate = [st[f"thr_{t}"]["mean_sign_flip_rate"] * 100 for t in thr]
fig, ax = plt.subplots(figsize=(7, 5.5))
ax.axhline(50, color=GREY, ls="--", lw=1.5, label="chance (50%)")
ax.plot(thr, rate, "o-", color=ORANGE, lw=2.5, ms=9, mec="black", mew=0.6, label="observed sign-flip rate")
for t, r in zip(thr, rate):
    ax.annotate(f"{r:.1f}%", (t, r), textcoords="offset points", xytext=(6, 8), fontsize=9)
ax.set_xlabel("magnitude stratum  min(|r$_A$|, |r$_B$|)", fontsize=12)
ax.set_ylabel("pairs reversing correlation sign (%)", fontsize=12)
ax.set_title("Sign is preserved: direct sign-flip rate is below chance and\nfalls to ~1% at meaningful correlation magnitude (2,052 comparisons)", fontsize=11)
ax.set_ylim(-3, 58); ax.legend(fontsize=11); ax.grid(alpha=0.25)
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
