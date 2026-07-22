from __future__ import annotations

from pathlib import Path
import math
import traceback
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(r"E:/GBM")
RESULTS = ROOT / "results"
INPUT = RESULTS / "e4_spot_gene_corr_ILC_total_prerank_combined_neuro_hits.csv"
OUT_CSV = RESULTS / "e4_spot_gene_corr_ILC_total_positive_neuro_top.csv"
OUT_JPG = RESULTS / "fig_e4_ilc_total_positive_neuro_pathways.jpg"
OUT_PDF = RESULTS / "fig_e4_ilc_total_positive_neuro_pathways.pdf"
OUT_ERR = RESULTS / "e4_ilc_total_positive_neuro_plot_error.txt"


def main() -> None:
    df = pd.read_csv(INPUT)
    df = df[(df["direction"] == "positive") & (df["neuro_related"] == True)].copy()
    if df.empty:
        raise SystemExit("No positive neuro-related pathways found.")

    df["FDR q-val"] = pd.to_numeric(df["FDR q-val"], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.sort_values(["FDR q-val", "NES"], ascending=[True, False])

    top_reactome = df[df["Gene_set"] == "Reactome_2022"].head(10)
    top_gobp = df[df["Gene_set"] == "GO_Biological_Process_2023"].head(10)
    top = pd.concat([top_reactome, top_gobp], ignore_index=True)
    top.to_csv(OUT_CSV, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    for ax, gene_set, title in zip(
        axes,
        ["Reactome_2022", "GO_Biological_Process_2023"],
        ["Reactome", "GO Biological Process"],
    ):
        sub = top[top["Gene_set"] == gene_set].copy()
        if sub.empty:
            ax.axis("off")
            continue
        sub = sub.sort_values("NES", ascending=True)
        sizes = sub["NES"].abs() * 160
        colors = -sub["FDR q-val"].clip(lower=1e-300).map(math.log10)
        sc = ax.scatter(sub["NES"], range(len(sub)), s=sizes, c=colors, cmap="Reds", alpha=0.9)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["Term"], fontsize=8)
        ax.set_xlabel("NES")
        ax.set_title(title, fontsize=12, weight="bold")
        ax.axvline(0, linestyle="--", linewidth=0.8, color="#777777")
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("-log10(FDR)")

    fig.suptitle("E4 ILC_total positive neuro-related pathways", fontsize=15, weight="bold")
    fig.savefig(OUT_JPG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)

    print(top[["Gene_set", "Term", "NES", "FDR q-val"]].to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        OUT_ERR.write_text(traceback.format_exc(), encoding="utf-8")
        raise
