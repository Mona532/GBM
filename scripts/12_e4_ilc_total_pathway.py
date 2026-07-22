from __future__ import annotations

from pathlib import Path
import re

import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(r"E:/GBM")
RESULTS = ROOT / "results"
INPUT = RESULTS / "e4_spot_gene_corr_ILC_total.csv"
LIB_DIR = RESULTS / "pathway_library_cache"

BAD_LABEL = re.compile(r"^(?:RPL|RPS|MT-|MTRNR|LINC|MALAT1$|NEAT1$)", re.I)
NEURO = re.compile(r"(?:neuro|neuron|synap|axon|dendrit|glutam|gaba|seroton|dopamin|cholin|nervous|transmission|vesicle)", re.I)

LIBRARIES = [
    ("Reactome_2022", "reactome"),
    ("GO_Biological_Process_2023", "gobp"),
]


def ensure_gmt(library_name: str) -> Path:
    LIB_DIR.mkdir(parents=True, exist_ok=True)
    gmt_path = LIB_DIR / f"{library_name}.gmt"
    if not gmt_path.exists():
        gp.get_library(name=library_name, organism="Human", save=str(gmt_path))
    return gmt_path


def load_ranking() -> pd.DataFrame:
    df = pd.read_csv(INPUT)
    df = df.loc[~df["gene"].astype(str).str.contains(BAD_LABEL, na=False)].copy()
    df = df[["gene", "rho", "pvalue", "fdr"]].dropna()
    df = df.sort_values(["rho", "pvalue"], ascending=[False, True])
    df = df.drop_duplicates("gene", keep="first")
    return df


def run_prerank(rnk: pd.DataFrame, library_name: str) -> pd.DataFrame:
    gmt_path = ensure_gmt(library_name)
    pre = gp.prerank(
        rnk=rnk[["gene", "rho"]],
        gene_sets=str(gmt_path),
        min_size=15,
        max_size=500,
        permutation_num=1000,
        outdir=None,
        no_plot=True,
        seed=1,
        threads=4,
        verbose=False,
    )
    res = pre.res2d.copy()
    if res.empty:
        return res

    res["Gene_set"] = library_name
    res["Term"] = res["Term"].astype(str)
    res["FDR q-val"] = pd.to_numeric(res["FDR q-val"], errors="coerce")
    res["NOM p-val"] = pd.to_numeric(res["NOM p-val"], errors="coerce")
    res["NES"] = pd.to_numeric(res["NES"], errors="coerce")
    res["ES"] = pd.to_numeric(res["ES"], errors="coerce")
    res["direction"] = res["NES"].map(lambda x: "positive" if x > 0 else "negative")
    res["neuro_related"] = res["Term"].str.contains(NEURO, na=False)
    return res


def save_tables(res: pd.DataFrame, short_name: str) -> None:
    res.to_csv(RESULTS / f"e4_spot_gene_corr_ILC_total_prerank_{short_name}.csv", index=False)
    res[res["neuro_related"]].to_csv(
        RESULTS / f"e4_spot_gene_corr_ILC_total_prerank_{short_name}_neuro_hits.csv",
        index=False,
    )


def plot_summary(combined: pd.DataFrame) -> None:
    plot_df = combined[combined["FDR q-val"] <= 0.25].copy()
    if plot_df.empty:
        plot_df = combined.nsmallest(12, "FDR q-val").copy()

    keep = []
    for direction in ["positive", "negative"]:
        for gene_set in plot_df["Gene_set"].drop_duplicates():
            sub = plot_df[(plot_df["direction"] == direction) & (plot_df["Gene_set"] == gene_set)]
            sub = sub.sort_values("FDR q-val").head(12)
            keep.append(sub)
    plot_df = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    if plot_df.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    gene_sets = list(plot_df["Gene_set"].drop_duplicates())
    for i, direction in enumerate(["positive", "negative"]):
        for j, gene_set in enumerate(gene_sets):
            ax = axes[i, j]
            sub = plot_df[(plot_df["direction"] == direction) & (plot_df["Gene_set"] == gene_set)].copy()
            if sub.empty:
                ax.axis("off")
                continue
            sub = sub.sort_values("NES")
            colors = ["#B2182B" if x else "#4D4D4D" for x in sub["neuro_related"]]
            sizes = sub["NES"].abs() * 140
            ax.scatter(sub["NES"], range(len(sub)), s=sizes, c=colors, alpha=0.85)
            ax.axvline(0, linestyle="--", linewidth=0.8, color="#777777")
            ax.set_yticks(range(len(sub)))
            ax.set_yticklabels(sub["Term"], fontsize=8)
            ax.set_xlabel("NES")
            ax.set_title(f"{direction} | {gene_set}", fontsize=11, weight="bold")
            ax.grid(axis="x", linestyle="--", alpha=0.25)

    fig.suptitle("E4 ILC_total prerank pathway enrichment", fontsize=14, weight="bold")
    fig.savefig(RESULTS / "fig_e4_spot_gene_corr_ILC_total_prerank_pathway.jpg", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "fig_e4_spot_gene_corr_ILC_total_prerank_pathway.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rnk = load_ranking()
    all_res = []
    for library_name, short_name in LIBRARIES:
        res = run_prerank(rnk, library_name)
        if res.empty:
            continue
        save_tables(res, short_name)
        all_res.append(res)

    if not all_res:
        raise SystemExit("No prerank results.")

    combined = pd.concat(all_res, ignore_index=True)
    combined = combined.sort_values(["FDR q-val", "NOM p-val"], ascending=[True, True])
    combined.to_csv(RESULTS / "e4_spot_gene_corr_ILC_total_prerank_combined.csv", index=False)
    combined[combined["neuro_related"]].to_csv(
        RESULTS / "e4_spot_gene_corr_ILC_total_prerank_combined_neuro_hits.csv",
        index=False,
    )
    plot_summary(combined)

    print(f"ranked_genes={len(rnk)}")
    print(f"terms={len(combined)}")
    print(
        combined.loc[
            combined["neuro_related"],
            ["Gene_set", "direction", "Term", "NES", "FDR q-val"],
        ]
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
