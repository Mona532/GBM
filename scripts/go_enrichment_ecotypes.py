"""
GO enrichment on NMF ecotype marker genes (top 200 DEGs per ecotype vs rest).
"""
import pandas as pd, numpy as np
from pathlib import Path
import gseapy as gp
import warnings; warnings.filterwarnings("ignore")

ROOT = Path(r"E:/GBM/results")
OUT = ROOT / "go_enrichment"; OUT.mkdir(parents=True, exist_ok=True)

ECO_NAMES = {
    "E1": "Lymphocyte TLS",
    "E2": "ILC-enriched TLS",
    "E3": "Myeloid-vascular TLS",
    "E4": "Glial-CD4 TLS",
}

all_results = []

for eco in ["E1", "E2", "E3", "E4"]:
    # Load full DEG results (all genes ranked, not just top 30)
    deg = pd.read_csv(ROOT / f"tls_compnmf_rank4_{eco}_markers_vs_rest.csv")

    # Take top 200 up-regulated genes (sorted by logFC or already ranked)
    if "logFC" in deg.columns:
        up_genes = deg[deg["logFC"] > 0].sort_values("logFC", ascending=False).head(200)
    elif "avg_log2FC" in deg.columns:
        up_genes = deg[deg["avg_log2FC"] > 0].sort_values("avg_log2FC", ascending=False).head(200)
    else:
        # Just take first 200 (already ranked)
        up_genes = deg.head(200)

    gene_list = up_genes["gene"].dropna().tolist()
    print(f"\n{eco} ({ECO_NAMES[eco]}): {len(gene_list)} genes for GO")

    if len(gene_list) < 10:
        print(f"  Skipping — too few genes")
        continue

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=["GO_Biological_Process_2023", "GO_Molecular_Function_2023",
                        "KEGG_2021_Human", "Reactome_2022"],
            organism="human",
            outdir=str(OUT / eco),
            no_plot=True,
            cutoff=0.05,
        )

        if enr.results is not None and len(enr.results) > 0:
            res = enr.results
            res["ecotype"] = eco
            res["ecotype_name"] = ECO_NAMES[eco]
            all_results.append(res)

            # Print top 5 GO:BP terms
            bp = res[res["Gene_set"] == "GO_Biological_Process_2023"].head(5)
            print(f"  Top GO:BP:")
            for _, row in bp.iterrows():
                print(f"    {row['Term']} (FDR={row['Adjusted P-value']:.1e}, genes={row['Overlap']})")
    except Exception as e:
        print(f"  GO failed: {e}")

if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(ROOT / "tls_nmf_ecotype_go_enrichment.csv", index=False)
    print(f"\nSaved: {ROOT}/tls_nmf_ecotype_go_enrichment.csv")
    print(f"Total enriched terms: {len(combined)}")
else:
    print("\nNo enrichment results.")
