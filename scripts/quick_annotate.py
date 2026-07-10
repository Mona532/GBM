import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

adata = sc.read("E:/GBM/GBM_DATA/5sample_final/GBM_5sample_final.h5ad")
markers_df = pd.read_csv("E:/GBM/GBM_DATA/5sample_final/markers_df.csv")
markers_df["cluster"] = markers_df["cluster"].astype(str)

clusters = adata.obs["leiden_r0.6"].cat.categories
for cl in sorted(clusters, key=lambda x: int(x)):
    top10 = markers_df[markers_df["cluster"] == cl]["names"].head(10).tolist()
    n = adata.obs["leiden_r0.6"].value_counts().get(cl, 0)
    print(f"Cluster {cl} ({n:,} cells): {', '.join(top10)}")

anno_map = {
    "0": "T/NK cytotoxic", "1": "TAM SPP1+", "2": "T/NK helper",
    "3": "Microglia homeo", "4": "TAM MHC-II+", "5": "Neutrophil",
    "6": "TAM APOE+", "7": "Tumor GSC/OPC",