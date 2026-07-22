"""
Extract all h5ad files to R-friendly format (robust to different file structures)
"""
import sys, os, gc, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import scanpy as sc
import scipy.sparse as sp
import scipy.io
import numpy as np

def safe_str(s):
    """Convert any string to ASCII-safe representation"""
    if isinstance(s, bytes):
        s = s.decode('utf-8', errors='replace')
    return str(s).encode('ascii', errors='replace').decode('ascii')

def write_lines(fp, lines):
    with open(fp, 'w', encoding='utf-8', errors='replace') as f:
        for line in lines:
            f.write(safe_str(line) + '\n')

DATA_DIR = "E:/GBM/spatial_data_visium/spatial_data_visium/anndata"
OUT_BASE = "E:/GBM/results/spalinker_input"
os.makedirs(OUT_BASE, exist_ok=True)

files = sorted([f for f in os.listdir(DATA_DIR) if f.endswith('.h5ad')])
print(f"Total files: {len(files)}")

# ---- Find first file with Cell state abundances ----
cell_names = []
gene_names_sample = []
ref_file = None
for fname in files:
    fp = os.path.join(DATA_DIR, fname)
    try:
        a = sc.read_h5ad(fp)
        if 'feature_types' in a.var.columns:
            cm = a.var['feature_types'] == 'Cell state abundances'
            if cm.sum() > 0:
                cell_names = list(a.var_names[cm])
                gene_names_sample = list(a.var_names[a.var['feature_types'] == 'Gene Expression'])
                ref_file = fname
                print(f"Reference file: {fname} ({len(cell_names)} cell states, {len(gene_names_sample)} genes)")
                break
    except:
        continue

if ref_file is None:
    print("ERROR: No file with Cell state abundances found!")
    sys.exit(1)

# Print cell states
print("\nCell state names:")
for i, n in enumerate(cell_names):
    print(f"  [{i}] {n}")

# Identify B/Plasma and T cell states
b_states = [n for n in cell_names if any(k in n.lower() for k in
    ['b cell','b_cell','plasma','bcell','b lymph','plasmablast','germinal center'])]
t_states = [n for n in cell_names if any(k in n.lower() for k in
    ['t cell','t_cell','cd4','cd8','tcell','t lymph','cd4+','cd8+'])]

print(f"\nB/Plasma: {len(b_states)} states")
for s in b_states: print(f"  {s}")
print(f"T cell: {len(t_states)} states")
for s in t_states: print(f"  {s}")

# Save metadata
with open(os.path.join(OUT_BASE, "cell_state_info.txt"), 'w', encoding='utf-8') as f:
    f.write("ALL_CELL_STATES:\n")
    for n in cell_names: f.write(f"  {n}\n")
    f.write("\nB_PLASMA_STATES:\n")
    for n in b_states: f.write(f"  {n}\n")
    f.write("\nT_CELL_STATES:\n")
    for n in t_states: f.write(f"  {n}\n")

# ---- Process all files ----
skipped = 0
done = 0
for fi, fname in enumerate(files):
    tag = fname.replace('.h5ad', '')
    out_dir = os.path.join(OUT_BASE, tag)
    os.makedirs(out_dir, exist_ok=True)

    if os.path.exists(os.path.join(out_dir, "done.txt")):
        done += 1
        continue

    print(f"[{fi+1}/{len(files)}] {tag} ...", end=' ', flush=True)

    try:
        adata = sc.read_h5ad(os.path.join(DATA_DIR, fname))

        # Check for feature_types or alternative column name
        ft_col = None
        for col_name in ['feature_types', 'feature_type', 'type']:
            if col_name in adata.var.columns:
                ft_col = col_name
                break

        if ft_col is None:
            print("SKIP (no feature_types column)")
            skipped += 1
            continue

        ft = adata.var[ft_col]

        # Gene expression
        gene_mask = ft == 'Gene Expression'
        if gene_mask.sum() == 0:
            # Try alternative: everything that's not cell abundance / other
            gene_mask = ~ft.isin(['Cell state abundances', 'Spatial niche abundances',
                                   'Histopath annotation overlap'])
        if gene_mask.sum() == 0:
            print("SKIP (no gene expression)")
            skipped += 1
            continue

        gene_mat = adata[:, gene_mask].X
        if sp.issparse(gene_mat):
            sp.io.mmwrite(os.path.join(out_dir, "counts.mtx"), gene_mat)
        gene_names = [safe_str(n) for n in adata.var_names[gene_mask]]
        write_lines(os.path.join(out_dir, "genes.tsv"), gene_names)

        # Cell state abundances
        cell_mask = ft == 'Cell state abundances'
        cell_mat = None
        if cell_mask.sum() > 0:
            cell_mat = adata[:, cell_mask].X
            if sp.issparse(cell_mat):
                cell_dense = cell_mat.toarray()
            else:
                cell_dense = np.array(cell_mat)
            np.savetxt(os.path.join(out_dir, "cell_props.tsv"), cell_dense, delimiter='\t')

        # Barcodes
        write_lines(os.path.join(out_dir, "barcodes.tsv"), [safe_str(b) for b in adata.obs_names])

        # Spatial coords
        np.savetxt(os.path.join(out_dir, "spatial.tsv"),
                   adata.obsm['spatial'], delimiter='\t', header='x\ty', comments='')

        # Mark done
        with open(os.path.join(out_dir, "done.txt"), 'w') as f:
            f.write('ok')
        n_cell = cell_mask.sum() if cell_mask.sum() > 0 else 0
        print(f"OK spots={adata.n_obs} genes={gene_mask.sum()} cells={n_cell}")
        gc.collect()

    except Exception as e:
        print(f"ERROR: {e}")
        skipped += 1

print(f"\nDone! Processed: {done+len(files)-skipped-done}, Skipped: {skipped}, Already done: {done}")
