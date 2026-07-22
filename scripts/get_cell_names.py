import h5py, sys
fp = 'E:/GBM/spatial_data_visium/spatial_data_visium/anndata/AT3-BRA5-FO-1_1.h5ad'
with h5py.File(fp, 'r') as f:
    var_names = f['var/_index'][:]
    var_ft = f['var/feature_type'][:]
    cell_idx = [i for i, ft in enumerate(var_ft) if ft.decode() == 'Cell state abundances']
    print(f'Cell states: {len(cell_idx)}')
    for i in cell_idx:
        print(var_names[i].decode())
