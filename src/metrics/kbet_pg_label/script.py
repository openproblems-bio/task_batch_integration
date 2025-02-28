from tqdm import tqdm
import sys
import numpy as np
import anndata as ad
import pegasus as pg
import pegasusio
from scipy.sparse import csr_matrix
from joblib import Parallel, delayed


def compute_kbet(mmdata, *args, **kwargs):
    stat_mean, pvalue_mean, accept_rate = pg.calc_kBET(mmdata, *args, **kwargs)
    return accept_rate


## VIASH START
par = {
    'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'name': 'foo',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

n_threads = -1  # use all available threads

print('Read input...', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns
adata.X = csr_matrix(adata.shape)
print(adata, flush=True)

print('Compute cell type averaged kBET...', flush=True)
cell_types = adata.obs['cell_type'].unique()
scores = []
mmdata_per_cell_type = []

# collect all adata subsets
for cell_type in cell_types:
    ad_sub = adata[adata.obs['cell_type'] == cell_type]
    if ad_sub.obs['batch'].nunique() <= 1:
        print(f'Skipping cell type {cell_type} because it\'s present in only one batch', flush=True)
        continue
    _mmdata = pegasusio.MultimodalData(ad_sub.copy())
    mmdata_per_cell_type.append(_mmdata)

# compute kBET scores
scores = Parallel(n_jobs=n_threads)(
    delayed(compute_kbet)(
        _mmdata,
        attr="batch",
        rep="emb",
        K=50,
        use_cache=False,
        n_jobs=1,
    ) for _mmdata in tqdm(
        mmdata_per_cell_type,
        desc=f'Compute per cell type with {n_threads} threads',
        miniters=1,
    )
)
score = np.nanmean(scores)
print('Cell type averaged kBET score:', score, flush=True)

print('Create output AnnData object', flush=True)
metric_name = meta['name']
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ metric_name ],
        'metric_values': [ score ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
