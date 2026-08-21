import sys
import anndata as ad
import pegasus as pg
import pegasusio
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_processed.h5ad',
    'input_solution': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/solution.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'name': 'ksim',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

n_threads = meta["cpus"] or -1

print('Read input...', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns
print(adata, flush=True)

print('Convert to pegasusio.MultimodalData...', flush=True)
adata.X = csr_matrix(adata.shape)
mmdata = pegasusio.MultimodalData(adata)

print('Compute kSIM...', flush=True)
# K and min_rate are pegasus' own defaults, spelled out here so the metric stays
# stable across pegasus releases. Note that pegasus caps K at sqrt(n_obs), so the
# effective neighbourhood is smaller for datasets with fewer than 625 cells.
ksim_mean, ksim_accept_rate = pg.calc_kSIM(
    mmdata,
    attr='cell_type',
    rep='emb',
    K=25,
    min_rate=0.9,
    n_jobs=n_threads,
)
print('kSIM mean:', ksim_mean, flush=True)
print('kSIM acceptance rate:', ksim_accept_rate, flush=True)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ 'ksim_mean', 'ksim_accept_rate' ],
        'metric_values': [ ksim_mean, ksim_accept_rate ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
