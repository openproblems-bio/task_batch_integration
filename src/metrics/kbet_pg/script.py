from tqdm import tqdm
import sys
import numpy as np
import anndata as ad
import pegasus as pg
import pegasusio
from scipy.sparse import csr_matrix


def compute_kbet(mmdata, *args, **kwargs):
    stat_mean, pvalue_mean, accept_rate = pg.calc_kBET(mmdata, *args, **kwargs)
    return accept_rate


## VIASH START
par = {
    'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'name': 'kbet_pg',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

n_threads = 4

print('Read input...', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns
print(adata, flush=True)

print('Convert to pegasusio.MultimodalData...', flush=True)
adata.X = csr_matrix(adata.shape)
mmdata = pegasusio.MultimodalData(adata)

print('Compute global kBET...', flush=True)
score = compute_kbet(
    mmdata,
    attr="batch",
    rep="emb",
    K=50,
    n_jobs=n_threads,
)
print('Global kBET score:', score, flush=True)

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
