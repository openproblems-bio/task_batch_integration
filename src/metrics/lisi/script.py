import sys
import numpy as np
import anndata as ad
from scib.metrics.lisi import lisi_graph_py

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


print('Read input', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsp='obsp', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

print('compute iLISI score...', flush=True)
ilisi_scores = lisi_graph_py(
    adata=adata,
    obs_key='batch',
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
)
ilisi = np.nanmedian(ilisi_scores)
ilisi = (ilisi - 1) / (adata.obs['batch'].nunique() - 1)

print('compute cLISI scores...', flush=True)
clisi_scores = lisi_graph_py(
    adata=adata,
    obs_key="cell_type",
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
)
clisi = np.nanmedian(clisi_scores)
nlabs = adata.obs["cell_type"].nunique()
clisi = (nlabs - clisi) / (nlabs - 1)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ 'ilisi', 'clisi' ],
        'metric_values': [ ilisi, clisi ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
