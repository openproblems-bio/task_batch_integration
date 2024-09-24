import sys
import scanpy as sc
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input': 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'name': 'foo',
    'config': 'bar'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)

print('Run Combat', flush=True)
adata.X = sc.pp.combat(adata, key='batch', inplace=False)

print("Store output", flush=True)
output = sc.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
    },
    layers={
        'corrected_counts': csr_matrix(adata.X),
    }
)

print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
