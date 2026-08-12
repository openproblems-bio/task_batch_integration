import sys
import anndata as ad
import scanorama
import numpy as np

## VIASH START
par = {
    'input': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
    'output': 'output.h5ad',
    'dimred': 100
}
meta = {
    'name': 'scanorama-integrate',
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

print('Run scanorama', flush=True)
split = []
batch_categories = adata.obs['batch'].cat.categories
for b in batch_categories:
    split.append(adata[adata.obs['batch'] == b].copy())
scanorama.integrate_scanpy(split, dimred=par["dimred"])

#From https://colab.research.google.com/drive/1CebA3Ow4jXITK0dW5el320KVTX_szhxG
result = np.zeros((adata.shape[0], split[0].obsm["X_scanorama"].shape[1]))
for i, b in enumerate(batch_categories):
    result[adata.obs['batch'] == b] = split[i].obsm["X_scanorama"]


print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
    },
    obsm={
        'X_emb': result
    },
    shape=adata.shape,
)

print("Write output to file", flush=True)
output.write(par['output'], compression='gzip')
