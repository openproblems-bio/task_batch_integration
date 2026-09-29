import sys

import anndata as ad
import scanorama

## VIASH START
par = {
    'input': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
    'output': 'output.h5ad',
    'dimred': 100
}
meta = {
    'name': 'scanorama',
    'resources_dir': 'src/utils'
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
split = [adata[adata.obs['batch'] == batch].copy() for batch in adata.obs['batch'].cat.categories]
corrected = scanorama.correct_scanpy(split, return_dimred=True, dimred=par['dimred'])

# scanorama returns one object per batch, with the genes sorted by name -- put the
# cells and genes back in the order of the input before storing the output
corrected = ad.concat(corrected)
corrected = corrected[adata.obs_names, adata.var_names]

print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={
        'corrected_counts': corrected.X,
    },
    obsm={
        'X_emb': corrected.obsm['X_scanorama'],
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
