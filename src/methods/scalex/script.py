import sys
import anndata as ad
import scalex

## VIASH START
par = {
    'input': 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'name' : 'foo',
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

print('Run SCALEX', flush=True)
adata = scalex.SCALEX(
    adata,
    batch_key="batch",
    ignore_umap=True,
    impute=adata.obs["batch"].cat.categories[0],
    processed=True,
    max_iteration=40,
    min_features=None,
    min_cells=None,
    n_top_features=0,
    outdir=None,
    gpu=0,
)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={
        'corrected_counts': adata.layers["impute"],
    },
    obsm={
        'X_emb': adata.obsm['latent'],
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
