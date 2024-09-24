import anndata as ad
import mnnpy

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

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])
adata.X = adata.layers['normalized']
del adata.layers['normalized']
del adata.layers['counts']

print('Run mnn', flush=True)
split = []
batch_categories = adata.obs['batch'].cat.categories
for i in batch_categories:
    split.append(adata[adata.obs['batch'] == i].copy())
corrected, _, _ = mnnpy.mnn_correct(
        *split,
        batch_key='batch',
        batch_categories=batch_categories,
        index_unique=None
    )

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
    },
    layers={
        'corrected_counts': corrected.X,
    }
)


print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
