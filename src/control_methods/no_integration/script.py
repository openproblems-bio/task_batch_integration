import anndata as ad

## VIASH START
par = {
    'input_dataset': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
    'output': 'output.h5ad',
}
## VIASH END

print('Read input_dataset', flush=True)
adata = ad.read_h5ad(par['input_dataset'])

print("Create output", flush=True)
output = ad.AnnData(
    layers={"corrected_counts": adata.layers["counts"]},
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={"X_emb": adata.obsm["X_pca"]},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"]
    }
)
# could also add knn

print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
