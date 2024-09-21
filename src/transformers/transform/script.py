import anndata as ad
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/integrated.h5ad',
    'ouput': 'output.h5ad'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

assert 'corrected_counts' in adata.layers.keys() or \
    'X_emb' in adata.obsm.keys() or \
    'neighbors' in adata.uns.keys(), \
    "Input data is missing necessary information. " \
    "Expected 'corrected_counts' in layers, " \
    "'X_emb' in obsm, or 'neighbors' in uns.\n" \
    f"Found: {adata}"

if 'X_emb' not in adata.obsm:
    print('Run PCA...', flush=True)
    adata.obsm['X_emb'] = sc.pp.pca(
        adata.layers['corrected_counts'],
        n_comps=50,
        use_highly_variable=False,
        svd_solver='arpack',
        return_info=False
    )

if 'neighbors' not in adata.uns:
    print('Run kNN...', flush=True)
    sc.pp.neighbors(adata, use_rep='X_emb')

print("Store outputs", flush=True)
adata.write_h5ad(par['output'], compression='gzip')