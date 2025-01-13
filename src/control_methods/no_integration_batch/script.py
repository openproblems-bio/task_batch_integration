import sys
import scanpy as sc
import numpy as np

## VIASH START

par = {
    'input_dataset': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality': 'foo',
    'config': 'bar'
}

## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input_dataset', flush=True)
adata = read_anndata(
    par['input_dataset'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)
adata.var["highly_variable"] = adata.var["hvg"]
print(adata, flush=True)

print("Process dataset", flush=True)
adata.obsm["X_emb"] = np.zeros((adata.shape[0], 50), dtype=float)
for batch in adata.obs["batch"].unique():
    batch_idx = adata.obs["batch"] == batch

    if np.sum(batch_idx) <= 50:
        n_comps = np.sum(batch_idx) - 1
        solver = "full"
        print(f"Batch '{batch}' has 50 or less cells. Using the 'full' solver with {n_comps} components.", flush=True)
    else:
        n_comps = 50
        solver = "arpack"

    adata.obsm["X_emb"][batch_idx, :n_comps] = sc.tl.pca(
        adata[batch_idx].copy(),
        n_comps=n_comps,
        use_highly_variable=True,
        svd_solver=solver,
        copy=True,
    ).obsm["X_pca"]

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['name']
adata.write_h5ad(par['output'], compression='gzip')
