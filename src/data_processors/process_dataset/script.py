import sys
import anndata as ad
import scanpy as sc
import openproblems as op
import numpy as np

## VIASH START
par = {
    "input": "resources_test/common/pancreas/dataset.h5ad",
    "hvgs": 2000,
    "obs_label": "cell_type",
    "obs_batch": "batch",
    "output": "output.h5ad"
}
meta = {
    "config": "target/nextflow/data_processors/process_dataset/.config.vsh.yaml",
    "resources_dir": "target/nextflow/data_processors/process_dataset"
}
## VIASH END

# import helper functions
sys.path.append(meta["resources_dir"])
from subset_h5ad_by_format import subset_h5ad_by_format

print(">> Load config", flush=True)
config = op.project.read_viash_config(meta["config"])

print("Read input", flush=True)
adata = ad.read_h5ad(par["input"])

def compute_batched_hvg(adata, n_hvgs):
    if n_hvgs > adata.n_vars or n_hvgs <= 0:
        hvg_list = adata.var_names.tolist()
    else:
        import scib
        scib_adata = adata.copy()
        del scib_adata.layers["counts"]
        scib_adata.X = scib_adata.layers["normalized"].copy()
        hvg_list = scib.pp.hvg_batch(
            scib_adata,
            batch_key="batch",
            target_genes=n_hvgs,
            adataOut=False
        )
    return hvg_list

n_hvgs = par["hvgs"]
n_components = adata.obsm["X_pca"].shape[1]

if adata.n_vars < n_hvgs:
    n_hvgs = adata.n_vars
    hvg_list = adata.var_names.tolist()
else:
    print(f"Select {par['hvgs']} highly variable genes", flush=True)
    hvg_list = compute_batched_hvg(adata, n_hvgs=n_hvgs)

print(">> Recompute HVG", flush=True)
out = sc.pp.highly_variable_genes(
  adata,
  layer="normalized",
  n_top_genes=n_hvgs,
  flavor="cell_ranger",
  inplace=False
)
adata.var["hvg"] = out["highly_variable"].values
adata.var["hvg_score"] = out["dispersions_norm"].values

print(">> Recompute PCA", flush=True)
hvg_mask = adata.var_names.isin(hvg_list)
X_pca, loadings, variance, variance_ratio = sc.pp.pca(
    adata.layers["normalized"],
    n_comps=n_components,
    mask_var=hvg_mask,
    return_info=True
)
adata.obsm["X_pca"] = X_pca
adata.varm["pca_loadings"] = np.zeros(shape=(adata.n_vars, n_components))
adata.varm["pca_loadings"][hvg_mask, :] = loadings.T
adata.uns["pca_variance"] = {
    "variance": variance,
    "variance_ratio": variance_ratio
}

print(">> Recompute neighbors", flush=True)
del adata.uns["knn"]
del adata.obsp["knn_connectivities"]
del adata.obsp["knn_distances"]
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, key_added="knn")

print(">> Create output object", flush=True)
output_dataset = subset_h5ad_by_format(
    adata[:, adata.var_names.isin(hvg_list)].copy(),
    config,
    "output_dataset"
)
output_solution = subset_h5ad_by_format(
    adata,
    config,
    "output_solution"
)

print("Writing adatas to file", flush=True)
output_dataset.write_h5ad(par["output_dataset"], compression="gzip")
output_solution.write_h5ad(par["output_solution"], compression="gzip")
