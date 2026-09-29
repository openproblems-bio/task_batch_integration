import sys

import anndata as ad
import numpy as np
import openproblems as op
import scanpy as sc
import scib

## VIASH START
par = {
    "input": "resources_test/common/pancreas/dataset.h5ad",
    "hvgs": 2000,
    "obs_label": "cell_type",
    "obs_batch": "batch",
    "output": "output.h5ad",
}
meta = {
    "config": "target/nextflow/data_processors/process_dataset/.config.vsh.yaml",
    "resources_dir": "target/nextflow/data_processors/process_dataset",
}
## VIASH END

# import helper functions
sys.path.append(meta["resources_dir"])
from subset_h5ad_by_format import subset_h5ad_by_format


# ----------- FUNCTIONS -----------
def compute_batched_hvg(adata, n_hvgs):
    if n_hvgs > adata.n_vars or n_hvgs <= 0:
        hvg_list = adata.var_names.tolist()
    else:
        scib_adata = adata.copy()
        del scib_adata.layers["counts"]
        scib_adata.X = scib_adata.layers["normalized"].copy()
        hvg_list = scib.pp.hvg_batch(
            scib_adata, batch_key="batch", target_genes=n_hvgs, adataOut=False
        )
    return hvg_list


# ----------- SCRIPT -----------
print(
    f"====== Process dataset (scanpy v{sc.__version__}, scib v{scib.__version__}) ======",
    flush=True,
)

print("\n>>> Loading config...", flush=True)
config = op.project.read_viash_config(meta["config"])

print("\n>>> Reading input files...", flush=True)
print(f"Input H5AD file: '{par['input']}'", flush=True)
adata = ad.read_h5ad(par["input"])
print(adata, flush=True)

n_hvgs = par["hvgs"]
if adata.n_vars < n_hvgs:
    print(f"\n>>> Using all {adata.n_vars} features as batch-aware HVGs...", flush=True)
    n_hvgs = adata.n_vars
    hvg_list = adata.var_names.tolist()
else:
    print(f"\n>>> Computing {n_hvgs} batch-aware HVGs...", flush=True)
    hvg_list = compute_batched_hvg(adata, n_hvgs=n_hvgs)

adata.var["batch_hvg"] = adata.var_names.isin(hvg_list)

n_components = adata.obsm["X_pca"].shape[1]
print(f"\n>>> Computing PCA with {n_components} components using HVGs...", flush=True)
X_pca, loadings, variance, variance_ratio = sc.pp.pca(
    adata.layers["normalized"],
    n_comps=n_components,
    mask_var=adata.var["batch_hvg"],
    return_info=True,
)
adata.obsm["X_pca"] = X_pca
adata.varm["pca_loadings"] = np.zeros(shape=(adata.n_vars, n_components))
adata.varm["pca_loadings"][adata.var["batch_hvg"], :] = loadings.T
adata.uns["pca_variance"] = {"variance": variance, "variance_ratio": variance_ratio}

print("\n>>> Computing neighbours using PCA...", flush=True)
del adata.uns["knn"]
del adata.obsp["knn_connectivities"]
del adata.obsp["knn_distances"]
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, key_added="knn")

print("\n>>> Recomputing standard HVG scores...", flush=True)
out = sc.pp.highly_variable_genes(
    adata, layer="normalized", n_top_genes=n_hvgs, flavor="cell_ranger", inplace=False
)
adata.var["hvg"] = out["highly_variable"].values
adata.var["hvg_score"] = out["dispersions_norm"].values

print("\n>>> Creating dataset output object...", flush=True)
output_dataset = subset_h5ad_by_format(
    adata[:, adata.var_names.isin(hvg_list)].copy(), config, "output_dataset"
)
print(output_dataset, flush=True)

print("\n>>> Creating solution output object...", flush=True)
output_solution = subset_h5ad_by_format(adata, config, "output_solution")
print(output_solution, flush=True)

print("\n>>> Writing output to files...", flush=True)
print(f"Dataset H5AD file: '{par['output_dataset']}'", flush=True)
output_dataset.write_h5ad(par["output_dataset"], compression="gzip")
print(f"Solution H5AD file: '{par['output_solution']}'", flush=True)
output_solution.write_h5ad(par["output_solution"], compression="gzip")

print("\n>>> Done!", flush=True)
