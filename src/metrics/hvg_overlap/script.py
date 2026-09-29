import sys

import anndata as ad
import numpy as np
import scanpy as sc
from scib.metrics import hvg_overlap
from scib.utils import split_batches

## VIASH START
par = {
    "input_integrated": "resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad",
    "input_solution": "resources_test/task_batch_integration/cxg_immune_cell_atlas/solution.h5ad",
    "output": "output.h5ad",
}
meta = {"name": "foo", "resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print("Read input", flush=True)
adata_solution = read_anndata(
    par["input_solution"], X="layers/normalized", obs="obs", var="var", uns="uns"
)
adata_integrated = read_anndata(
    par["input_integrated"],
    X="layers/corrected_counts",
    obs="obs",
    var="var",
    uns="uns",
)

print("Copy batch information", flush=True)
adata_integrated.obs["batch"] = adata_solution.obs["batch"]

print("Remove batches with insufficient genes", flush=True)
adata_list = split_batches(adata_solution, "batch", hvg=adata_integrated.var_names)
skip_batches = []
for adata_batch in adata_list:
    sc.pp.filter_genes(adata_batch, min_cells=1)
    n_hvg_tmp = np.minimum(500, int(0.5 * adata_batch.n_vars))
    if n_hvg_tmp < 500:
        batch = adata_batch.obs["batch"][0]
        print(
            f"Batch '{batch}' has insufficient genes (0.5 * {adata_batch.n_vars} < 500) and will be skipped",
            flush=True,
        )
        skip_batches.append(batch)

adata_solution = adata_solution[~adata_solution.obs["batch"].isin(skip_batches)]
adata_integrated = adata_integrated[~adata_integrated.obs["batch"].isin(skip_batches)]

print("Compute score", flush=True)
score = hvg_overlap(
    adata_solution[:, adata_solution.var_names.isin(adata_integrated.var_names)],
    adata_integrated,
    batch_key="batch"
)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": adata_solution.uns["dataset_id"],
        "normalization_id": adata_solution.uns["normalization_id"],
        "method_id": adata_integrated.uns["method_id"],
        "metric_ids": [meta["name"]],
        "metric_values": [score],
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
