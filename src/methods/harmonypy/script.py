import sys
import anndata as ad
import numpy as np
import harmonypy as hm

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad",
    "output": "output.h5ad"
}
meta = {
    "name": "harmonypy",
    "resources_dir": "src/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print(">> Read input", flush=True)
adata = read_anndata(
    par["input"],
    obs="obs",
    obsm="obsm",
    var="var",
    uns="uns"
)

print(">> Run harmonypy", flush=True)
out = hm.run_harmony(
  adata.obsm["X_pca"],
  adata.obs,
  "batch"
)

print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": out.Z_corr.transpose()
    },
    shape=adata.shape,
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"],
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
