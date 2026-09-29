import anndata as ad
import pandas as pd

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
    "clusterings": ["output.h5ad", "output2.h5ad"],
    "output": "output3.h5ad",
}
## VIASH END

print("Read clusterings", flush=True)
clusterings = []
for clus_file in par["clusterings"]:
    adata = ad.read_h5ad(clus_file)
    obs_filt = adata.obs.filter(regex='leiden_[0-9.]+')
    clusterings.append(obs_filt)

print("Merge clusterings", flush=True)
merged = pd.concat(clusterings, axis=1)

print("Read input", flush=True)
input = ad.read_h5ad(par["input"])

input.obsm["clustering"] = merged

print("Store outputs", flush=True)
input.write_h5ad(par["output"], compression="gzip")
