import sys
import anndata as ad

## VIASH START
par = {
    'input_dataset': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
    'output': 'output.h5ad',
}
meta = {
    "resources_dir": "src/tasks/batch_integration/control_methods/"
}
## VIASH END

# add helper scripts to path
sys.path.append(meta["resources_dir"])
from utils import _randomize_features

print('Read input_dataset', flush=True)
adata = ad.read_h5ad(par['input_dataset'])

print("Randomise", flush=True)
corrected_counts = _randomize_features(
    adata.layers["normalized"],
    partition=adata.obs["cell_type"],
)

output = ad.AnnData(
    layers={"corrected_counts": corrected_counts},
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"]
    }
)

print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
