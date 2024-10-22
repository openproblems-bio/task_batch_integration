import sys
import anndata as ad
from scib.metrics import hvg_overlap

## VIASH START
par = {
    'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
    'input_solution': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/solution.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'name': 'foo',
    "resources_dir": "src/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print('Read input', flush=True)
adata_solution = read_anndata(
    par['input_solution'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)
adata_integrated = read_anndata(
    par['input_integrated'],
    X='layers/corrected_counts',
    obs='obs',
    var='var',
    uns='uns'
)

print("Copy batch information", flush=True)
adata_integrated.obs['batch'] = adata_solution.obs['batch']

print('Compute score', flush=True)
score = hvg_overlap(
    adata_solution,
    adata_integrated,
    batch_key="batch"
)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": adata_solution.uns['dataset_id'],
        'normalization_id': adata_solution.uns['normalization_id'],
        "method_id": adata_integrated.uns['method_id'],
        "metric_ids": [meta['name']],
        "metric_values": [score]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
