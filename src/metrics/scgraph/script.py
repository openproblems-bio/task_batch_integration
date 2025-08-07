import anndata as ad
import sys
from scgraph import scGraph
import pandas
import numpy

## VIASH START
par = {
  'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
  'output': 'output.h5ad',
  "trim_rate": 0.05,
  "thres_batch": 100,
  "thres_celltype": 10,
  "only_umap": True
}
meta = {
  'name': 'scgraph'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Reading input files', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', obsp='obsp', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns
print('cell type:')
adata.obs['cell_type']

print('Compute metrics', flush=True)
uns_metric_ids = [ 'scgraph' ]
uns_metric_values = [ 0.5 ]

from scgraph import scGraph

# Initialize the graph analyzer
scgraph = scGraph(
    adata_path=par["input_integrated"],
    batch_key="batch",
    label_key="cell_type", # why doesn't this work!!
    trim_rate=par["trim_rate"],
    thres_batch=par["thres_batch"],
    thres_celltype=par["thres_celltype"],
    only_umap=par["only_umap"],
)

results = scgraph.main()

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['name'] ],
        'metric_values': [ results ]
    }
)

print("Write output AnnData to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')