import sys
import anndata as ad
from scib.metrics import cell_cycle
import numpy as np

## VIASH START
par = {
    'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
    'output': 'output.h5ad'
}

meta = {
    'name': 'foo'
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
    obs='obs',
    obsm='obsm',
    uns='uns'
)

print("Copy batch information", flush=True)
adata_integrated.obs['batch'] = adata_solution.obs['batch']

print('Use gene symbols for features', flush=True)
adata_solution.var_names = adata_solution.var['feature_name']

dataset_organism = adata_solution.uns['dataset_organism']
translator = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
}
assert dataset_organism in translator, \
    f"dataset organism ({dataset_organism}) not in translator!"

print('Compute score', flush=True)
score = cell_cycle(
    adata_solution,
    adata_integrated,
    batch_key='batch',
    embed='X_emb',
    organism=translator[dataset_organism],
)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata_solution.uns['dataset_id'],
        'normalization_id': adata_solution.uns['normalization_id'],
        'method_id': adata_integrated.uns['method_id'],
        'metric_ids': [ meta['name'] ],
        'metric_values': [ score ]
    }
)


print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
