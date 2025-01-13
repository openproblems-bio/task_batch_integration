import sys
import pandas as pd
import anndata as ad
from scib.metrics import isolated_labels_f1

## VIASH START
par = {
    'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
    'output': 'output.h5ad',
    "resolutions": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
}

meta = {
    'name': 'foo',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', obsp='obsp', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

# get existing clusters
cluster_df = adata.obsm.get('clustering', pd.DataFrame(index=adata.obs_names))
adata.obs = pd.concat([adata.obs, cluster_df], axis=1)

print('compute score', flush=True)
score = isolated_labels_f1(
    adata,
    label_key="cell_type",
    batch_key="batch",
    cluster_key="leiden",
    resolutions=par["resolutions"],
    embed=None,
    iso_threshold=None,
    verbose=True,
)
print(score, flush=True)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['name'] ],
        'metric_values': [ score ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
