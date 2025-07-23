import sys
import pandas as pd
import anndata as ad
import scanpy as sc
from scib.metrics.clustering import cluster_optimal_resolution
from scib.metrics import ari, nmi

## VIASH START
par = {
    'adata_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
    'output': 'output.h5ad',
    "resolutions": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
}

meta = {
    'name': 'foo'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', obsp='obsp', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

print('Run optimal Leiden clustering', flush=True)
cluster_key = "leiden"

# get existing clusters
cluster_df = adata.obsm.get('clustering', pd.DataFrame(index=adata.obs_names))
adata.obs = pd.concat([adata.obs, cluster_df], axis=1)

cluster_optimal_resolution(
    adata=adata,
    label_key="cell_type",
    cluster_key=cluster_key,
    cluster_function=sc.tl.leiden,
    resolutions=par["resolutions"],
)

print('Compute ARI score', flush=True)
ari_score = ari(adata, cluster_key=cluster_key, label_key="cell_type")

print('Compute NMI score', flush=True)
nmi_score = nmi(adata, cluster_key=cluster_key, label_key="cell_type")

print('Compute ARI score with batches', flush=True)
ari_b_score = 1 - ari(adata, cluster_key=cluster_key, label_key="batch")

print('Compute NMI score with batches', flush=True)
nmi_b_score = 1 - nmi(adata, cluster_key=cluster_key, label_key="batch")

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        "method_id": adata.uns['method_id'],
        "metric_ids": [ "ari", "nmi", "ari-batch", "nmi-batch" ],
        "metric_values": [ ari_score, nmi_score, ari_b_score, nmi_b_score ]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
