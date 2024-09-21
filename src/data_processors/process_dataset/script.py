import sys
import anndata as ad
import openproblems as op

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'hvgs': 2000,
    'obs_label': 'cell_type',
    'obs_batch': 'batch',
    'subset_hvg': False,
    'output': 'output.h5ad'
}
meta = {
    "config": "target/nextflow/batch_integration/process_dataset/.config.vsh.yaml",
    "resources_dir": "src/common/helper_functions"
}
## VIASH END

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_h5ad_by_format import subset_h5ad_by_format

print(">> Load config", flush=True)
config = op.project.read_viash_config(meta["config"])

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])

def compute_batched_hvg(adata, n_hvgs):
    if n_hvgs > adata.n_vars or n_hvgs <= 0:
        hvg_list = adata.var_names.tolist()
    else:
        import scib
        scib_adata = adata.copy()
        scib_adata.X = scib_adata.layers['normalized'].copy()
        hvg_list = scib.pp.hvg_batch(
            scib_adata,
            batch_key='batch',
            target_genes=n_hvgs,
            adataOut=False
        )
    adata.var['hvg'] = adata.var_names.isin(hvg_list)
    return adata

print(f'Select {par["hvgs"]} highly variable genes', flush=True)
adata_with_hvg = compute_batched_hvg(input, n_hvgs=par['hvgs'])

if par['subset_hvg']:
    print('Subsetting to HVG dimensions', flush=True)
    adata_with_hvg = adata_with_hvg[:, adata_with_hvg.var['hvg']].copy()

print(">> Figuring out which data needs to be copied to which output file", flush=True)
# use par arguments to look for label and batch value in different slots
slot_mapping = {
    "obs": {
        "label": par["obs_label"],
        "batch": par["obs_batch"],
    }
}
print(">> Create output object", flush=True)
output_dataset = subset_h5ad_by_format(
    adata_with_hvg,
    config,
    "output_dataset",
    slot_mapping
)
output_solution = subset_h5ad_by_format(
    adata_with_hvg,
    config,
    "output_solution",
    slot_mapping
)

print('Writing adatas to file', flush=True)
output_dataset.write_h5ad(par['output_dataset'], compression='gzip')
output_solution.write_h5ad(par['output_solution'], compression='gzip')
