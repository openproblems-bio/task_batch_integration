import anndata as ad
import scanpy as sc
import drvi
from drvi.model import DRVI
from drvi.utils.misc import hvg_batch
import pandas as pd
import numpy as np
import warnings
import sys
import scipy.sparse

## VIASH START
par = {
  'input': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
  'output': 'output.h5ad',
  'n_hvg': 2000,
  'n_epochs': 400
}
meta = {
  'name': 'drvi'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print('Reading input files', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/counts',
    obs='obs',
    var='var',
    uns='uns'
)

if par["n_hvg"]: 
     print(f"Select top {par['n_hvg']} high variable genes", flush=True) 
     idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]] 
     adata = adata[:, idx].copy() 

print('Train model with DRVI', flush=True)

DRVI.setup_anndata(
    adata,
    categorical_covariate_keys=["batch"],
    is_count_data=False,
)

model = DRVI(
    adata,
    categorical_covariates=["batch"],
    n_latent=128,
    encoder_dims=[128, 128],
    decoder_dims=[128, 128],
)
model

model.train(
    max_epochs=par["n_epochs"],
    early_stopping=False,
    early_stopping_patience=20,
    plan_kwargs={
        "n_epochs_kl_warmup": par["n_epochs"],
    },
)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs.copy(),
    var=adata.var.copy(),
    obsm={
        "X_emb": model.get_latent_representation(),
    },
    uns={
        "dataset_id": adata.uns.get("dataset_id", "unknown"),
        "normalization_id": adata.uns.get("normalization_id", "unknown"),
        "method_id": meta["name"],
    },
)

print("Write output AnnData to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
