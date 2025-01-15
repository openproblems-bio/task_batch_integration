import os
import sys
import tarfile
import tempfile
import zipfile

import anndata as ad
import gdown
import scgpt
import torch

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
    "input": "resources_test/.../input.h5ad",
    "output": "output.h5ad",
    "model_name": "scGPT_human",
    "model": "scGPT_human",
    "n_hvg": 3000,
}
meta = {"name": "scgpt"}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata
from exit_codes import exit_non_applicable

print(f"====== scGPT version {scgpt.__version__} ======", flush=True)

print("\n>>> Reading input files...", flush=True)
print(f"Input H5AD file: '{par['input']}'", flush=True)
adata = read_anndata(par["input"], X="layers/counts", obs="obs", var="var", uns="uns")

if adata.uns["dataset_organism"] != "homo_sapiens":
    exit_non_applicable(
        f"scGPT can only be used with human data "
        f"(dataset_organism == \"{adata.uns['dataset_organism']}\")"
    )

print(adata, flush=True)

print("\n>>> Preprocessing data...", flush=True)
if par["n_hvg"]:
    print(f"Selecting top {par['n_hvg']} highly variable genes", flush=True)
    idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][: par["n_hvg"]]
    adata = adata[:, idx].copy()

print(adata, flush=True)

if par["model"] is None:
    print(f"\n>>> Downloading '{par['model_name']}' model...", flush=True)
    model_drive_ids = {
        "scGPT_human": "1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
        "scGPT_CP": "1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB",
    }
    drive_path = (
        f"https://drive.google.com/drive/folders/{model_drive_ids[par['model_name']]}"
    )
    model_temp = tempfile.TemporaryDirectory()
    model_dir = model_temp.name
    print(f"Downloading from '{drive_path}'", flush=True)
    gdown.download_folder(drive_path, output=model_dir, quiet=True)
else:
    if os.path.isdir(par["model"]):
        print(f"\n>>> Using model directory...", flush=True)
        model_temp = None
        model_dir = par["model"]
    else:
        model_temp = tempfile.TemporaryDirectory()
        model_dir = model_temp.name

        if zipfile.is_zipfile(par["model"]):
            print(f"\n>>> Extracting model from .zip...", flush=True)
            print(f".zip path: '{par['model']}'", flush=True)
            with zipfile.ZipFile(par["model"], "r") as zip_file:
                zip_file.extractall(model_dir)
        elif tarfile.is_tarfile(par["model"]) and par["model"].endswith(
            ".tar.gz"
        ):
            print(f"\n>>> Extracting model from .tar.gz...", flush=True)
            print(f".tar.gz path: '{par['model']}'", flush=True)
            with tarfile.open(par["model"], "r:gz") as tar_file:
                tar_file.extractall(model_dir)
                model_dir = os.path.join(model_dir, os.listdir(model_dir)[0])
        else:
            raise ValueError(
                f"The 'model' argument should be a directory a .zip file or a .tar.gz file"
            )

print(f"Model directory: '{model_dir}'", flush=True)

print("\n>>> Embedding data...", flush=True)
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Device: '{device}'", flush=True)
embedded = scgpt.tasks.embed_data(
    adata,
    model_dir,
    gene_col="feature_name",
    batch_size=64,
    use_fast_transformer=False,  # Disable fast-attn as not installed
    device=device,
    return_new_adata=True,
)

print("\n>>> Storing output...", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": embedded.X,
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"],
    },
)
print(output)

print("\n>>> Writing output to file...", flush=True)
print(f"Output H5AD file: '{par['output']}'", flush=True)
output.write_h5ad(par["output"], compression="gzip")

if model_temp is not None:
    print("\n>>> Cleaning up temporary directories...", flush=True)
    model_temp.cleanup()

print("\n>>> Done!", flush=True)
