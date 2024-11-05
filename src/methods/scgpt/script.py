import anndata as ad
import scgpt
import sys
import gdown
import tempfile

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input': 'resources_test/.../input.h5ad',
  'output': 'output.h5ad',
  "model" : "scGPT_human",
  "n_hvg": 3000
}
meta = {
  'name': 'scgpt'
}
## VIASH END

print(f"====== scGPT version {scgpt.__version__} ======", flush=True)

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print("\n>>> Reading input files...", flush=True)
print(f"Input H5AD file: '{par['input']}'", flush=True)
adata = read_anndata(par["input"], X="layers/counts", obs="obs", var="var", uns="uns")

if adata.uns["dataset_organism"] != "homo_sapiens":
    raise ValueError(
        f"scGPT can only be used with human data "
        f"(dataset_organism == \"{adata.uns['dataset_organism']}\")"
    )

print(adata, flush=True)

print("\n>>> Preprocessing data...", flush=True)
if par["n_hvg"]:
    print(f"Selecting top {par['n_hvg']} highly variable genes", flush=True)
    idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    adata = adata[:, idx].copy()

print(adata, flush=True)

print("\n>>> Downloading model...", flush=True)
model_drive_ids = {
    "scGPT_human" : "1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
    "scGPT_CP" : "1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB"
}
drive_path = f"https://drive.google.com/drive/folders/{model_drive_ids[par['model']]}"
model_dir = tempfile.TemporaryDirectory()
print(f"Downloading from '{drive_path}'", flush=True)
gdown.download_folder(drive_path, output=model_dir.name, quiet = True)
print(f"Model directory: '{model_dir.name}'", flush=True)

print('Preprocess data', flush=True)
# ... preprocessing ...

print('Train model', flush=True)
# ... train model ...

print('Generate predictions', flush=True)
# ... generate predictions ...

print("Write output AnnData to file", flush=True)
output = ad.AnnData(

)
output.write_h5ad(par['output'], compression='gzip')

print("\n>>> Cleaning up temporary directories...", flush=True)
model_dir.cleanup()

print("\n>>> Done!", flush=True)
