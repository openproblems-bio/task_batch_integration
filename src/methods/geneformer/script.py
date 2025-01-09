import os
import sys
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd
from geneformer import EmbExtractor, TranscriptomeTokenizer
from huggingface_hub import hf_hub_download

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
    "input": "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
    "output": "output.h5ad",
    "model": "gf-12L-95M-i4096",
}
meta = {"name": "geneformer"}
## VIASH END

n_processors = os.cpu_count()

print(">>> Reading input...", flush=True)
sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

adata = read_anndata(par["input"], X="layers/counts", obs="obs", var="var", uns="uns")

if adata.uns["dataset_organism"] != "homo_sapiens":
    raise ValueError(
        f"Geneformer can only be used with human data "
        f"(dataset_organism == '{adata.uns['dataset_organism']}')"
    )

# Set adata.var_names to gene IDs
adata.var_names = adata.var["feature_id"]
is_ensembl = all(var_name.startswith("ENSG") for var_name in adata.var_names)
if not is_ensembl:
    raise ValueError(f"Geneformer requires adata.var_names to contain ENSEMBL gene ids")

print(f">>> Getting settings for model '{par['model']}'...", flush=True)
model_split = par["model"].split("-")
model_details = {
    "layers": model_split[1],
    "dataset": model_split[2],
    "input_size": int(model_split[3][1:]),
}
print(model_details, flush=True)

print(">>> Getting model dictionary files...", flush=True)
if model_details["dataset"] == "95M":
    dictionaries_subfolder = "geneformer"
elif model_details["dataset"] == "30M":
    dictionaries_subfolder = "geneformer/gene_dictionaries_30m"
else:
    raise ValueError(f"Invalid model dataset: {model_details['dataset']}")
print(f"Dictionaries subfolder: '{dictionaries_subfolder}'")

dictionary_files = {
    "ensembl_mapping": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=dictionaries_subfolder,
        filename=f"ensembl_mapping_dict_gc{model_details['dataset']}.pkl",
    ),
    "gene_median": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=dictionaries_subfolder,
        filename=f"gene_median_dictionary_gc{model_details['dataset']}.pkl",
    ),
    "gene_name_id": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=dictionaries_subfolder,
        filename=f"gene_name_id_dict_gc{model_details['dataset']}.pkl",
    ),
    "token": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=dictionaries_subfolder,
        filename=f"token_dictionary_gc{model_details['dataset']}.pkl",
    ),
}

print(">>> Creating working directory...", flush=True)
work_dir = TemporaryDirectory()
input_dir = os.path.join(work_dir.name, "input")
os.makedirs(input_dir)
tokenized_dir = os.path.join(work_dir.name, "tokenized")
os.makedirs(tokenized_dir)
embedding_dir = os.path.join(work_dir.name, "embedding")
os.makedirs(embedding_dir)
print(f"Working directory: '{work_dir.name}'", flush=True)

print(">>> Preparing data...", flush=True)
adata.var["ensembl_id"] = adata.var_names
adata.obs["n_counts"] = np.ravel(adata.X.sum(axis=1))
adata.write_h5ad(os.path.join(input_dir, "input.h5ad"))
print(adata)

print(">>> Tokenizing data...", flush=True)
special_token = model_details["dataset"] == "95M"
print(f"Input size: {model_details['input_size']}, Special token: {special_token}")
tokenizer = TranscriptomeTokenizer(
    nproc=n_processors,
    model_input_size=model_details["input_size"],
    special_token=special_token,
    gene_median_file=dictionary_files["gene_median"],
    token_dictionary_file=dictionary_files["token"],
    gene_mapping_file=dictionary_files["ensembl_mapping"],
)
tokenizer.tokenize_data(input_dir, tokenized_dir, "tokenized", file_format="h5ad")

print(f">>> Getting model files for model '{par['model']}'...", flush=True)
model_files = {
    "model": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=par["model"],
        filename="model.safetensors",
    ),
    "config": hf_hub_download(
        repo_id="ctheodoris/Geneformer",
        subfolder=par["model"],
        filename="config.json",
    ),
}
model_dir = os.path.dirname(model_files["model"])

print(">>> Extracting embeddings...", flush=True)
embedder = EmbExtractor(
    emb_mode="cell", max_ncells=None, token_dictionary_file=dictionary_files["token"]
)
embedder.extract_embs(
    model_dir,
    os.path.join(tokenized_dir, "tokenized.dataset"),
    embedding_dir,
    "embedding",
)
embedding = pd.read_csv(os.path.join(embedding_dir, "embedding.csv")).to_numpy()

print(">>> Storing outputs...", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": embedding,
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"],
    },
)
print(output)

print(">>> Writing output AnnData to file...", flush=True)
output.write_h5ad(par["output"], compression="gzip")
print(">>> Done!")
