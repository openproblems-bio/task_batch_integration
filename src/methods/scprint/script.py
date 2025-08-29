import os
import sys

import anndata as ad
import scprint
import torch
from huggingface_hub import hf_hub_download
from scdataloader import Preprocessor
from scdataloader.utils import load_genes
from scprint import scPrint
from scprint.tasks import Embedder

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
    "output": "output.h5ad",
    "model_name": "v2-medium",
    "model": None,
}
meta = {"name": "scprint"}
## VIASH END

sys.path.append(meta["resources_dir"])
from exit_codes import exit_non_applicable
from read_anndata_partial import read_anndata

print(f"====== scPRINT version {scprint.__version__} ======", flush=True)

# Set suggested PyTorch environment variable
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

print("\n>>> Reading input data...", flush=True)
input = read_anndata(par["input"], X="layers/counts", obs="obs", var="var", uns="uns")
if (
    "organism_ontology_term_id" not in input.obs.columns
    and "dataset_organism" in input.uns
):
    if input.uns["dataset_organism"] == "homo_sapiens":
        input.obs["organism_ontology_term_id"] = "NCBITaxon:9606"
    elif input.uns["dataset_organism"] == "mus_musculus":
        input.obs["organism_ontology_term_id"] = "NCBITaxon:10090"
    else:
        exit_non_applicable(
            f"scPRINT requires human or mouse data, not '{input.uns['dataset_organism']}'"
        )
adata = input.copy()

print("\n>>> Preprocessing data...", flush=True)
preprocessor = Preprocessor(
    min_valid_genes_id=min(0.9 * adata.n_vars, 10000),  # 90% of features up to 10,000
    # Turn off cell filtering to return results for all cells
    filter_cell_by_counts=False,
    min_nnz_genes=False,
    do_postp=False,
    # Skip ontology checks
    skip_validate=True,
)
adata = preprocessor(adata)

model_checkpoint_file = par["model"]
if model_checkpoint_file is None:
    print(f"\n>>> Downloading '{par['model_name']}' model...", flush=True)
    model_checkpoint_file = hf_hub_download(
        repo_id="jkobject/scPRINT", filename=f"{par['model_name']}.ckpt"
    )

print(f"Model checkpoint file: '{model_checkpoint_file}'", flush=True)

if torch.cuda.is_available():
    print("CUDA is available, using GPU", flush=True)
    transformer = "flash"
else:
    print("CUDA is not available, using CPU", flush=True)
    transformer = "normal"

# make sure that you check if you have a GPU with flashattention or not (see README)
try:
    m = torch.load(model_checkpoint_file)
# if not use this instead since the model weights are by default mapped to GPU types
except RuntimeError:
    m = torch.load(model_checkpoint_file, map_location=torch.device("cpu"))

# both are for compatibility issues with different versions of the pretrained model, so we need to load it with the correct transformer
if "prenorm" in m["hyper_parameters"]:
    m["hyper_parameters"].pop("prenorm")
    torch.save(m, model_checkpoint_file)
if "label_counts" in m["hyper_parameters"]:
    # you need to set precpt_gene_emb=None otherwise the model will look for its precomputed gene embeddings files although they were already converted into model weights, so you don't need this file for a pretrained model
    model = scPrint.load_from_checkpoint(
        model_checkpoint_file,
        precpt_gene_emb=None,
        classes=m["hyper_parameters"]["label_counts"],
        transformer=transformer,
    )
else:
    model = scPrint.load_from_checkpoint(
        model_checkpoint_file, precpt_gene_emb=None, transformer=transformer
    )
del m
# this might happen if you have a model that was trained with a different set of genes than the one you are using in the ontology (e.g. newer ontologies), While having genes in the onlogy not in the model is fine. the opposite is not, so we need to remove the genes that are in the model but not in the ontology
missing = set(model.genes) - set(load_genes(model.organisms).index)
if len(missing) > 0:
    print(
        "Warning: some genes missmatch exist between model and ontology: solving...",
    )
    model._rm_genes(missing)

# again if not on GPU you need to convert the model to float32
if not torch.cuda.is_available():
    model = model.to(torch.float32)

# you can perform your inference on float16 if you have a GPU, otherwise use float64
dtype = torch.float16 if torch.cuda.is_available() else torch.float32

# the models are often loaded with some parts still displayed as "cuda" and some as "cpu", so we need to make sure that the model is fully on the right device
model = model.to("cuda" if torch.cuda.is_available() else "cpu")

print("\n>>> Embedding data...", flush=True)
n_cores = min(len(os.sched_getaffinity(0)), 24)
print(f"Using {n_cores} worker cores")
embedder = Embedder(
    how="random expr",
    batch_size=par["batch_size"],
    max_len=par["max_len"],
    add_zero_genes=0,
    num_workers=n_cores,
    doclass=False,
    doplot=False,
    pred_embedding=["cell_type_ontology_term_id"],
    keep_all_cls_pred=False,
    output_expression="none",
    save_every=30_000,
    dtype=dtype,
)
embedded, _ = embedder(model, adata, cache=False)

print("\n>>> Storing output...", flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    var=input.var[[]],
    obsm={
        "X_emb": embedded.obsm["scprint_emb"],
    },
    uns={
        "dataset_id": input.uns["dataset_id"],
        "normalization_id": input.uns["normalization_id"],
        "method_id": meta["name"],
    },
)
print(output)

print("\n>>> Writing output AnnData to file...", flush=True)
output.write_h5ad(par["output"], compression="gzip")

print("\n>>> Done!", flush=True)
