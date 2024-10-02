import sys
import anndata as ad
import scimilarity

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad",
    "output": "output.h5ad",
    "model": "model_v1.1",
}
meta = {
    "name": "scvi",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print("Read input", flush=True)
adata = read_anndata(par["input"], X="layers/counts", obs="obs", var="var", uns="uns")

if adata.uns["dataset_organism"] != "homo_sapiens":
    raise ValueError(
        f"SCimilarity can only be used with human data "
        f"(dataset_organism == \"{adata.uns['dataset_organism']}\")"
    )

print("Load SCimilarity model", flush=True)
scimilarity_embedding = scimilarity.cell_embedding.CellEmbedding(
    model_path=par["model"]
)
print("SCimilarity version:", scimilarity.__version__)

print("Create input data", flush=True)
# Some of the functions modify the adata so make sure we have a copy
input = ad.AnnData(X=adata.X.copy(), layers={"counts": adata.X.copy()})
# Set input.var_names to gene symbols
input.var_names = adata.var["feature_name"]

print("Align datasets", flush=True)

# Check the number of genes in the dataset and reduce the overlap threshold if
# necessary (mostly for subsampled test datasets)
gene_overlap_threshold = 5000
if 0.8 * input.n_vars < gene_overlap_threshold:
    from warnings import warn

    warn(
        f"The number of genes in the dataset ({input.n_vars}) "
        f"is less than or close to {gene_overlap_threshold}. "
        f"Setting gene_overlap_threshold to 0.8 * n_var ({int(0.8 * input.n_vars)})."
    )
    gene_overlap_threshold = int(0.8 * input.n_vars)

input = scimilarity.utils.align_dataset(
    input,
    scimilarity_embedding.gene_order,
    gene_overlap_threshold=gene_overlap_threshold,
)
input = scimilarity.utils.consolidate_duplicate_symbols(input)

print("Normalizing dataset", flush=True)
input = scimilarity.utils.lognorm_counts(input)

print("Get cell embeddings", flush=True)
cell_embeddings = scimilarity_embedding.get_embeddings(input.X)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": cell_embeddings,
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"],
    },
)
print(output)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
