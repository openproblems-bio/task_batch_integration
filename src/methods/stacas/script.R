requireNamespace("anndata", quietly = TRUE)
suppressPackageStartupMessages({
  library(STACAS)
  library(Matrix)
  library(SeuratObject)
  library(Seurat)
})

## VIASH START
par <- list(
  input = "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "stacas"
)
## VIASH END

cat("Reading input file\n")
adata <- anndata::read_h5ad(par[["input"]])

cat("Create Seurat object\n")
# Only loading normalized values, as raw counts are not needed

# Transpose because Seurat expects genes in rows, cells in columns
normalized <- Matrix::t(adata$layers[["normalized"]])
# Convert to a regular sparse matrix first and then to dgCMatrix
normalized <- as(as(normalized, "CsparseMatrix"), "dgCMatrix")

# Create Seurat object
seurat_obj <- Seurat::CreateSeuratObject(counts = normalized,
                                         meta.data = adata$obs)
# Manually assign pre-normalized values to the "data" slot
seurat_obj@assays$RNA$data   <- normalized
seurat_obj@assays$RNA$counts <- NULL # remove counts


# Obtain anchor features from the preprocessing pipeline
anchor.features <- head(adata$var[order(adata$var$hvg_score, decreasing = T), "feature_id"], 2000)

cat("Run STACAS\n")
object_integrated <- seurat_obj |>
      Seurat::SplitObject(split.by = "batch") |>
      STACAS::Run.STACAS(anchor.features = anchor.features) 


cat("Store outputs\n")
output <- anndata::AnnData(
    uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  ),
  obs = adata$obs,
  var = adata$var,
  obsm = list(
    X_emb = object_integrated@reductions$pca@cell.embeddings
  )
)

cat("Write output AnnData to file\n")
output$write_h5ad(par[["output"]], compression = "gzip")
