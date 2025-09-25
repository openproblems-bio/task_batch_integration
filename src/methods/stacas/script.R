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
# Transpose because Seurat expects genes in rows, cells in columns
counts_r <- Matrix::t(adata$layers[["counts"]])
normalized_r <- Matrix::t(adata$layers[["normalized"]])
# Convert to a regular sparse matrix first and then to dgCMatrix
counts_c <- as(as(counts_r, "CsparseMatrix"), "dgCMatrix")
normalized_c <- as(as(normalized_r, "CsparseMatrix"), "dgCMatrix")

# Create Seurat object with raw counts, these are needed to compute Variable Genes
seurat_obj <- Seurat::CreateSeuratObject(counts = counts_c,
                                         meta.data = adata$obs)
# Manually assign pre-normalized values to the "data" slot
seurat_obj@assays$RNA$data <- normalized_c

cat("Run STACAS\n")
object_integrated <- seurat_obj |>
      Seurat::SplitObject(split.by = "batch") |>
      STACAS::Run.STACAS() 

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
