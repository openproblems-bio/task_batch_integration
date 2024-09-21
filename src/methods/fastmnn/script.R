cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  library(Matrix, warn.conflicts = FALSE)
  requireNamespace("batchelor", quietly = TRUE)
  library(SingleCellExperiment, warn.conflicts = FALSE)
})

## VIASH START
par <- list(
  input = 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad',
  output = 'output.h5ad'
)
meta <- list(
  name = "mnn_correct_feature"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

# TODO: pass output of 'multiBatchNorm' to fastMNN

cat("Run mnn\n")
out <- suppressWarnings(batchelor::fastMNN(
  t(adata$layers[["normalized"]]),
  batch = adata$obs[["batch"]]
))

layer <- as(SummarizedExperiment::assay(out, "reconstructed"), "sparseMatrix")
obsm <- SingleCellExperiment::reducedDim(out, "corrected")

cat("Reformat output\n")
output <- anndata::AnnData(
  layers = list(
    corrected_counts = t(layer)
  ),
  obsm = list(
    X_emb = obsm
  ),
  shape = adata$shape,
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  )
)

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
