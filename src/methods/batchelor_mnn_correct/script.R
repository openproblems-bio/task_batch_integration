cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  requireNamespace("Matrix", quietly = TRUE)
  requireNamespace("batchelor", quietly = TRUE)
  library(SingleCellExperiment, warn.conflicts = FALSE)
})
## VIASH START
par <- list(
  input = 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
  output = 'output.h5ad'
)
meta <- list(
  name = "mnn_correct_feature"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

cat("Run mnnCorrect\n")
out <- suppressWarnings(batchelor::mnnCorrect(
  Matrix::t(adata$layers[["normalized"]]),
  batch = adata$obs[["batch"]]
))

cat("Reformat output\n")
layer <- SummarizedExperiment::assay(out, "corrected")

cat("Store outputs\n")
output <- anndata::AnnData(
  layers = list(
    corrected_counts = out |>
      SummarizedExperiment::assay("corrected") |>
      t() |>
      as("sparseMatrix")
  ),
  obs = adata$obs[, c()],
  var = adata$var[, c()],
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  )
)

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
