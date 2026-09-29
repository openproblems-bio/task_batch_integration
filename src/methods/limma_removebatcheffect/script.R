cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  library(Matrix, warn.conflicts = FALSE)
  library(limma, warn.conflicts = FALSE)
})

## VIASH START
par <- list(
  input = 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
  output = 'output.h5ad'
)
meta <- list(
  name = "limma_removebatcheffect"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

cat("Extract data and metadata\n")
# Extract normalized data
expr_data <- t(adata$layers[["normalized"]])
batch_info <- adata$obs[["batch"]]
obs <- adata$obs
var <- adata$var

# Convert to dgCMatrix if needed
if (inherits(expr_data, "dgRMatrix")) {
  dense_temp <- as.matrix(expr_data)
  expr_data <- as(dense_temp, "dgCMatrix")
}

cat("Apply limma removeBatchEffect\n")
# Create design matrix (intercept only, as we want to preserve all biological variation)
design <- matrix(1, nrow = ncol(expr_data), ncol = 1)
colnames(design) <- "Intercept"
rownames(design) <- colnames(expr_data)

# Apply batch correction using limma's removeBatchEffect
corrected_data <- limma::removeBatchEffect(
  x = expr_data,
  batch = batch_info,
  design = design
)

cat("Prepare output\n")
# Create output AnnData object with corrected feature matrix
output <- anndata::AnnData(
  obs = obs[, c()],
  var = var[, c()],
  layers = list(
    corrected_counts = t(corrected_data)
  ),
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  )
)

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")

cat("Finished\n")
