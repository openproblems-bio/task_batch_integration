cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("harmony", quietly = TRUE)

## VIASH START
par <- list(
  input = 'resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad',
  output = 'output.h5ad'
)
meta <- list(
  name = "harmony"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

cat("Run harmony\n")
out <- harmony::RunHarmony(
  data_mat = adata$obsm[["X_pca"]],
  meta_data = adata$obs[["batch"]]
)

cat("Store outputs\n")
output <- anndata::AnnData(
  obs = adata$obs[, c()],
  var = adata$var[, c()],
  obsm = list(
    X_emb = out
  ),
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  )
)

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
