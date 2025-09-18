library(anndata)
library(scIntegrationMetrics)

## VIASH START
par <- list(
  input_integrated = "resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_processed.h5ad",
  input_solution = "resources_test/task_batch_integration/cxg_immune_cell_atlas/solution.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "cilisi"
)
## VIASH END

cat("Reading input files\n")
adata <- anndata::read_h5ad(par[["input_integrated"]])
solution <- anndata::read_h5ad(par[["input_solution"]])
embeddings <- adata$obsm[["X_emb"]]
metadata <- solution$obs

cat("Compute CiLISI metrics...\n")
lisisplit <-
  scIntegrationMetrics::compute_lisi_splitBy(
    X = embeddings,
    meta_data = metadata,
    label_colnames = "batch",
    perplexity = 30,
    split_by_colname = "cell_type",
    normalize = TRUE,
    min.cells.split = 10,
    min.vars.label = 2
)
# average CiLISI
cilisi <- mean(unlist(lisisplit))
# Mean per cell type
cilisi_means <- mean(sapply(lisisplit, function(x) mean(x[, 1])))

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = c(1,2),
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = adata$uns[["method_id"]],
    metric_ids = c("cilisi", "cilisi_means"),
    metric_values = list(cilisi, cilisi_means)
  )
)
output$write_h5ad(par[["output"]], compression = "gzip")
