library(anndata)
library(scMerge)
library(Matrix)
library(stats)

## VIASH START
par <- list(
  input = "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "semisupervised_scmerge2"
)
## VIASH END

cat("Reading input files\n")
adata <- anndata::read_h5ad(par$input)

anndataToSemiSupervisedScMerge2 <- function(adata, top_n = 1000, verbose = TRUE) {
  counts <- t(as.matrix(adata$layers[["counts"]]))
  rownames(counts) <- as.character(adata$var_names)
  colnames(counts) <- as.character(adata$obs_names)

  seg_df <- scSEGIndex(exprs_mat = counts)
  seg_df <- seg_df[order(seg_df$segIdx, decreasing = TRUE), , drop = FALSE]
  ctl <- rownames(seg_df)[seq_len(min(top_n, nrow(seg_df)))]

  exprsMat <- t(as.matrix(adata$layers[["normalized"]]))
  rownames(exprsMat) <- as.character(adata$var_names)
  colnames(exprsMat) <- as.character(adata$obs_names)

  batch     <- as.character(adata$obs$batch)
  cellTypes <- as.character(adata$obs$cell_type)

  scMerge2_res <- scMerge2(
    exprsMat = exprsMat,
    batch = batch,
    cellTypes = cellTypes,
    ctl = ctl,
    verbose = verbose
  )

  return(scMerge2_res)
}


cat("Run semi-supervised scMerge2\n")

scMerge2_res <- anndataToSemiSupervisedScMerge2(adata, top_n = 1000, verbose = TRUE)


cat("Store output\n")
corrected_mat <- scMerge2_res$newY
embedding <- prcomp(t(corrected_mat))$x[, 1:10, drop = FALSE]
rownames(embedding) <- colnames(corrected_mat)

output <- anndata::AnnData(
  X = NULL,
  obs = adata$obs[, c()],
  var = NULL,
  obsm = list(
    X_emb = embedding[as.character(adata$obs_names), , drop = FALSE]  # match input cells
  ),
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name 
  ),
  shape = adata$shape
)

cat("Write output AnnData to file\n")
output$write_h5ad(par[["output"]], compression = "gzip")
