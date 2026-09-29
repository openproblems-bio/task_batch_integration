cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)
requireNamespace("scMerge", quietly = TRUE)
requireNamespace("BiocParallel", quietly = TRUE)
requireNamespace("BiocSingular", quietly = TRUE)
requireNamespace("methods", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
  output = "output.h5ad",
  cell_type_aware = FALSE,
  n_control_genes = 1000L,
  n_dim = 50L
)
meta <- list(
  name = "scmerge2",
  cpus = 1L
)
## VIASH END

n_cpus <- if (is.null(meta$cpus)) 1L else meta$cpus
bpparam <- BiocParallel::MulticoreParam(workers = n_cpus)

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

# both scSEGIndex and scMerge2 want genes in the rows. scMerge2 only coerces a *dense* matrix to
# CsparseMatrix, so a row-compressed one would slip past that check and break later on
exprs_mat <- methods::as(Matrix::t(adata$layers[["normalized"]]), "CsparseMatrix")
rownames(exprs_mat) <- as.character(adata$var_names)
colnames(exprs_mat) <- as.character(adata$obs_names)

cat("Select stably expressed genes\n")
seg_df <- scMerge::scSEGIndex(exprs_mat = exprs_mat, BPPARAM = bpparam)
seg_df <- seg_df[order(seg_df$segIdx, decreasing = TRUE), , drop = FALSE]
ctl <- rownames(seg_df)[seq_len(min(par$n_control_genes, nrow(seg_df)))]

cat("Run scMerge2\n")
out <- scMerge::scMerge2(
  exprsMat = exprs_mat,
  batch = as.character(adata$obs$batch),
  cellTypes = if (par$cell_type_aware) as.character(adata$obs$cell_type) else NULL,
  ctl = ctl,
  use_bpparam = bpparam,
  use_bsparam = BiocSingular::RandomParam(),
  verbose = TRUE
)

cat("Compute embedding\n")
newY <- out$newY
stopifnot(ncol(newY) == adata$n_obs)
if (!is.null(colnames(newY))) {
  newY <- newY[, adata$obs_names, drop = FALSE]
}
# subtracting the estimated unwanted variation makes the corrected matrix dense, so let
# BiocSingular stream it rather than materialising it in one go
corrected <- t(newY)
n_dim <- min(par$n_dim, min(dim(corrected)) - 1L)
embedding <- BiocSingular::runPCA(
  corrected,
  rank = n_dim,
  center = TRUE,
  scale = FALSE,
  BSPARAM = BiocSingular::RandomParam(),
  BPPARAM = bpparam
)$x
rownames(embedding) <- adata$obs_names

cat("Store output\n")
output <- anndata::AnnData(
  obs = adata$obs[, c()],
  var = adata$var[, c()],
  obsm = list(
    X_emb = embedding
  ),
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  )
)

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
