cat(">> Load dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(scMerge)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

## VIASH START
par <- list(
  input = "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "unsupervised_scmerge2"
)
## VIASH END

cat("Reading input files\n")
adata <- anndata::read_h5ad(par$input)
adata$obs["batch"] <- sub("\\+", "plus", adata$obs[["batch"]]) # Replace "+"" characters in batch names

anndataToScMerge2 <- function(adata, seg_list, layer = "normalized", verbose = FALSE) {
  exprsMat_all <- t(as.matrix(adata$layers[[layer]]))
  batch_all <- as.character(adata$obs$batch)

  valid_cells <- !is.na(batch_all)
  exprsMat <- exprsMat_all[, valid_cells, drop = FALSE]
  batch <- batch_all[valid_cells]
  
  ctl_flat <- unlist(seg_list, recursive = FALSE)
  ctl_matches <- ctl_flat[grepl("scSEG$", names(ctl_flat))]
  if (length(ctl_matches) == 0) {
    stop("No stably expressed gene (scSEG) list found in the provided seg_list.")
  }
  ctl <- ctl_matches[[1]]

  scMerge2_res <- scMerge2(
    exprsMat = exprsMat,
    batch = batch,
    ctl = ctl,
    verbose = verbose
  )

  return(scMerge2_res)
}

data("segList_ensemblGeneID") # only for human and mouse- is that okay?

cat("Run scMerge2\n")

scMerge2_res <- anndataToScMerge2(
  adata = adata,
  seg_list = segList_ensemblGeneID,
  layer = "normalized",
  verbose = TRUE
)


cat("Store output\n")
corrected_mat <- scMerge2_res$newY

# PCA as embedding - is this right?
embedding <- prcomp(t(corrected_mat))$x[, 1:10]

rownames(embedding) <- colnames(corrected_mat)

output <- anndata::AnnData(
  X = NULL,
  obs = adata$obs[, c()],
  var = NULL,
  obsm = list(
    X_emb = embedding[rownames(adata), , drop = FALSE]  # match input cells
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
