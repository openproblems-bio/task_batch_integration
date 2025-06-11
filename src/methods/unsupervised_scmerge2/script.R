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

  # Check overlap with human/mouse scSEG lists
  gene_ids <- rownames(exprsMat)
  species <- NULL
  best_match <- 0

  for (organism in names(seg_list)) {
    scseg_name <- paste0(organism, "_scSEG")
    seg_genes <- seg_list[[organism]][[scseg_name]]
    overlap <- length(intersect(gene_ids, seg_genes))

    if (overlap > best_match) {
      best_match <- overlap
      species <- organism
    }
  }

  if (is.null(species) || best_match == 0) {
    stop("No match found between gene IDs in exprsMat and scSEG lists for human or mouse. ",
         "Please ensure you're using Ensembl IDs for human or mouse, or provide a custom SEG list.")
  }

  message("Detected species: ", species, " (matched ", best_match, " genes)")

  ctl <- seg_list[[species]][[paste0(species, "_scSEG")]]

  scMerge2_res <- scMerge2(
    exprsMat = exprsMat,
    batch = batch,
    ctl = ctl,
    verbose = verbose
  )

  return(scMerge2_res)
}

data("segList_ensemblGeneID")

cat("Run scMerge2\n")

scMerge2_res <- anndataToScMerge2(
  adata = adata,
  seg_list = segList_ensemblGeneID,
  layer = "normalized",
  verbose = TRUE
)


cat("Store output\n")
corrected_mat <- scMerge2_res$newY

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
