suppressPackageStartupMessages({
  library(anndataR)
  library(Seurat)
})

## VIASH START
par <- list(
  input = "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
  output = "output.h5ad",
  dims = NULL,
  k_anchor = NULL,
  k_filter = NULL,
  k_score = NULL
)
meta <- list(
  name = "seurat_rpca"
)
## VIASH END

cat("Reading input file\n")
adata <- read_h5ad(par[["input"]])

cat("Create Seurat object\n")
# The benchmark's log_cp10k normalization is Seurat's LogNormalize, so map it to
# the "data" layer instead of calling NormalizeData().
seurat_obj <- adata$as_Seurat(
  x_mapping = NULL,
  layers_mapping = c(data = "normalized"),
  assay_metadata_mapping = FALSE,
  reduction_mapping = FALSE,
  graph_mapping = FALSE,
  misc_mapping = FALSE
)
# Use the benchmark's HVGs instead of FindVariableFeatures() so that feature
# selection is the same across methods.
VariableFeatures(seurat_obj) <- adata$var_names[adata$var$hvg]

cat("Split layers by batch, scale and run PCA\n")
# Seurat v5 integration workflow, see
# https://satijalab.org/seurat/articles/seurat5_integration
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$batch)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
pca_args <- list(object = seurat_obj, verbose = FALSE)
if (!is.null(par$dims)) pca_args$npcs <- par$dims
seurat_obj <- do.call(RunPCA, pca_args)

cat("Run RPCAIntegration\n")
# Only forward arguments that were set so Seurat's own defaults apply otherwise
integrate_args <- list(
  object = seurat_obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)
if (!is.null(par$dims)) integrate_args$dims <- seq_len(par$dims)
if (!is.null(par$k_anchor)) integrate_args$k.anchor <- par$k_anchor
if (!is.null(par$k_filter)) integrate_args$k.filter <- par$k_filter
if (!is.null(par$k_score)) integrate_args$k.score <- par$k_score
seurat_obj <- do.call(IntegrateLayers, integrate_args)

cat("Store outputs\n")
output <- AnnData(
  obs = adata$obs[, character(0)],
  var = adata$var[, character(0)],
  obsm = list(
    X_emb = Embeddings(seurat_obj, reduction = "integrated.rpca")
  ),
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  )
)

cat("Write output AnnData to file\n")
output$write_h5ad(par[["output"]])
