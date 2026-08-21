requireNamespace("anndata", quietly = TRUE)
suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
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
  name = "seurat_cca"
)
## VIASH END

cat("Reading input file\n")
adata <- anndata::read_h5ad(par[["input"]])

cat("Create Seurat object\n")
# Seurat expects genes in rows, cells in columns, as a dgCMatrix
normalized <- Matrix::t(adata$layers[["normalized"]])
normalized <- as(as(normalized, "CsparseMatrix"), "dgCMatrix")

seurat_obj <- Seurat::CreateSeuratObject(counts = normalized, meta.data = adata$obs)
# The benchmark's log_cp10k normalization is Seurat's LogNormalize, so assign it to
# the "data" layer instead of calling NormalizeData().
seurat_obj[["RNA"]]$data <- normalized
seurat_obj[["RNA"]]$counts <- NULL

# Use the benchmark's HVGs instead of FindVariableFeatures() so that feature
# selection is the same across methods.
VariableFeatures(seurat_obj) <- rownames(adata$var)[adata$var$hvg]

cat("Split layers by batch, scale and run PCA\n")
# Seurat v5 integration workflow, see
# https://satijalab.org/seurat/articles/seurat5_integration
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$batch)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
pca_args <- list(object = seurat_obj, verbose = FALSE)
if (!is.null(par$dims)) pca_args$npcs <- par$dims
seurat_obj <- do.call(RunPCA, pca_args)

cat("Run CCAIntegration\n")
# Only forward arguments that were set so Seurat's own defaults apply otherwise
integrate_args <- list(
  object = seurat_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
if (!is.null(par$dims)) integrate_args$dims <- seq_len(par$dims)
if (!is.null(par$k_anchor)) integrate_args$k.anchor <- par$k_anchor
if (!is.null(par$k_filter)) integrate_args$k.filter <- par$k_filter
if (!is.null(par$k_score)) integrate_args$k.score <- par$k_score
seurat_obj <- do.call(IntegrateLayers, integrate_args)

cat("Store outputs\n")
output <- anndata::AnnData(
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$name
  ),
  obs = adata$obs,
  var = adata$var,
  obsm = list(
    X_emb = Embeddings(seurat_obj, reduction = "integrated.cca")
  )
)

cat("Write output AnnData to file\n")
output$write_h5ad(par[["output"]], compression = "gzip")
