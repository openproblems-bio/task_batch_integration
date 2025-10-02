cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  library(Matrix, warn.conflicts = FALSE)
  library(Seurat, warn.conflicts = FALSE)
  library(SeuratObject, warn.conflicts = FALSE)
})

## VIASH START
par <- list(
  input = 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
  output = 'output.h5ad',
  dims = 30L,
  k_anchor = 5L,
  k_filter = 200L,
  k_score = 30L
)
meta <- list(
  name = "seurat_cca"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

cat("Create Seurat object using precomputed data\n")
# Extract preprocessed data
norm_data <- t(adata$layers[["normalized"]])
obs <- adata$obs
var <- adata$var

# Convert to dgCMatrix if needed (Seurat v5 compatibility)
if (inherits(norm_data, "dgRMatrix")) {
  dense_temp <- as.matrix(norm_data)
  norm_data <- as(dense_temp, "dgCMatrix")
}

# Ensure proper dimnames for other matrix types
rownames(norm_data) <- rownames(var)
colnames(norm_data) <- rownames(obs)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = norm_data,
  meta.data = obs,
  assay = "RNA"
)

# In Seurat v5, we need to set the data layer for normalized data
seurat_obj[["RNA"]]$data <- norm_data

cat("Set highly variable genes from input\n")
hvg_genes <- rownames(adata$var)[adata$var$hvg]
cat("Using", length(hvg_genes), "HVGs from input dataset\n")
VariableFeatures(seurat_obj) <- hvg_genes

cat("Split by batch and perform CCA integration\n")
# Split the object by batch
seurat_list <- SplitObject(seurat_obj, split.by = "batch")

# Find integration anchors using CCA
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = hvg_genes,
  dims = seq_len(par$dims),
  k.anchor = par$k_anchor,
  k.filter = par$k_filter,
  k.score = par$k_score,
  verbose = FALSE
)

# Integrate the data
integrated <- IntegrateData(
  anchorset = anchors,
  dims = seq_len(par$dims),
  verbose = FALSE
)

cat("Scale integrated data and run PCA\n")
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = par$dims, verbose = FALSE)

cat("Generate UMAP embedding\n")
integrated <- RunUMAP(
  integrated,
  reduction = "pca",
  dims = seq_len(par$dims),
  verbose = FALSE
)

cat("Extract embedding\n")
embedding <- Embeddings(integrated, reduction = "umap")

cat("Store outputs\n")
output <- anndata::AnnData(
  obs = adata$obs,
  var = adata$var,
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

cat("Finished\n")
