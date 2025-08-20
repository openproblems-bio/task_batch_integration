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
  name = "seurat_rpca"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

cat("Create Seurat object\n")
# Extract preprocessed data
counts_data <- t(adata$layers[["counts"]])
norm_data <- t(adata$layers[["normalized"]])
obs <- adata$obs
var <- adata$var

# Convert to dgCMatrix if needed (Seurat v5 compatibility)
if (inherits(counts_data, "dgRMatrix")) {
  dense_temp <- as.matrix(counts_data)
  counts_data <- as(dense_temp, "dgCMatrix")
}

if (inherits(norm_data, "dgRMatrix")) {
  dense_temp <- as.matrix(norm_data)
  norm_data <- as(dense_temp, "dgCMatrix")
}

# Ensure proper dimnames
rownames(counts_data) <- rownames(var)
colnames(counts_data) <- rownames(obs)
rownames(norm_data) <- rownames(var)
colnames(norm_data) <- rownames(obs)

# Create Seurat object from counts
seurat_obj <- CreateSeuratObject(
  counts = counts_data,
  meta.data = obs,
  assay = "RNA"
)

# Add normalized data layer
LayerData(seurat_obj, layer = "data") <- norm_data

# Use existing HVGs from the dataset
hvg_genes <- rownames(adata$var)[adata$var$hvg]
cat("Using", length(hvg_genes), "HVGs from input dataset\n")
VariableFeatures(seurat_obj) <- hvg_genes

# Use existing PCA from input dataset
pca_embeddings <- adata$obsm[["X_pca"]]
rownames(pca_embeddings) <- colnames(seurat_obj)
colnames(pca_embeddings) <- paste0("PC_", seq_len(ncol(pca_embeddings)))

seurat_obj[["pca"]] <- CreateDimReducObject(
  embeddings = pca_embeddings,
  key = "PC_",
  assay = "RNA"
)

cat("Split object by batch\n")
# Split the object by batch for traditional integration
seurat_list <- SplitObject(seurat_obj, split.by = "batch")

# For RPCA, we need to scale and run PCA on each dataset
cat("Scale data and run PCA for RPCA integration\n")
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(x, features = hvg_genes, verbose = FALSE)
  x <- RunPCA(x, features = hvg_genes, npcs = par$dims, verbose = FALSE)
  return(x)
})

cat("Perform RPCA integration\n")
# Find integration anchors using RPCA
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = hvg_genes,
  reduction = "rpca",
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

cat("Scale and run PCA on integrated data\n")
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = par$dims, verbose = FALSE)

cat("Generate UMAP embedding\n")
integrated <- RunUMAP(integrated, reduction = "pca", dims = seq_len(par$dims), verbose = FALSE)

cat("Extract embedding\n")
embedding <- Embeddings(integrated, reduction = "umap")

cat("Store outputs\n")
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

cat("Finished\n")
