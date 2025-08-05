requireNamespace("anndata", quietly = TRUE)
library(CellMixS)
library(zellkonverter)
library(SingleCellExperiment)

## VIASH START
par <- list(
  input_integrated = "resources_test/task_batch_integration/
                      cxg_immune_cell_atlas/integrated_full.h5ad",
  output = "output.h5ad",
  k=70
)
meta <- list(
  name = "cellmixs"
)
## VIASH END

cat("Reading input .h5ad file into SingleCellExperiment\n")
adata <- anndata::read_h5ad(par$input_integrated)
sce <- zellkonverter::readH5AD(par$input_integrated)
sce$batch <- as.character(adata$obs[["batch"]])

cat("Available columns in colData(sce):\n")
print(colnames(colData(sce)))

cat("Run CellMixS\n")
sce <- CellMixS::cms(
  sce,
  k = par[["k"]],
  group = "batch"
)

cat("Convert CellMixS results into AnnData\n")

obs_df <- as.data.frame(colData(sce))

var_df <- as.data.frame(rowData(sce))

X <- as.matrix(assay(sce))

obsm_list <- list()
for (rd_name in reducedDimNames(sce)) {
  obsm_list[[rd_name]] <- reducedDim(sce, rd_name)
}

output <- anndata::AnnData(
  X = X,
  obs = obs_df,
  var = var_df,
  obsm = obsm_list
)

cat("Write output AnnData to file\n")
output$write_h5ad(par[["output"]], compression = "gzip")

cat("Done.\n")
