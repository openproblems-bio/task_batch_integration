import anndata as ad
import scanpy as sc

## VIASH START
par = {
    "input_integrated": "resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/integrated.h5ad",
    "input_dataset": "resources_test/task_batch_integration/cxg_mouse_pancreas_atlas/dataset.h5ad",
    "expected_method_types": ["feature"],
    "ouput": "output.h5ad"
}
## VIASH END

print("Read input", flush=True)
integrated = ad.read_h5ad(par["input_integrated"])
dataset = ad.read_h5ad(par["input_dataset"])

print(f"Format of integrated: {integrated}", flush=True)
print(f"Format of dataset: {dataset}", flush=True)

print("Checking shapes", flush=True)
assert integrated.shape[0] == dataset.shape[0], "Number of cells do not match"
assert integrated.shape[1] == dataset.shape[1], "Number of genes do not match"

print("Checking index", flush=True)
if not integrated.obs.index.equals(dataset.obs.index):
    assert integrated.obs.index.sort_values().equals(dataset.obs.index.sort_values()), "Cell names do not match"
    print("Reordering cells", flush=True)
    integrated = integrated[dataset.obs.index]

if "corrected_counts" in integrated.layers.keys() and \
    not integrated.var.index.equals(dataset.var.index):
    assert integrated.var.index.sort_values().equals(dataset.var.index.sort_values()), "Gene names do not match"
    print("Reordering genes", flush=True)
    integrated = integrated[:, dataset.var.index]

print("Checking method output based on type", flush=True)
if "feature" in par["expected_method_types"]:
    assert "corrected_counts" in integrated.layers.keys(), \
        "Method output is missing \"corrected_counts\" layer."

if "embedding" in par["expected_method_types"]:
    assert "X_emb" in integrated.obsm.keys(), \
        "Method output is missing \"X_emb\" obsm"

if "neighbors" in par["expected_method_types"]:
    assert "neighbors" in integrated.uns.keys(), \
        "Method output is missing \"neighbors\" uns"

if "corrected_counts" in integrated.layers and "X_emb" not in integrated.obsm:
    print("Run PCA...", flush=True)
    integrated.obsm["X_emb"] = sc.pp.pca(
        integrated.layers["corrected_counts"],
        n_comps=50,
        use_highly_variable=False,
        svd_solver="arpack",
        return_info=False
    )

if "X_emb" in integrated.obsm and "neighbors" not in integrated.uns:
    print("Run kNN...", flush=True)
    sc.pp.neighbors(integrated, use_rep="X_emb")

print("Store outputs", flush=True)
integrated.write_h5ad(par["output"], compression="gzip")
