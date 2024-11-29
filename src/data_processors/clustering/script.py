import anndata as ad
import scanpy as sc
USE_GPU = False
try:
    import subprocess
    assert subprocess.run('nvidia-smi', shell=True, stdout=subprocess.DEVNULL).returncode == 0
    from rapids_singlecell.tl import leiden, louvain
    USE_GPU = True
except Exception as e:
    from scanpy.tools import leiden, louvain

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
    "output": "output.h5ad",
    "algorithm": "leiden",
    "resolution": 0.8,
    "level": 1
}
## VIASH END

n_cell_cpu = 300_000

print("Read input", flush=True)
input = ad.read_h5ad(par["input"])

key = f'{par["algorithm"]}_r{par["resolution"]}_l{par["level"]}'

algorithm_map = {
    'louvain': louvain,
    'leiden': leiden,
}

cluster_fun = algorithm_map.get(par["algorithm"], KeyError(f"Algorithm {par['algorithm']} not supported"))

cluster_fun(
    input,
    resolution=par["resolution"],
    neighbors_key="knn",
    key_added=key,
)

clust = input.obs[key]

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=input.obs[[key]],
)
output.write_h5ad(par["output"], compression="gzip")
