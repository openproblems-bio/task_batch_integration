import sys
import anndata as ad

# check if we can use GPU
USE_GPU = False
try:
    import subprocess
    assert subprocess.run('nvidia-smi', shell=True, stdout=subprocess.DEVNULL).returncode == 0
    from rapids_singlecell.tl import leiden
    USE_GPU = True
except Exception as e:
    from scanpy.tl import leiden

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
    "output": "output.h5ad",
    "resolution": 0.8,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

n_cell_cpu = 300_000

print("Read input", flush=True)
input = read_anndata(par["input"], obs='obs', obsp='obsp', uns='uns')

key = f'leiden_{par["resolution"]}'
kwargs = dict()
if not USE_GPU:
    kwargs |= dict(
        flavor='igraph',
        n_iterations=2,
    )

print(f"Run Leiden clustering with {kwargs}", flush=True)
leiden(
    input,
    resolution=par["resolution"],
    key_added=key,
    **kwargs,
)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=input.obs[[key]],
)
output.write_h5ad(par["output"], compression="gzip")
