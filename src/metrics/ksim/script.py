import anndata as ad
import sys
import pegasus as pg
import pegasusio
from scipy.sparse import csr_matrix


sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print('Reading input files', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns
print(adata)

print('Convert to pegasusio.MultimodalData...', flush=True)
adata.X = csr_matrix(adata.shape)
mmdata = pegasusio.MultimodalData(adata)

print('Compute metrics', flush=True)
score = pg.calc_kSIM(
    mmdata,
    attr='cell_type',
    rep='emb',
    K=par["K"],
    min_rate=par["min_rate"],
    n_jobs=par["n_jobs"],
    random_state=par["random_state"],
    use_cache=par["use_cache"]
)
print("score:", score)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['name'] ],
        'metric_values': [ score ]
    }
)

print("Write output AnnData to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
