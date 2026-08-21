import sys

import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors

## VIASH START
par = {
    'input': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_pca_components': 100,
    'n_neighbors_per_batch': 10,
    'total_k_neighbors': 50
}
meta = {
    'name': 'era_combat_bbknn',
    'resources_dir': 'src/utils'
}
## VIASH END

sys.path.append(meta['resources_dir'])
from read_anndata_partial import read_anndata


def batch_aware_knn(X, batches, n_neighbors_per_batch, total_k_neighbors):
    """Construct a batch-balanced nearest neighbour graph.

    Neighbours are searched within each batch separately, so that every batch is
    represented among the candidates of a cell. The candidates are then merged and
    the `total_k_neighbors` closest ones are retained.

    Returns a symmetric distance matrix and its binary connectivities matrix.
    """
    # fit one nearest neighbour model per batch, capping k at the size of the batch
    batch_indices = {batch: np.where(batches == batch)[0] for batch in np.unique(batches)}
    models = {}
    for batch, indices in batch_indices.items():
        k = min(n_neighbors_per_batch, len(indices) - 1)
        if k >= 1:
            models[batch] = (NearestNeighbors(n_neighbors=k, metric='euclidean').fit(X[indices]), k)

    # query every batch for every cell and keep the closest candidates overall
    rows = []
    cols = []
    data = []
    for query_indices in batch_indices.values():
        candidate_cols = []
        candidate_data = []
        for batch, (model, k) in models.items():
            distances, indices = model.kneighbors(X[query_indices], n_neighbors=k)
            candidate_cols.append(batch_indices[batch][indices])
            candidate_data.append(distances)
        candidate_cols = np.hstack(candidate_cols)
        candidate_data = np.hstack(candidate_data)

        # a cell is never a neighbour of itself
        candidate_data[candidate_cols == query_indices[:, None]] = np.inf

        selection = np.argsort(candidate_data, axis=1)[:, :total_k_neighbors]
        selected_data = np.take_along_axis(candidate_data, selection, axis=1)
        selected_cols = np.take_along_axis(candidate_cols, selection, axis=1)
        keep = np.isfinite(selected_data)

        rows.append(np.repeat(query_indices, keep.sum(axis=1)))
        cols.append(selected_cols[keep])
        data.append(selected_data[keep])

    # symmetrise, keeping the larger distance of a reciprocal pair
    n_obs = X.shape[0]
    distances = csr_matrix(
        (np.hstack(data), (np.hstack(rows), np.hstack(cols))),
        shape=(n_obs, n_obs)
    )
    distances = distances.maximum(distances.T)
    distances.eliminate_zeros()

    connectivities = distances.copy()
    connectivities.data[:] = 1.0

    return distances, connectivities


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/counts',
    obs='obs',
    var='var',
    uns='uns'
)

print('Normalize data', flush=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

print('Run ComBat', flush=True)
sc.pp.combat(adata, key='batch')

print('Run PCA', flush=True)
n_comps = min(par['n_pca_components'], adata.n_vars, adata.n_obs - 1)
sc.tl.pca(adata, n_comps=n_comps, svd_solver='arpack')

print('Build batch-aware kNN graph', flush=True)
distances, connectivities = batch_aware_knn(
    adata.obsm['X_pca'],
    batches=adata.obs['batch'].values,
    n_neighbors_per_batch=par['n_neighbors_per_batch'],
    total_k_neighbors=par['total_k_neighbors']
)

print('Store output', flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        'X_emb': adata.obsm['X_pca'],
    },
    obsp={
        'connectivities': connectivities,
        'distances': distances,
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['name'],
        'neighbors': {
            'connectivities_key': 'connectivities',
            'distances_key': 'distances',
            'params': {
                'n_neighbors': par['total_k_neighbors'],
                'n_neighbors_per_batch': par['n_neighbors_per_batch'],
                'n_pcs': n_comps,
                'method': 'era_combat_bbknn',
                'metric': 'euclidean'
            }
        }
    }
)

print('Write output to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
