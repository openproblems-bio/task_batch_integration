import anndata as ad
import sys
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances


def silhouette_batch_custom(
        adata,
        batch_key,
        group_key,
        embed,
        metric='euclidean',
        return_all=False,
        scale=True,
        verbose=True,
        between_cluster_distances='nearest'
):
    """
    Modification of silhouette_batch from scib package (custom silhouette_samples function) to prevent confounding by nested batch effects.
    Absolute silhouette score of batch labels subsetted for each group. Groups are usually cell types in this context.
    between_cluster_distances='nearest' is equivalent to scib original implementation

    :param batch_key: batches to be compared against
    :param group_key: group labels to be subsetted by e.g. cell type
    :param embed: name of column in adata.obsm
    :param metric: see sklearn silhouette score
    :param scale: if True, scale between 0 and 1
    :param return_all: if True, return all silhouette scores and label means
        default False: return average width silhouette (ASW)
    :param between_cluster_distances: one out of 'mean_other', 'furthest', 'nearest'
    :param verbose:
    :return:
        average width silhouette ASW
        mean silhouette per group in pd.DataFrame
        Absolute silhouette scores per group label
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')

    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])

    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue
        
        #Modified
        sil_per_group = silhouette_samples_custom(
            adata_group.obsm[embed],
            adata_group.obs[batch_key],
            metric=metric,
            between_cluster_distances=between_cluster_distances,
        )

        # take only absolute value
        sil_per_group = [abs(i) for i in sil_per_group]

        if scale:
            # scale s.t. highest number is optimal
            sil_per_group = [1 - i for i in sil_per_group]

        #sil_all = sil_all.append(
        #    pd.DataFrame({
        #        'group': [group] * len(sil_per_group),
        #        'silhouette_score': sil_per_group
        #    })
        #)
        
        sil_all = pd.concat([sil_all, pd.DataFrame({
                'group': [group] * len(sil_per_group),
                'silhouette_score': sil_per_group
            })], ignore_index=True)

    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    asw = sil_means['silhouette_score'].mean()

    if verbose:
        print(f'mean silhouette per cell: {sil_means}')

    if return_all:
        return asw, sil_means, sil_all

    return asw

def silhouette_samples_custom(X, labels, metric="euclidean", between_cluster_distances="nearest"):
    """
    Compute the (modified) silhouette scores for the dataset X with the given labels, considering the chosen between_cluster_distances.
    Nearest corresponds to original silhouette definition. 

    For a fast and efficient implementation of the BRAS metric consider the implementation we provide as part of the scib-metrics package. https://github.com/yoseflab/scib-metrics 

    Parameters:
    X : array-like, shape (n_samples, n_features)
        Feature array.
    labels : array-like, shape (n_samples,)
        Labels of each point.
        
    metric : metric for distance calculation, default:"euclidean", alternatives, e.g., "cosine"
    
    between_cluster_distances: one out of "mean_other", "furthest", "nearest"


    Returns:
    score : float
        The (modified) silhouette score.
    """

    # Number of clusters
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    # If there's only one cluster or no clusters, return 0 as silhouette score cannot be computed
    if n_clusters == 1 or n_clusters == 0:
        return 0

    # Initialize silhouette scores
    silhouette_scores = np.zeros(len(X))

    # Calculate pairwise distance matrix
    #distance_matrix = np.linalg.norm(X[:, np.newaxis] - X, axis=2)
    distance_matrix = pairwise_distances(X, metric=metric)
    
    for i in range(len(X)):
        # Points in the same cluster
        same_cluster = labels == labels[i]
        other_clusters = labels != labels[i]
        # Exclude the current point for intra-cluster distance
        same_cluster[i] = False

        # a: Mean distance from i to all other points in the same cluster
        if np.sum(same_cluster) == 0:
            silhouette_scores[i] = 0
            continue
        
        a = np.mean(distance_matrix[i, same_cluster])

        # b: Mean distance from i to all points in the furthest different cluster
        if between_cluster_distances == "furthest":
            b = np.max([
                np.mean(distance_matrix[i, labels == label]) 
                for label in unique_labels if label != labels[i]
            ])
        
        # b: Mean distance from i to all points in any other cluster
        elif between_cluster_distances == "mean_other":
            b = np.mean(distance_matrix[i, other_clusters]) 
            
        # b: Mean distance from i to all points in the nearest different cluster
        else:
            b = np.min([
                np.mean(distance_matrix[i, labels == label]) 
                for label in unique_labels if label != labels[i]
            ])

        # Silhouette score for point i
        silhouette_scores[i] = (b - a) / max(a, b)
    
    return silhouette_scores
## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input_integrated': 'resources_test/.../integrated.h5ad',
  'input_solution': 'resources_test/.../solution.h5ad',
  'output': 'output.h5ad',
  'input_integrated': 'resources_test/task_batch_integration/cxg_immune_cell_atlas/integrated_full.h5ad',
  'output': 'output.h5ad',
}
meta = {
  'name': 'bras'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print('Reading input files', flush=True)
adata = read_anndata(par['input_integrated'], obs='obs', obsm='obsm', uns='uns')
adata.obs = read_anndata(par['input_solution'], obs='obs').obs
adata.uns |= read_anndata(par['input_solution'], uns='uns').uns

print('Compute metrics', flush=True)
score = silhouette_batch_custom(
    adata,
    batch_key='batch',
    group_key="cell_type",
    embed='X_emb',
)

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
