"""
Common utilities for MLflow-based methods.
"""
import os
import tempfile

import anndata as ad
import pandas as pd
import sklearn.neighbors


def create_temp_h5ad(
    adata, layers=None, obs=None, var=None, obsm=None, varm=None, uns=None
):
    """
    Create a temporary H5AD file with specified data from an AnnData object.

    Args:
        adata: Input AnnData object
        layers: List of layer names to include (e.g., ["counts"])
        obs: List of obs column names to include (e.g., ["batch"])
        var: Dict mapping var column names to new names (e.g., {"feature_id": "ensembl_id"})
        obsm: List of obsm keys to include
        varm: List of varm keys to include
        uns: List of uns keys to include

    Returns:
        tuple: (h5ad_file, input_adata) where h5ad_file is the NamedTemporaryFile and
               input_adata is the created AnnData object
    """
    # Extract X from layers or use X directly
    if layers and len(layers) > 0:
        X = adata.layers[layers[0]].copy()
    else:
        X = adata.X.copy()

    # Create new AnnData
    input_adata = ad.AnnData(X=X)

    # Set var_names
    input_adata.var_names = adata.var_names

    # Add obs columns
    if obs:
        for obs_key in obs:
            if obs_key in adata.obs:
                input_adata.obs[obs_key] = adata.obs[obs_key].values

    # Add var columns (with optional renaming)
    if var:
        for old_name, new_name in var.items():
            if old_name in adata.var:
                input_adata.var[new_name] = adata.var[old_name].values

    # Add obsm
    if obsm:
        for obsm_key in obsm:
            if obsm_key in adata.obsm:
                input_adata.obsm[obsm_key] = adata.obsm[obsm_key].copy()

    # Add varm
    if varm:
        for varm_key in varm:
            if varm_key in adata.varm:
                input_adata.varm[varm_key] = adata.varm[varm_key].copy()

    # Add uns
    if uns:
        for uns_key in uns:
            if uns_key in adata.uns:
                input_adata.uns[uns_key] = adata.uns[uns_key]

    # Write to temp file
    h5ad_file = tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False)
    input_adata.write(h5ad_file.name)

    return h5ad_file, input_adata


def embed(adata, model, layers=None, obs=None, var=None, model_params=None, process_adata=None):
    """
    Embed data using an MLflow model.

    Args:
        adata: Input AnnData object to embed
        model: Loaded MLflow model
        layers: List of layer names to include (e.g., ["counts"])
        obs: List of obs column names to include (e.g., ["batch"])
        var: Dict mapping var column names to new names (e.g., {"feature_id": "ensembl_id"})
        model_params: Optional dict of parameters to pass to model.predict()
        process_adata: Optional function to process input_adata before writing (e.g., to add defaults)

    Returns:
        np.ndarray: Embeddings for the input data
    """
    print("Writing temporary input H5AD file...", flush=True)
    h5ad_file, input_adata = create_temp_h5ad(adata, layers=layers, obs=obs, var=var)
    
    # Apply any post-processing to input_adata
    if process_adata:
        process_adata(input_adata)
    
    print(f"Temporary H5AD file: '{h5ad_file.name}'", flush=True)
    print(input_adata, flush=True)

    # Re-write the file after processing
    input_adata.write(h5ad_file.name)

    print("Running model...", flush=True)
    input_df = pd.DataFrame({"input_uri": [h5ad_file.name]})
    if model_params:
        embedding = model.predict(input_df, params=model_params)
    else:
        embedding = model.predict(input_df)

    # Clean up
    h5ad_file.close()
    os.unlink(h5ad_file.name)

    return embedding


def embed_and_classify(
    train_adata,
    test_adata,
    model,
    layers=None,
    obs=None,
    var=None,
    model_params=None,
    process_adata=None,
    n_neighbors=5,
):
    """
    Generic pipeline for embedding data and training a kNN classifier.

    Args:
        train_adata: Training AnnData object with labels
        test_adata: Test AnnData object to predict
        model: Loaded MLflow model
        layers: List of layer names to include (e.g., ["counts"])
        obs: List of obs column names to include (e.g., ["batch"])
        var: Dict mapping var column names to new names (e.g., {"feature_id": "ensembl_id"})
        model_params: Optional dict of parameters to pass to model.predict()
        process_adata: Optional function to process input_adata before writing (e.g., to add defaults)
        n_neighbors: Number of neighbors for kNN classifier

    Returns:
        np.ndarray: Predicted labels for test data
    """
    # Embed training data
    print("\n>>> Embedding training data...", flush=True)
    embedding_train = embed(
        train_adata, model, layers=layers, obs=obs, var=var,
        model_params=model_params, process_adata=process_adata
    )

    # Train kNN classifier
    print("\n>>> Training kNN classifier...", flush=True)
    classifier = sklearn.neighbors.KNeighborsClassifier(n_neighbors=n_neighbors)
    classifier.fit(embedding_train, train_adata.obs["label"].astype(str))

    # Embed test data
    print("\n>>> Embedding test data...", flush=True)
    embedding_test = embed(
        test_adata, model, layers=layers, obs=obs, var=var,
        model_params=model_params, process_adata=process_adata
    )

    # Classify
    print("\n>>> Classifying test data...", flush=True)
    predictions = classifier.predict(embedding_test)

    return predictions
