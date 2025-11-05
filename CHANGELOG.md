# task_batch_integration devel

## New functionality

* Added `metrics/kbet_pg` and `metrics/kbet_pg_label` components (PR #52).
* Added `methods/stacas` new method (PR #58).
    - Add non-supervised version of STACAS tool for integration of single-cell transcriptomics data. This functionality enables correction of batch effects while preserving biological variability without requiring prior cell type annotations.
* Added `method/drvi` component (PR #61).
* Added `method/fadvi` component.
    - Add FActor Disentangled Variantional Inference (FADVI) for dimentionality reduction
* Added `ARI_batch` and `NMI_batch` to `metrics/clustering_overlap` (PR #68).

* Added `metrics/cilisi` new metric component (PR #57).
    - ciLISI measures batch mixing in a cell type-aware manner by computing iLISI within each cell type and normalizing
        the scores between 0 and 1. Unlike iLISI, ciLISI preserves sensitivity to biological variance and avoids favoring
        overcorrected datasets with removed cell type signals.
        We propose adding this metric to substitute iLISI.

## Minor changes

* Un-pin the scPRINT version and update parameters (PR #51)
* Update scPRINT to better handle large datasets, including a new default model (PR #54)

## Bug fixes

* Update scPRINT to use latest stable version (PR #70)
* Fix kbet dependencies to numpy<2 and scipy<=1.13 (PR #78).

# task_batch_integration 2.0.0

A major update to the OpenProblems framework, switching from a Python-based framework to a Viash + Nextflow-based framework. This update features the same concepts as the previous version, but with a new implementation that is more flexible, scalable, and maintainable.

## Migration

* Added expected input/output interfaces in `src/api` and document them in `README.md`.

* Store common resources used across tasks in a git submodule `common`.

* Methods, metrics, workflows and other components are implemented as Viash components with a per-component Docker image.

## New functionality

* Switched to larger datasets derived from CELLxGENE.

* Added SCimilarity (PR #3).

* Added Geneformer (PR #6).

* Added UCE method (PR #7).

* Added scGPT zero shot (PR #8, #16).

* Added scPRINT (PR #13).

* Added scGPT fine-tuned (PR #17).


## Major changes

* Prefilter batches in HVG overlap metric (PR #9).

* Precompute clustering for some metrics (PR #18).


## Minor changes

* Add arguments for filtering methods in benchmarking workflow (PR #4).

* Update compute environment (PR #5).

* Adjust resources (PR #10).

* Update dependency components (PR #10).

* Update API formats (PR #21, #28, #31).

* Add support for zebrafish and C. elegans (PR #22).

* Bump scIB to v1.1.7 (PR #30).

* Update common submodule (PR #29).

## Bug fixes

* Multiple fixes prior to release (PR #24, #25, #26, #27, #32, #34, #36, #37, #39, #41, #42, #43, #44).

## Documentation

* Update documentation (PR #45).


# task_batch_integration 1.0.0

This version can be found [here](https://github.com/openproblems-bio/openproblems/tree/v1.0.0/openproblems/tasks/_batch_integration).
