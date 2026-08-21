# task_batch_integration devel

## New functionality

* Added `metrics/kbet_pg` and `metrics/kbet_pg_label` components (PR #52).
* Add `methods/ss_stacas` new method (PR #59).
    - Add semi-supervised version of STACAS tool for integration of single-cell transcriptomics data. This functionality leverages partial or imperfect knowledge of cell identity to improve integration quality by preserving biological variation while correcting for batch effects.
* Added `methods/stacas` new method (PR #58).
    - Add non-supervised version of STACAS tool for integration of single-cell transcriptomics data. This functionality enables correction of batch effects while preserving biological variability without requiring prior cell type annotations.
* Added `method/drvi` component (PR #61).
* Added `method/sca` component (PR #89).
    - Add Surprisal Causal Analysis (SCA) for dimensionality reduction
* Added `method/fadvi` component.
    - Add FActor Disentangled Variantional Inference (FADVI) for dimentionality reduction
* Added `ARI_batch` and `NMI_batch` to `metrics/clustering_overlap` (PR #68).
* Added `methods/condo` new method (PR #83).
    - ConDo (Confounded Domain Adaptation) is a feature-space batch correction
        method that fits a linear transform to match the conditional distribution
        of each batch's features given cell type. Batches are integrated by walking
        a compatibility graph (batches sharing at least one cell type are connected)
        and iteratively merging each best-scoring neighbour into a growing target
        pool. Affine and location-scale variants are provided.

* Added `metrics/cilisi` new metric component (PR #57).
    - ciLISI measures batch mixing in a cell type-aware manner by computing iLISI within each cell type and normalizing
        the scores between 0 and 1. Unlike iLISI, ciLISI preserves sensitivity to biological variance and avoids favoring
        overcorrected datasets with removed cell type signals.
        We propose adding this metric to substitute iLISI.

* Added `method/limma_removebatcheffect` component (PR #79).

* Added ComBat-Seq method (PR #55).
* Added `methods/scmerge2` component (PR #63).

## Minor changes

* Un-pin the scPRINT version and update parameters (PR #51)
* Update scPRINT to better handle large datasets, including a new default model (PR #54)
* Credit contributors missing from the authors list, and fix Martin Kim's orcid (PR #100).

## Bug fixes

* Update scPRINT to use latest stable version (PR #70)
* Fix kbet dependencies to numpy<2 and scipy<=1.13 (PR #78).
* Fix `render_readme` crashing on `comp_process_integration.yaml`'s absolute `__merge__` paths.
* Fix `metrics/kbet_pg` and `metrics/kbet_pg_label` failing to build: drop the zarr, pandas and numpy bounds pegasuspy
    doesn't ask for, and pin `setuptools<81` so `pkg_resources` is still there.
* Fix `methods/harmonypy` writing a transposed embedding: harmonypy 2.0 returns `Z_corr` as cells x pcs.
* Fix `methods/stacas` erroring on `GetAssayData(slot = )`, defunct since SeuratObject 5.0.0. Tracks STACAS master
    until upstream cuts a release containing the fix.
* Fix `methods/scalex` failing to build: its numpy and torch bounds contradicted what scalex asks for. Keeps
    `numpy<2`, since scalex still calls `np.Inf`.
* Fix `methods/pyliger` failing to build: louvain has no python 3.12 wheel and needs cmake to build igraph from source.
* Bump `methods/cellplm`, `methods/condo`, `methods/drvi` and `metrics/bras` from base image `:1.0.0` to `:1`, so their
    `openproblems` is new enough for the component tests in `common`.
* Fix `methods/geneformer` failing to build: pip's `--filter=blob:none` clone of the huggingface repo no longer works,
    so clone it ourselves. Also needs `transformers<5`, which still has `SpecialTokensMixin`.
* Fix `methods/cellplm` failing to build: drop the pytorch base image's broken `/usr/local/bin/cmake` shim, which
    shadowed the apt cmake that louvain needs.
* Fix `methods/scgpt_zeroshot` and `methods/scgpt_finetuned` failing to build: stop installing `flash-attn`. Both
    scripts pass `use_fast_transformer=False`, so it was never used. Behind it sat three more pins with no python 3.12
    wheels, now `numpy<2`, `torchtext==0.17.2` and `transformers==4.36.2`.
* Fix `methods/scprint` erroring with `Triton only support CUDA 10.0 or higher`: point triton at the ptxas it bundles,
    since its version table stops at CUDA 12 and the base image now ships CUDA 13.
* Fix `methods/scgpt_finetuned` erroring on `'csr_matrix' object has no attribute 'A'`: `.A` was removed from
    scipy sparse matrices in scipy 1.14.
* Fix `methods/geneformer` erroring with a 404 on every dictionary and model file: upstream reorganised the repo for
    Geneformer V2 and dropped the gc95M assets from `main`, so pin the downloads to the last revision that has them.
* Give `methods/fadvi` a `gpu` label. Without one it was scheduled on the CPU partition with no GPU attached, so it
    trained on the CPU and hit its walltime on every dataset.
* Give `methods/scalex` a `midcpu` label instead of `lowcpu`. It is CPU-bound only because its engine has no CUDA.

* Split Scanorama into two methods/scores
    - Split Scanorama into embedding (integrate) and count-correction (correct) modes, instead of running both together. 
        This makes clear what the reported score(s) are describing, and also corrects the misleadingly low score that 
        the combined method receives. The scores for each component  are in line with their scores from v1, where the modes 
        were separated.  

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

* Added Density-Adaptive BBSG method.


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
