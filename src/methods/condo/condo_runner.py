"""Shared runner for the ConDo batch-integration method.

Runs the agglomerative batch integrator: pick the seed batch with the
highest per-batch pre-integration silhouette of cell_type on X_pca, then
iteratively merge each next-best compatible (cell-type-overlapping)
neighbour by fitting a ConDo adapter conditioned on cell_type. Batches
in disconnected components of the compatibility graph are passed through
untouched.

Supports two orthogonal axes (set via the ``par`` dict in the viash
script):

* ``transform_type`` : 'location-scale' | 'affine'
* ``rep``            : 'features' | 'pca'
    - 'features': fit on normalized expression, write corrected_counts.
    - 'pca':      fit on obsm['X_pca'], write obsm['X_emb'] (embedding).
* ``hvg_only`` (features only): restrict fit to var['hvg'] columns; pass
  through non-HVGs at source values.
"""
from __future__ import annotations

import sys
from typing import Any

import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix, issparse


def _to_dense(x: Any) -> np.ndarray:
    return x.toarray() if issparse(x) else np.asarray(x)


def _build_adapter(par: dict[str, Any]):
    from condo import ConDoAdapterMMD

    return ConDoAdapterMMD(
        transform_type=par["transform_type"],
        bootstrap_fraction=float(par.get("bootstrap_fraction", 1.0)),
        n_epochs=int(par.get("n_epochs", 5)),
        learning_rate=float(par.get("learning_rate", 1e-3)),
        mmd_size=int(par.get("mmd_size", 40)),
        batch_size=int(par.get("batch_size", 8)),
        weight_decay=float(par.get("weight_decay", 1e-4)),
        patience=int(par.get("patience", 3)),
        random_state=int(par.get("random_state", 42)),
        verbose=0,
        device=par.get("device", "cpu"),
    )


def _read_input(par: dict, meta: dict) -> ad.AnnData:
    sys.path.append(meta["resources_dir"])
    from read_anndata_partial import read_anndata

    rep = par.get("rep", "features")
    # The agglomerative seed selection needs obsm['X_pca'] for per-batch
    # pre-integration silhouette; pull it in both rep modes.
    if rep == "features":
        return read_anndata(
            par["input"], X="layers/normalized",
            obs="obs", obsm="obsm", var="var", uns="uns",
        )
    if rep == "pca":
        return read_anndata(
            par["input"], obs="obs", obsm="obsm", var="var", uns="uns"
        )
    raise ValueError(f"Unknown rep: {rep!r}")


def _pick_target_by_pre_asw(
    adata: ad.AnnData, batches: np.ndarray, cell_types: np.ndarray
) -> tuple[str, dict]:
    """Compute per-batch pre-integration silhouette of cell_type on
    ``obsm['X_pca']``. Returns ``(argmax_batch_label, per_batch_asw_table)``.
    Used as the seed/scoring criterion for the agglomerative integrator."""
    from scib.metrics import silhouette

    if "X_pca" not in adata.obsm:
        raise ValueError(
            "agglomerative seed selection requires obsm['X_pca']; "
            f"obsm keys present: {list(adata.obsm.keys())}"
        )
    batch_labels = np.unique(batches)
    per_batch: dict[str, float] = {}
    for b in batch_labels:
        mask = batches == b
        if mask.sum() < 4:
            per_batch[b] = float("-inf")
            continue
        sub = adata[mask].copy()
        cts = sub.obs["cell_type"].astype(str)
        keep_cts = cts.value_counts()[cts.value_counts() >= 2].index
        keep_mask = cts.isin(keep_cts).values
        if keep_mask.sum() < 4 or len(keep_cts) < 2:
            per_batch[b] = float("-inf")
            continue
        sub2 = sub[keep_mask].copy()
        try:
            s = float(silhouette(sub2, label_key="cell_type", embed="X_pca"))
        except Exception:
            s = float("-inf")
        per_batch[b] = s
    best = max(per_batch.items(), key=lambda kv: kv[1])
    return best[0], per_batch


def _select_feature_columns(adata: ad.AnnData, hvg_only: bool) -> np.ndarray | None:
    """Return a boolean column mask if HVG-only is requested, else None."""
    if not hvg_only:
        return None
    if "hvg" not in adata.var.columns:
        raise ValueError(
            "--hvg_only requires a boolean 'hvg' column in var; "
            f"available columns: {list(adata.var.columns)}"
        )
    mask = adata.var["hvg"].astype(bool).values
    if mask.sum() == 0:
        raise ValueError("--hvg_only requested but no var['hvg'] entries are True")
    return mask


def run_condo(par: dict, meta: dict) -> None:
    rep = par.get("rep", "features")
    hvg_only = bool(par.get("hvg_only", False)) and rep == "features"

    print(f">> Read input (rep={rep}, hvg_only={hvg_only})", flush=True)
    adata = _read_input(par, meta)

    # condo's product_prior dispatches on dtype.kind ∈ {'U','S'} for discrete
    # confounders; pandas .astype(str).values returns object dtype which
    # slips past that check.
    batches = np.asarray(adata.obs["batch"].astype(str).values, dtype="U")
    cell_types = np.asarray(adata.obs["cell_type"].astype(str).values, dtype="U")

    # Per-batch pre-integration silhouette of cell_type on obsm['X_pca'].
    # Used by agglomerative_integrate as the seed criterion (argmax) and
    # the neighbor-ranking score at each merge step.
    _, per_batch_asw = _pick_target_by_pre_asw(adata, batches, cell_types)
    print(">> Agglomerative seed = argmax pre_asw", flush=True)
    for b, s in sorted(per_batch_asw.items(), key=lambda kv: -kv[1]):
        print(f"    pre_asw[{b}] = {s:.4f}", flush=True)

    if rep == "features":
        Y_full = _to_dense(adata.X).astype(np.float64)
        hvg_mask = _select_feature_columns(adata, hvg_only)
        if hvg_mask is not None:
            print(
                f">> HVG-only: fitting on {int(hvg_mask.sum())}/{Y_full.shape[1]} genes",
                flush=True,
            )
            Y = Y_full[:, hvg_mask]
        else:
            Y = Y_full
    else:  # pca
        Y_full = None
        hvg_mask = None
        Y = np.asarray(adata.obsm["X_pca"], dtype=np.float64)

    from condo.batch_integration import agglomerative_integrate

    def _adapter_factory():
        return _build_adapter(par)

    result = agglomerative_integrate(
        Y=Y,
        batches=batches,
        confounders=cell_types,
        batch_score=per_batch_asw,
        adapter_factory=_adapter_factory,
        verbose=True,
    )
    Y_out = result.Y_out

    print(">> Build output", flush=True)
    if rep == "features":
        if hvg_mask is None:
            corrected = Y_out
        else:
            # Start from the original full-feature matrix; overwrite the HVG
            # columns with the adapted values. Non-HVG columns are pass-through.
            corrected = Y_full.copy()
            corrected[:, hvg_mask] = Y_out
        output = ad.AnnData(
            obs=adata.obs[[]],
            var=adata.var[[]],
            layers={"corrected_counts": csr_matrix(corrected)},
            uns={
                "dataset_id": adata.uns["dataset_id"],
                "normalization_id": adata.uns["normalization_id"],
                "method_id": meta["name"],
            },
        )
    else:  # pca -> embedding output
        output = ad.AnnData(
            obs=adata.obs[[]],
            var=adata.var[[]],
            obsm={"X_emb": Y_out},
            uns={
                "dataset_id": adata.uns["dataset_id"],
                "normalization_id": adata.uns["normalization_id"],
                "method_id": meta["name"],
            },
        )

    print(">> Write output", flush=True)
    output.write_h5ad(par["output"], compression="gzip")
