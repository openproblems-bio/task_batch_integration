import sys
import anndata as ad
import numpy as np
import pandas as pd
import torch

# Memory-efficient OpenProblems BERNN adapter uses BERNN internals below.

## VIASH START
par = {
    "input": "resources_test/task_batch_integration/cxg_immune_cell_atlas/dataset.h5ad",
    "output": "output.h5ad",
    "n_hvg": 2000,
    "n_epochs": 200,
    "warmup": 50,
    "early_stop": 5,
    "early_warmup_stop": 5,
    "batch_size": 1024,
    "n_layers": 2,
    "layer1": 512,
    "dloss": "inverseTriplet",
}
meta = {"name": "bernn"}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

import gc
from sklearn.preprocessing import StandardScaler
from bernn.openproblems import (
    TrainingConfig,
    _OpenProblemsFinalEpochTrainer,
    _dense_frame,
    _select_feature_indices,
)
import bernn.dl.train.train_ae_classifier_holdout as _bernn_holdout_module

_ORIGINAL_BERNN_SCALE_DATA = _bernn_holdout_module.scale_data


def _memory_efficient_standard_scale_data(scale, data, device="cpu"):
    # OpenProblems BERNN uses standard scaling. Scale the full train matrix once
    # in float32 and reuse it as "all" instead of materializing train/valid/test/all.
    if scale != "standard":
        return _ORIGINAL_BERNN_SCALE_DATA(scale, data, device)

    train = data["inputs"]["train"]
    scaler = StandardScaler()
    scaler.fit(np.asarray(train, dtype=np.float32))

    for group in ("train", "valid", "test"):
        current = data["inputs"].get(group)
        if current is None or current.shape[0] == 0:
            continue
        columns = current.columns
        index = current.index
        values = np.asarray(current, dtype=np.float32)
        values = scaler.transform(values).astype(np.float32, copy=False)
        data["inputs"][group] = pd.DataFrame(values, columns=columns, index=index)

    data["inputs"]["all"] = data["inputs"]["train"]
    return data, scaler


class _MemoryEfficientOpenProblemsTrainer(_OpenProblemsFinalEpochTrainer):
    def _alias_all_to_train(self):
        if self.data is None:
            return
        for key in ("inputs", "names", "labels", "cats", "batches", "orders", "sets"):
            if key in self.data and "train" in self.data[key]:
                self.data[key]["all"] = self.data[key]["train"]

    def _prepare_data(self, *args, **kwargs):
        super()._prepare_data(*args, **kwargs)
        if getattr(self, "_no_internal_validation", False):
            self._alias_all_to_train()
            gc.collect()

    def make_samples_weights(self):
        super().make_samples_weights()
        if getattr(self, "_no_internal_validation", False):
            self._alias_all_to_train()
            gc.collect()


def _fit_openproblems_memory_efficient_once(
    adata,
    *,
    layer="normalized",
    batch_key="batch",
    label_key="cell_type",
    n_hvg=2000,
    dloss="inverseTriplet",
    n_epochs=200,
    warmup=50,
    early_stop=5,
    early_warmup_stop=5,
    batch_size=256,
    n_layers=2,
    layer1=256,
    scaler="standard",
    learning_rate=1e-3,
    weight_decay=1e-5,
    dropout=0.1,
    margin=1.0,
    smoothing=0.1,
    nu=1.0,
    random_state=1,
    num_workers=0,
    device=None,
    return_trainer=False,
):
    if batch_key not in adata.obs:
        raise KeyError(f"AnnData obs is missing required batch key {batch_key!r}")
    if label_key not in adata.obs:
        raise KeyError(f"AnnData obs is missing required label key {label_key!r}")

    feature_idx = _select_feature_indices(adata, n_hvg)
    frame = _dense_frame(adata, feature_idx, layer)
    labels = adata.obs[label_key].astype(str).to_numpy()
    batches = adata.obs[batch_key].astype(str).to_numpy()

    # BERNN training loaders use drop_last=True. When the batch size exceeds
    # the dataset size, train_classifier returns None because there are no batches.
    # Preserve the requested batch size on large datasets and use two full
    # batches on smaller datasets, including the 600-cell Viash fixture.
    requested_bs = int(batch_size)
    batch_size = min(requested_bs, max(2, len(frame) // 2))
    if batch_size != requested_bs:
        print(f"BERNN: batch size {requested_bs} -> {batch_size} for {len(frame)} cells", flush=True)

    # The dense float32 frame is now the training source of truth. Release the
    # input expression matrix before BERNN creates training tensors.
    if layer in ("X", None):
        try:
            adata.X = None
        except Exception:
            pass
    gc.collect()

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    config = TrainingConfig(
        optimize_hyperparams=False, dloss=dloss, class_triplet=False,
        variational=False, kan=False, n_layers=int(n_layers), layer1=int(layer1),
        tied_weights=False, use_mapping=True, rec_loss="l1", scaler=scaler,
        log1p=False, use_l1=True, prune_network=False, update_grid=False,
        n_epochs=int(n_epochs), warmup=int(warmup), n_repeats=1,
        bs=int(batch_size), num_workers=int(num_workers), groupkfold=True,
        device=device, dataset=str(getattr(adata, "uns", {}).get("dataset_id", "openproblems")),
        exp_id="bernn_openproblems_exact",
        early_stop=int(early_stop),
        early_warmup_stop=int(early_warmup_stop),
    )
    config.lr = float(learning_rate)
    config.wd = float(weight_decay)
    config.dropout = float(dropout)
    config.margin = float(margin)
    config.smoothing = float(smoothing)
    config.nu = float(nu)
    config.thres = 0.0
    config.gamma = 0.1
    config.beta = 0.0

    trainer = _MemoryEfficientOpenProblemsTrainer(
        config=config, groupkfold=True, pools=False, keep_models=False,
        log_inputs=False, log_plots=False, log_tb=False, log_mlflow=False,
        log_dvclive=False,
    )
    params = {
        "lr": float(learning_rate), "dropout": float(dropout),
        "wd": float(weight_decay), "margin": float(margin),
        "smoothing": float(smoothing), "scaler": scaler, "gamma": 0.1,
        "beta": 0.0, "nu": float(nu), "thres": 0.0,
        "prune_threshold": 0.0, "warmup": int(warmup), "l1": 0.0,
        "reg_entropy": 0.0, "layer1": int(layer1),
        "n_layers": int(n_layers),
    }
    for idx in range(2, int(n_layers) + 1):
        params[f"layer{idx}"] = max(16, int(layer1) // (2 ** (idx - 1)))

    # Small deterministic resubstitution monitors satisfy legacy bookkeeping
    # without holding full-size valid/test copies. They remain part of train.
    monitor_n = min(len(frame), max(256, min(1024, int(batch_size))))
    monitor_idx = np.unique(np.linspace(0, len(frame) - 1, monitor_n, dtype=np.int64))
    monitor_frame = frame.iloc[monitor_idx]
    monitor_labels = labels[monitor_idx]
    monitor_batches = batches[monitor_idx]

    _bernn_holdout_module.scale_data = _memory_efficient_standard_scale_data
    try:
        trainer.fit(
            frame,
            labels,
            groups_train=batches,
            X_valid=monitor_frame,
            y_valid=monitor_labels,
            groups_valid=monitor_batches,
            X_test=monitor_frame,
            y_test=monitor_labels,
            groups_test=monitor_batches,
            params=params,
            internal_validation=False,
        )
    finally:
        _bernn_holdout_module.scale_data = _ORIGINAL_BERNN_SCALE_DATA

    trainer.restore_openproblems_final_state()

    # Encode in chunks so inference never creates another full dense copy plus
    # a second full torch tensor.
    encoded_parts = []
    chunk_size = max(2048, int(batch_size) * 4)
    for start in range(0, len(frame), chunk_size):
        stop = min(len(frame), start + chunk_size)
        encoded = trainer.transform(
            frame.iloc[start:stop],
            groups_test=batches[start:stop],
        )
        encoded_parts.append(np.asarray(encoded, dtype=np.float32))
    embedding = (
        np.concatenate(encoded_parts, axis=0)
        if encoded_parts
        else np.empty((0, 0), dtype=np.float32)
    )

    if return_trainer:
        return embedding, trainer
    return embedding


print("Read input", flush=True)
adata = read_anndata(
    par["input"],
    X="layers/normalized",
    obs="obs",
    var="var",
    uns="uns",
)

device = "cuda:0" if torch.cuda.is_available() else "cpu"
print(
    f"Run BERNN OpenProblems full-dataset fit on {adata.shape[0]} cells x {adata.shape[1]} genes "
    f"using {device}; no train/valid/test holdout and final epoch output",
    flush=True,
)

embedding, trainer = _fit_openproblems_memory_efficient_once(
    adata,
    layer="X",
    batch_key="batch",
    label_key="cell_type",
    n_hvg=int(par["n_hvg"]),
    dloss=str(par["dloss"]),
    n_epochs=int(par["n_epochs"]),
    warmup=int(par["warmup"]),
    early_stop=int(par["early_stop"]),
    early_warmup_stop=int(par["early_warmup_stop"]),
    batch_size=int(par["batch_size"]),
    n_layers=int(par["n_layers"]),
    layer1=int(par["layer1"]),
    scaler="standard",
    learning_rate=1e-3,
    weight_decay=1e-5,
    dropout=0.0,
    margin=1.0,
    smoothing=0.1,
    nu=1.0,
    num_workers=0,
    device=device,
    return_trainer=True,
)

print(
    f"Extracted final-epoch BERNN embedding after epoch {trainer._openproblems_last_epoch}",
    flush=True,
)
print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=pd.DataFrame(index=[]),
    obsm={"X_emb": np.asarray(embedding, dtype=np.float32)},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["name"],
    },
)
output.write_h5ad(par["output"], compression="gzip")
