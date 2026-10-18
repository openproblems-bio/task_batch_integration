"""Exercise production parameter wiring without importing scGPT or PyTorch."""

import ast
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock

import yaml


COMPONENT = Path(__file__).resolve().parent
SCRIPT = COMPONENT / "script.py"
TREE = ast.parse(SCRIPT.read_text(encoding="utf-8"))
CONFIG = yaml.safe_load((COMPONENT / "config.vsh.yaml").read_text(encoding="utf-8"))


def assignment(name):
    return next(
        node
        for node in TREE.body
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == name for t in node.targets)
    )


def execute(nodes, namespace):
    # Compile the original nodes unchanged, retaining production line numbers.
    exec(
        compile(ast.Module(body=nodes, type_ignores=[]), str(SCRIPT), "exec"), namespace
    )


def run_training(par):
    """Run the actual epoch loop and embedding call with only ML work stubbed."""
    namespace = dict.fromkeys(
        [
            "tokenized_train",
            "tokenized_valid",
            "train_batch_labels",
            "valid_batch_labels",
            "model",
            "scaler",
            "optimizer",
            "criterion",
            "criterion_dab",
            "device",
            "adata",
        ]
    )
    namespace.update(
        par=par,
        vocab=["gene"],
        model_settings={"per_seq_batch_sample": True},
        prepare_data=Mock(return_value=("train_data", "valid_data")),
        prepare_dataloader=Mock(side_effect=lambda data, **kwargs: (data, kwargs)),
        train=Mock(),
        evaluate=Mock(return_value=(1.0, 0.1)),
        scheduler=Mock(),
        time=SimpleNamespace(time=lambda: 0),
        copy=SimpleNamespace(deepcopy=lambda model: model),
        scgpt=SimpleNamespace(logger=Mock(), tasks=SimpleNamespace(embed_data=Mock())),
        best_val_loss=float("inf"),
        best_model_dir=SimpleNamespace(name="unused_model_path"),
    )
    epoch_loop = next(
        node
        for node in TREE.body
        if isinstance(node, ast.For)
        and isinstance(node.target, ast.Name)
        and node.target.id == "epoch"
    )
    execute(
        [assignment("hyperparameters"), epoch_loop, assignment("embedded")], namespace
    )
    return namespace


class TrainingParametersTest(unittest.TestCase):
    def test_config_exposes_integer_defaults(self):
        arguments = {arg["name"]: arg for arg in CONFIG["arguments"]}
        for name, default in (("epochs", 15), ("batch_size", 64)):
            with self.subTest(parameter=name):
                self.assertIn("--" + name, arguments)
                argument = arguments["--" + name]
                self.assertEqual(argument["type"], "integer")
                self.assertEqual(argument["default"], default)
                self.assertEqual(argument["min"], 1)
                self.assertEqual(
                    argument["info"]["test_default"],
                    2 if name == "epochs" else 4,
                )

    def test_committed_par_matches_config_defaults(self):
        namespace = {}
        execute([assignment("par")], namespace)
        for name, default in (("epochs", 15), ("batch_size", 64)):
            with self.subTest(parameter=name):
                self.assertEqual(namespace["par"].get(name), default)

    def check_training(self, epochs, batch_size):
        namespace = run_training({"epochs": epochs, "batch_size": batch_size})
        self.assertEqual(namespace["train"].call_count, epochs)
        self.assertEqual(namespace["evaluate"].call_count, epochs)
        self.assertEqual(namespace["scheduler"].step.call_count, epochs)
        for epoch, call in enumerate(namespace["train"].call_args_list, start=1):
            self.assertEqual(call.args[-1], epoch)
            self.assertEqual(call.args[1][0], "train_data")
            self.assertEqual(call.args[1][1]["batch_size"], batch_size)
            self.assertEqual(call.args[8]["epochs"], epochs)
            self.assertEqual(call.args[8]["batch_size"], batch_size)
        for call in namespace["evaluate"].call_args_list:
            self.assertEqual(call.args[1][0], "valid_data")
            self.assertEqual(call.args[1][1]["batch_size"], batch_size)

    def test_default_training(self):
        self.check_training(15, 64)

    def test_overridden_training(self):
        self.check_training(2, 8)

    def test_overridden_epochs(self):
        self.check_training(2, 64)

    def test_overridden_batch_size(self):
        self.check_training(15, 8)

    def test_embedding_batch_size(self):
        for batch_size in (64, 8):
            with self.subTest(batch_size=batch_size):
                namespace = run_training({"epochs": 2, "batch_size": batch_size})
                kwargs = namespace["scgpt"].tasks.embed_data.call_args.kwargs
                self.assertEqual(kwargs["batch_size"], batch_size)
                self.assertFalse(kwargs["use_fast_transformer"])

    def test_other_hyperparameters_unchanged(self):
        namespace = run_training({"epochs": 2, "batch_size": 8})
        hyperparameters = namespace["hyperparameters"].copy()
        hyperparameters.pop("epochs")
        hyperparameters.pop("batch_size")
        self.assertEqual(
            hyperparameters,
            {
                "n_tokens": 1,
                "GEPC": True,
                "ecs_thres": 0.8,
                "dab_weight": 1.0,
                "mask_ratio": 0.4,
                "lr": 1e-4,
                "dropout": 0.2,
                "schedule_ratio": 0.9,
                "log_interval": 100,
                "fast_transformer": False,
                "pre_norm": False,
                "amp": True,
            },
        )


if __name__ == "__main__":
    unittest.main()
