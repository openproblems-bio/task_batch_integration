#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/task_batch_integration

mkdir -p $DATASET_DIR

# process dataset
viash run src/process_dataset/config.vsh.yaml -- \
  --input "$RAW_DATA/cxg_mouse_pancreas_atlas/dataset.h5ad" \
  --output_dataset "$DATASET_DIR/cxg_mouse_pancreas_atlas/dataset.h5ad" \
  --output_solution "$DATASET_DIR/cxg_mouse_pancreas_atlas/solution.h5ad"

# run one method
viash run src/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/cxg_mouse_pancreas_atlas/dataset.h5ad \
  --output $DATASET_DIR/cxg_mouse_pancreas_atlas/integrated.h5ad

# run transformer
viash run src/transformers/transform/config.vsh.yaml -- \
    --input $DATASET_DIR/cxg_mouse_pancreas_atlas/integrated.h5ad \
    --output $DATASET_DIR/cxg_mouse_pancreas_atlas/integrated_full.h5ad

# run one metric
viash run src/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/cxg_mouse_pancreas_atlas/integrated.h5ad \
    --input_solution $DATASET_DIR/cxg_mouse_pancreas_atlas/solution.h5ad \
    --output $DATASET_DIR/cxg_mouse_pancreas_atlas/score.h5ad

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "resources_test/task_batch_integration" \
  s3://openproblems-data/resources_test/task_batch_integration \
  --delete --dryrun
