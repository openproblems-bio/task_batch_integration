#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# remove this when you have implemented the script
echo "TODO: replace the commands in this script with the sequence of components that you need to run to generate test_resources."
echo "  Inside this script, you will need to place commands to generate example files for each of the 'src/api/file_*.yaml' files."
exit 1

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/task_template

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  --publish_dir "$DATASET_DIR" \
  --id "pancreas" \
  --input "$RAW_DATA/pancreas/dataset.h5ad" \
  --output_train '$id/train.h5ad' \
  --output_test '$id/test.h5ad' \
  --output_solution '$id/solution.h5ad' \
  --output_state '$id/state.yaml'

# run one method
viash run src/methods/knn/config.vsh.yaml -- \
    --input_train $DATASET_DIR/pancreas/train.h5ad \
    --input_test $DATASET_DIR/pancreas/test.h5ad \
    --output $DATASET_DIR/pancreas/prediction.h5ad

# run one metric
viash run src/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/pancreas/prediction.h5ad \
    --input_solution $DATASET_DIR/pancreas/solution.h5ad \
    --output $DATASET_DIR/pancreas/score.h5ad

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "$DATASET_DIR" s3://openproblems-data/resources_test/task_template \
  --delete --dryrun
