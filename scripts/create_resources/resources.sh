#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# remove this when you have implemented the script
echo "TODO: once the 'process_datasets' workflow is implemented, update this script to use it."
echo "  Step 1: replace 'task_batch_integration' with the name of the task in the following command."
echo "  Step 2: replace the rename keys parameters to fit your process_dataset inputs"
echo "  Step 3: replace the settings parameter to fit your process_dataset outputs"
echo "  Step 4: remove this message"
exit 1

cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources/datasets/**/state.yaml
rename_keys: 'input:output_dataset'
output_state: '$id/state.yaml'
settings: '{"output_train": "$id/train.h5ad", "output_test": "$id/test.h5ad"}'
publish_dir: s3://openproblems-data/resources/task_batch_integration/datasets/
HERE

tw launch https://github.com/openproblems-bio/task_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_batch_integration,process_datasets
