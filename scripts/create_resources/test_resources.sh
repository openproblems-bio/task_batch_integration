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
nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  --input "$RAW_DATA/cxg_immune_cell_atlas/dataset.h5ad" \
  --publish_dir "$DATASET_DIR/cxg_immune_cell_atlas" \
  --output_dataset dataset.h5ad \
  --output_solution solution.h5ad \
  --output_state state.yaml \
  -c common/nextflow_helpers/labels_ci.config \

# run one method
viash run src/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/cxg_immune_cell_atlas/dataset.h5ad \
  --output $DATASET_DIR/cxg_immune_cell_atlas/integrated.h5ad

# process integration
nextflow run . \
  -main-script target/nextflow/data_processors/process_integration/main.nf \
  -profile docker \
  --input_dataset "$RAW_DATA/cxg_immune_cell_atlas/dataset.h5ad" \
  --input_integrated "$DATASET_DIR/cxg_immune_cell_atlas/integrated.h5ad" \
  --expected_method_types feature \
  --publish_dir "$DATASET_DIR/cxg_immune_cell_atlas" \
  --output integrated_processed.h5ad \
  -c common/nextflow_helpers/labels_ci.config \

# run one metric
viash run src/metrics/graph_connectivity/config.vsh.yaml -- \
    --input_integrated $DATASET_DIR/cxg_immune_cell_atlas/integrated_processed.h5ad \
    --input_solution $DATASET_DIR/cxg_immune_cell_atlas/solution.h5ad \
    --output $DATASET_DIR/cxg_immune_cell_atlas/score.h5ad

# write the state file
cat > $DATASET_DIR/cxg_immune_cell_atlas/state.yaml << HERE
id: cxg_immune_cell_atlas
output_dataset: !file dataset.h5ad
output_solution: !file solution.h5ad
output_integrated: !file integrated.h5ad
output_integrated_full: !file integrated_full.h5ad
output_score: !file score_mod1.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "resources_test/task_batch_integration" \
  s3://openproblems-data/resources_test/task_batch_integration \
  --delete --dryrun
