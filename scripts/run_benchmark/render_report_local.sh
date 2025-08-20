#!/bin/bash

# fail on error
set -e

# ensure we're in the root of the repo
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# set input and output directories
TASK=task_batch_integration
BASE_DIR="s3://openproblems-data/resources/$TASK/results"
OUTPUT_DIR="output/report"

# find subdir in bucket with latest date
DATE=$(aws s3 ls $BASE_DIR --recursive | awk '{print $4}' | grep 'task_info.yaml' | sort -r | head -n 1 | sed 's#.*/run_\(.*\)/[^/]*$#\1#')

INPUT_DIR="$BASE_DIR/run_$DATE"
TASK_STRIP_PREFIX=$(echo $TASK | sed 's/task_//')

echo "Processing $DATE -> $OUTPUT_DIR"

# start the run
extra_filters=()
# extra_filters=(
#   --datasets_exclude "cellxgene_census/hypomap;cellxgene_census/mouse_pancreas_atlas"
#   --metrics_exclude "hvg_overlap"
# )

nextflow run openproblems-bio/openproblems \
  -r build/main \
  -main-script target/nextflow/reporting/process_task_results/main.nf \
  -profile docker \
  -resume \
  -latest \
  -c common/nextflow_helpers/labels_ci.config \
  --id "$TASK/run_$DATE" \
  --input_scores "$INPUT_DIR/score_uns.yaml" \
  --input_dataset_info "$INPUT_DIR/dataset_uns.yaml" \
  --input_method_configs "$INPUT_DIR/method_configs.yaml" \
  --input_metric_configs "$INPUT_DIR/metric_configs.yaml" \
  --input_trace "$INPUT_DIR/trace.txt" \
  --input_task_info "$INPUT_DIR/task_info.yaml" \
  --output_state '$id/state.yaml' \
  --output_combined '$id/combined_output.json' \
  --output_report '$id/report.html' \
  --output_dataset_info '$id/dataset_info.json' \
  --output_method_info '$id/method_info.json' \
  --output_metric_info '$id/metric_info.json' \
  --output_results '$id/results.json' \
  --output_quality_control '$id/quality_control.json' \
  --publish_dir "$OUTPUT_DIR" \
  "${extra_filters[@]}"

