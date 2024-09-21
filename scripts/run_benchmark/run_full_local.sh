#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# NOTE: depending on the the datasets and components, you may need to launch this workflow
# on a different compute platform (e.g. a HPC, AWS Cloud, Azure Cloud, Google Cloud).
# please refer to the nextflow information for more details:
# https://www.nextflow.io/docs/latest/

# remove this when you have implemented the script
echo "TODO: once the 'run_benchmark' workflow has been implemented, update this script to use it."
echo "  Step 1: replace 'task_batch_integration' with the name of the task in the following command."
echo "  Step 2: replace the rename keys parameters to fit your run_benchmark inputs"
echo "  Step 3: replace the settings parameter to fit your run_benchmark outputs"
echo "  Step 4: remove this message"
exit 1

set -e

echo "Running benchmark on test data"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="resources/results/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
input_states: resources/datasets/**/state.yaml
rename_keys: 'input_train:output_train;input_test:output_test;input_solution:output_solution'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# run the benchmark
nextflow run openproblems-bio/task_batch_integration \
  --revision build/main \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file /tmp/params.yaml
