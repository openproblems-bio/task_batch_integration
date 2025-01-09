include { checkItemAllowed } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

// construct list of methods and control methods
methods = [
  // embed_cell_types,
  // embed_cell_types_jittered,
  no_integration,
  // no_integration_batch,
  // shuffle_integration,
  // shuffle_integration_by_batch,
  // shuffle_integration_by_cell_type,
  // batchelor_fastmnn,
  // batchelor_mnn_correct,
  // bbknn,
  combat,
  // geneformer,
  // harmony,
  // harmonypy,
  // liger,
  // mnnpy,
  // pyliger,
  // scalex,
  // scanorama,
  // scanvi,
  // scgpt,
  // scimilarity.run(
  //   args: [model: file("s3://openproblems-work/cache/scimilarity-model_v1.1.tar.gz")]
  // ),
  // scprint,
  // scvi,
  // uce.run(
  //   args: [model: file("s3://openproblems-work/cache/uce-model-v5.zip")]
  // )
]

// construct list of metrics
metrics = [
  asw_batch,
  // asw_label,
  // cell_cycle_conservation,
  clustering_overlap,
  // graph_connectivity,
  // hvg_overlap,
  // isolated_label_asw,
  // isolated_label_f1,
  // kbet,
  // lisi,
  // pcr
]

workflow run_wf {
  take:
  input_ch

  main:

  /****************************
   * EXTRACT DATASET METADATA *
   ****************************/
  dataset_ch = input_ch
    // store join id
    | map{ id, state ->
      [id, state + ["_meta": [join_id: id]]]
    }

    // extract the dataset metadata
    | extract_uns_metadata.run(
      fromState: [input: "input_solution"],
      toState: { id, output, state ->
        state + [
          dataset_uns: readYaml(output.output).uns
        ]
      }
    )

  /***************
   * RUN METHODS *
   ***************/

  score_ch = dataset_ch

    // run all methods
    | runEach(
      components: methods,

      // use the 'filter' argument to only run a method on the normalisation the component is asking for
      filter: { id, state, comp ->
        def norm = state.dataset_uns.normalization_id
        def pref = comp.config.info.preferred_normalization
        // if the preferred normalisation is none at all,
        // we can pass whichever dataset we want
        def norm_check = (norm == "log_cp10k" && pref == "counts") || norm == pref
        def method_check = checkItemAllowed(
          comp.config.name,
          state.methods_include,
          state.methods_exclude,
          "methods_include",
          "methods_exclude"
        )

        method_check && norm_check
      },

      // define a new 'id' by appending the method name to the dataset id
      id: { id, state, comp ->
        id + "." + comp.config.name
      },

      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: { id, state, comp ->
        def new_args = [:]
        if (comp.config.info.type == "method") {
          new_args.input = state.input_dataset
        } else if (comp.config.info.type == "control_method") {
          new_args.input_dataset = state.input_dataset
          new_args.input_solution = state.input_solution
        }
        new_args
      },

      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          method_id: comp.config.name,
          method_types: comp.config.info.method_types,
          method_output: output.output
        ]
      }
    )

  /******************
   * PROCESS OUTPUT *
   ******************/

    | process_integration.run(
      fromState: [
        input_integrated: "method_output",
        input_dataset: "input_dataset",
        expected_method_types: "method_types"
      ],
      toState: { id, output, state ->
        // Add method types to the state
        // This is done here because state can't be passed from the processing subworkflow
        def method_types_cleaned = []
        if ("feature" in state.method_types) {
          method_types_cleaned += ["feature", "embedding", "graph"]
        } else if ("embedding" in state.method_types) {
          method_types_cleaned += ["embedding", "graph"]
        } else if ("graph" in state.method_types) {
          method_types_cleaned += ["graph"]
        }

        def new_state = state + [
          method_output_processed: output.output,
          method_types_cleaned: method_types_cleaned
        ]

        new_state
      }
    )

  /***************
   * RUN METRICS *
   ***************/

    // run all metrics
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "." + comp.config.name
      },
      filter: { id, state, comp ->
        comp.config.info.metric_type in state.method_types_cleaned
      },
      // use 'fromState' to fetch the arguments the component requires from the overall state
      fromState: [
        input_solution: "input_solution",
        input_integrated: "method_output_processed"
      ],
      // use 'toState' to publish that component's outputs to the overall state
      toState: { id, output, state, comp ->
        state + [
          metric_id: comp.config.name,
          metric_output: output.output
        ]
      }
    )

    // extract the scores
    | extract_uns_metadata.run(
      key: "extract_scores",
      fromState: [input: "metric_output"],
      toState: { id, output, state ->
        state + [
          score_uns: readYaml(output.output).uns
        ]
      }
    )

    | joinStates { ids, states ->
      // store the scores in a file
      def score_uns = states.collect{it.score_uns}
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)

      ["output", [output_scores: score_uns_file]]
    }

  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  // TODO: can we store everything below in a separate helper function?

  // extract the dataset metadata
  meta_ch = dataset_ch
    // only keep one of the normalization methods
    | filter{ id, state ->
      state.dataset_uns.normalization_id == "log_cp10k"
    }
    | joinStates { ids, states ->
      // store the dataset metadata in a file
      def dataset_uns = states.collect{state ->
        def uns = state.dataset_uns.clone()
        uns.remove("normalization_id")
        uns
      }
      def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(dataset_uns_yaml_blob)

      // store the method configs in a file
      def method_configs = methods.collect{it.config}
      def method_configs_yaml_blob = toYamlBlob(method_configs)
      def method_configs_file = tempFile("method_configs.yaml")
      method_configs_file.write(method_configs_yaml_blob)

      // store the metric configs in a file
      def metric_configs = metrics.collect{it.config}
      def metric_configs_yaml_blob = toYamlBlob(metric_configs)
      def metric_configs_file = tempFile("metric_configs.yaml")
      metric_configs_file.write(metric_configs_yaml_blob)

      // store the task info in a file
      def viash_file = meta.resources_dir.resolve("_viash.yaml")

      // create output state
      def new_state = [
        output_dataset_info: dataset_uns_file,
        output_method_configs: method_configs_file,
        output_metric_configs: metric_configs_file,
        output_task_info: viash_file,
        _meta: states[0]._meta
      ]

      ["output", new_state]
    }

  // merge all of the output data
  output_ch = score_ch
    | mix(meta_ch)
    | joinStates{ ids, states ->
      def mergedStates = states.inject([:]) { acc, m -> acc + m }
      [ids[0], mergedStates]
    }

  emit:
  output_ch
}
