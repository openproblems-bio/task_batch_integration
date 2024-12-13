include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    | check_dataset_with_schema.run(
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "input")
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.input,
          "schema": schemaYaml
        ]
      },
      toState: { id, output, state ->
        // read the output to see if dataset passed the qc
        def checks = readYaml(output.output)
        state + [
          "dataset": checks["exit_code"] == 0 ? state.input : null,
        ]
      }
    )

    // remove datasets which didn't pass the schema check
    | filter { id, state ->
      state.dataset != null
    }

    // precompute clustering of the input data at various resolutions
    | flatMap { id, state ->
      state.resolutions.collect { resolution ->
        def newId = "${id}_r${resolution}"
        def newState = state + [
          "resolution": resolution,
          "prevId": id
        ]
        [newId, newState]
      }
    }
    
    // precompute clustering at this resolution
    | precompute_clustering_run.run(
      fromState: ["input": "dataset", "resolution": "resolution"],
      toState: ["output_clustering": "output"]
    )

    // group by original dataset id
    | map{id, state ->
      [state.prevId, state]
    }
    | groupTuple()

    // merge the clustering results into one state
    | map{ id, states -> 
      if (states.size() == 0) {
        throw new RuntimeException("Expected at least one state, but got ${states.size()}")
      }
      if (states.size() != states[0].resolutions.size()) {
        throw new RuntimeException("Expected ${states[0].resolutions.size()} states, but got ${states.size()}")
      }

      def clusterings = states.collect { it.output_clustering }
      def newState = states[0] + ["clusterings": clusterings]
      
      [id, newState]
    }

    // merge clustering results into dataset h5ad
    | precompute_clustering_merge.run(
      fromState: ["input": "dataset", "clusterings": "clusterings"],
      toState: ["dataset": "output"]
    )

    // process the dataset
    | process_dataset.run(
      fromState: [ input: "dataset" ],
      toState: [
        output_dataset: "output_dataset",
        output_solution: "output_solution"
      ]
    )

    // only output the files for which an output file was specified
    | setState(["output_dataset", "output_solution"])

  emit:
  output_ch
}
