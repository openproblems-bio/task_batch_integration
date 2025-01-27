workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // transform output formats
    | transform.run(
      fromState: [
        input_integrated: "input_integrated",
        input_dataset: "input_dataset",
        expected_method_types: "expected_method_types"
      ],
      toState: { id, output, state ->
        def new_state = state + [
          method_output_cleaned: output.output,
        ]

        new_state
      }
    )

    // collect clustering resolutions
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

    // precompute clustering at one resolution
    | precompute_clustering_run.run(
      fromState: [
        input: "method_output_cleaned",
        resolution: "resolution"
      ],
      toState: ["output_clustering": "output"]
    )

    // group by original dataset id
    | map{id, state ->
      // groupKey() allows us to set a size for each group based on the state
      // which means each group can continue once it is complete
      [groupKey(state.prevId, state.resolutions.size()), state]
    }
    // Group and sort by resolution to ensure the order is consistent
    | groupTuple(sort: { res1, res2 -> res1.resolution <=> res2.resolution })

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

      [id.toString(), newState]
    }

    // merge clustering results into dataset h5ad
    | precompute_clustering_merge.run(
      fromState: [
        input: "method_output_cleaned",
        clusterings: "clusterings"
      ],
      toState: [output : "output"]
    )

    // only output what is defined in config
    | setState(["output"])

  emit:
  output_ch
}
