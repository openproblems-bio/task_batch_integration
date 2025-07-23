#!/bin/bash

# template for adding new datasets
#   - id: cellxgene_census/
#     species:
#     census_version: "2023-07-25"
#     obs_value_filter: "dataset_id == ''"
#     obs_batch:
#     dataset_name:
#     dataset_summary:
#     dataset_description:
#     dataset_url:
#     dataset_reference:
#     dataset_organism:

# not sure which dataset ids to use
#   - id: cellxgene_census/human_brain_atlas
#     species: homo_sapiens
#     census_version: "2023-07-25"
#     obs_value_filter: "dataset_id == ''" # <--- ?
#     obs_batch: donor_id
#     dataset_name:  Human Brain Atlas
#     dataset_summary: Single-Cell DNA Methylation and 3D Genome Human Brain Atlas
#     dataset_description: Delineating the gene regulatory programs underlying complex cell types is fundamental for understanding brain functions in health and disease. Here, we comprehensively examine human brain cell epigenomes by probing DNA methylation and chromatin conformation at single-cell resolution in over 500,000 cells from 46 brain regions. We identified 188 cell types and characterized their molecular signatures. Integrative analyses revealed concordant changes in DNA methylation, chromatin accessibility, chromatin organization, and gene expression across cell types, cortical areas, and basal ganglia structures. With these resources, we developed scMCodes that reliably predict brain cell types using their methylation status at select genomic sites. This multimodal epigenomic brain cell atlas provides new insights into the complexity of cell type-specific gene regulation in the adult human brain.
#     dataset_url: https://cellxgene.cziscience.com/collections/fdebfda9-bb9a-4b4b-97e5-651097ea07b0
#     dataset_reference: tian2023singlecell
#     dataset_organism: homo_sapiens

cat > "/tmp/params.yaml" << 'HERE'
param_list:
 - id: cellxgene_census/hnoca
    species: homo_sapiens
    census_version: "2024-11-20"
    obs_value_filter: "dataset_id == '3a805b8d-c9d8-4c9b-881f-0eefa635974d'"
    obs_batch: donor_id
    dataset_name: Human Neural Organoids Cell Atlas
    dataset_summary: An integrated transcriptomic cell atlas of human neural organoids
    dataset_description: PLEASE NOTE: 1. The metadata field `cell_type` corresponds to a manual mapping of the original author annotations (metadata field `cell_type_original`) to the Cell Ontology. For the harmonised cell type, region, and neurotransmitter-transporter annotations, please refer to the metadata fields starting with `annot_` in the Author Categories. 2. For the HNOCA extended, you can find the harmonised cell type annotation covering all cells (including the extension datasets) in the `annot_level_2_extended` metadata field. 3. The metadata field `tissue` corresponds to the target tissue of the employed organoid differentiation protocol. For example, a cell originating from a sample generated using a cortical differentiation protocol will be annotated as `cerebral cortex`. 4. The data deposited here contains slightly fewer cells than in the data associated with the original publication. This is due to the removal of some cells with identical expression values as required by the CellxGene schema. You can find the full object at the Zenodo Data Source link to the right. PUBLICATION ABSTRACT: Neural tissues generated from human pluripotent stem cells in vitro (known as neural organoids) are becoming useful tools to study human brain development, evolution and disease. The characterization of neural organoids using single-cell genomic methods has revealed a large diversity of neural cell types with molecular signatures similar to those observed in primary human brain tissue. However, it is unclear which domains of the human nervous system are covered by existing protocols. It is also difficult to quantitatively assess variation between protocols and the specific cell states in organoids as compared to primary counterparts. Single-cell transcriptome data from primary tissue and neural organoids derived with guided or un-guided approaches and under diverse conditions combined with large-scale integrative analyses make it now possible to address these challenges. Recent advances in computational methodology enable the generation of integrated atlases across many data sets. Here, we integrated 36 single-cell transcriptomics data sets spanning 26 protocols into one integrated human neural organoid cell atlas (HNOCA) totaling over 1.7 million cells. We harmonize cell type annotations by incorporating reference data sets from the developing human brain. By mapping to the developing human brain reference, we reveal which primary cell states have been generated in vitro, and which are under-represented. We further compare transcriptomic profiles of neuronal populations in organoids to their counterparts in the developing human brain. To support rapid organoid phenotyping and quantitative assessment of new protocols, we provide a programmatic interface to browse the atlas and query new data sets, and showcase the power of the atlas to annotate new query data sets and evaluate new organoid protocols. Taken together, the HNOCA will be useful to assess the fidelity of organoids, characterize perturbed and diseased states and facilitate protocol development in the future.
    dataset_url: https://cellxgene.cziscience.com/collections/de379e5f-52d0-498c-9801-0f850823c847
    dataset_reference: He2024integrated
    dataset_organism: homo_sapiens

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: s3://openproblems-data/resources/datasets/task_batch_integration
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withLabel: highmem {
    memory = '350GB'
  }
  withName: '.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/scrnaseq/process_cellxgene_census/main.nf \
  --workspace 53907369739130 \  <-- might need to be changed
  --compute-env 7gRyww9YNGb0c6BUBtLhDP \  <-- might need to be changed
  --params-file "/tmp/params.yaml" \  <-- should be adjusted to whatever files we are interested in
  --config /tmp/nextflow.config \ <-- might need to be changed according to computational requirements
  --labels cellxgene_census,dataset_loader