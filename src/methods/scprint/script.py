import anndata as ad
from scdataloader import Preprocessor
import sys

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input': 'resources_test/.../input.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'name': 'scprint'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata

print(">>> Reading input data...", flush=True)
input = read_anndata(par["input"], X="layers/counts", obs="obs", var="var", uns="uns")

print(">>> Setting ontology term IDs...", flush=True)
# For now, set all ontology term IDs to 'unknown' but these could be used for
# cellxgene datasets that have this information
print("NOTE: All ontology term IDs except organism are set to 'unknown'", flush=True)
if input.uns["dataset_organism"] == "homo_sapiens":
    input.obs["organism_ontology_term_id"] = "NCBITaxon:9606"
elif input.uns["dataset_organism"] == "mus_musculus":
    input.obs["organism_ontology_term_id"] = "NCBITaxon:10090"
else:
    raise ValueError(f"scPRINT requires human or mouse data, not '{input.uns['dataset_organism']}'")
input.obs["self_reported_ethnicity_ontology_term_id"] = "unknown"
input.obs["disease_ontology_term_id"] = "unknown"
input.obs["cell_type_ontology_term_id"] = "unknown"
input.obs["development_stage_ontology_term_id"] = "unknown"
input.obs["tissue_ontology_term_id"] = "unknown"
input.obs["assay_ontology_term_id"] = "unknown"
input.obs["sex_ontology_term_id"] = "unknown"

print('\n>>> Preprocessing data...', flush=True)
preprocessor = Preprocessor(
    # Lower this threshold for test datasets
    min_valid_genes_id = 1000 if input.n_vars < 2000 else 10000,
    # Turn off cell filtering to return results for all cells
    filter_cell_by_counts = False,
    min_nnz_genes = False,
    do_postp=False
)
processed = preprocessor(input)

print('Train model', flush=True)
# ... train model ...

print('Generate predictions', flush=True)
# ... generate predictions ...

print("Write output AnnData to file", flush=True)
output = ad.AnnData(

)
output.write_h5ad(par['output'], compression='gzip')
