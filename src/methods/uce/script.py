import anndata as ad

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input': 'resources_test/.../input.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'name': 'my_python_method'
}
## VIASH END

print('Reading input files', flush=True)
input = ad.read_h5ad(par['input'])

print('Preprocess data', flush=True)
# ... preprocessing ...

print('Train model', flush=True)
# ... train model ...

print('Generate predictions', flush=True)
# ... generate predictions ...

print("Write output AnnData to file", flush=True)
output = ad.AnnData(

)
output.write_h5ad(par['output'], compression='gzip')
