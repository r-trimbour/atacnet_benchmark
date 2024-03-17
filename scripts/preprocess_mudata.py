import mudata as md
import scanpy as sc

mudata = md.read_h5mu(snakemake.input['mudata'])
print(mudata)

atac = mudata['atac']

sc.pp.subsample(atac, n_obs=500, random_state=0)
atac = atac[:, :10000]
atac = atac[atac.X.sum(1) > 0, :]
atac = atac[:, atac.X.sum(0) > 0]

atac.var_names = atac.var_names.str.replace(':', '_', regex=False)
atac.var_names = atac.var_names.str.replace('-', '_', regex=False)

print(atac.var_names)
print(atac.obs_names)
print(atac)
print(atac.X.sum(0).min(), atac.X.sum(0).min())

atac.write_h5ad(snakemake.output['anndata'])

atac.to_df().to_csv(snakemake.output['csv'], sep='\t')
