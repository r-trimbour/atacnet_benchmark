import pandas as pd
import anndata as ad
import scanpy as sc
import atacnet as an
import atacnet.metacells

def extract_links(
    adata,  # AnnData object
    key=None,  # key from adata.varp
    columns=['row', 'col', 'weight']  # output col names (e.g. 'TF', 'gene', 'score')
    ):
    if key is None:  # if only one key (I guess often), no need to precise key
                     # maybe replace by a default one later if one is prvided in AnnData package
        if len(list(adata.varp)) == 1:
            key = list(adata.varp)[0]
        else:
            raise "Several keys were found in adata.varp: {}, please precise which keyword use (arg 'key')".format(list(adata.varp))

    return pd.DataFrame(
        [a for a in zip(
            [adata.var_names[i] for i in adata.varp[key].row],
            [adata.var_names[i] for i in adata.varp[key].col],
            adata.varp[key].data)],
        columns=columns
        ).sort_values(by=columns[2], ascending=False)

# Parameters
distance_threshold = snakemake.params['distance_threshold']

# Load data
atac = ad.read_csv(snakemake.input["atac"], delimiter='\t')
print(atac.var_names)


# Add region infos
an.add_region_infos(atac, sep=('_', '_'))
print(atac)
print(atac.X)
print(atac.X.sum(0).min(), atac.X.sum(0).min())

#atac = an.metacells.compute_metacells(atac,
 #                              k = snakemake.params["number_cells_per_clusters"])


atac = ad.read_csv(snakemake.input["cicero_cds"], delimiter=' ').transpose()
print(atac)
print(atac.X)
print(atac.X.sum(0).min(), atac.X.sum(0).min())
# Add region infos
an.add_region_infos(atac, sep=('-', '-'))


# Compute network
an.compute_atac_network(
    atac,
    window_size=distance_threshold,
    unit_distance=1000,
    distance_constraint=distance_threshold/2,
    n_samples=100,
    n_samples_maxtry=100,
    max_alpha_iteration=100,
)
#    n_jobs=1
#)

# Save results
peak_layer = extract_links(
    atac,
    columns=['peak1', 'peak2', 'score'])
peak_layer.to_csv(snakemake.output[0], sep='\t', index=False)
