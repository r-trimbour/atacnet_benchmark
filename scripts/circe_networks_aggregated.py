import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import circe as ci
import circe.metacells
from rich.progress import track
import scipy as sp


# Parameters
distance_threshold = snakemake.params['distance_threshold']

atac = ad.read_csv(snakemake.input["cicero_cds"], delimiter=' ').transpose()
print(atac)
# Add region infos
atac = ci.add_region_infos(atac, sep=('-', '-'))


# Compute network
ci.compute_atac_network(
    atac,
    window_size=distance_threshold,
    unit_distance=1000,
    distance_constraint=distance_threshold/2,
    n_samples=100,
    n_samples_maxtry=500,
    max_alpha_iteration=100,
)

# Save results
peak_layer = ci.extract_atac_links(
    atac,
    columns=['peak1', 'peak2', 'score'])
peak_layer.to_csv(snakemake.output[0], sep='\t', index=False)
