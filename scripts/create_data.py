import warnings
warnings.filterwarnings("ignore", category=FutureWarning, message=r".*Reordering categories will always return a new Categorical object.*")
warnings.filterwarnings("ignore", category=FutureWarning, message=r".*is_categorical_dtype is deprecated and will be removed in a future version.*")

import numpy as np
import pandas as pd
import anndata as ad

import scipy as sp
from scipy import linalg
from sklearn.datasets import make_sparse_spd_matrix
#import matplotlib.pyplot as plt

import atacnet
from atacnet.atacnet import get_distances_regions

print("load params")
# parameters for fake data
nb_cells = snakemake.params['nb_cells']
nb_chr = snakemake.params['nb_chr']
nb_regions_per_chr = snakemake.params['nb_regions_per_chr']
between_reg = snakemake.params['between_reg']
size_reg = snakemake.params['size_reg']
sep = snakemake.params['sep']
distance_threshold = snakemake.params['distance_threshold']

print("Initiate data")
# Inititate structure of fake single-cell atac-seq data 
counts = []
for chr in range(nb_chr):
    counts.append(pd.DataFrame(np.zeros((nb_cells, nb_regions_per_chr)),
                        index=['Cell_'+j for j in map(str, range(nb_cells))],
                        columns=['chr'+str(chr)+sep+str(i)+sep+str(i+size_reg)
                                 for i in
                                 range(1, nb_regions_per_chr*between_reg+1,
                                       between_reg)])
                )
# Create AnnData object
atac = ad.AnnData(pd.concat(counts, axis=1))

# Add region infos as per-variables columns
atacnet.add_region_infos(atac)

# Generate fake correlations between regions
n_samples, n_features = nb_cells, nb_regions_per_chr
prng = np.random.RandomState(1)
prec = make_sparse_spd_matrix(
    n_features, alpha=0.99, smallest_coef=0.4, largest_coef=0.7, random_state=prng
)
cov = linalg.inv(prec)

# cov with only potential connections
possible_co = sp.sparse.csr_matrix(get_distances_regions(atac)<distance_threshold/2)[:cov.shape[0],:cov.shape[1]]
possible_co = sp.sparse.coo_matrix(possible_co).toarray() + sp.sparse.coo_matrix(possible_co).toarray().T 
cov = np.eye(len(cov))*np.diag(cov) + possible_co*cov 
d = np.sqrt(np.diag(cov))
cov /= d
cov /= d[:, np.newaxis]
prec *= d
prec *= d[:, np.newaxis]
X = prng.multivariate_normal(np.zeros(n_features), cov, size=n_samples)
X -= X.mean(axis=0)
X /= X.std(axis=0)

X_ = np.concatenate([X]*nb_chr, axis=1)
atac.X = X_

# Save data as h5ad
atac.write_h5ad(snakemake.output[0])
# Save data as csv.gz
df_atac = pd.DataFrame(atac.X, index=atac.obs_names, columns=atac.var_names)
df_atac.to_csv(snakemake.output[1], sep='\t')
print(atac)
cov_df = pd.DataFrame(cov)
cov_df.index = atac.var_names[atac.var['chromosome'] == 'chr0']
cov_df.columns = atac.var_names[atac.var['chromosome'] == 'chr0']
cov_df.replace(0, np.nan, inplace=True)
cov_df = cov_df.unstack()
cov_df.to_csv(snakemake.output[2], sep='\t')
