import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import atacnet as an
import atacnet.metacells
from rich.progress import track
import scipy as sp

def compute_atac_network(
    AnnData,
    window_size=500000,
    unit_distance=1000,
    distance_constraint=250000,
    s=0.75,
    max_alpha_iteration=100,
    distance_parameter_convergence=1e-22,
    max_elements=200,
    n_samples=100,
    n_samples_maxtry=500
):
    """
    Compute co-accessibility scores between regions in a sparse matrix, stored
    in the varp slot of the passed AnnData object.
    Scores are computed using 'sliding_graphical_lasso'.

    1. First, the function calculates the optimal penalty coefficient alpha.
        Alpha is calculated by averaging alpha values from 'n_samples' windows,
        such as there's less than 5% of possible long-range edges
        (> distance_constraint) and less than 20% co-accessible regions
        (regardless of distance constraint) in each window.

    2. Then, it will calculate co-accessibility scores between regions in a
    sliding window of size 'window_size' and step 'window_size/2'.
        Results should be very similar to Cicero's results. There is a strong
        correlation between Cicero's co-accessibility scores and the ones
        calculated by this function. However, absolute values are not the same,
        because Cicero uses a different method to apply Graphical Lasso.

    3. Finally, it will average co-accessibility scores across windows.

    Parameters
    ----------
    AnnData : AnnData object
        AnnData object with var_names as region names.
    window_size : int, optional
        Size of sliding window, in which co-accessible regions can be found.
        The default is 500000.
    unit_distance : int, optional
        Distance between two regions in the matrix, in base pairs.
        The default is 1000.
    distance_constraint : int, optional
        Distance threshold for defining long-range edges.
        It is used to fit the penalty coefficient alpha.
        The default is 250000.
    s : float, optional
        Parameter for penalizing long-range edges. The default is 0.75 and
        should probably not be changed, unless you know what you are doing.
    max_alpha_iteration : int, optional
        Maximum number of iterations to calculate optimal penalty coefficient.
        The default is 100.
    distance_parameter_convergence : float, optional
        Convergence parameter for alpha (penalty) coefficiant calculation.
        The default is 1e-22.
    max_elements : int, optional
        Maximum number of regions in a window. The default is 200.
    n_samples : int, optional
        Number of windows used to calculate optimal penalty coefficient alpha.
        The default is 100.
    n_samples_maxtry : int, optional
        Maximum number of windows to try to calculate optimal penalty
        coefficient alpha. Should be higher than n_samples. The default is 500.

    Returns
    -------
    None.

    """

    AnnData.varp['atac_network'] = sliding_graphical_lasso(
        AnnData=AnnData,
        window_size=window_size,
        unit_distance=unit_distance,
        distance_constraint=distance_constraint,
        s=s,
        max_alpha_iteration=max_alpha_iteration,
        distance_parameter_convergence=distance_parameter_convergence,
        max_elements=max_elements,
        n_samples=n_samples,
        n_samples_maxtry=n_samples_maxtry
    )

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

def sliding_graphical_lasso(
    AnnData,
    window_size: int = 500_000,
    unit_distance=1_000,
    distance_constraint=250_000,
    s=0.75,
    max_alpha_iteration=100,
    distance_parameter_convergence=1e-22,
    max_elements=200,
    n_samples=500,
    n_samples_maxtry=00,
):
    """
    Extract sliding submatrix from a sparse correlation matrix.

    WARNING: might look generalised for many overlaps but is not yet at the
    end, that's why 'start_sliding' is hard coded as list of 2 values.

    Parameters
    ----------
    AnnData : AnnData object
        AnnData object with var_names as region names.
    window_size : int, optional
        Size of the sliding window, where co-accessible regions can be found.
        The default is 500000.
    unit_distance : int, optional
        Distance between two regions in the matrix, in base pairs.
        The default is 1000.
    distance_constraint : int, optional
        Distance threshold for defining long-range edges.
        It is used to fit the penalty coefficient alpha.
        The default is 250000.
    s : float, optional
        Parameter for penalizing long-range edges. The default is 0.75 and
        should probably not be changed, unless you know what you are doing.
    max_alpha_iteration : int, optional
        Maximum number of iterations to calculate optimal penalty coefficient.
        The default is 100.
    distance_parameter_convergence : float, optional
        Convergence parameter for alpha (penalty) coefficiant calculation.
        The default is 1e-22.
    max_elements : int, optional
        Maximum number of regions in a window. The default is 200.
    n_samples : int, optional
        Number of windows used to calculate optimal penalty coefficient alpha.
        The default is 100.
    n_samples_maxtry : int, optional
        Maximum number of windows to try to calculate optimal penalty
        coefficient alpha. Should be higher than n_samples. The default is 500.
    """

    # print("Calculating penalty coefficient alpha...")
    alpha = an.atacnet.average_alpha(
        AnnData,
        window_size=window_size,
        unit_distance=unit_distance,
        n_samples=n_samples,
        n_samples_maxtry=n_samples_maxtry,
        max_alpha_iteration=max_alpha_iteration,
        s=s,
        distance_constraint=distance_constraint,
        distance_parameter_convergence=distance_parameter_convergence,
        max_elements=max_elements,
    )


    start_slidings = [0, int(window_size / 2)]

    results = {}
    regions_list = AnnData.var_names
    # Get global indices of regions
    map_indices = {regions_list[i]: i for i in range(len(regions_list))}

    for k in start_slidings:
        slide_results = {}
        slide_results["scores"] = np.array([])
        slide_results["idx"] = np.array([])
        slide_results["idy"] = np.array([])
        if k == 0:
            print("Starting to process chromosomes : {}".format(
                AnnData.var["chromosome"].unique()))
        else:
            print("Finishing to process chromosomes : {}".format(
                AnnData.var["chromosome"].unique()))
        for chromosome in track(
            AnnData.var["chromosome"].unique(),
            description="Calculating co-accessibility: {}/2".format(
                1 if k == 0 else 2),):
            # Get start positions of windows
            window_starts = [
                i
                for i in range(
                    k,
                    AnnData.var["end"][
                        AnnData.var["chromosome"] == chromosome].max(),
                    window_size,
                )
            ]

            for start in window_starts:
                end = start + window_size
                # Get global indices of regions in the window
                idx = np.where(
                    ((AnnData.var["chromosome"] == chromosome)
                     & (AnnData.var["start"] >= start)
                     & (AnnData.var["start"] <= end)))[0]

                # already global ?
                # Get global indices of regions in the window
                # idx = [map_indices[i] for i in regions_list[idx]]

                if idx is None or len(idx) <= 1:
                    # print("Less than two regions in window")
                    continue

                # Get submatrix
                if sp.sparse.issparse(AnnData.X):
                    window_accessibility = AnnData.X[:, idx].toarray()
                    window_scores = np.cov(window_accessibility, rowvar=False)
                    window_scores = window_scores + 1e-4 * np.eye(
                        len(window_scores))

                else:
                    window_accessibility = AnnData.X[:, idx].copy()
                    window_scores = np.cov(window_accessibility, rowvar=False)
                    window_scores = window_scores + 1e-4 * np.eye(
                        len(window_scores))

                distance = an.atacnet.get_distances_regions(AnnData[:, idx])

                # Test if distance is negative
                if np.any(distance < 0):
                    raise ValueError(
                        """
                        Distance between regions should be
                        positive. You might have overlapping
                        regions.
                        """
                    )

                window_penalties = an.atacnet.calc_penalty(
                    alpha,
                    distance=distance,
                    unit_distance=unit_distance)

                # Initiating graphical lasso
                graph_lasso_model = an.quic_graph_lasso.QuicGraphicalLasso(
                    init_method="precomputed",
                    lam=window_penalties,
                    tol=1e-4,
                    max_iter=1e4,
                    auto_scale=False,
                )

                # Fit graphical lasso
                graph_lasso_model.fit(window_scores)
                

                # Names of regions in the window
                window_region_names = AnnData.var_names[idx].copy()

                # Transform to correlation matrix
                scores = an.atacnet.cov_to_corr(graph_lasso_model.covariance_)
                # remove diagonal
                scores = scores - np.eye(scores.shape[0])

                # convert to sparse matrix the results
                corrected_scores = sp.sparse.coo_matrix(
                    scores)

                # Convert corrected_scores column
                # and row indices to global indices
                idx = [
                    map_indices[name]
                    for name in window_region_names[corrected_scores.row]
                ]
                idy = [
                    map_indices[name]
                    for name in window_region_names[corrected_scores.col]
                ]

                # Add the "sub" resuls to the global sparse matrix
                slide_results["scores"] = np.concatenate(
                    [slide_results["scores"], corrected_scores.data]
                )
                slide_results["idx"] = np.concatenate(
                    [slide_results["idx"], idx]
                    )
                slide_results["idy"] = np.concatenate(
                    [slide_results["idy"], idy]
                    )

            # Create sparse matrix
            results["window_" + str(k)] = sp.sparse.coo_matrix(
                (slide_results["scores"],
                 (slide_results["idx"], slide_results["idy"])),
                shape=(AnnData.X.shape[1], AnnData.X.shape[1]),
            )
    results = an.atacnet.reconcile(results)

    print("Done !")
    return results

# Parameters
distance_threshold = snakemake.params['distance_threshold']

# Load data
#atac = ad.read_csv(snakemake.input["atac"], delimiter='\t')
#print(atac.var_names)
#
#
## Add region infos
#an.add_region_infos(atac, sep=('_', '_'))
#print(atac)
#print(atac.X)
#print(atac.X.sum(0).min(), atac.X.sum(0).min())
#
#atac = an.metacells.compute_metacells(atac,
 #                              k = snakemake.params["number_cells_per_clusters"])


atac = ad.read_csv(snakemake.input["cicero_cds"], delimiter=' ').transpose()
print(atac)
# Add region infos
an.add_region_infos(atac, sep=('-', '-'))


# Compute network
an.compute_atac_network(
    atac,
    window_size=distance_threshold,
    unit_distance=1000,
    distance_constraint=distance_threshold/2,
    n_samples=100,
    n_samples_maxtry=500,
    max_alpha_iteration=300,
)

# Save results
peak_layer = extract_links(
    atac,
    columns=['peak1', 'peak2', 'score'])
peak_layer.to_csv(snakemake.output[0], sep='\t', index=False)
