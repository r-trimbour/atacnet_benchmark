rule create_fake_data:
    output:
        "data/test1.h5ad",
        "data/test1.csv.gz",
        "data/test1_ground_truth.tsv"
    params:
        nb_cells = 300,
        nb_chr = 10,
        nb_regions_per_chr = 200,
        between_reg = 20000,
        size_reg = 500,
        sep = '_',
        distance_threshold = 500000
    conda:
        "../envs/atac_networks.yaml"
    script:
        "../scripts/create_data.py"

rule download_test2_data:
    output:
        mudata = "data/test2.h5mu",
        #"data/{sample}_ground_truth.tsv"
    params:
        url = "https://figshare.com/ndownloader/files/44639476?private_link=b0840d90e42e37fa165f"
    conda:
        "../envs/atac_networks.yaml"
    shell:
        """
        wget {params.url} -O {output.mudata}
        """

rule preprocess_test2_data:
    input:
        mudata = "data/test2.h5mu"
    output:
        anndata = "data/test2.h5ad",
        csv = "data/test2.csv.gz"
    conda:
        "../envs/atac_networks.yaml"
    script:
        "../scripts/preprocess_mudata.py"