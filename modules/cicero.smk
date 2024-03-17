rule run_cicero:
    input:
        atac_path = "data/{sample}.csv.gz"
    params:
        distance_threshold = 500000,
        number_cells_per_clusters = 10
    singularity:
        "envs/hummus2.sif"
    output:
        network = "results/{sample}/cicero/cicero_results.csv",
        cicero_cds = "results/{sample}/cicero/cicero_cds.tsv"
    script:
        "../scripts/cicero_save_cds.R"