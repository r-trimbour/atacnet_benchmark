#rule run_cicero:
#    input:
#        atac_path = "data/{sample}.csv.gz"
#    params:
#        distance_threshold = 500000,
#        number_cells_per_clusters = 10,
#        distance_parameter = "data/{sample}_distance_parameter.csv",
#        cicero_cds_rds = "results/{sample}/cicero/cicero_cds.RDS"
#    singularity:
#        "envs/hummus2.sif"
#    output:
#        network = "results/{sample}/cicero/cicero_results.csv",
#        cicero_cds = "results/{sample}/cicero/cicero_cds.tsv"
#    script:
#        "../scripts/cicero_fixed_params.R"


rule save_pseudocells_cicero:
    input:
        atac_path = "data/{sample}.csv.gz"
    output:
        cicero_cds_tsv = "results/{sample}/cicero/cicero_cds.tsv",
        cicero_cds_rds = "results/{sample}/cicero/cicero_cds.RDS",
    params:
        number_cells_per_clusters = 10,
        distance_parameter = "data/{sample}_distance_parameter.csv",
    singularity:
        "envs/hummus2.sif"
    script:
        "../scripts/cicero_save_cds.R"


# just redo for time
rule time_cicero:
    input:
        cicero_cds_rds = "results/{sample}/cicero/cicero_cds.RDS"
    output:
        network = "results/{sample}/cicero/cicero_results.csv",
    params:
        distance_threshold = 500_000,
        sample_num = 1_500
    benchmark:
        "results/{sample}/cicero_time.txt"
    singularity:
        "envs/hummus2.sif"
    script:
        "../scripts/timed_cicero.R"

# just redo for time
#rule time_cicero:
#    input:
#        cicero_cds_rds = "results/{sample}/cicero/cicero_cds.RDS"
#    output:
#        network = "results/{sample}/cicero/cicero_results.csv",
#    params:
#        distance_threshold = 500_000,
#        sample_num = 100
#    benchmark:
#        "results/{sample}/cicero_time.txt"
#    singularity:
#        "envs/hummus2.sif"
#    script:
#        "../scripts/timed_cicero2.R"