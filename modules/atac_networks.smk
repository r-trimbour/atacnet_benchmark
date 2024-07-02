rule run_atac_networks:
    input:
        cicero_cds = "results/{sample}/cicero/cicero_cds.tsv",
        #atac = "data/{sample}.csv.gz"
    params:
        distance_threshold = 500_000,
        number_cells_per_clusters = 10
    benchmark:
        "results/{sample}/atacnet_time.txt"
    conda:
        "../envs/atac_networks.yaml"
    output:
        "results/{sample}/atac_networks/atac_networks_results.csv"
    script:
        "../scripts/atac_networks_aggregated.py"

        
#rule run_atac_networks:
#    input:
#        "data/{sample}.h5ad"
#    params:
#        distance_threshold = 500000
#    conda:
#        "../envs/atac_networks.yaml"
#    output:
#        "results/{sample}/atac_networks/atac_networks_results.csv"
