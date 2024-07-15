# Benchmarking  Cicero / atac_networks (python)
from snakemake.utils import min_version
min_version("7.0.0")

configfile: "config/data_config.yaml"
samples = ["test2"]

# Target rule to define the desired final output
rule all:
    input:
        #expand("results/{sample}/correlation_matrix.png", sample=samples),
        expand("results/{sample}/circe/circe_networks_results.csv", sample=samples),
        #expand("results/{sample}/cicero/cicero_results.csv", sample=samples)

# 4. Compare Cicero and atac_networks
rule compare:
    input:
        "results/{sample}/cicero/cicero_results.csv",
        "results/{sample}/atac_networks/atac_networks_results.csv"
    output:
        "results/{sample}/correlation_matrix.png"
    script:
        "scripts/compare.py"

# 3. Module to run Cicero
module cicero:
    snakefile: "modules/cicero.smk"
    config: config
use rule * from cicero

# 2. Module to run circe
module circe:
    snakefile: "modules/circe.smk"
    config: config
use rule * from circe

# 1. Generate/Preprocess the fake scATAC-seq data
module data_generation:
    snakefile: "modules/data_processing.smk"
    config: config
use rule * from data_generation
