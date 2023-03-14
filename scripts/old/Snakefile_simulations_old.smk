"""
Workflow to run simulations for 'Optimal strategies for learning multi-ancestry polygenic scores vary across traits'
"""

# snakemake --snakefile Snakefile_simulations.smk --profile profile/ --cluster-config cluster_config.yaml --use-conda -j1 test

from snakemake.utils import min_version
from snakemake.utils import Paramspace
import pandas as pd
min_version("5.26")

paramspace = Paramspace(pd.read_csv("data/sim_params.csv"))
paramspace_cv = Paramspace(paramspace.loc[(paramspace['prop_afr'] == 0.1) & (paramspace['beta_cor'] == 0.7) & (paramspace['iter'].between(1,5))])

configfile: "config.yaml"
ruleorder: run_simul > run_simul_cv

onsuccess: print("finished successfully") # insert more useful code here, e.g. send email to yourself
onerror: print("finished with errors") # insert more useful code here, e.g. send email to yourself

#---#
# Just include necessary genetic map as part of the repo
rule download_genetic_maps:
    input: "Snakefile_simulations.smk"
    output:
        dir = directory("data/stdpopsim_cache")
    conda: "envs/stdpopsim.yml" 
    shell: "stdpopsim -c output.dir download-genetic-maps "
           "HomSap HapMapII_GRCh37"

#---#

rule generate_treesequence:
    #input: rules.download_genetic_maps.output.dir
    input: 
    output: "data/ooa.trees"
    params: NSAMPLES
    conda: "envs/stdpopsim.yml"
    shell: "stdpopsim HomSap {params} {params} {params} "
           "-d OutOfAfrica_3G09 "
           "-c chr20 "
           "-g HapMapII_GRCh37 "
           "-s 42 "
           "-o {output}"

#---#

rule process_treesequence:
    input: "data/ooa.trees" 
    output: "data/ooa.vcf.gz"
    conda: "envs/stdpopsim.yml"
    script: "scripts/simulations/02_process_ts.py"

#---#

rule convert_treesequence:
    input: "data/ooa.vcf.gz"
    output: "data/ooa.pgen"
    conda: "envs/plink2.yml"
    params:
        cores = 4
    shell: "plink2 --vcf {input} --threads {params.cores} --memory 40000 "
           "--make-pgen vzs psam-cols=fid --set-all-var-ids @:# --out data/ooa"


rule process_pgen:
    input: "data/ooa.pgen"
    output: "data/ooa.pgen"
    conda: "envs/plink2.yml"
    params:
        cores = 4
    shell: "plink2 --pfile data/ooa --threads {params.cores} --memory 40000 "
           "--set-all-var-ids @:#"

# --- #
rule simul_wrapper:
    input:
        expand(
            "output/simulations/{params}.csv",
            params=paramspace.instance_patterns
        )

rule simul_cv_wrapper:
    input:
        expand(
            "output/simulations/cv/{params}.csv",
            params=paramspace_cv.instance_patterns
        )

rule run_simul:
    input: "data/ooa.pgen"
    output: f"output/simulations/{paramspace.wildcard_pattern}.csv"
    conda: "envs/snpnet.yml"
    script: "scripts/simulations/04a_run.R"

rule run_simul_cv:
    input: "data/ooa.pgen"
    output: f"output/simulations/cv/{paramspace.wildcard_pattern}.csv"
    conda: "envs/snpnet.yml"
    script: "scripts/simulations/04b_run_cv.R"

rule calculate_h2:
    input: "data/sim_params.csv"
    output: "output/simulations/h2.csv"
    conda: "envs/snpnet.yml"
    script: "scripts/simulations/04c_calculate_h2.R"