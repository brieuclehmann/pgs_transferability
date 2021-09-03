"""
Workflow to run UKB analysis for PGS transferability manuscript
inputs: Snakefile_simulations.smk (dummy input)
output: data/ooa.trees
"""

from snakemake.utils import min_version
from snakemake.utils import Paramspace
import pandas as pd
import numpy as np
min_version("5.26")

configfile: "config.yaml"

NSAMPLES = config["NSAMPLES"]

onsuccess: print("finished successfully") # insert more useful code here, e.g. send email to yourself
onerror:
    print("finished with errors")

### Set up Paramspaces ###
ancestries = ["AFR"]#, "CSA", "AMR", "EAS", "MID"]
chromosomes = [i + 1 for i in range(22)] + ["X"]
maf_params = [(x, y) for x in ancestries for y in chromosomes]
ps_maf = Paramspace(pd.DataFrame(maf_params, columns = ["ancestry", "chrom"]))

basic_df = pd.read_csv("data/ukb_basic_params.csv")
chrom_df = pd.DataFrame({"chrom":chromosomes})
pow_df = pd.DataFrame({"prop_min":np.repeat(0.1, 6), "pow":[0, 0.2, 0.4, 0.6, 0.8, 1]})
basic_lasso_df = basic_df.merge(pow_df, how = 'outer').merge(chrom_df, how = 'cross')
basic_lasso_df['pow'] = basic_lasso_df['pow'].fillna(0) 

ps_basic = Paramspace(basic_df)
ps_height = Paramspace(ps_basic.loc[ps_basic['pheno'] == "A4080"])
ps_basic_lasso = Paramspace(basic_lasso_df)
ps_height_lasso = Paramspace(basic_lasso_df.loc[basic_lasso_df['pheno'] == "A4080"])
ps_height_lasso.geno_pattern = 'ancestry~{min_ancestry}/chrom~{chrom}'
ps_height_lasso.train_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/prop_min~{prop_min}/fold~{fold}'

ps_pred = Paramspace(ps_height_lasso.drop("chrom", 1).drop_duplicates())
ps_pred.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'
ps_score = Paramspace(ps_height.drop(["fold", "prop_min"], 1).drop_duplicates())
#---#
# Extract and preprocess phenotypes
rule ukb_funpack:
    output: "data/all_ukb_vars.tsv"
    conda: "envs/funpack.yml" 
    shell: "funpack -ow -q -vi 0 -ex data/w12788_20210809.csv "
           "-v 50 -v 93 -v 4080 -v 21001 -v 5983 "
           "-v 30020 -v 30040 -v 30070 -v 30090 -v 30100 "
           "-v 30110 -v 30120 -v 30130 -v 30140 -v 30210 -v 30300 "
           "-v 20002 -v 22146 -v 6152 -v 3741 -v 6159 -v 41270 "
           "-v 31 -v 34 -v 52 -v 21000 -v 22000 "
           "data/all_ukb_vars.tsv "
           "data/ukbb12788/ukbb_download_42801/ukb42801.csv "
           "data/ukbb12788/ukbb_download_2010607/ukb44775.csv "
           "data/ukbb12788/ukbb_download_35062/ukb35062.csv "
           "data/ukbb12788/ukbb_download_7749/ukb7749.csv "
           "data/ukbb12788/ukbb_download_27864/ukb27864.csv "
           "data/ukbb12788/ukbb_download_37088/ukb37088.csv "

rule combine_phenotypes:
    input: "data/all_ukb_vars.tsv", "scripts/ukb_02_preprocess.R"
    output: "data/all_vars.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_02_preprocess.R"

#---#

rule basic_train_wrapper:
    input:
        expand("data/train_ids/{params}.txt", params=ps_basic.instance_patterns)

rule basic_train_split:
    input: "data/all_vars.tsv"
    output: f"data/train_ids/{ps_basic.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_03_training_sets.R"

rule ukb_maf_wrapper:
    input:
        expand("data/snp_ids/{params}.txt", params=ps_maf.instance_patterns)

rule ukb_maf:
    input: "data/full_variant_qc_metrics.txt"
    output: f"data/snp_ids/{ps_maf.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_04_get_ancestry_maf.R"

rule convert_pgen_wrapper:
    input:
        expand("data/genotypes/{params}.pgen", params=ps_maf.instance_patterns)

rule convert_pgen:
    input: f"data/snp_ids/{ps_maf.wildcard_pattern}.txt"
    output:
        multiext(f"data/genotypes/{ps_maf.wildcard_pattern}", ".pgen", ".pvar.zst", ".psam")
    conda: "envs/plink2.yml"
    params:
        imp_dir = "/well/ukbb-wtchg/v3/imputation/",
        cores = 6,
        memory = 12000 * 6
    shell: "plink2 --bgen {params.imp_dir}ukb_imp_chr{wildcards.chrom}_v3.bgen "
           " ref-first "
           " --sample data/links/ukb12788_imp_chr{wildcards.chrom}_v3.sample "
           " --extract {input} "
           " --threads {params.cores} --memory {params.memory} "
           " --make-pgen 'vzs' "
           " --set-all-var-ids @:#_\$r_\$a "
           "--new-id-max-allele-len 700"
           " --out data/genotypes/ancestry~{wildcards.ancestry}/chrom~{wildcards.chrom}"

rule basic_lasso_wrapper:
    input:
        expand("output/ukb/{params}.RDS", params=ps_height_lasso.instance_patterns)

rule basic_lasso:
    input:
        f"data/genotypes/{ps_height_lasso.geno_pattern}.pgen",
        f"data/train_ids/{ps_height_lasso.train_pattern}.txt"
    output: f"output/ukb/{ps_height_lasso.wildcard_pattern}.RDS"
    conda: "envs/snpnet.yml"
    params:
        ncores = lambda wildcards: 6 if wildcards["prop_min"] == '0.0' else 1
    script: "scripts/ukb_05_fit_lasso.R"

rule basic_pred_wrapper:
    input:
        expand("output/ukb/{params}/pred.tsv", params=ps_pred.instance_patterns)

rule basic_pred:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_pred.wildcard_pattern)
    output: f"output/ukb/{ps_pred.wildcard_pattern}/pred.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06_predict.R"

rule basic_scores_wrapper:
    input:
        expand("output/ukb/{params}/scores.tsv", params=ps_score.instance_patterns)

rule basic_scores:
    input:
        rules.basic_pred_wrapper.input
    output: f"output/ukb/{ps_score.wildcard_pattern}/scores.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_07_performance.R"
