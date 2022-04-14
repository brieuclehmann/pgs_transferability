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

ruleorder: basic_lasso > full_lasso
ruleorder: basic_scores > full_scores
ruleorder: basic_beta > full_beta
ruleorder: basic_pred > full_pred
ruleorder: basic_train_split > full_train_split
#ruleorder: full_get_genotype > full_get_maf
#ruleorder: genotype_convert_pgen > full_convert_pgen

##########################
### Set up Paramspaces ###
##########################
n_fold = 5
all_ancestries = ["CSA", "AMR", "EAS", "MID"] # AFR handled separately
chromosomes = [i + 1 for i in range(22)] + ["X"]
pow_range = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"]
pow_range = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
pow_prop_mix = [(0.1, fold + 1, pow) for pow in pow_range for fold in range(n_fold)]
pow_prop_min = [(prop, fold + 1, 0) for fold in range(n_fold) for prop in [1]]
pow_prop_maj = [(prop, fold + 1, 0) for fold in range(n_fold) for prop in [0]]
pow_prop_all_txt = ["prop_min~" + f'{prop:.1f}' + "/fold~" + str(fold) + "/pow~" + f'{power:.1f}' for (prop, fold, power) in pow_prop_mix + pow_prop_min + pow_prop_maj]
pow_prop_sub_txt = ["prop_min~" + f'{prop:.1f}' + "/fold~" + str(fold) + "/pow~" + f'{power:.1f}' for (prop, fold, power) in pow_prop_mix + pow_prop_min]
#print(pow_prop_txt)

cts_pheno = ["A21001", "A30100", "A30040", "A30090", "A30300", "A30130", "A30070", "A30210", "A30120", "A50"]
bin_pheno = ["NC1226", "NC1111", "I48", "K57", "N81"]
all_pheno = cts_pheno + bin_pheno

maf_params = [(x, y) for x in all_ancestries for y in chromosomes]
ps_maf = Paramspace(pd.DataFrame(maf_params, columns = ["ancestry", "chrom"]))
ps_geno = Paramspace(pd.DataFrame(["AFR"] + all_ancestries, columns = ["ancestry"]))

afr_params = [("AFR", y) for y in chromosomes]
ps_afr = Paramspace(pd.DataFrame(afr_params, columns = ["ancestry", "chrom"]))

chrom_df = pd.DataFrame({"chrom":chromosomes})
pow_df = pd.DataFrame(
    {"prop_min":np.repeat(0.1, 6), "pow":pow_range}
    )

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
################################
### AFR-only 'basic' analysis ##
################################

# Set up basic parametes
basic_df = pd.read_csv("data/ukb_params.csv")
basic_df = basic_df.loc[basic_df['pheno'].isin(all_pheno)]
basic_df = basic_df.loc[basic_df['min_ancestry'] == "AFR"]
#basic_df = basic_df.loc[basic_df['pheno'] == "A50"]  
basic_df = basic_df.merge(pow_df, how = 'outer').merge(chrom_df, how = 'cross')
basic_df['pow'] = basic_df['pow'].fillna(0)

ps_basic = Paramspace(basic_df)
ps_basic_lasso = Paramspace(basic_df)
ps_basic_lasso.geno_pattern = 'ancestry~{min_ancestry}/chrom~{chrom}'
ps_basic_lasso.pred_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/prop_min~{prop_min}/fold~{fold}/pow~{pow}'
ps_basic_lasso.full_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/prop_min~{prop_min}/fold~{fold}/pow~{pow}/chrom~{chrom}'
ps_basic_lasso.train_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/prop_min~{prop_min}/fold~{fold}'
ps_basic_train = Paramspace(basic_df.drop(["chrom", "pow"], 1).drop_duplicates()) 

ps_pred = Paramspace(basic_df.drop("chrom", 1).drop_duplicates())
ps_pred.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'
ps_score = Paramspace(basic_df.drop(["fold", "prop_min", "chrom", "pow"], 1).drop_duplicates())
ps_score.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'

rule basic_train_split:
    input: "data/all_vars.tsv"
    output: f"data/train_ids/{ps_basic_train.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_03_training_sets.R"

rule basic_lasso:
    input:
        f"data/genotypes/{ps_basic_lasso.geno_pattern}.pgen",
        f"data/train_ids/{ps_basic_lasso.train_pattern}.txt"
    output: f"output/ukb/{ps_basic_lasso.wildcard_pattern}.RDS"
    conda: "envs/snpnet.yml"
    wildcard_constraints: chrom="\w{1,2}"
    params:
        ncores = lambda wildcards: 6 if wildcards["prop_min"] == '0.0' else 2
    script: "scripts/ukb_05_fit_lasso.R"

rule basic_pred_wrapper:
    input:
        expand("output/ukb/{params}/pred.tsv", params=ps_pred.instance_patterns)

rule basic_pred:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_pred.wildcard_pattern)
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/pred.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06a_predict.R"


rule basic_beta:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_pred.wildcard_pattern)
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/beta.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06a_effects.R"

rule basic_decomposition:
    input: f"output/ukb/{ps_pred.wildcard_pattern}/beta.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/decomp.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_09a_decomposition.R"

rule basic_fst:
    input: expand("output/ukb/{pattern}/{pow_prop}/beta.tsv", pattern = ps_score.wildcard_pattern, pow_prop = pow_prop_all_txt)
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_score.wildcard_pattern}/fst.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_09_fst.R"

#   " --make-pgen --out " f"output/ukb/{ps_pred.wildcard_pattern}/chrom~{wildcards.chrom}" 
#    shell: "./plink2 --pmerge-list mergelist.txt pfile-vzs "
#       " --pmerge-list-dir data/genotypes/ancestry~{wildcards.ancestry}/ "
#        " --pheno data/all_vars_plink.tsv "
#        " --extract {output.nonzero} "
#        " --fst pop 'base=EUR' "
#        " --out" f"output/ukb/{ps_pred.wildcard_pattern}/plink2"
# --pfile vzs data/genotypes/ancestry~{wildcards.ancestry}/chrom~{wildcards.chrom} "

rule basic_nonzero:
    input: f"output/ukb/{ps_pred.wildcard_pattern}/beta.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/nonzero.tsv"
    shell:  "cat {input} | awk '{{print $1}}' | tail -n +2 | sed 's|\(.*\)_.*|\1|' > {output} "

rule basic_fst_extract:
    input: f"output/ukb/{ps_pred.wildcard_pattern}/nonzero.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/plink2.yml"
    output: multiext(f"output/ukb/{ps_pred.wildcard_pattern}/pgs_geno", ".pgen", ".pvar.zst", ".psam")
    params: prefix=f"output/ukb/{ps_pred.wildcard_pattern}/pgs_geno"
    shell: "./plink2 --pmerge-list mergelist.txt pfile-vzs "
           " --pmerge-list-dir data/genotypes/ancestry~{wildcards.min_ancestry}/ "
           " --extract {input} --threads 4 "
           " --make-pgen --out {params.prefix} "

rule basic_fst_extract2:
    input: f"output/ukb/{ps_basic_lasso.pred_pattern}/nonzero.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/plink2.yml"
    output: expand("output/ukb/{pattern}_geno.{ext}", pattern = ps_basic_lasso.wildcard_pattern, ext = ["pgen", "pvar.zst", "psam"])
    params:
      prefix=f"output/ukb/{ps_basic_lasso.wildcard_pattern}_geno", geno=f"data/genotypes/{ps_basic_lasso.geno_pattern}"
    shell: "plink2 --pfile vzs {params.geno} "
           " --extract {input} --threads 4 "
           " --make-pgen --out {params.prefix} "

rule basic_scores:
    input:
        expand("output/ukb/{pattern}/{pow_prop}/pred.tsv", pattern = ps_score.wildcard_pattern, pow_prop = pow_prop_all_txt)
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_score.wildcard_pattern}/scores.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_07b_performance.R"


rule basic_wrapper:
    input:
        expand("output/ukb/{params}/scores.tsv", params=ps_score.instance_patterns),
        expand("output/ukb/{params}/decomp.tsv", params=ps_pred.instance_patterns),
        expand("output/ukb/{params}/fst.tsv", params=ps_score.instance_patterns)
        #expand("output/ukb/{pattern}_geno.pgen", pattern = ps_basic_lasso.instance_patterns)
        #expand("output/ukb/{params}/geno.pgen", params=ps_pred.instance_patterns)


#######################
## Enhanced analysis ##
#######################

# Set up full
full_pheno = ["A50", "A30040", "A30090", "NC1111", "N81"]

full_df = pd.read_csv("data/ukb_params.csv")
full_df = full_df.loc[full_df['pheno'].isin(full_pheno)]
full_df = full_df.loc[full_df['min_ancestry'].isin(all_ancestries + ["AFR"])]
full_df = full_df.loc[~(full_df['min_ancestry'].isin(["EAS", "MID", "AMR"]) & (full_df['pheno'] == "N81"))]
full_df = full_df.merge(pow_df, how = 'outer').merge(chrom_df, how = 'cross')
full_df['pow'] = full_df['pow'].fillna(0)
full_df2 = full_df.copy()
full_df = full_df.loc[~(full_df['min_ancestry'].isin(all_ancestries) & (full_df['prop_min'] == "0.0"))]

ps_full = Paramspace(full_df)
ps_full_lasso = Paramspace(full_df)
ps_full_lasso.geno_pattern = 'ancestry~{min_ancestry}/chrom~{chrom}'
ps_full_lasso.train_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/prop_min~{prop_min}/fold~{fold}'
ps_full_train = Paramspace(full_df.drop(["chrom", "pow"], 1).drop_duplicates()) 

ps_full_pred = Paramspace(full_df2.drop("chrom", 1).drop_duplicates())
ps_full_pred.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'
ps_full_score = Paramspace(full_df.drop(["fold", "prop_min", "chrom", "pow"], 1).drop_duplicates())
ps_full_score.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'

rule full_train_split:
    input: "data/all_vars.tsv"
    output: f"data/train_ids/{ps_full_train.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_03_training_sets.R"

rule full_get_maf:
    input: "data/full_variant_qc_metrics.txt"
    output: f"data/snp_ids/{ps_maf.wildcard_pattern}.txt"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_04_get_ancestry_maf.R"

rule full_convert_pgen:
    input: f"data/snp_ids/{ps_maf.wildcard_pattern}.txt"
    output:
        multiext(f"data/genotypes/{ps_maf.wildcard_pattern}", ".pgen", ".pvar.zst", ".psam")
    conda: "envs/plink2.yml"
    wildcard_constraints: chrom="\w{1,2}" 
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

rule full_lasso:
    input:
        f"data/genotypes/{ps_full_lasso.geno_pattern}.pgen",
        f"data/train_ids/{ps_full_lasso.train_pattern}.txt"
    output: f"output/ukb/{ps_full_lasso.wildcard_pattern}.RDS"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    params:
        ncores = lambda wildcards: 6 if wildcards["prop_min"] == '0.0' else 2
    script: "scripts/ukb_05_fit_lasso.R"

rule full_pred:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_full_pred.wildcard_pattern)
    output: f"output/ukb/{ps_full_pred.wildcard_pattern}/pred.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06a_predict.R"


rule full_beta:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_full_pred.wildcard_pattern)
    output: f"output/ukb/{ps_full_pred.wildcard_pattern}/beta.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06a_effects.R"


rule full_scores:
    input:
        expand("output/ukb/{pattern}/{pow_prop}/pred.tsv", pattern = ps_full_score.wildcard_pattern, pow_prop = pow_prop_sub_txt)
    output: f"output/ukb/{ps_full_score.wildcard_pattern}/scores.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_07b_performance.R"


rule full_wrapper:
    input:
        expand("output/ukb/{params}/scores.tsv", params=ps_full_score.instance_patterns),
        expand("output/ukb/{params}/decomp.tsv", params=ps_full_pred.instance_patterns),
        expand("output/ukb/{params}/fst.tsv", params=ps_full_score.instance_patterns) 

#######################
## Genotype analysis ##
#######################

rule full_get_genotype:
    input: "data/full_variant_qc_metrics.txt"
    output: f"data/snp_ids/{ps_maf.wildcard_pattern}_genotyped.txt"
    wildcard_constraints:
        ancestry="\w{3}",
        chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_04b_get_ancestry_maf_genotype.R"

rule genotype_convert_pgen:
    input: f"data/snp_ids/{ps_maf.wildcard_pattern}_genotyped.txt"
    output:
        multiext(f"data/genotypes/{ps_maf.wildcard_pattern}_genotyped", ".pgen", ".pvar.zst", ".psam")
    wildcard_constraints: ancestry="\w{3}"
    conda: "envs/plink2.yml"
    params:
        imp_dir = "/well/ukbb-wtchg/v3/imputation/",
        ncores = 6,
        memory = 12000 * 6
    shell: "plink2 --bgen {params.imp_dir}ukb_imp_chr{wildcards.chrom}_v3.bgen "
           " ref-first "
           " --sample data/links/ukb12788_imp_chr{wildcards.chrom}_v3.sample "
           " --extract {input} "
           " --threads {params.ncores} --memory {params.memory} "
           " --make-pgen 'vzs' "
           " --set-all-var-ids @:#_\$r_\$a "
           "--new-id-max-allele-len 700"
           " --out data/genotypes/ancestry~{wildcards.ancestry}/chrom~{wildcards.chrom}_genotyped"

rule genotype_lasso:
    input:
        f"data/genotypes/{ps_basic_lasso.geno_pattern}_genotyped.pgen",
        f"data/train_ids/{ps_basic_lasso.train_pattern}.txt"
    output: f"output/ukb/{ps_basic_lasso.wildcard_pattern}_genotyped.RDS"
    conda: "envs/snpnet.yml"
    params:
        ncores = lambda wildcards: 6 if wildcards["prop_min"] == '0.0' else 2
    script: "scripts/ukb_05_fit_lasso.R"

rule genotype_pred:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}_genotyped.RDS", chroms=chromosomes, pattern = ps_full_pred.wildcard_pattern)
    output: f"output/ukb/{ps_full_pred.wildcard_pattern}/pred_genotyped.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06b_genotype_predict.R"


rule genotype_beta:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}_genotyped.RDS", chroms=chromosomes, pattern = ps_full_pred.wildcard_pattern)
    output: f"output/ukb/{ps_full_pred.wildcard_pattern}/beta_genotyped.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06b_genotype_effects.R"

def score_input(wildcards):
    return expand("output/ukb/{pattern}/{pow_prop}/pred_genotyped.tsv", pattern = ps_full_score.wildcard_pattern, pow_prop = pow_prop_all_txt if wildcards.min_ancestry == "AFR" else pow_prop_sub_txt)

rule genotype_scores:
    input: score_input
    output: f"output/ukb/{ps_full_score.wildcard_pattern}/scores_genotyped.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_07b_genotype_performance.R"


rule genotype_wrapper:
    input:
        expand("output/ukb/{params}/scores_genotyped.tsv", params=ps_full_score.instance_patterns),
        expand("output/ukb/{params}/beta_genotyped.tsv", params=ps_full_pred.instance_patterns) 

##########################
## Sample size analysis ##
##########################

frac_range = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
frac_df = pd.DataFrame({"frac":frac_range})
type_df = pd.DataFrame({"type":["min", "maj"]})


# Set up sample size analyses
ss_df = pd.read_csv("data/ukb_params.csv")
ss_df = ss_df.drop("prop_min", 1).drop_duplicates()
ss_df = ss_df.loc[((ss_df['pheno'].isin(all_pheno)) & (ss_df['min_ancestry'] == "AFR")) | ((ss_df['pheno'].isin(full_pheno)) & (ss_df['min_ancestry'].isin(all_ancestries)))]
ss_df = ss_df.loc[ss_df['min_ancestry'].isin(["AFR"])]
ss_df = ss_df.merge(type_df, how = 'cross').merge(frac_df, how = 'cross').merge(chrom_df, how = 'cross')

ps_ss = Paramspace(ss_df)
ps_ss_lasso = Paramspace(ss_df)
ps_ss_lasso.geno_pattern = 'ancestry~{min_ancestry}/chrom~{chrom}'
ps_ss_lasso.train_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/fold~{fold}/type~{type}/frac~{frac}'
ps_ss_train = Paramspace(ss_df.drop(["chrom"], 1).drop_duplicates()) 

ps_ss_pred = Paramspace(ss_df.drop("chrom", 1).drop_duplicates())
ps_ss_pred.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'
ps_ss_score = Paramspace(ss_df.drop(["fold", "type", "chrom", "frac"], 1).drop_duplicates())
ps_ss_score.score_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}'

frac_all_txt = ["fold~" + str(fold + 1) + "/type~" + frac_type + "/frac~" + f'{frac:.1f}' for fold in range(n_fold) for frac_type in ['min', 'maj'] for frac in frac_range]

rule samplesize_train_split:
    input: "data/all_vars.tsv"
    output: f"data/train_ids/{ps_ss_train.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_03b_training_sets.R"

rule samplesize_lasso:
    input:
        f"data/genotypes/{ps_ss_lasso.geno_pattern}.pgen",
        f"data/train_ids/{ps_ss_lasso.train_pattern}.txt"
    output: f"output/ukb/{ps_ss_lasso.wildcard_pattern}.RDS"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    params:
        ncores = 2
    script: "scripts/ukb_05b_fit_lasso.R"

rule samplesize_pred:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_ss_pred.wildcard_pattern)
    output: f"output/ukb/{ps_ss_pred.wildcard_pattern}/pred.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06b_predict.R"


rule samplesize_beta:
    input:
        expand("output/ukb/{pattern}/chrom~{chroms}.RDS", chroms=chromosomes, pattern = ps_ss_pred.wildcard_pattern)
    output: f"output/ukb/{ps_ss_pred.wildcard_pattern}/beta.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_06b_effects.R"

rule samplesize_scores:
    input:
        expand("output/ukb/{pattern}/{frac_pattern}/pred.tsv", pattern = ps_ss_score.wildcard_pattern, frac_pattern = frac_all_txt)
    output: f"output/ukb/{ps_ss_score.wildcard_pattern}/samplesize_scores.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb_07c_performance.R"

rule samplesize_wrapper:
    input:
        expand("output/ukb/{params}/samplesize_scores.tsv", params=ps_ss_score.instance_patterns),
        #expand("output/ukb/{params}/beta.tsv", params=ps_ss_pred.instance_patterns) 

#######################
## Old wrapper rules ##
#######################

rule basic_train_wrapper:
    input:
        expand("data/train_ids/{params}.txt", params=ps_basic.instance_patterns)

rule basic_maf_wrapper:
    input:
        expand("data/snp_ids/{params}.txt", params=ps_afr.instance_patterns)

rule convert_pgen_wrapper:
    input:
        expand("data/genotypes/{params}.pgen", params=ps_maf.instance_patterns)

rule basic_lasso_wrapper:
    input:
        expand("output/ukb/{params}.RDS", params=ps_basic_lasso.instance_patterns)

rule basic_scores_wrapper:
    input:
        expand("output/ukb/{params}/scores.tsv", params=ps_score.instance_patterns)

rule basic_beta_wrapper:
    input:
        expand("output/ukb/{params}/beta.tsv", params=ps_pred.instance_patterns)

rule full_pred_wrapper:
    input:
        expand("output/ukb/{params}/pred.tsv", params=ps_full_pred.instance_patterns)
#ps_height = Paramspace(ps_basic.loc[ps_basic['pheno'].isin(cts_pheno)])
#ps_height_lasso = Paramspace(basic_df.loc[basic_df['pheno'].isin(cts_pheno)])
#ps_height_lasso.geno_pattern = 'ancestry~{min_ancestry}/chrom~{chrom}'
#ps_height_lasso.train_pattern = 'pheno~{pheno}/min_ancestry~{min_ancestry}/prop_min~{prop_min}/fold~{fold}'
