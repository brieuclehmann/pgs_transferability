"""
Workflow to run UKB analysis for 'Optimal strategies for learning multi-ancestry polygenic scores vary across traits'
"""

# Import relevant libraries
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

# snakemake --snakefile Snakefile_ukb.smk --profile profile/ --cluster-config cluster_config.yaml --use-conda -j1 test

##########################
### Set up Paramspaces ###
##########################

# Global parameters
n_fold = 5
min_ancestries = ["AFR", "CSA", "AMR", "EAS", "MID"]
chromosomes = [i + 1 for i in range(22)] + ["X"]
pow_range = [0, 0.2, 0.4, 0.6, 0.8, 1.0] # gamma
pow_ext_range = [1.2, 1.4]

# Phenotype UKB codes
pheno_all = ["A50", "A30040", "N81", "NC1111"]
pheno_afr = [
    "A21001", "A30100", "A30090", "A30300", "A30130",
    "A30070", "A30210", "A30120", "NC1226", "I48", "K57"
    ]

# For setting up file names
pow_prop_multi = [(0.1, fold + 1, pow) for pow in pow_range for fold in range(n_fold)]
pow_prop_multi_ext = [(0.1, fold + 1, pow) for pow in pow_ext_range for fold in range(n_fold)]
pow_prop_min = [(prop, fold + 1, 0) for fold in range(n_fold) for prop in [1]]
pow_prop_maj = [(prop, fold + 1, 0) for fold in range(n_fold) for prop in [0]]
pow_prop_dual = [(-1, fold + 1, pow) for pow in pow_range for fold in range(n_fold)]
pow_prop_all = pow_prop_min + pow_prop_multi + pow_prop_maj

pow_prop_all_txt = [
    "prop_min~" + f'{prop:.1f}' + "/fold~" + str(fold) + "/pow~" + f'{power:.1f}' 
    for (prop, fold, power) 
    in pow_prop_all
]

pow_prop_sub_txt = [
    "prop_min~" + f'{prop:.1f}' + "/fold~" + str(fold) + "/pow~" + f'{power:.1f}'
    for (prop, fold, power)
    in pow_prop_multi + pow_prop_min
]

pow_prop_dual_txt = [
    "prop_min~" + f'{prop:.1f}' + "/fold~" + str(fold) + "/pow~" + f'{power:.1f}' 
    for (prop, fold, power) 
    in pow_prop_dual
]


pow_prop_ext = pow_prop_all + pow_prop_multi_ext
pow_prop_ext_txt = [
    "prop_min~" + f'{prop:.1f}' + "/fold~" + str(fold) + "/pow~" + f'{power:.1f}'
    for (prop, fold, power)
    in pow_prop_multi_ext
]


pow_prop_df = pd.DataFrame(pow_prop_all, columns = ['prop_min', 'fold', 'pow'])
pow_prop_dual_df = pd.DataFrame(pow_prop_dual, columns = ['prop_min', 'fold', 'pow'])
pow_prop_ext_df = pd.DataFrame(pow_prop_ext, columns = ['prop_min', 'fold', 'pow'])
chrom_df = pd.DataFrame({"chrom":chromosomes})
ancestry_df = pd.DataFrame({"min_ancestry":min_ancestries})
variant_df = pd.DataFrame({"v":["tagged", "imputed"]})
variant_df = pd.DataFrame({"v":["imputed"]})
pheno_df = pd.DataFrame({"pheno":pheno_all})


# AFR-only analyses #
afr_df = pd.DataFrame({
    "v": "imputed",
    "pheno":pheno_afr,
    "min_ancestry":"AFR"
})
afr_df = afr_df.merge(pow_prop_df, how = 'cross')

dual_df = pd.DataFrame({
    "v": "imputed",
    "pheno":pheno_all,
    "min_ancestry":"AFR"
})
dual_df = dual_df.merge(pow_prop_dual_df, how = 'cross')

# Remaining analyses #
all_df = variant_df.merge(pheno_df, how = 'cross')
all_df = all_df.merge(ancestry_df, how = 'cross')
all_df = all_df.merge(pow_prop_ext_df, how = 'cross')
all_df = all_df.loc[ # Insufficient FGP cases in these ancestry groups
    ~(all_df['min_ancestry'].isin(["EAS", "MID", "AMR"]) & 
    (all_df['pheno'] == "N81"))
] 
all_df = all_df.loc[ # Don't retrain EUR-only for other minority ancestries
    ~ ((all_df['min_ancestry'] != "AFR") & (all_df['prop_min'] == 0.0))
]
all_df = all_df.loc[ # Remove due to convergence issues for NC1111 (asthma)
    ~ (
    (all_df['min_ancestry'] == "MID") & 
    (all_df['pow'].isin([1.2, 1.4])) &
    (all_df['pheno'] == "NC1111")
    )
]



full_df = pd.concat([afr_df, dual_df, all_df])
full_df = full_df.merge(chrom_df, how = 'cross')

full_pred_df = full_df.copy()


# Paramspaces
ps_full = Paramspace(full_df)
ps_geno = Paramspace(full_df[['v', 'min_ancestry', 'chrom']].drop_duplicates())
ps_train = Paramspace(full_df.drop(["v", "chrom", "pow"], 1).drop_duplicates())
ps_pred = Paramspace(full_df.drop("chrom", 1).drop_duplicates())
ps_score = Paramspace(full_df[['v', 'pheno', 'min_ancestry']].drop_duplicates())

ps_full.geno_pattern = ps_geno.wildcard_pattern
ps_full.pred_pattern = ps_pred.wildcard_pattern
ps_full.train_pattern = ps_train.wildcard_pattern
ps_pred.score_pattern = ps_score.wildcard_pattern


################################
########### Workflow ###########
################################

#-#-#-#-# Prepare data #-#-#-#-#

### Phenotypes ###

# https://git.fmrib.ox.ac.uk/fsl/funpack
# Note that the data files are specific to UKBB Application 12788
rule ukb_funpack:
    output: "data/all_ukb_vars.tsv"
    conda: "envs/funpack.yml" 
    shell: "funpack -ow -q -vi 0 -ex data/w12788_20220222.csv "
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
    input: "data/all_ukb_vars.tsv", "scripts/ukb/02_preprocess.R"
    output: "data/all_vars.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/02_preprocess.R"

### Genotypes ###

rule filter_variants:
    input: "data/full_variant_qc_metrics.txt"
    output: f"data/snp_ids/{ps_geno.wildcard_pattern}.txt"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/03_filter_variants.R"

rule convert_pgen:
    input: f"data/snp_ids/{ps_geno.wildcard_pattern}.txt"
    output:
        multiext(f"data/genotypes/{ps_geno.wildcard_pattern}", ".pgen", ".pvar.zst", ".psam")
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
           " --out data/genotypes/v~{wildcards.v}/min_ancestry~{wildcards.min_ancestry}/chrom~{wildcards.chrom}"

rule pgen_wrapper:
    input:
        expand(
            "data/genotypes/{params}.pgen",
            params = ps_geno.instance_patterns
        )

#-#-#-#-# Estimate PGS #-#-#-#-#

rule test_train_wrapper:
    input:
        expand(
            "data/train_ids/{params}.txt",
            params = ps_train.instance_patterns
        )

rule test_train_split:
    input: "data/all_vars.tsv"
    output: f"data/train_ids/{ps_train.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/04a_training_sets.R"

rule fit_lasso:
    input:
        f"data/genotypes/{ps_full.geno_pattern}.pgen",
        f"data/train_ids/{ps_full.train_pattern}.txt"
    output: f"output/ukb/{ps_full.wildcard_pattern}.RDS"
    conda: "envs/snpnet.yml"
    wildcard_constraints: chrom="\w{1,2}"
    params:
        ncores = lambda wildcards: 16 if wildcards["prop_min"] in ['0.0', '-1.0'] else 2
    script: "scripts/ukb/05a_fit_lasso.R"


#-#-#-#-# Process output #-#-#-#-#

rule extract_effects:
    input:
        expand(
            "output/ukb/{pattern}/chrom~{chroms}.RDS",
            chroms = chromosomes,
            pattern = ps_pred.wildcard_pattern
        )
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/beta.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/06_effects.R"


rule predict_traits:
    input:
        expand(
            "output/ukb/{pattern}/chrom~{chroms}.RDS", 
            chroms=chromosomes, 
            pattern = ps_pred.wildcard_pattern
        )
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/pred.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/07a_predict.R"


def eval_input(wildcards): # Only predict for EUR-only when min_ancestry == "AFR"
    if wildcards.min_ancestry == "AFR":
        pow_prop = pow_prop_all_txt
        if wildcards.pheno in pheno_all:
            pow_prop = pow_prop + pow_prop_dual_txt + pow_prop_ext_txt
    else:
        pow_prop = pow_prop_sub_txt
        if wildcards.pheno in pheno_all:
            pow_prop = pow_prop + pow_prop_ext_txt
        if (wildcards.pheno == "NC1111") and (wildcards.min_ancestry == "MID"):
            pow_prop = pow_prop_sub_txt
    return expand(
        "output/ukb/{pattern}/{pow_prop}/pred.tsv", 
        pattern = ps_score.wildcard_pattern, 
        pow_prop = pow_prop
    )

rule evaluate_pgs:
    input: eval_input #expand("output/ukb/{pattern}/{pow_prop}/pred.tsv", pattern = ps_score.wildcard_pattern, pow_prop = pow_prop_all_txt)
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_score.wildcard_pattern}/scores.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/08a_performance.R"

rule var_decomposition:
    input: f"output/ukb/{ps_pred.wildcard_pattern}/beta.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    output: f"output/ukb/{ps_pred.wildcard_pattern}/decomp.tsv"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/09_decomposition.R"


#-#-#-#-# Wrapper rule #-#-#-#-#

rule ukb_wrapper:
    input:
        expand("output/ukb/{params}/scores.tsv", params=ps_score.instance_patterns),
        expand("output/ukb/{params}/decomp.tsv", params=ps_pred.instance_patterns)


##########################
## Sample size analysis ##
##########################

frac_range = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
frac_df = pd.DataFrame({"frac":frac_range})
type_df = pd.DataFrame({"type":["min", "maj"]})

# Set up sample size analyses
ss_df = full_df.loc[(full_df['min_ancestry'] == "AFR") & (full_df['v'] == "imputed")]
ss_df = ss_df.drop(["prop_min", "pow"], 1).drop_duplicates()
ss_df = ss_df.merge(type_df, how = 'cross')
ss_df = ss_df.merge(frac_df, how = 'cross')
ss_df = ss_df[[c for c in ss_df if c not in ['chrom']] + ['chrom']]

ps_ss = Paramspace(ss_df)
ps_ss_train = Paramspace(ss_df.drop(["v", "chrom"], 1).drop_duplicates())
ps_ss_pred = Paramspace(ss_df.drop("chrom", 1).drop_duplicates())
ps_ss_score = Paramspace(ss_df[['v', 'pheno', 'min_ancestry']].drop_duplicates())

ps_ss.geno_pattern = ps_geno.wildcard_pattern
ps_ss.pred_pattern = ps_ss_pred.wildcard_pattern
ps_ss.train_pattern = ps_ss_train.wildcard_pattern
ps_ss_pred.score_pattern = ps_ss_score.wildcard_pattern


frac_all_txt = [
    "fold~" + str(fold + 1) + "/type~" + frac_type + "/frac~" + f'{frac:.1f}'
    for fold in range(n_fold)
    for frac_type in ['min', 'maj']
    for frac in frac_range
]

rule ss_test_train_split:
    input: "data/all_vars.tsv"
    output: f"data/train_ids/{ps_ss_train.wildcard_pattern}.txt"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/04b_training_sets.R"

rule ss_fit_lasso:
    input:
        f"data/genotypes/{ps_ss.geno_pattern}.pgen",
        f"data/train_ids/{ps_ss.train_pattern}.txt"
    output: f"output/ukb/{ps_ss.wildcard_pattern}.RDS"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    params:
        ncores = 2
    script: "scripts/ukb/05b_fit_lasso.R"

rule ss_predict_traits:
    input:
        expand(
            "output/ukb/{pattern}/chrom~{chroms}.RDS",
            chroms=chromosomes,
            pattern = ps_ss_pred.wildcard_pattern
        )
    output: f"output/ukb/{ps_ss_pred.wildcard_pattern}/pred.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/07b_predict.R"

rule ss_evaluate_pgs:
    input:
        expand(
            "output/ukb/{pattern}/{frac_pattern}/pred.tsv",
            pattern = ps_ss_score.wildcard_pattern,
            frac_pattern = frac_all_txt
        )
    output: f"output/ukb/{ps_ss_score.wildcard_pattern}/samplesize_scores.tsv"
    wildcard_constraints: chrom="\w{1,2}"
    conda: "envs/snpnet.yml"
    script: "scripts/ukb/08b_performance.R"

rule ss_wrapper:
    input:
        expand(
            "output/ukb/{params}/samplesize_scores.tsv",
            params=ps_ss_score.instance_patterns
        )
