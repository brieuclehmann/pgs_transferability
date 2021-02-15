# This script uses the fmrib unpack tool to extract variables from the main
# csv file. See the R script preprocess.R for preprocessing post extraction.

### Body measures:

# 48 - Waist circumference
# 50 - Standing height
# 93 - Systolic BP
# 94 - Diastolic BP
# 4079 - Diastolic BP (auto)
# 4080 - Systolic BP (auto)
# 21001 - BMI

INDIR=data/ukbb12788/ukbb_download_42801/

funpack -ow -q \
  -v 48 -v 50 -v 93 -v 94 -v 4079 -v 4080 -v 21001 \
  -vi 0 \
  data/body_phenotypes_raw.tsv ${INDIR}ukb42801.csv

### Blood measurement

# 30000 - WBC
# 30010 - RBC
# 30020 - haemoglobin
# 30030 - haematocrit
# 30040 - MCV
# 30050 - MCH
# 30080 - platelet
# 30120 - lymphocyte
# 30130 - monocyte
# 30710 - C reactive protein
INDIR=data/ukbb12788/ukbb_download_35062/

funpack -ow -q \
  -v 30000 -v 30010 -v 30020 -v 30030 -v 30040 -v 30050 -v 30080 \
  -v 30120 -v 30130  -v 30710 \
  -vi 0 \
  data/blood_phenotypes_raw.tsv ${INDIR}ukb35062.csv

### General covariates

# 31 - Sex
# 34 - Year of birth
# 52 - Month of birth
# 21000 - Self-reported ancestry
# 22000 - Genotype batch
# 22009 - Genetic PCs

INDIR=data/ukbb12788/ukbb_download_37088/

funpack -ow -q \
  -v 31 -v 34 -v 52 -v 21000 -v 22000 -v 22009 \
  -vi 0 \
  data/covariates_raw.tsv ${INDIR}ukb37088.csv
