# This script uses the fmrib unpack tool to extract variables from the main
# csv file. See the R script preprocess.R for preprocessing post extraction.

### Body measures:

# 50 - Standing height
# 93 - Systolic BP
# 4080 - Systolic BP (auto)
# 21001 - BMI

INDIR=data/ukbb12788/ukbb_download_42801/
IND_EXCLUDE=data/w12788_20210809.csv

funpack -ow -q \
  -v 50 -v 93 -v 4080 -v 21001 \
  -vi 0 -ex ${IND_EXCLUDE} \
  data/body_phenotypes_raw.tsv ${INDIR}ukb42801.csv 

# 5983 ECG heart rate

INDIR=data/ukbb12788/ukbb_download_2010607/

funpack -ow -q \
  -v 5983 \
  -vi 0 -ex ${IND_EXCLUDE} \
  data/ecg_raw.tsv ${INDIR}ukb44775.csv

### Blood measurement

# 30020 - haemoglobin
# 30040 - MCV
# 30070 - RBC distribution
# 30090 - platelet crit
# 30100 - mean platelet volume
# 30110 - platelet distribution width
# 30120 - lymphocyte count
# 30130 - monocyte count
# 30140 - neutrophill count
# 30210 - eosinophill percentage
# 30300 - reticulocyte count
INDIR=data/ukbb12788/ukbb_download_35062/

funpack -ow -q \
  -v 30020 -v 30040 -v 30070 -v 30090 -v 30100 \
  -v 30110 -v 30120 -v 30130 -v 30140 -v 30210 -v 30300 \
  -vi 0 -ex ${IND_EXCLUDE} \
  data/blood_phenotypes_raw.tsv ${INDIR}ukb35062.csv


### Miscellaneous disease phenotypes

# 20002 - Non-cancer illness codes (9 total)
# 22146 - Age hayfever diagnosed
# 6152 - Blood clot (6152_5)

INDIR=data/ukbb12788/ukbb_download_7749/

funpack -ow -q \
  -v 20002 -v 22146 -v 6152 \
  -vi 0 -ex ${IND_EXCLUDE} \
  data/misc_phenotypes1_raw.tsv ${INDIR}ukb7749.csv

# 3741 - Stomach pain
# 6159 - No pain experienced (6159_100)
INDIR=data/ukbb12788/ukbb_download_27864/

funpack -ow -q \
  -v 3741 -v 6159 \
  -vi 0 -ex ${IND_EXCLUDE} \
  data/misc_phenotypes2_raw.tsv ${INDIR}ukb27864.csv

### ICD-10 codes
INDIR=/well/mcvean/ukbb12788/ukbb_download_27864/

funpack -ow -q \
  -v 41270 \
  -vi 0 -ex ${IND_EXCLUDE} \
  data/icd10_raw.tsv ${INDIR}ukb27864.csv

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
  -vi 0 -ex ${IND_EXCLUDE} \
  data/covariates_raw.tsv ${INDIR}ukb37088.csv
