#!/bin/bash
# Script to download UKBB genetic correlation estimates from 
# https://ukbb-rg.hail.is and UKBB heritability estimates from
# https://nealelab.github.io/UKBB_ldsc/

mkdir data

# Download genetic correlation estimates
wget https://www.dropbox.com/sh/qvuil7op8bw68fm/AACBkXBnj4rPUrlRq3Hm514ha?dl=1 -O data/genetic_correlation_matrices.zip

# Unzip and extract
unzip genetic_correlation_matrices.zip -d genetic_correlation_matrices
gunzip genetic_correlation_matrices/geno_correlation.r2.gz


# Download heritability estimates
wget https://www.dropbox.com/s/8vca84rsslgbsua/ukb31063_h2_topline.02Oct2019.tsv.gz?dl=1 -O data/ukb31063_h2_topline.02Oct2019.tsv.gz

# Extract
gunzip data/ukb31063_h2_topline.02Oct2019.tsv.gz 