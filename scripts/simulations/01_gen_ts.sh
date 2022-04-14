#!/bin/bash
# Script to simulate tree sequence data from an Out-of-Africa model using 
# stdpopsim

NSAMPLES=400000

mkdir data

stdpopsim HomSap $NSAMPLES $NSAMPLES $NSAMPLES \
  -d OutOfAfrica_3G09 \
  -c chr20 \
  -g HapMapII_GRCh37 \
  -s 52 \
  -o data/ooa.trees
