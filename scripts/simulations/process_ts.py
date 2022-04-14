# Script to process tree sequence and convert to vcf for plink

import tskit
import numpy as np
import gzip

ts = tskit.load("data/ooa.trees")

# Compute minor allele frequency
allele_freq = np.zeros((ts.num_populations, ts.num_sites))
for pop_id in range(ts.num_populations):
    ind_sub = ts.samples(population_id = pop_id)
    ts_pop = ts.simplify(ind_sub, filter_sites = False)
    this_freq = np.array([g.genotypes.mean() for g in ts_pop.variants()])
    allele_freq[pop_id,:] = this_freq

np.savetxt("data/maf.txt", allele_freq.T)
maf_mat = np.fmin(allele_freq, 1 - allele_freq)
max_maf = np.amax(maf_mat, 0)

# Apply MAF threshold
maf_thresh = 0.01
ind_del = np.where(max_maf < maf_thresh)
ts_keep = ts.delete_sites(ind_del)

# Save output as VCF
n_dip_indv = int(ts_keep.num_samples / 2)
indv_names = [f"tsk{str(i)}indv" for i in range(n_dip_indv)]
with gzip.open("data/ooa.vcf.gz", "wt") as vcf_file:
    ts_keep.write_vcf(vcf_file, ploidy=2, individual_names=indv_names,
                      position_transform = "legacy")