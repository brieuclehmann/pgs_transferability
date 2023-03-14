# Script to truncate tree sequence to first 10% of chromosome

tskit <- reticulate::import("tskit")
ts <- tskit$load("data/ooa.trees")

samp_keep <- c(1:4000, 200001:204000) - 1L
ts = ts$simplify(samp_keep)

ind_keep <- seq(0, ts$num_sites*0.1)
ind_del <- setdiff(seq(0, ts$num_sites - 1), ind_keep)

ts10 <- ts$delete_sites(ind_del)
ts10$dump("data/ooa_chr20_10.trees")
