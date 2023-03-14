library(tibble)

all_phenos <- c(
    "A21001", "A30100", "A30040", "A30090", "A30300", "A30130",
    "A30070", "A30210", "A30120", "A50", "NC1226", "NC1111", "I48", "K57", "N81"
)
min_ancestries <- c("AFR", "EAS", "MID", "CSA", "AMR")

fracs <- seq(0, 1, 0.2)
nfolds <- 5

pow_range <- c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4")
chroms <- c(1:22, "X")

time_df <- tibble()

for (vtype in c("imputed", "tagged")) {
    for (pheno in all_phenos) {
        for (min_ancestry in min_ancestries) {
            for (prop_min in c("-1.0", "0.0", "0.1", "1.0")) {
                for (fold in seq_len(nfolds)) {
                    for (pow in pow_range) {
                        for (chrom in chroms) {

                            outfile <- file.path(
                                "output/ukb/", paste0("v~", vtype),
                                paste0("pheno~", pheno),
                                paste0("min_ancestry~", min_ancestry),
                                paste0("prop_min~", prop_min),
                                paste0("fold~", fold),
                                paste0("pow~", pow),
                                paste0("chrom~", chrom, ".RDS")
                            )

                            if (file.exists(outfile)) {
                                out <- readRDS(outfile)
                                this_time <- out$time

                                this_df <- tibble(vtype = vtype, pheno = pheno, min_ancestry = min_ancestry,
                                                  prop_min = prop_min, fold = fold, pow = pow, chrom = chrom,
                                                  user_self = this_time[1], sys_self = this_time[2], elapsed = this_time[3],
                                                  user_child = this_time[4], sys_child = this_time[5])

                                time_df <- rbind(time_df, this_df)
                            }
                        }
                    }
                }
            }
        }
    }
}

write_tsv(time_df, "output/times.tsv")
