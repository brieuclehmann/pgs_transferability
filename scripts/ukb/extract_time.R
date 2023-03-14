library(readr)
library(tibble)

outfile <- snakemake@output[[1]]

this_dir <- dirname(outfile)
chroms <- c(1:22, "X")

time_df <- tibble()
for (chrom in chroms) {
    infile <- paste0(this_dir, "/chrom~", chrom, ".RDS")
    out <- readRDS(infile)
    this_time <- out$time

    time_df <- rbind(time_df, this_time)
}
names(time_df) <- c("user_self", "sys_self", "elapsed", "user_child", "sys_child")
time_df$chrom <- chroms

write_tsv(time_df, outfile)
