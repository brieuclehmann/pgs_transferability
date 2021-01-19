# Combine output for a simulation

library(foreach)
library(dplyr)
source("sim_utils.R")
n_rep <- 50

beta_cor <- 0.5
prop_afr <- 0.1

this_sim <- paste0("beta_cor", beta_cor, "-prop_afr", prop_afr)

out_df <- foreach(iter = seq_len(n_rep), .combine = rbind) %do% {
	
	pheno_df <- readr::read_tsv(paste0("temp/", this_sim, iter, "_pheno.tsv"))
	pred_df <- readr::read_tsv(paste0("temp/", this_sim, iter, "_pred.tsv")) %>%
		inner_join(pheno_df, by = c("n", "pop")) 
		
	out <- pred_df %>%
		group_by(r, pow, pop) %>%
		summarise(bias = mean(pred) - mean(risk),
				  r2 = compute_r2(risk, pred),
				  r2_y = compute_r2(pheno, pred),
				  mse = mean((pred - risk)^2),
				  cor = cor(pred, risk), 
				  .groups = "drop")
	}
	
	out$iter <- iter
	
	out
}

dir.create("output", showWarnings=FALSE)
out_file <- paste0("output/", this_sim, ".csv")
readr::write_csv(out_df, out_file)
