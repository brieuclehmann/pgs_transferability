rsync -arv \
	rxa753@cluster1.bmrc.ox.ac.uk:/well/holmes/users/rxa753/pgs_transferability/output/simulations/ \
	output/simulations
	
rsync -narv --include='*samplesize_scores.tsv' --include='*/' --exclude='*' \
	rxa753@cluster1.bmrc.ox.ac.uk:/well/holmes/users/rxa753/pgs_transferability/output/ukb/ \
	output/ukb
	
scp rxa753@cluster1.bmrc.ox.ac.uk:/well/holmes/users/rxa753/pgs_transferability/output/ntrain.tsv \
	output
	
