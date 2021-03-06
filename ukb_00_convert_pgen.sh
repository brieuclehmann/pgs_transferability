# This script merges and converts UKB bed to pgen for use with snpnet

# Path to ukb genotype data
INDIR=data/ukb/v2/calls/

# Path to .fam file associated with ukb genotype calls
FAM=data/ukb12788_cal_chr1_v2_s488264.fam
for chr in {1..22}
do
  echo ${INDIR}ukb_cal_chr${chr}_v2.bed ${INDIR}ukb_snp_chr${chr}_v2.bim $FAM >> data/ukb_merge.txt
done

mdkir -p data/ukb_cal_pgen

plink --bed ${INDIR}ukb_cal_chr1_v2.bed \
  --bim ${INDIR}ukb_snp_chr1_v2.bim \
  --fam data/ukb12788_cal_chr1_v2_s488264.fam \
  --merge-list data/ukb_merge.txt \
  --make-bed \
  --out data/ukb_cal_pgen/ukb_cal_tmp

plink2 --bfile data/ukb_cal_pgen/ukb_cal_tmp \
  --make-pgen 'vzs' \
  --out data/ukb_cal_pgen/ukb_cal_all
