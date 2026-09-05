## desc: Filter SNPs from imputed data and convert to bed files

library(data.table)
setDTthreads(4)
library(dplyr)

## a recent version of Plink is needed so --exclude-palindromic-snps flag is available
plink2 <- "/path/to/plink2"



## filter imputed genotypes using the same filters as Stamp et al.
## - INFO score >= 0.8
## - remove structural variants
## - remove strand-ambiguous SNPs
## - missingness = 0%
## - MAF >= 1%
## - HWE p-value > 1e-6
for (chr in 1:22) {
    system(paste(plink2,
              ' --pgen ../data/imputed-genotypes/ukb_imp_chr', chr, '_v3.pgen ',
              '--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr', chr, '_v3.pvar ',
              '--psam ../data/sample-ids/ukb-103076-imp-auto-s486989.psam ',
              '--keep ../data/sample-ids/filtered/wb_all-ids.tab ',
              '--extract-if-info "IIS>=0.8" ',  # imputation quality score >= 0.8
              '--snps-only ',  # no structural variants
              '--exclude-palindromic-snps ',  # no strand-ambiguous SNPs (AT or CG)
              '--geno 0.00 ',  # no missing genotypes
              '--maf 0.01 ',  # MAF >= 0.1
              '--hwe 0.000001 0 ',  # HWE p-value >= 1e-6
              '--make-bed ',
              '--threads 4 ',
              '--memory 60000 ',
              '--out ../results/05-new-app/sme/genotype-calls/ukb_imputed_filtered_chr', chr, sep = ""))
}
