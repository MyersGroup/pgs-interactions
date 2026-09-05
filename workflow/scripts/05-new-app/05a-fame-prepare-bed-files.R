## desc: Prepare bed/bim/fam files for running FAME

library(data.table)
setDTthreads(4)
library(dplyr)

plink2 <- "/path/to/plink2"  # Plink 2.0
plink1 <- "/path/to/plink"   # Plink 1.9


## filter array genotypes using same filters as Fu et al.
## - missingness <= 1%
## - MAF >= 1%
## - HWE p-value > 1e-7
for (chr in 1:22) {
    system(paste0(plink2,
              " --bed ../data/genotype-calls/ukb_cal_chr", chr, "_v2.bed ",
              "--bim ../data/genotype-calls/ukb_snp_chr", chr, "_v2.bim ",
              "--fam ../data/sample-ids/ukb22418_c1_b0_v2_s487957.fam ",
              "--keep ../data/sample-ids/filtered/wb_all-ids.tab ",
              "--geno 0.01 ",
              "--maf 0.01 ",
              "--hwe 1e-7 0 ",
              "--make-bed ",
              "--threads 4 ",
              "--memory 60000 ",
              "--out ../results/05-new-app/fame/genotype-calls/ukb_array_filtered_chr", chr))
}


## prepare bed files with 1000 PGS SNPs + 22 LHS SNPs in simulated traits + 144 SNP*PGS hit SNPs

## load coefficients from null simulation with 1000 SNPs and no clustering/transformation
null_sim_f <- paste0("../data/phenotypes/simulations-prep/coeff/sim-coeff-1k_rg-chr", seq(1, 22), ".tab")
null_sim   <- lapply(null_sim_f, function(x) fread(x, data.table = FALSE))
null_sim   <- do.call("rbind", null_sim)

## load 22 LHS SNPs
inter <- fread("../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab", data.table = FALSE)
inter_lhs <- unique(inter$ID.1)

## load SNP*PGS interaction hits
hits <- fread("../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab", data.table = FALSE)
## sort by physical position
hits <- hits %>% arrange(CHR, POS)
## keep only unique hits
hits <- hits[, c("CHR", "POS", "ID", "rsID")]
hits <- hits[!duplicated(hits),]

## aggregate all SNPs to add
add_df <- data.frame(CHR = NA,
                     ID  = c(null_sim$ID, inter_lhs, hits$ID))
add_df <- add_df[!duplicated(add_df),]                       
add_df$CHR <- as.numeric(sapply(strsplit(add_df$ID, split = ":"), `[[`, 1))
add_df$POS <- as.numeric(sapply(strsplit(add_df$ID, split = ":"), `[[`, 2))
add_df$CHR_POS <- paste0(add_df$CHR, "_", add_df$POS)

## load SNPs in filtered bed files
bim_f <- paste0("../results/05-new-app/fame/genotype-calls/ukb_array_filtered_chr", seq(1, 22), ".bim")
bim <- lapply(bim_f, function(x) fread(x, col.names = c("CHR", "rsID", "cM", "POS", "REF", "ALT"), data.table = FALSE))
bim <- do.call("rbind", bim)
bim$CHR_POS <- paste0(bim$CHR, "_", bim$POS)

## remove those that are already in bed files
## (checked that it makes no difference to remove by position only or by full ID)
add_df <- add_df[!(add_df$CHR_POS %in% bim$CHR_POS),]

for (chr in 1:22) {

    ## write filtered set of variants to add from imputed data
    add_df_chr <- add_df[add_df$CHR == chr,]
    fwrite(list(add_df_chr$ID),
           file = paste0("../results/05-new-app/fame/genotype-calls/ukb_imputed_to_add_chr", chr, ".tab"),
           row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')

    ## make bed file
    system(paste0(plink2,
                  " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                  "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                  "--psam ../data/sample-ids/ukb-103076-imp-auto-s486989.psam ",
                  "--keep ../data/sample-ids/filtered/wb_all-ids.tab ",
                  "--extract ../results/05-new-app/fame/genotype-calls/ukb_imputed_to_add_chr", chr, ".tab ",
                  "--make-bed ",
                  "--threads 4 ",
                  "--memory 60000 ",
                  "--out ../results/05-new-app/fame/genotype-calls/ukb_imputed_to_add_chr", chr))
}


## concatenate files with filtered and to-be-added imputed variants by chromosome using Plink 1.9
for (chr in 1:22) {
    system(paste0(plink1,
                  " --bfile ../results/05-new-app/fame/genotype-calls/ukb_array_filtered_chr", chr, " ",
                  "--bmerge ../results/05-new-app/fame/genotype-calls/ukb_imputed_to_add_chr", chr, " ",
                  "--threads 4 ",
                  "--memory 60000 ",
                  "--make-bed ",
                  "--out ../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_chr", chr))
}


## concatenate resulting files across chromosomes
## make list of files to concatenate
bed_files <- paste0("../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_chr", seq(1, 22))
fwrite(list(bed_files),
       file = paste0("../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_all_files.tab"),
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')

## merge
system(paste0(plink2,
              " --pmerge-list ../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_all_files.tab bfile ",
              "--delete-pmerge-result ",
              "--make-bed ",
              "--threads 4 ",
              "--memory 60000 ",
              "--out ../results/05-new-app/fame/genotype-calls/ukb_all_chr"))

## delete temporary files
system("rm ../results/05-new-app/fame/genotype-calls/ukb_array_filtered_chr*")
system("rm ../results/05-new-app/fame/genotype-calls/ukb_imputed_to_add_chr*")
system("rm ../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_chr*")
system("rm ../results/05-new-app/fame/genotype-calls/ukb_array_imputed_concat_all_files.tab")


## prepare table listing phenotype-target SNP pairs for which we should run FAME
traits_sim <- fread("../data/phenotypes/clean/code-desc-map-sim.tab", data.table = FALSE)
traits_sim_int <- traits_sim$field_id[grepl("_qn", traits_sim$field_id) &
                                      grepl("_int_", traits_sim$field_id) &
                                      grepl("_all_", traits_sim$field_id)]

to_run <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(to_run) <- c("phen", "snp")

for (phen in traits_sim_int) {
    to_run <- rbind(to_run,
                    data.frame(phen = phen,
                               snp  = gsub(":", "_", inter_lhs)))
}

hits <- fread("../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab", data.table = FALSE)
for (i in 1:nrow(hits)) {
    to_run <- rbind(to_run,
                    data.frame(phen = hits$field_id[i],
                               snp  = gsub(":", "_", hits$ID[i])))
}

fwrite(to_run,
       file = "../results/05-new-app/fame/phen-snp-key.tab",
       col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')
