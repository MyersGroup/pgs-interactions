## desc: Prepare bed/bim/fam files for running SME

library(argparser)
library(data.table)
setDTthreads(4)
library(dplyr)

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen", help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--snp",  help = "Target SNP",     nargs = 1)
argv <- parse_args(p)

phen <- argv$phen
snp_ <- argv$snp
snp  <- gsub("_", ":", snp_)

samples <- "wb_all"

plink2 <- "/path/to/plink2"  # Plink 2.0
plink1 <- "/path/to/plink"   # Plink 1.9



## prepare bed files with 1000 PGS SNPs + 22 LHS SNPs + SNPs in inferred PGSs in simulated traits
if (startsWith(phen, "f.sim")) {
    ## load coefficients from null simulation with 1000 SNPs and no clustering/transformation
    null_sim_f <- paste0("../data/phenotypes/simulations-prep/coeff/sim-coeff-1k_rg-chr", seq(1, 22), ".tab")
    null_sim   <- lapply(null_sim_f, function(x) fread(x, data.table = FALSE))
    null_sim   <- do.call("rbind", null_sim)

    ## load 22 LHS SNPs
    inter <- fread("../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab", data.table = FALSE)
}

## load coefficients of trait's iterative PGS
r2 <- 0.9; kb <- 500
coeff <- list()
for (chr in 1:22) {
    coeff_f <- paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", samples,
                      "/coeff/final/", phen, "-coeff-chr", chr, ".tab")
    if (file.exists(coeff_f)) {
        coeff[[chr]] <- fread(coeff_f, data.table = FALSE)
    }
}
coeff <- do.call("rbind", coeff)

## aggregate all SNPs to add
if (startsWith(phen, "f.sim")) {
    add_df <- data.frame(CHR = NA,
                         ID  = unique(c(null_sim$ID, inter$ID.1, coeff$ID, snp)))
} else {
    add_df <- data.frame(CHR = NA,
                         ID  = unique(c(coeff$ID, snp)))
}
add_df <- add_df[!duplicated(add_df),]
add_df$CHR <- as.numeric(sapply(strsplit(add_df$ID, split = ":"), `[[`, 1))

## load SNPs in filtered bed files
bim_f <- paste0("../results/05-new-app/sme/genotype-calls/ukb_imputed_filtered_chr", seq(1, 22), ".bim")
bim <- lapply(bim_f, function(x) fread(x, col.names = c("CHR", "ID", "cM", "POS", "REF", "ALT"), data.table = FALSE))
bim <- do.call("rbind", bim)

## remove those that are already in bed files
add_df <- add_df[!(add_df$ID %in% bim$ID),]

for (chr in 1:22) {

    ## write filtered set of variants to add from imputed data
    add_df_chr <- add_df[add_df$CHR == chr,]

    ## if no new variants to add, copy original filtered bed file for concatenation below
    if (nrow(add_df_chr) == 0) {
        for (ext in c("bed", "bim", "fam")) {
            system(paste0("cp ../results/05-new-app/sme/genotype-calls/ukb_imputed_filtered_chr", chr, ".", ext, " ",
                          "../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_chr", chr, ".", ext))
        }
        next
    }
    
    fwrite(list(add_df_chr$ID),
           file = paste0("../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_to_add_chr", chr, ".tab"),
           row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')

    ## make bed file with variants to add
    system(paste0(plink2,
                  " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                  "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                  "--psam ../data/sample-ids/ukb-103076-imp-auto-s486989.psam ",
                  "--keep ../data/sample-ids/filtered/wb_all-ids.tab ",
                  "--extract ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_to_add_chr", chr, ".tab ",
                  "--make-bed ",
                  "--threads 4 ",
                  "--memory 60000 ",
                  "--out ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_to_add_chr", chr))

    ## concatenate files with filtered and to-be-added imputed variants
    system(paste0(plink1,
                  " --bfile ../results/05-new-app/sme/genotype-calls/ukb_imputed_filtered_chr", chr, " ",
                  "--bmerge ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_to_add_chr", chr, " ",
                  "--threads 4 ",
                  "--memory 60000 ",
                  "--make-bed ",
                  "--out ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_chr", chr))
}


## concatenate resulting files across chromosomes
## make list of files to concatenate
bed_files <- paste0("../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_chr", seq(1, 22))
fwrite(list(bed_files),
       file = paste0("../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_all_files.tab"),
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')

## merge
system(paste0(plink2,
              " --pmerge-list ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_all_files.tab bfile ",
              "--delete-pmerge-result ",
              "--make-bed ",
              "--threads 4 ",
              "--memory 60000 ",
              "--out ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_all_chr"))

## delete temporary files
system(paste0("rm ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_to_add_chr*.*"))
system(paste0("rm ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_chr*.*"))
system(paste0("rm ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_imputed_concat_all_files.tab"))
