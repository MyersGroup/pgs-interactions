## Aggregate independent hits obtained through stepwise regression

library(data.table)
setDTthreads(4)
library(dplyr)
library(hash)
library(stringr)
library(tidyr)

options(warn=2)  # turn warnings into errors

source("scripts/misc/fn-load_geno.R")

samples <- "wb_all"



## load list of phenotypes + short descriptions
names <- fread("../data/phenotypes/clean/code-desc-map-real.tab", data.table = FALSE)
names_hash <- hash(names$field_id, paste0(names$desc_full, " (QN)"))

phen_codes_qn <- names$field_id[1:97]

## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)



## LOCO hits
hits_loco <- list()
for (phen in phen_codes_qn) {

    ## start with QN results
    hits <- FALSE
    
    if (file.size(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-loco-indep-hits.tab")) > 0) {
        hits <- TRUE

        ## load independent hits
        main <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                             phen, "-", samples, "-loco-indep-hits-stats.tab"),
                      data.table = FALSE)

        ## retrieve coeffs and p-values from Plink output
        main$snp_BETA     <- NA
        main$snp_LOG10_P  <- NA
        main$pgs_BETA     <- NA
        main$pgs_LOG10_P  <- NA
        main$int_BETA     <- NA
        main$int_LOG10_P  <- NA
        main$orig_BETA    <- NA
        main$orig_LOG10_P <- NA
        
        for (chr in unique(main$CHR)) {

            snps_chr <- main$ID[main$CHR == chr]

            ## interaction GWAS        
            gwas_chr <- fread(paste0("../results/03-interaction-gwas/plink-output/", phen, "/snp-pgs/",
                                     phen, ".", samples, ".loco.chr", chr, ".", phen, ".res.cov.pgs_full_mc.ac.ctr.glm.linear.gz"),
                              data.table = FALSE)
            gwas_chr <- gwas_chr[gwas_chr$ID %in% snps_chr,]
            for (snp in snps_chr) {
                ## SNP
                main$snp_BETA[main$ID == snp]    <- gwas_chr$BETA[gwas_chr$ID == snp & gwas_chr$TEST == "ADD"]
                main$snp_LOG10_P[main$ID == snp] <- gwas_chr$LOG10_P[gwas_chr$ID == snp & gwas_chr$TEST == "ADD"]
                ## PGS
                main$pgs_BETA[main$ID == snp]    <- gwas_chr$BETA[gwas_chr$ID == snp & gwas_chr$TEST == paste0("pgs_loco", chr, "_mc")]
                main$pgs_LOG10_P[main$ID == snp] <- gwas_chr$LOG10_P[gwas_chr$ID == snp & gwas_chr$TEST == paste0("pgs_loco", chr, "_mc")]
                ## interaction term
                main$int_BETA[main$ID == snp]    <- gwas_chr$BETA[gwas_chr$ID == snp & gwas_chr$TEST == paste0("ADDxpgs_loco", chr, "_mc")]
                main$int_LOG10_P[main$ID == snp] <- gwas_chr$LOG10_P[gwas_chr$ID == snp & gwas_chr$TEST == paste0("ADDxpgs_loco", chr, "_mc")]
            }
            rm(gwas_chr)

            ## original GWAS
            orig_chr <- fread(paste0("../results/01-gwas/plink-output/", phen, "/", samples, ".chr", chr, ".", phen, ".res.cov.glm.linear.gz"),
                              data.table = FALSE)
            orig_chr <- orig_chr[orig_chr$ID %in% snps_chr,]
            for (snp in snps_chr) {
                main$orig_BETA[main$ID == snp]    <- orig_chr$BETA[orig_chr$ID == snp]
                main$orig_LOG10_P[main$ID == snp] <- orig_chr$LOG10_P[orig_chr$ID == snp]
            }
            rm(orig_chr)
        }

        ## add trait info
        main$field_id <- phen
        main$short_desc <- values(names_hash, phen)

    }

    if (hits) {
        hits_loco[[phen]] <- main
    }
}
hits_loco <- do.call("rbind", hits_loco)

## rearrange and export
hits_loco$PGS <- "LOCO"
hits_loco <- hits_loco %>%
    select(field_id, short_desc,
           CHR, POS, ID, rsID, REF, ALT, A1, PGS,
           snp_BETA, snp_LOG10_P,
           pgs_BETA, pgs_LOG10_P,
           int_BETA, int_LOG10_P,
           orig_BETA, orig_LOG10_P)
colnames(hits_loco)[2] <- "desc_clean"
fwrite(hits_loco,
       file = "../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab",
       sep = "\t", na = "NA", quote = FALSE)



## full
hits_full <- list()
for (phen in phen_codes_qn) {

    ## start with QN results
    hits <- FALSE
    
    if (file.size(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-full-indep-hits.tab")) > 0) {
        hits <- TRUE

        ## load independent hits
        main <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                             phen, "-", samples, "-full-indep-hits-stats.tab"),
                      data.table = FALSE)

        ## retrieve coeffs and p-values from Plink output
        main$snp_BETA     <- NA
        main$snp_LOG10_P  <- NA
        main$pgs_BETA     <- NA
        main$pgs_LOG10_P  <- NA
        main$int_BETA     <- NA
        main$int_LOG10_P  <- NA
        main$orig_BETA    <- NA
        main$orig_LOG10_P <- NA
        
        for (chr in unique(main$CHR)) {

            snps_chr <- main$ID[main$CHR == chr]

            ## interaction GWAS        
            gwas_chr <- fread(paste0("../results/03-interaction-gwas/plink-output/", phen, "/snp-pgs/",
                                     phen, ".", samples, ".full.chr", chr, ".", phen, ".res.cov.pgs_full_mc.ac.ctr.glm.linear.gz"),
                              data.table = FALSE)
            gwas_chr <- gwas_chr[gwas_chr$ID %in% snps_chr,]
            for (snp in snps_chr) {
                ## SNP
                main$snp_BETA[main$ID == snp]    <- gwas_chr$BETA[gwas_chr$ID == snp & gwas_chr$TEST == "ADD"]
                main$snp_LOG10_P[main$ID == snp] <- gwas_chr$LOG10_P[gwas_chr$ID == snp & gwas_chr$TEST == "ADD"]
                ## PGS
                main$pgs_BETA[main$ID == snp]    <- gwas_chr$BETA[gwas_chr$ID == snp & gwas_chr$TEST == paste0("pgs_full_mc")]
                main$pgs_LOG10_P[main$ID == snp] <- gwas_chr$LOG10_P[gwas_chr$ID == snp & gwas_chr$TEST == paste0("pgs_full_mc")]
                ## interaction term
                main$int_BETA[main$ID == snp]    <- gwas_chr$BETA[gwas_chr$ID == snp & gwas_chr$TEST == paste0("ADDxpgs_full_mc")]
                main$int_LOG10_P[main$ID == snp] <- gwas_chr$LOG10_P[gwas_chr$ID == snp & gwas_chr$TEST == paste0("ADDxpgs_full_mc")]
            }
            rm(gwas_chr)

            ## original GWAS
            orig_chr <- fread(paste0("../results/01-gwas/plink-output/", phen, "/", samples, ".chr", chr, ".", phen, ".res.cov.glm.linear.gz"),
                              data.table = FALSE)
            orig_chr <- orig_chr[orig_chr$ID %in% snps_chr,]
            for (snp in snps_chr) {
                main$orig_BETA[main$ID == snp]    <- orig_chr$BETA[orig_chr$ID == snp]
                main$orig_LOG10_P[main$ID == snp] <- orig_chr$LOG10_P[orig_chr$ID == snp]
            }
            rm(orig_chr)
        }

        ## add trait info
        main$field_id <- phen
        main$short_desc <- values(names_hash, phen)

    }

    if (hits) {
        hits_full[[phen]] <- main
    }
}
hits_full <- do.call("rbind", hits_full)

## rearrange
hits_full$PGS <- "Full"
hits_full <- hits_full %>%
    select(field_id, short_desc,
           CHR, POS, ID, rsID, REF, ALT, A1, PGS,
           snp_BETA, snp_LOG10_P,
           pgs_BETA, pgs_LOG10_P,
           int_BETA, int_LOG10_P,
           orig_BETA, orig_LOG10_P)
colnames(hits_full)[2] <- "desc_clean"
fwrite(hits_full,
       file = "../results/03-interaction-gwas/indep-hits/aggregate/hits-full-qn.tab",
       sep = "\t", na = "NA", quote = FALSE)
