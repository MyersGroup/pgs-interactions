## Aggregate all independent pairwise hits

library(data.table)
setDTthreads(8)
library(dplyr)
library(pgenlibr)
library(stringr)
library(tidyr)

options(warn=2)  # turn warnings into errors

samples <- "wb_all"

source("scripts/misc/fn-load_geno.R")



## load sample file
sample_f <- fread("../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam", data.table = FALSE)

## load list of phenotypes + short descriptions
names <- fread("../data/phenotypes/clean/code-desc-map-real.tab", data.table = FALSE)
phen_codes_qn <- names$field_id[1:97]



## structure of data.frame with independent pairwise hits
indep_hits_empty <- data.frame(ID_1         = character(),
                               A1_1         = character(),
                               MAF_1        = numeric(),
                               BETA_1       = numeric(),
                               LOG10_P_1    = numeric(),
                               BETA_sq_1    = numeric(),
                               LOG10_P_sq_1 = numeric(),
                               ID_2         = character(),
                               A1_2         = character(),
                               MAF_2        = numeric(),
                               BETA_2       = numeric(),
                               LOG10_P_2    = numeric(),
                               BETA_sq_2    = numeric(),
                               LOG10_P_sq_2 = numeric(),
                               BETA_int     = numeric(),
                               LOG10_P_int  = numeric(),
                               h2           = numeric(),
                               LD           = numeric())



indep_hits <- data.frame(field_id   = character(),
                         short_desc = character(),
                         PGS        = character(),
                         indep_hits_empty)

for (phen in phen_codes_qn) {

    phen_short_desc <- names$short_desc[names$field_id == phen]

    ## LOCO
    indep_loco_f <- paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                           phen, "-", samples, "-loco-indep-hits.tab")

    if (file.size(indep_loco_f) > 0) {

        indep_loco <- fread(indep_loco_f, data.table = FALSE)

        for (snp in indep_loco$ID) {

            snp_ <- gsub(":", "_", snp)
            indep_snp_f <- paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                                  phen, "-", samples, "-loco-", snp_, "-pairwise-hits.tab")

            if (file.exists(indep_snp_f)) {
                if (file.size(indep_snp_f) > 0) {

                    indep_snp <- fread(indep_snp_f, data.table = FALSE)
                    indep_snp_add <- data.frame(field_id   = phen,
                                                short_desc = phen_short_desc,
                                                PGS = rep("LOCO", nrow(indep_snp)),
                                                indep_snp)
                    indep_hits <- rbind(indep_hits, indep_snp_add)

                }
            }
        }
    }


    ## Full
    indep_full_f <- paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                           phen, "-", samples, "-full-indep-hits.tab")

    if (file.size(indep_full_f) > 0) {

        indep_full <- fread(indep_full_f, data.table = FALSE)

        for (snp in indep_full$ID) {

            snp_ <- gsub(":", "_", snp)
            indep_snp_f <- paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                                  phen, "-", samples, "-full-", snp_, "-pairwise-hits.tab")

            if (file.exists(indep_snp_f)) {
                if (file.size(indep_snp_f) > 0) {

                    indep_snp <- fread(indep_snp_f, data.table = FALSE)
                    indep_snp_add <- data.frame(field_id   = phen,
                                                short_desc = phen_short_desc,
                                                PGS = rep("Full", nrow(indep_snp)),
                                                indep_snp)
                    indep_hits <- rbind(indep_hits, indep_snp_add)

                }
            }
        }
    }
}



## add chromosome fields
indep_hits$CHR_1 <- sapply(strsplit(indep_hits$ID_1, split = ":"), "[[", 1)
indep_hits$CHR_2 <- sapply(strsplit(indep_hits$ID_2, split = ":"), "[[", 1)



## add original GWAS summary stats
indep_hits$BETA_1_orig    <- NA
indep_hits$LOG10_P_1_orig <- NA
indep_hits$BETA_2_orig    <- NA
indep_hits$LOG10_P_2_orig <- NA

for (phen in unique(indep_hits$field_id)) {

    ## get SNP IDs
    snps <- data.frame(ID = c(indep_hits$ID_1[indep_hits$field_id == phen],
                              indep_hits$ID_2[indep_hits$field_id == phen]),
                       side = c(rep("LHS", length(indep_hits$ID_1[indep_hits$field_id == phen])),
                                rep("RHS", length(indep_hits$ID_2[indep_hits$field_id == phen]))))
    ## remove duplicates
    snps <- snps[!duplicated(snps),]
    ## sort by chromosome
    snps$CHR <- as.numeric(sapply(strsplit(snps$ID, split = ":"), "[[", 1))
    snps <- snps %>% arrange(CHR)
    
    ## load and add summary stats
    for (chr in unique(snps$CHR)) {

        orig <- fread(paste0("../results/01-gwas/plink-output/", phen, "/",
                             samples, ".chr", chr, ".", phen, ".res.cov.glm.linear.gz"),
                      data.table = FALSE)

        snps_chr <- snps[snps$CHR == chr,]

        for (i in 1:nrow(snps_chr)) {

            snp <- snps_chr$ID[i]
            snp_side <- snps_chr$side[i]

            if (snp_side == "LHS") {

                indep_hits$BETA_1_orig[indep_hits$field_id == phen &
                                       indep_hits$ID_1 == snp] <- orig$BETA[orig$ID == snp]
                indep_hits$LOG10_P_1_orig[indep_hits$field_id == phen &
                                          indep_hits$ID_1 == snp] <- orig$LOG10_P[orig$ID == snp]

            } else {

                indep_hits$BETA_2_orig[indep_hits$field_id == phen &
                                       indep_hits$ID_2 == snp] <- orig$BETA[orig$ID == snp]
                indep_hits$LOG10_P_2_orig[indep_hits$field_id == phen &
                                          indep_hits$ID_2 == snp] <- orig$LOG10_P[orig$ID == snp]

            }
        }
    }
}



## add count of samples with non-NA phenotypes and with interaction
indep_hits$n_sp_phen       <- NA
indep_hits$n_sp_int        <- NA
indep_hits$freq_int        <- NA
indep_hits$sum_int_al      <- NA
indep_hits$sum_prod_int_al <- NA

for (phen in unique(indep_hits$field_id)) {

    ## load phenotypes (residuals from basic GWAS)
    resid_train <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/",
                                phen, "-", samples, "-resid-covars.tab"),
                     data.table = FALSE)
    resid_train <- resid_train[!is.na(resid_train[[paste0(phen, ".res.cov")]]),]

    ## add number of non-NA samples
    indep_hits$n_sp_phen[indep_hits$field_id == phen] <- nrow(resid_train)


    ## count samples with interaction
    for (lhs_snp in unique(indep_hits$ID_1[indep_hits$field_id == phen])) {

        ## load genotypes
        geno_lhs <- load_geno(lhs_snp,
                              unique(indep_hits$A1_1[indep_hits$ID_1 == lhs_snp]),
                              resid_train$IID)
        
        for (rhs_snp in indep_hits$ID_2[indep_hits$field_id == phen &
                                        indep_hits$ID_1 == lhs_snp]) {

            ## load genotypes
            geno_rhs <- load_geno(rhs_snp,
                                  indep_hits$A1_2[indep_hits$field_id == phen &
                                                  indep_hits$ID_1 == lhs_snp &
                                                  indep_hits$ID_2 == rhs_snp],
                                  resid_train$IID)

            ## count of samples with non-zero genotypes at both positions
            geno_lhs_rhs <- left_join(geno_lhs, geno_rhs, by = "IID")
            geno_lhs_rhs <- geno_lhs_rhs[complete.cases(geno_lhs_rhs),]
            indep_hits$n_sp_int[indep_hits$field_id == phen &
                                indep_hits$ID_1 == lhs_snp &
                                indep_hits$ID_2 == rhs_snp] <- sum(geno_lhs_rhs[[lhs_snp]] > 0 &
                                                                   geno_lhs_rhs[[rhs_snp]] > 0)

            ## frequency computed as mean of the product of the two genotypes
            indep_hits$freq_int[indep_hits$field_id == phen &
                                indep_hits$ID_1 == lhs_snp &
                                indep_hits$ID_2 == rhs_snp] <- mean(geno_lhs_rhs[[lhs_snp]] * geno_lhs_rhs[[rhs_snp]])

            ## sum of genotypes at both positions
            geno_lhs_rhs <- geno_lhs_rhs[geno_lhs_rhs[[lhs_snp]] > 0 &
                                         geno_lhs_rhs[[rhs_snp]] > 0,]
            indep_hits$sum_int_al[indep_hits$field_id == phen &
                                  indep_hits$ID_1 == lhs_snp &
                                  indep_hits$ID_2 == rhs_snp] <- sum(colSums(geno_lhs_rhs[, 2:3]))
            indep_hits$sum_prod_int_al[indep_hits$field_id == phen &
                                       indep_hits$ID_1 == lhs_snp &
                                       indep_hits$ID_2 == rhs_snp] <- sum(geno_lhs_rhs[, 2] * geno_lhs_rhs[, 3])

        }
    }
}



## count duplicates: SNP_A interacts with SNP_B, SNP_B interacts with SNP_A
## create temporary column to identify duplicates
indep_hits <- indep_hits %>%
    rowwise() %>%
    mutate(test_dup = paste0(sort(c(ID_1, ID_2)), collapse = "-")) %>%
    mutate(test_dup = paste(field_id, test_dup, sep = "-"))
## count duplicates
indep_hits <- indep_hits %>%
    group_by(test_dup) %>%
    mutate(count = n()) %>%
    ungroup()
## make sure there are no triplicates
if (max(indep_hits$count > 2)) {
    stop("There are triplicate field_id-ID_1-ID_2 strings.")
}
## turn duplicate count into new indicator column
indep_hits <- indep_hits %>%
    mutate(reciprocal = case_when(count == 2 ~ TRUE,
                                  count == 1 ~ FALSE,
                                  TRUE ~ NA)) %>%
    select(-test_dup, -count) %>%
    as.data.frame()



## add rsIDs
map_posid_rsid_f <- paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", seq(1, 22), ".tab")
map_posid_rsid <- lapply(map_posid_rsid_f, function(x) fread(x, data.table = FALSE))
map_posid_rsid <- do.call("rbind", map_posid_rsid)
## LHS 
map_posid_rsid_add <- map_posid_rsid
colnames(map_posid_rsid_add) <- c("ID_1", "rsID_1")
indep_hits <- left_join(indep_hits, map_posid_rsid_add, by = "ID_1")
## RHS 
colnames(map_posid_rsid_add) <- c("ID_2", "rsID_2")
indep_hits <- left_join(indep_hits, map_posid_rsid_add, by = "ID_2")



## add external annotations
## VEP: (most severe) variant consequence
vep_f <- paste0("../data/annotations/vep/anno/vep-chr", seq(1, 22), ".most_severe.tab.gz")
vep <- lapply(vep_f, function(x) fread(x, skip = "#Uploaded_variation", data.table = FALSE))
vep <- do.call("rbind", vep)
vep <- vep[, c("#Uploaded_variation", "Consequence")]
colnames(vep) <- c("ID", "VEP_Consequence")
## LHS
vep_add <- vep
colnames(vep_add) <- c("ID_1", "VEP_Consequence_1")
indep_hits <- left_join(indep_hits, vep_add, by = "ID_1")
## RHS
colnames(vep_add) <- c("ID_2", "VEP_Consequence_2")
indep_hits <- left_join(indep_hits, vep_add, by = "ID_2")


## Annovar: affected/closest gene(s) from Gencode and refSeq
## Gencode
annovar_var_fn_f <- paste0("../data/annotations/annovar/gene-anno/annovar-chr", seq(1, 22), ".ensGene.variant_function.gz")
annovar_var_fn <- lapply(annovar_var_fn_f,
                         function(x) fread(x,
                                           header = FALSE,
                                           sep = "\t",
                                           col.names = c("ensGene_Func", "ensGene_Gene", "Comments"),
                                           data.table = FALSE))
annovar_var_fn <- do.call("rbind", annovar_var_fn)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$Comments, split = "comments: "), `[[`, 2)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$ID, split = ","), `[[`, 1)
## LHS
annovar_add <- annovar_var_fn[, c("ID", "ensGene_Func", "ensGene_Gene")]
colnames(annovar_add) <- c("ID_1", "ensGene_Func_1", "ensGene_Gene_1")
indep_hits <- left_join(indep_hits, annovar_add[, c("ID_1", "ensGene_Func_1", "ensGene_Gene_1")], by = "ID_1")
## RHS
colnames(annovar_add) <- c("ID_2", "ensGene_Func_2", "ensGene_Gene_2")
indep_hits <- left_join(indep_hits, annovar_add[, c("ID_2", "ensGene_Func_2", "ensGene_Gene_2")], by = "ID_2")
rm(annovar_var_fn, annovar_add)

## RefSeq
annovar_var_fn_f <- paste0("../data/annotations/annovar/gene-anno/annovar-chr", seq(1, 22), ".refGene.variant_function.gz")
annovar_var_fn <- lapply(annovar_var_fn_f,
                         function(x) fread(x,
                                           header = FALSE,
                                           sep = "\t",
                                           col.names = c("refGene_Func", "refGene_Gene", "Comments"),
                                           data.table = FALSE))
annovar_var_fn <- do.call("rbind", annovar_var_fn)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$Comments, split = "comments: "), `[[`, 2)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$ID, split = ","), `[[`, 1)
## LHS
annovar_add <- annovar_var_fn[, c("ID", "refGene_Func", "refGene_Gene")]
colnames(annovar_add) <- c("ID_1", "refGene_Func_1", "refGene_Gene_1")
indep_hits <- left_join(indep_hits, annovar_add[, c("ID_1", "refGene_Func_1", "refGene_Gene_1")], by = "ID_1")
## RHS
colnames(annovar_add) <- c("ID_2", "refGene_Func_2", "refGene_Gene_2")
indep_hits <- left_join(indep_hits, annovar_add[, c("ID_2", "refGene_Func_2", "refGene_Gene_2")], by = "ID_2")
rm(annovar_var_fn, annovar_add)



## round frequencies and sums of interaction alleles
indep_hits$freq_int        <- round(indep_hits$freq_int, digits = 8)
indep_hits$sum_int_al      <- round(indep_hits$sum_int_al, digits = 4)
indep_hits$sum_prod_int_al <- round(indep_hits$sum_prod_int_al, digits = 4)



## expected vs. observed number of samples with interaction
indep_hits$sum_prod_int_al_exp        <- 4 * indep_hits$MAF_1 * indep_hits$MAF_2 * indep_hits$n_sp_phen
indep_hits$sum_prod_int_al_rt_obs_exp <- indep_hits$sum_prod_int_al / indep_hits$sum_prod_int_al_exp
## round
indep_hits$sum_prod_int_al_exp        <- round(indep_hits$sum_prod_int_al_exp, digits = 4)
indep_hits$sum_prod_int_al_rt_obs_exp <- round(indep_hits$sum_prod_int_al_rt_obs_exp, digits = 8)



## export
indep_hits <- indep_hits %>%
    select(field_id, short_desc, n_sp_phen,
           CHR_1, ID_1, rsID_1, A1_1, MAF_1, PGS, VEP_Consequence_1, ensGene_Func_1, ensGene_Gene_1, refGene_Func_1, refGene_Gene_1,
           CHR_2, ID_2, rsID_2, A1_2, MAF_2,      VEP_Consequence_2, ensGene_Func_2, ensGene_Gene_2, refGene_Func_2, refGene_Gene_2,
           LD, 
           BETA_1_orig, LOG10_P_1_orig, BETA_1, LOG10_P_1, BETA_sq_1, LOG10_P_sq_1,
           BETA_2_orig, LOG10_P_2_orig, BETA_2, LOG10_P_2, BETA_sq_2, LOG10_P_sq_2,
           BETA_int, LOG10_P_int,
           h2, reciprocal, n_sp_int, freq_int, sum_int_al, sum_prod_int_al,
           sum_prod_int_al_exp, sum_prod_int_al_rt_obs_exp)
fwrite(indep_hits,
       "../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-all.tab",
       sep = "\t", na = "NA", quote = FALSE)
## save data.frame in RData file
saveRDS(indep_hits,
        file = "../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-all.RData")



## replace short description with longer description
names_clean_qn <- fread("../data/phenotypes/clean/code-desc-map-real", data.table = FALSE)
colnames(names_clean_qn)[4] <- "desc_clean"
names_clean_qn <- names_clean_qn[, c("field_id", "desc_clean")]
names_clean_reg <- names_clean_qn
names_clean_qn$desc_clean <- paste0(names_clean_qn$desc_clean, " (QN)")
names_clean_reg$field_id <- substr(names_clean_reg$field_id, 1, nchar(names_clean_reg$field_id) - 3)
names_clean <- rbind(names_clean_qn, names_clean_reg)
indep_hits <- left_join(indep_hits, names_clean[, c("field_id", "desc_clean")], by = "field_id")

## further filtering
indep_hits <- indep_hits %>%
    select(field_id, desc_clean,
           CHR_1, ID_1, rsID_1, A1_1, MAF_1, PGS, VEP_Consequence_1, ensGene_Func_1, ensGene_Gene_1, refGene_Func_1, refGene_Gene_1,
           CHR_2, ID_2, rsID_2, A1_2, MAF_2,      VEP_Consequence_2, ensGene_Func_2, ensGene_Gene_2, refGene_Func_2, refGene_Gene_2,
           BETA_1_orig, LOG10_P_1_orig, BETA_1, LOG10_P_1, BETA_sq_1, LOG10_P_sq_1,
           BETA_2_orig, LOG10_P_2_orig, BETA_2, LOG10_P_2, BETA_sq_2, LOG10_P_sq_2,
           BETA_int, LOG10_P_int,
           freq_int)

## "high-confidence set": GWS & LOCO & off-chr & both MAF>=0.01
indep_hits_hc <- indep_hits[indep_hits$LOG10_P_int >= -log10(5e-8) &
                            indep_hits$PGS == "LOCO" &
                            (indep_hits$CHR_1 != indep_hits$CHR_2) &
                            (indep_hits$MAF_1 >= 0.01 & indep_hits$MAF_2 >= 0.01),]
## > dim(indep_hits_hc)
## [1] 38 38
fwrite(indep_hits_hc,
       "../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-high_conf.tab",
       sep = "\t", na = "NA", quote = FALSE)

## others: GWS & LOCO & off-chr & interaction frequency >= 0.001 (& not in high-confidence set)
indep_hits_others <- indep_hits[indep_hits$LOG10_P_int >= -log10(5e-8) &
                                indep_hits$PGS == "LOCO" &
                                (indep_hits$CHR_1 != indep_hits$CHR_2) &
                                (indep_hits$freq_int >= 0.001 &
                                (indep_hits$MAF_1 < 0.01 | indep_hits$MAF_2 < 0.01)),]
## > dim(indep_hits_others)
## [1] 75 38
fwrite(indep_hits_others,
       "../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-others.tab",
       sep = "\t", na = "NA", quote = FALSE)
