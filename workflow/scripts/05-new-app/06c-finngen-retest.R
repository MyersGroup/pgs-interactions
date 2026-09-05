library(data.table)
setDTthreads(2)
library(dplyr)
library(pgenlibr)


## load sample file
sample_f <- fread("../data/sample-ids/ukb-103076-imp-auto-s486989.psam", data.table = FALSE)

## function for loading genotypes
load_geno <- function(var_ids, a1, sample_ids_filtered) {

    ## extract chromosome
    chr <- as.numeric(unique(sapply(strsplit(var_ids, split = ":"), `[[`, 1)))
    if (length(chr) > 1) { stop("Cannot have SNPs from more than one chromosome.") }

    ## PGEN/PVAR files
    f.pgen <- paste0("../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen")
    f.pvar <- paste0("../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar")
    pvar <- pgenlibr::NewPvar(f.pvar)
    pgen <- pgenlibr::NewPgen(f.pgen, pvar=pvar)

    ## load dosages from PGEN
    var_num <- rep(NA_real_, length(var_ids))
    for (i in 1:length(var_ids)) {
        var_num[i] <- pgenlibr::GetVariantsById(pvar = pvar, id = var_ids[i])
    }
    geno <- as.matrix(ReadList(pgen, var_num))

    ## subset to filtered set of sample IIDs
    geno <- geno %>% data.frame() %>% setNames(var_ids)
    geno$IID <- sample_f$IID
    geno <- geno[geno$IID %in% sample_ids_filtered,]

    ## flip allele where needed to match A1
    var_ids_ref <- sapply(strsplit(var_ids, split = ":"), `[[`, 3)
    ref_snps <- (a1 == var_ids_ref)
    if (sum(ref_snps) > 0) {
        ref_snps_cols <- seq(1, ncol(geno) - 1)[ref_snps]
        geno[, ref_snps_cols] <- 2 - geno[, ref_snps_cols]
    }

    ## export
    geno <- geno[, c(ncol(geno), 1:(ncol(geno) - 1))]
    return(geno)
}



## load list of WB sample IIDs
wb_ids <- fread("../data/sample-ids/filtered/wb_all-ids.tab", data.table = FALSE)

## load table of SNP*PGS hits
hits_loco <- fread("../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab", data.table = FALSE)

## load list of FinnGen traits to analyse
lab_traits <- fread("../data/sum-stats/finngen/ukb-finngen-lab-values-match.tab", data.table = FALSE)

## subset to FinnGen traits
main <- hits_loco[hits_loco$field_id %in% lab_traits$field_id,]



## add performance of our own PGS (mean corrected)
main$pgs_R2_train_mc <- NA
main$pgs_R2_test_mc  <- NA

for (i in 1:nrow(main)) {

    phen <- main$field_id[i]

    ## training set
    ## load phenotype
    phen_df <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-wb_all-resid-covars.tab"), data.table = FALSE)
    colnames(phen_df)[3] <- "phen"

    ## load PGS
    pgs_df_mc <- fread(paste0("../results/05-new-app/mean-corrected/", phen, "/r2_0.9-kb_500/", phen, "-wb_all-pgs-full_mc.tab"), data.table = FALSE)

    ## merge 
    phen_pgs_df <- left_join(phen_df[, c("IID", "phen")], pgs_df_mc[, c("IID","pgs_full_mc")], by = "IID")

    ## run regressions
    reg_mc <- lm(phen ~ pgs_full_mc, data = phen_pgs_df)
    main$pgs_R2_train_mc[i] <- summary(reg_mc)$r.squared

    ## test set
    ## load phenotype
    phen_df <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-other_gb_test_all-resid-covars.tab"), data.table = FALSE)
    colnames(phen_df)[3] <- "phen"

    ## load PGS
    pgs_df_mc <- fread(paste0("../results/05-new-app/mean-corrected/", phen, "/r2_0.9-kb_500/", phen, "-other_gb_test_all-pgs-full_mc.tab"), data.table = FALSE)

    ## merge 
    phen_pgs_df <- left_join(phen_df[, c("IID", "phen")], pgs_df_mc[, c("IID","pgs_full_mc")], by = "IID")

    ## run regressions
    reg_mc <- lm(phen ~ pgs_full_mc, data = phen_pgs_df)
    main$pgs_R2_test_mc[i] <- summary(reg_mc)$r.squared
}



## add performance of FinnGen PGS (not mean corrected)
main$fg_pgs_R2_train <- NA
main$fg_pgs_R2_test  <- NA

for (i in 1:nrow(main)) {

    phen <- main$field_id[i]
    omopid <- lab_traits$OMOPID[lab_traits$field_id == phen]

    ## training set
    ## load phenotype
    phen_df <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-wb_all-resid-covars.tab"), data.table = FALSE)
    colnames(phen_df)[3] <- "phen"

    ## load PGS
    pgs_df_fg_full <- fread(paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-all-pgs-full.tab"), data.table = FALSE)

    ## merge 
    phen_pgs_df <- left_join(phen_df[, c("IID", "phen")], pgs_df_fg_full[, c("IID","pgs_full")], by = "IID")

    ## run regressions
    reg <- lm(phen ~ pgs_full, data = phen_pgs_df)
    main$fg_pgs_R2_train[i] <- summary(reg)$r.squared

    ## test set
    ## load phenotype
    phen_df <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-other_gb_test_all-resid-covars.tab"), data.table = FALSE)
    colnames(phen_df)[3] <- "phen"

    ## merge 
    phen_pgs_df <- left_join(phen_df[, c("IID", "phen")], pgs_df_fg_full[, c("IID","pgs_full")], by = "IID")

    ## run regressions
    reg <- lm(phen ~ pgs_full, data = phen_pgs_df)
    main$fg_pgs_R2_test[i] <- summary(reg)$r.squared
}



## retest interactions using FinnGen PGS
main$fg_snp_BETA       <- NA
main$fg_snp_LOG10_P    <- NA
main$fg_snp_sq_BETA    <- NA
main$fg_snp_sq_LOG10_P <- NA
main$fg_pgs_BETA       <- NA
main$fg_pgs_LOG10_P    <- NA
main$fg_int_BETA       <- NA
main$fg_int_LOG10_P    <- NA

for (i in 1:nrow(main)) {

    phen <- main$field_id[i]
    omopid <- lab_traits$OMOPID[lab_traits$field_id == phen]
    snp_id <- main$ID[i]
    snp_a1 <- main$A1[i]
    snp_chr <- main$CHR[i]

    ## load simple INT phenotype
    phen_df <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-wb_all-resid-covars.tab"),
                     data.table = FALSE)
    colnames(phen_df)[3] <- "phen"

    ## load FinnGen PGS
    pgs_df_fg_loco <- fread(paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-all-pgs-loco.tab"), data.table = FALSE)
    pgs_df_fg_loco <- pgs_df_fg_loco[, c("IID", paste0("pgs_loco", snp_chr))]
    colnames(pgs_df_fg_loco)[2] <- "pgs_loco"

    ## join
    reg_df <- left_join(phen_df[, c("IID", "phen")], pgs_df_fg_loco, by = "IID")

    ## load genotypes
    geno <- load_geno(snp_id, snp_a1, wb_ids$IID)
    colnames(geno)[2] <- "snp"
    ## check whether correlation between genotype and its square is < 0.999
    if (!(cor(geno[, 2], geno[, 2]^2) < 0.999)) {
        stop("Correlation between SNP genotype and its square too high.")
    }

    ## join
    reg_df <- left_join(reg_df, geno, by = "IID")

    ## run regression
    reg <- lm(phen ~ snp*pgs_loco + I(snp^2), data = reg_df)
    coeffs <- summary(reg)$coeff
    
    ## add to table
    ## SNP
    main$fg_snp_BETA[i]       <- coeffs["snp", 1]
    main$fg_snp_LOG10_P[i]    <- - (pt(abs(coeffs["snp", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
    ## SNP^2
    main$fg_snp_sq_BETA[i]    <- coeffs["I(snp^2)", 1]
    main$fg_snp_sq_LOG10_P[i] <- - (pt(abs(coeffs["I(snp^2)", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
    ## PGS
    main$fg_pgs_BETA[i]       <- coeffs["pgs_loco", 1]
    main$fg_pgs_LOG10_P[i]    <- - (pt(abs(coeffs["pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
    ## interaction term
    main$fg_int_BETA[i]       <- coeffs["snp:pgs_loco", 1]
    main$fg_int_LOG10_P[i]    <- - (pt(abs(coeffs["snp:pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
}

## save
saveRDS(main,
        file = "../results/05-new-app/finngen/finngen-retest-interactions.rds")
