## Do stepwise regression to find independent SNPs among SNP*PGS interaction hits

library(argparser)
library(data.table)
library(dplyr)
library(pgenlibr)
library(RcppArmadillo)
library(stringr)

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples", help = "Sample code",    nargs = 1)
p <- add_argument(p, "--chr",     help = "Chromosome",     nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples
p1      <- 5e-8



## structure of main code below:
##
## LOCO:
##
## if (nrow(main_loco_filt) > 0) {  # if there are any GWS LOCO SNPs
##     if (ncol(geno_loco) == 1) {  # if there's only one GWS SNP
##     } else {  # if there's more than one hit, run stepwise regression
##     }
## }
##
## Full:
##
## if (nrow(main_full_filt) > 0) {  # if there are any GWS full SNPs
##     if (length(indep_hits_loco_id[[chr]]) > 0) {  # if there are LOCO hits (from above)
##         if (nrow(main_full_filt) > 0) {  # if there are full hits which are *not* independent LOCO hits
##             if (ncol(geno_full) == 1) {  # if there is only one full hit, check if it's significant after including LOCO hits
##             } else {  # if there's more than one hit, run stepwise regression
##             }
##         }
##     } else {
##         if (ncol(geno_full) == 1) {  # if there's only one GWS SNP
##         } else {  # if there's more than one hit, run stepwise regression
##         } 
##     }
## }



## load phenotype
phen_df <- fread(paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_0.9-kb_500/",
                        phen, "-", samples, "-resid-covars-pgs-full_mc-ac-centre.tab"),
                 data.table = FALSE)

## load sample file
sample_f <- fread("../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam", data.table = FALSE)

## load PGS
f_pgs_loco <- paste0("../results/02-pgs/mean-corrected/", phen, "/r2_0.9-kb_500/",
                phen, "-", samples, "-pgs-loco", seq(1, 22), "_mc.tab")
pgs_df_loco <- lapply(f_pgs_loco, function(x) fread(x, data.table = FALSE))

pgs_df_full <- fread(paste0("../results/02-pgs/mean-corrected/", phen, "/r2_0.9-kb_500/",
                            phen, "-", samples, "-pgs-full_mc.tab"),
                     data.table = FALSE)

## confirm same order
phen_df <- phen_df[order(match(phen_df$IID, sample_f$IID)),]
if (!all.equal(sample_f$IID[sample_f$IID %in% phen_df$IID], phen_df$IID)) {
    stop("Samples not in same order in phenotype file and sample file.")
}

for (chr in 1:22) {
    pgs_df_loco[[chr]] <- pgs_df_loco[[chr]][order(match(pgs_df_loco[[chr]]$IID, sample_f$IID)),]
    if (!all.equal(sample_f$IID[sample_f$IID %in% pgs_df_loco[[chr]]$IID], pgs_df_loco[[chr]]$IID)) {
        stop("Samples not in same order in PGS file and sample file.")
    }
}

pgs_df_full <- pgs_df_full[order(match(pgs_df_full$IID, sample_f$IID)),]
if (!all.equal(sample_f$IID[sample_f$IID %in% pgs_df_full$IID], pgs_df_full$IID)) {
    stop("Samples not in same order in PGS file and sample file.")
}


indep_hits_loco_id  <- vector("list", length = 22)
indep_hits_loco_ind <- vector("list", length = 22)
indep_hits_full_id  <- vector("list", length = 22)
indep_hits_full_ind <- vector("list", length = 22)

## lists to store additive effect + dominance + SNP*PGS test results
indep_hits_loco_add_sign_LOG10_P <- vector("list", length = 22)
indep_hits_full_add_sign_LOG10_P <- vector("list", length = 22)
indep_hits_loco_dom_sign_LOG10_P <- vector("list", length = 22)
indep_hits_full_dom_sign_LOG10_P <- vector("list", length = 22)
indep_hits_loco_snp_pgs_sign_LOG10_P <- vector("list", length = 22)
indep_hits_full_snp_pgs_sign_LOG10_P <- vector("list", length = 22)

## create directories for sample ID lists
if (!dir.exists(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/loco/"))) {
    dir.create(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/loco/"))
}
if (!dir.exists(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/full/"))) {
    dir.create(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/full/"))
}

main_loco <- list()
main_full <- list()

for (chr in 1:22) {

    ## LOCO ---------------------------------------------------------------------
    
    ## load PLINK2 SNP*PGS interaction GWAS output
    main_loco[[chr]] <- fread(paste0("../results/03-interaction-gwas/plink-output/", phen, "/snp-pgs/",
                                     phen, ".", samples, ".loco.chr", chr, ".", phen, ".res.cov.pgs_full_mc.ac.ctr.glm.linear.gz"),
                              data.table = FALSE)
    main_loco[[chr]] <- main_loco[[chr]][main_loco[[chr]]$TEST == paste0("ADDxpgs_loco", chr, "_mc"),]
    colnames(main_loco[[chr]])[1] <- "CHR"

    ## filter by p-value
    main_loco_filt <- main_loco[[chr]] %>%
        filter(LOG10_P >= -log10(p1)) %>%
        arrange(-LOG10_P)

    ## imputed genotypes
    f.pgen = paste0("../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen")
    f.pvar = paste0("../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar")

    pvar <- pgenlibr::NewPvar(f.pvar)
    pgen <- pgenlibr::NewPgen(f.pgen, pvar=pvar)


    if (nrow(main_loco_filt) > 0) {

        ## load dosages from PGEN
        var_ids <- main_loco_filt$ID
        var_num <- rep(NA_real_, length(var_ids))
        for (i in 1:length(var_ids)) {
            var_num[i] <- pgenlibr::GetVariantsById(pvar = pvar, id = var_ids[i])
        }

        geno_loco <- as.matrix(ReadList(pgen, var_num))


        if (ncol(geno_loco) == 1) {
            ## even with only on SNP, must check if it remains significant once we add the dominance term

            ## subset to sample IIDs with non-missing phenotype (residual)
            geno_loco <- geno_loco[sample_f$IID %in% phen_df$IID, , drop = FALSE]
            geno_loco_samples <- geno_loco

            geno_loco <- geno_loco[!is.na(phen_df[, 3]), , drop = FALSE]
            
            ## flip allele where needed to match A1
            if (main_loco_filt$A1 == main_loco_filt$REF) {
                geno_loco[, 1] <- 2 - geno_loco[, 1]
                geno_loco_samples[, 1] <- 2 - geno_loco_samples[, 1]
            }


            ## phenotype vector
            phen_vec <- phen_df[!is.na(phen_df[, 3]), 3]

            ## PGS vector
            pgs_vec_loco <- pgs_df_loco[[chr]][!is.na(phen_df[, 3]), 3]

            ## initial covars matrix
            X <- as.matrix(cbind(rep(1, nrow(geno_loco)),  # intercept
                                 pgs_vec_loco,
                                 geno_loco[, 1],  # most significant SNP
                                 geno_loco[, 1]^2,  # dominance deviation
                                 pgs_vec_loco * geno_loco[, 1]))  # most significant SNP * PGS

            regression <- fastLm(X,
                                 phen_vec)

            if (summary(regression)$coeff[nrow(summary(regression)$coeff), 4] <= p1) {
                ## add to list of independent hits
                indep_hits_loco_id[[chr]]  <- var_ids[1]
                indep_hits_loco_ind[[chr]] <- 1
            }

        } else {  # run stepwise regression
            
            ## subset to sample IIDs with non-missing phenotype (residual)
            geno_loco <- geno_loco[sample_f$IID %in% phen_df$IID,]
            geno_loco_samples <- geno_loco  ## save version with all samples belonging to sample code, useful to get samples with 0/1/2 mutations below

            geno_loco <- geno_loco[!is.na(phen_df[, 3]), , drop = FALSE]

            ## flip allele where needed to match A1
            ref_snps <- (main_loco_filt$A1 == main_loco_filt$REF)
            if (sum(ref_snps) > 0) {
                ref_snps_cols <- seq(1, ncol(geno_loco))[ref_snps]
                geno_loco[, ref_snps_cols] <- 2 - geno_loco[, ref_snps_cols]
                geno_loco_samples[, ref_snps_cols] <- 2 - geno_loco_samples[, ref_snps_cols]
            }


            ## phenotype vector
            phen_vec <- phen_df[!is.na(phen_df[, 3]), 3]

            ## PGS vector
            pgs_vec_loco <- pgs_df_loco[[chr]][!is.na(phen_df[, 3]), 3]

            ## initial covars matrix
            X <- as.matrix(cbind(rep(1, nrow(geno_loco)),  # intercept
                                 pgs_vec_loco))

            ## iterate
            for (i in 1:ncol(geno_loco)) {
                if (length(indep_hits_loco_id[[chr]]) > 0) {
                    ## check squared correlation with existing variables
                    cor_prev <- cor(geno_loco[, i], geno_loco[, indep_hits_loco_ind[[chr]]])
                    if (max(cor_prev^2) > 0.9) {
                        next
                    }
                }
                
                X_add <- cbind(X,
                               geno_loco[, i],  # most significant SNP
                               geno_loco[, i]^2,  # dominance deviation
                               pgs_vec_loco * geno_loco[, i])  # most significant SNP * PGS
                
                regression <- fastLm(X_add,
                                     phen_vec)

                if (summary(regression)$coeff[nrow(summary(regression)$coeff), 4] <= p1) {
                    X <- X_add
                    indep_hits_loco_id[[chr]]  <- c(indep_hits_loco_id[[chr]], var_ids[i])
                    indep_hits_loco_ind[[chr]] <- c(indep_hits_loco_ind[[chr]], i)
                }
            }

            if (length(indep_hits_loco_id[[chr]]) > 0) {

                ## restore positional order
                indep_hits_loco_df <- data.frame(ID  = indep_hits_loco_id[[chr]],
                                                 ind = indep_hits_loco_ind[[chr]],
                                                 stringsAsFactors = FALSE)
                indep_hits_loco_df <- indep_hits_loco_df[order(match(indep_hits_loco_df$ID, main_loco[[chr]]$ID)),]

                ## add to list of independent hits
                indep_hits_loco_id[[chr]]  <- indep_hits_loco_df$ID
                indep_hits_loco_ind[[chr]] <- indep_hits_loco_df$ind

            }
        }

        ## compute signed log p-values for additive effect + dominance deviation + SNP*PGS all together
        if (length(indep_hits_loco_id[[chr]]) > 0) {
            for (i in 1:length(indep_hits_loco_id[[chr]])) {
                snp_id <- indep_hits_loco_id[[chr]][i]
                snp_ind <- indep_hits_loco_ind[[chr]][i]

                X <- as.matrix(cbind(rep(1, nrow(geno_loco)),  # intercept
                                     pgs_vec_loco,
                                     geno_loco[, snp_ind],  # independent hit
                                     geno_loco[, snp_ind]^2,  # dominance deviation
                                     pgs_vec_loco * geno_loco[, snp_ind]))  # SNP * PGS
                
                regression <- fastLm(X,
                                     phen_vec)

                indep_hits_loco_add_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[3, 4]) * sign(summary(regression)$coeff[3, 1])
                indep_hits_loco_dom_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[4, 4]) * sign(summary(regression)$coeff[4, 1])
                indep_hits_loco_snp_pgs_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[5, 4]) * sign(summary(regression)$coeff[5, 1])
            }
        }
    }





 


    ## full ---------------------------------------------------------------------

    ## load PLINK2 SNP*PGS interaction GWAS output
    main_full[[chr]] <- fread(paste0("../results/03-interaction-gwas/plink-output/", phen, "/snp-pgs/",
                                     phen, ".", samples, ".full.chr", chr, ".", phen, ".res.cov.pgs_full_mc.ac.ctr.glm.linear.gz"),
                              data.table = FALSE)
    main_full[[chr]] <- main_full[[chr]][main_full[[chr]]$TEST == "ADDxpgs_full_mc",]
    colnames(main_full[[chr]])[1] <- "CHR"

    ## filter by p-value
    main_full_filt <- main_full[[chr]] %>%
        filter(LOG10_P >= -log10(p1)) %>%
        arrange(-LOG10_P)


    if (nrow(main_full_filt) > 0) {

        if (length(indep_hits_loco_id[[chr]]) > 0) {  # if there are LOCO hits

            main_full_filt <- main_full_filt %>%
                filter(!(ID %in% indep_hits_loco_id[[chr]]))
            
            if (nrow(main_full_filt) > 0) {  # if there are full hits which are not independent LOCO hits

                ## load dosages from PGEN
                var_ids <- main_full_filt$ID

                var_num <- rep(NA_real_, length(var_ids))
                for (i in 1:length(var_ids)) {
                    var_num[i] <- pgenlibr::GetVariantsById(pvar = pvar, id = var_ids[i])
                }

                geno_full <- as.matrix(ReadList(pgen, var_num))


                if (ncol(geno_full) == 1) {  # if there is only one full hit

                    ## subset to sample IIDs with non-missing phenotype (residual)
                    geno_full <- geno_full[sample_f$IID %in% phen_df$IID, , drop = FALSE]
                    geno_full_samples <- geno_full  ## save version with all samples belonging to sample code, useful to get samples with 0/1/2 mutations below

                    geno_full <- geno_full[!is.na(phen_df[, 3]), , drop = FALSE]

                    ## flip allele where needed to match A1
                    if (main_full_filt$A1 == main_full_filt$REF) {
                        geno_full[, 1] <- 2 - geno_full[, 1]
                        geno_full_samples[, 1] <- 2 - geno_full_samples[, 1]
                    }


                    ## phenotype vector
                    phen_vec <- phen_df[!is.na(phen_df[, 3]), 3]

                    ## PGS vector
                    pgs_vec_loco <- pgs_df_loco[[chr]][!is.na(phen_df[, 3]), 3]
                    pgs_vec_full <- pgs_df_full[!is.na(phen_df[, 3]), 3]

                    ## initial covars matrix
                    X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                         pgs_vec_full,
                                         geno_loco[, indep_hits_loco_ind[[chr]]],  # LOCO independent SNPs
                                         geno_loco[, indep_hits_loco_ind[[chr]]]^2,  # dominance deviation
                                         pgs_vec_full * geno_loco[, indep_hits_loco_ind[[chr]]]))


                    ## check squared correlation with existing variables
                    cor_prev <- cor(geno_full[, 1],
                                    cbind(geno_loco[, indep_hits_loco_ind[[chr]]]))
                    if (max(cor_prev^2) > 0.9) {
                        next
                    }


                    X_add <- cbind(X,
                                   geno_full[, 1],  # most significant SNP
                                   geno_full[, 1]^2,  # dominance deviation
                                   pgs_vec_full * geno_full[, 1])  # most significant SNP * PGS
                    
                    regression <- fastLm(X_add,
                                         phen_vec)

                    if (summary(regression)$coeff[nrow(summary(regression)$coeff), 4] <= p1) {

                        ## add to list of independent hits
                        indep_hits_full_id[[chr]]  <- c(indep_hits_full_id[[chr]], var_ids[1])
                        indep_hits_full_ind[[chr]] <- c(indep_hits_full_ind[[chr]], 1)

                        ## restore positional order
                        indep_hits_full_df <- data.frame(ID  = indep_hits_full_id[[chr]],
                                                         ind = indep_hits_full_ind[[chr]],
                                                         stringsAsFactors = FALSE)
                        indep_hits_full_df <- indep_hits_full_df[order(match(indep_hits_full_df$ID, main_full[[chr]]$ID)),]

                        indep_hits_full_id[[chr]]  <- indep_hits_full_df$ID
                        indep_hits_full_ind[[chr]] <- indep_hits_full_df$ind

                        ## compute signed log p-values for additive effect + dominance deviation + SNP*PGS all together
                        snp_id <- indep_hits_full_id[[chr]]
                        snp_ind <- indep_hits_full_ind[[chr]]

                        X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                             pgs_vec_full,
                                             geno_full[, snp_ind],  # independent hit
                                             geno_full[, snp_ind]^2,  # dominance deviation
                                             pgs_vec_full * geno_full[, snp_ind]))  # SNP * PGS
                        
                        regression <- fastLm(X,
                                             phen_vec)

                        indep_hits_full_add_sign_LOG10_P[[chr]] <- -log10(summary(regression)$coeff[3, 4]) * sign(summary(regression)$coeff[3, 1])
                        indep_hits_full_dom_sign_LOG10_P[[chr]] <- -log10(summary(regression)$coeff[4, 4]) * sign(summary(regression)$coeff[4, 1])
                        indep_hits_full_snp_pgs_sign_LOG10_P[[chr]] <- -log10(summary(regression)$coeff[5, 4]) * sign(summary(regression)$coeff[5, 1])

                    }

                } else {  # run stepwise regression
                    
                    ## subset to sample IIDs with non-missing phenotype (residual)
                    geno_full <- geno_full[sample_f$IID %in% phen_df$IID,]
                    geno_full_samples <- geno_full  ## save version with all samples belonging to sample code, useful to get samples with 0/1/2 mutations below

                    geno_full <- geno_full[!is.na(phen_df[, 3]), , drop = FALSE]

                    ## flip allele where needed to match A1
                    ref_snps <- (main_full_filt$A1 == main_full_filt$REF)
                    if (sum(ref_snps) > 0) {
                        ref_snps_cols <- seq(1, ncol(geno_full))[ref_snps]
                        geno_full[, ref_snps_cols] <- 2 - geno_full[, ref_snps_cols]
                        geno_full_samples[, ref_snps_cols] <- 2 - geno_full_samples[, ref_snps_cols]
                    }


                    ## phenotype vector
                    phen_vec <- phen_df[!is.na(phen_df[, 3]), 3]

                    ## PGS vector
                    pgs_vec_loco <- pgs_df_loco[[chr]][!is.na(phen_df[, 3]), 3]
                    pgs_vec_full <- pgs_df_full[!is.na(phen_df[, 3]), 3]

                    ## initial covars matrix
                    X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                         pgs_vec_full,
                                         geno_loco[, indep_hits_loco_ind[[chr]]],  # LOCO independent SNPs
                                         geno_loco[, indep_hits_loco_ind[[chr]]]^2,  # dominance deviation
                                         pgs_vec_full * geno_loco[, indep_hits_loco_ind[[chr]]]))

                    ## iterate
                    for (i in 1:ncol(geno_full)) {
                        ## (we should check for correlation from the first variable because there are LOCO hits to account for)
                        ## check squared correlation with existing variables
                        cor_prev <- cor(geno_full[, i],
                                        cbind(geno_loco[, indep_hits_loco_ind[[chr]]],
                                              geno_full[, indep_hits_full_ind[[chr]]]))
                        if (max(cor_prev^2) > 0.9) {
                            next
                        }

                        X_add <- cbind(X,
                                       geno_full[, i],  # most significant SNP
                                       geno_full[, i]^2,  # dominance deviation
                                       pgs_vec_full * geno_full[, i])  # most significant SNP * PGS
                        
                        regression <- fastLm(X_add,
                                             phen_vec)

                        if (summary(regression)$coeff[nrow(summary(regression)$coeff), 4] <= p1) {
                            X <- X_add

                            ## add to list of independent hits
                            indep_hits_full_id[[chr]]  <- c(indep_hits_full_id[[chr]], var_ids[i])
                            indep_hits_full_ind[[chr]] <- c(indep_hits_full_ind[[chr]], i)
                        }
                    }

                    if (length(indep_hits_full_id[[chr]]) > 0) {

                        ## restore positional order
                        indep_hits_full_df <- data.frame(ID  = indep_hits_full_id[[chr]],
                                                         ind = indep_hits_full_ind[[chr]],
                                                         stringsAsFactors = FALSE)
                        indep_hits_full_df <- indep_hits_full_df[order(match(indep_hits_full_df$ID, main_full[[chr]]$ID)),]

                        indep_hits_full_id[[chr]]  <- indep_hits_full_df$ID
                        indep_hits_full_ind[[chr]] <- indep_hits_full_df$ind

                        ## compute signed log p-values for additive effect + dominance deviation + SNP*PGS all together
                        for (i in 1:length(indep_hits_full_id[[chr]])) {
                            snp_id <- indep_hits_full_id[[chr]][i]
                            snp_ind <- indep_hits_full_ind[[chr]][i]

                            X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                                 pgs_vec_full,
                                                 geno_full[, snp_ind],  # independent hit
                                                 geno_full[, snp_ind]^2,  # dominance deviation
                                                 pgs_vec_full * geno_full[, snp_ind]))  # SNP * PGS
                            
                            regression <- fastLm(X,
                                                 phen_vec)

                            indep_hits_full_add_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[3, 4]) * sign(summary(regression)$coeff[3, 1])
                            indep_hits_full_dom_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[4, 4]) * sign(summary(regression)$coeff[4, 1])
                            indep_hits_full_snp_pgs_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[5, 4]) * sign(summary(regression)$coeff[5, 1])
                        }
                    }
                }
            }

        } else {  # if there are no LOCO hits

            ## load dosages from PGEN
            var_ids <- main_full_filt$ID
            var_num <- rep(NA_real_, length(var_ids))
            for (i in 1:length(var_ids)) {
                var_num[i] <- pgenlibr::GetVariantsById(pvar = pvar, id = var_ids[i])
            }

            geno_full <- as.matrix(ReadList(pgen, var_num))


            if (ncol(geno_full) == 1) {
                ## even with only on SNP, must check if it remains significant once we add the dominance term

                ## subset to sample IIDs with non-missing phenotype (residual)
                geno_full <- geno_full[sample_f$IID %in% phen_df$IID, , drop = FALSE]
                geno_full_samples <- geno_full

                geno_full <- geno_full[!is.na(phen_df[, 3]), , drop = FALSE]

                ## flip allele where needed to match A1
                if (main_full_filt$A1 == main_full_filt$REF) {
                    geno_full[, 1] <- 2 - geno_full[, 1]
                    geno_full_samples[, 1] <- 2 - geno_full_samples[, 1]
                }


                ## phenotype vector
                phen_vec <- phen_df[!is.na(phen_df[, 3]), 3]

                ## PGS vector
                pgs_vec_full <- pgs_df_full[!is.na(phen_df[, 3]), 3]

                ## initial covars matrix
                X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                     pgs_vec_full,
                                     geno_full[, 1],  # most significant SNP
                                     geno_full[, 1]^2,  # dominance deviation
                                     pgs_vec_full * geno_full[, 1]))  # most significant SNP * PGS

                regression <- fastLm(X,
                                     phen_vec)

                if (summary(regression)$coeff[nrow(summary(regression)$coeff), 4] <= p1) {
                    ## add to list of independent hits
                    indep_hits_full_id[[chr]]  <- var_ids[1]
                    indep_hits_full_ind[[chr]] <- 1
                }
                
            } else {  # run stepwise regression
                
                ## subset to sample IIDs with non-missing phenotype (residual)
                geno_full <- geno_full[sample_f$IID %in% phen_df$IID,]
                geno_full_samples <- geno_full  ## save version with all samples belonging to sample code, useful to get samples with 0/1/2 mutations below

                geno_full <- geno_full[!is.na(phen_df[, 3]), , drop = FALSE]

                ## flip allele where needed to match A1
                ref_snps <- (main_full_filt$A1 == main_full_filt$REF)
                if (sum(ref_snps) > 0) {
                    ref_snps_cols <- seq(1, ncol(geno_full))[ref_snps]
                    geno_full[, ref_snps_cols] <- 2 - geno_full[, ref_snps_cols]
                    geno_full_samples[, ref_snps_cols] <- 2 - geno_full_samples[, ref_snps_cols]
                }


                ## phenotype vector
                phen_vec <- phen_df[!is.na(phen_df[, 3]), 3]

                ## PGS vector
                pgs_vec_full <- pgs_df_full[!is.na(phen_df[, 3]), 3]

                ## initial covars matrix
                X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                     pgs_vec_full))
                
                ## iterate
                for (i in 1:ncol(geno_full)) {
                    if (length(indep_hits_full_id[[chr]]) > 0) {
                        ## check squared correlation with existing variables
                        cor_prev <- cor(geno_full[, i], geno_full[, indep_hits_full_ind[[chr]]])
                        if (max(cor_prev^2) > 0.9) {
                            next
                        }
                    }

                    X_add <- cbind(X,
                                   geno_full[, i],  # most significant SNP
                                   geno_full[, i]^2,  # dominance deviation
                                   pgs_vec_full * geno_full[, i])  # most significant SNP * PGS
                    
                    regression <- fastLm(X_add,
                                         phen_vec)

                    if (summary(regression)$coeff[nrow(summary(regression)$coeff), 4] <= p1) {
                        X <- X_add

                        ## add to list of independent hits
                        indep_hits_full_id[[chr]]  <- c(indep_hits_full_id[[chr]], var_ids[i])
                        indep_hits_full_ind[[chr]] <- c(indep_hits_full_ind[[chr]], i)
                    }
                }

                if (length(indep_hits_full_id[[chr]]) > 0) {

                    ## restore positional order
                    indep_hits_full_df <- data.frame(ID  = indep_hits_full_id[[chr]],
                                                     ind = indep_hits_full_ind[[chr]],
                                                     stringsAsFactors = FALSE)
                    indep_hits_full_df <- indep_hits_full_df[order(match(indep_hits_full_df$ID, main_full[[chr]]$ID)),]

                    indep_hits_full_id[[chr]]  <- indep_hits_full_df$ID
                    indep_hits_full_ind[[chr]] <- indep_hits_full_df$ind

                }
            }

            ## compute signed log p-values for additive effect + dominance deviation + SNP*PGS all together
            if (length(indep_hits_full_id[[chr]]) > 0) {
                for (i in 1:length(indep_hits_full_id[[chr]])) {
                    snp_id <- indep_hits_full_id[[chr]][i]
                    snp_ind <- indep_hits_full_ind[[chr]][i]

                    X <- as.matrix(cbind(rep(1, nrow(geno_full)),  # intercept
                                         pgs_vec_full,
                                         geno_full[, snp_ind],  # independent hit
                                         geno_full[, snp_ind]^2,  # dominance deviation
                                         pgs_vec_full * geno_full[, snp_ind]))  # SNP * PGS
                    
                    regression <- fastLm(X,
                                         phen_vec)

                    indep_hits_full_add_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[3, 4]) * sign(summary(regression)$coeff[3, 1])
                    indep_hits_full_dom_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[4, 4]) * sign(summary(regression)$coeff[4, 1])
                    indep_hits_full_snp_pgs_sign_LOG10_P[[chr]][i] <- -log10(summary(regression)$coeff[5, 4]) * sign(summary(regression)$coeff[5, 1])
                }
            }
        }
    }



    if (length(indep_hits_full_id[[chr]]) == 0 & length(indep_hits_loco_id[[chr]]) == 0) {
        next
    }



    ## export list of samples with 0/1/2 genotypes for each independent hit
    if (length(indep_hits_loco_id[[chr]]) > 0) {
        for (i in 1:length(indep_hits_loco_id[[chr]])) {
            sp_0 <- phen_df$IID[round(geno_loco_samples[, indep_hits_loco_ind[[chr]][i]]) == 0]
            sp_1 <- phen_df$IID[round(geno_loco_samples[, indep_hits_loco_ind[[chr]][i]]) == 1]
            sp_2 <- phen_df$IID[round(geno_loco_samples[, indep_hits_loco_ind[[chr]][i]]) == 2]

            fwrite(data.frame(`#FID` = sp_0, IID = sp_0, check.names = FALSE),
                   file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/loco/", gsub(":", "_", indep_hits_loco_id[[chr]][i]), "-gen0.tab"),
                   sep = "\t", na = "NA", quote = FALSE)
            fwrite(data.frame(`#FID` = sp_1, IID = sp_1, check.names = FALSE),
                   file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/loco/", gsub(":", "_", indep_hits_loco_id[[chr]][i]), "-gen1.tab"),
                   sep = "\t", na = "NA", quote = FALSE)
            fwrite(data.frame(`#FID` = sp_2, IID = sp_2, check.names = FALSE),
                   file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/loco/", gsub(":", "_", indep_hits_loco_id[[chr]][i]), "-gen2.tab"),
                   sep = "\t", na = "NA", quote = FALSE)
        }
    }

    if (length(indep_hits_full_id[[chr]]) > 0){
        for (i in 1:length(indep_hits_full_id[[chr]])) {
            sp_0 <- phen_df$IID[round(geno_full_samples[, indep_hits_full_ind[[chr]][i]]) == 0]
            sp_1 <- phen_df$IID[round(geno_full_samples[, indep_hits_full_ind[[chr]][i]]) == 1]
            sp_2 <- phen_df$IID[round(geno_full_samples[, indep_hits_full_ind[[chr]][i]]) == 2]

            fwrite(data.frame(`#FID` = sp_0, IID = sp_0, check.names = FALSE),
                   file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/full/", gsub(":", "_", indep_hits_full_id[[chr]][i]), "-gen0.tab"),
                   sep = "\t", na = "NA", quote = FALSE)
            fwrite(data.frame(`#FID` = sp_1, IID = sp_1, check.names = FALSE),
                   file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/full/", gsub(":", "_", indep_hits_full_id[[chr]][i]), "-gen1.tab"),
                   sep = "\t", na = "NA", quote = FALSE)
            fwrite(data.frame(`#FID` = sp_2, IID = sp_2, check.names = FALSE),
                   file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/full/", gsub(":", "_", indep_hits_full_id[[chr]][i]), "-gen2.tab"),
                   sep = "\t", na = "NA", quote = FALSE)
        }
    }
}



# load minor alleles
minal_f <- paste0("../data/imputed-wb_all-stats/minor-alleles/chr", seq(1, 22), ".tab")
minal <- lapply(minal_f, function(x) fread(x, data.table = FALSE))
minal <- do.call("rbind", minal)

# load rsIDs
map_rsid_f <- paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", seq(1, 22), ".tab")
map_rsid <- lapply(map_rsid_f, function(x) fread(x, header = TRUE, data.table = FALSE))
map_rsid <- do.call("rbind", map_rsid)

## original GWAS
orig_f <- paste0('../results/01-gwas/plink-output/', phen, '/', samples, '.chr', seq(1, 22), '.', phen, '.res.cov.glm.linear.gz')
orig <- lapply(orig_f, function(x) fread(x, data.table = FALSE))
orig <- do.call("rbind", orig)
colnames(orig)[1] <- 'CHR'
## compute signed p-values
orig$orig_sign_LOG10_P <- orig$LOG10_P * sign(orig$BETA)

## full PGS
main_full <- do.call("rbind", main_full)
## compute signed p-values
main_full$snp_pgs_full_sign_LOG10_P <- main_full$LOG10_P * sign(main_full$BETA)

## LOCO PGS
main_loco <- do.call("rbind", main_loco)
## compute signed p-values
main_loco$snp_pgs_loco_sign_LOG10_P <- main_loco$LOG10_P * sign(main_loco$BETA)

## abs. res. GWAS
res_abs_f <- paste0("../results/03-interaction-gwas/plink-output/", phen, "/res_abs/",
                    phen, ".", samples, ".full.chr", seq(1, 22), ".", phen, ".res_abs.cov.pgs_full_mc.ac.ctr.glm.linear.gz")
res_abs <- lapply(res_abs_f, function(x) fread(x, data.table = FALSE))
res_abs <- do.call("rbind", res_abs)
## compute signed p-values
res_abs$res_abs_full_sign_LOG10_P <- res_abs$LOG10_P * sign(res_abs$BETA)



## export list of independent hits for all chromosomes
if (length(unlist(indep_hits_loco_id)) > 0) {
    indep_hits_loco_id_df <- data.frame(ID = unlist(indep_hits_loco_id),
                                        add_loco_sign_LOG10_P = unlist(indep_hits_loco_add_sign_LOG10_P),
                                        dom_loco_sign_LOG10_P = unlist(indep_hits_loco_dom_sign_LOG10_P),
                                        snp_pgs_wdom_loco_sign_LOG10_P = unlist(indep_hits_loco_snp_pgs_sign_LOG10_P),
                                        stringsAsFactors = FALSE)

    ## add position and A1 allele
    indep_hits_loco_id_df$CHR <- as.numeric(word(indep_hits_loco_id_df$ID, 1, sep = ":"))
    indep_hits_loco_id_df <- merge(indep_hits_loco_id_df,
                                   minal[, c("ID", "POS", "REF", "ALT", "A1")],
                                   by = "ID")

    ## add rsID
    indep_hits_loco_id_df <- merge(indep_hits_loco_id_df,
                                   map_rsid,
                                   by.x = "ID", by.y = "POSID")
                                       
    ## reorder and export version without sumstats
    indep_hits_loco_id_df <- indep_hits_loco_id_df[order(match(indep_hits_loco_id_df$ID, minal$ID)),]
    indep_hits_loco_id_df_stats <- indep_hits_loco_id_df
    indep_hits_loco_id_df <- indep_hits_loco_id_df[, c("CHR", "POS", "ID", "rsID", "REF", "ALT", "A1")]
    fwrite(indep_hits_loco_id_df,
           file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-loco-indep-hits.tab"),
           sep = "\t", na = "NA", quote = FALSE)

    ## add signed log p-values
    indep_hits_loco_id_df_stats <- merge(indep_hits_loco_id_df_stats,
                                         main_full[, c("ID", "snp_pgs_full_sign_LOG10_P")],
                                         by = "ID")
    indep_hits_loco_id_df_stats <- merge(indep_hits_loco_id_df_stats,
                                         main_loco[, c("ID", "snp_pgs_loco_sign_LOG10_P")],
                                         by = "ID")
    indep_hits_loco_id_df_stats <- merge(indep_hits_loco_id_df_stats,
                                         res_abs[, c("ID", "res_abs_full_sign_LOG10_P")],
                                         by = "ID")
    indep_hits_loco_id_df_stats <- merge(indep_hits_loco_id_df_stats,
                                         orig[, c("ID", "orig_sign_LOG10_P")],
                                         by = "ID")

    indep_hits_loco_id_df_stats <- indep_hits_loco_id_df_stats[order(match(indep_hits_loco_id_df_stats$ID, minal$ID)),]
    indep_hits_loco_id_df_stats <- indep_hits_loco_id_df_stats[, c("CHR", "POS", "ID", "rsID", "REF", "ALT", "A1",
                                                                   "snp_pgs_full_sign_LOG10_P",
                                                                   "snp_pgs_loco_sign_LOG10_P",
                                                                   "res_abs_full_sign_LOG10_P",
                                                                   "orig_sign_LOG10_P",
                                                                   "add_loco_sign_LOG10_P",
                                                                   "dom_loco_sign_LOG10_P",
                                                                   "snp_pgs_wdom_loco_sign_LOG10_P")]
    fwrite(indep_hits_loco_id_df_stats,
           file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-loco-indep-hits-stats.tab"),
           sep = "\t", na = "NA", quote = FALSE)
} else {
    file.create(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-loco-indep-hits.tab"))
}

if (length(unlist(indep_hits_full_id)) > 0) {
    indep_hits_full_id_df <- data.frame(ID = unlist(indep_hits_full_id),
                                        add_full_sign_LOG10_P = unlist(indep_hits_full_add_sign_LOG10_P),
                                        dom_full_sign_LOG10_P = unlist(indep_hits_full_dom_sign_LOG10_P),
                                        snp_pgs_wdom_full_sign_LOG10_P = unlist(indep_hits_full_snp_pgs_sign_LOG10_P),
                                        stringsAsFactors = FALSE)

    ## add position and A1 allele
    indep_hits_full_id_df$CHR <- as.numeric(word(indep_hits_full_id_df$ID, 1, sep = ":"))
    indep_hits_full_id_df <- merge(indep_hits_full_id_df,
                                   minal[, c("ID", "POS", "REF", "ALT", "A1")],
                                   by = "ID")

    ## add rsID
    indep_hits_full_id_df <- merge(indep_hits_full_id_df,
                                   map_rsid,
                                   by.x = "ID", by.y = "POSID")
                                       
    ## reorder and export version without sumstats
    indep_hits_full_id_df <- indep_hits_full_id_df[order(match(indep_hits_full_id_df$ID, minal$ID)),]
    indep_hits_full_id_df_stats <- indep_hits_full_id_df
    indep_hits_full_id_df <- indep_hits_full_id_df[, c("CHR", "POS", "ID", "rsID", "REF", "ALT", "A1")]
    fwrite(indep_hits_full_id_df,
           file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-full-indep-hits.tab"),
           sep = "\t", na = "NA", quote = FALSE)

    ## add signed log p-values
    indep_hits_full_id_df_stats <- merge(indep_hits_full_id_df_stats,
                                         main_full[, c("ID", "snp_pgs_full_sign_LOG10_P")],
                                         by = "ID")
    indep_hits_full_id_df_stats <- merge(indep_hits_full_id_df_stats,
                                         main_loco[, c("ID", "snp_pgs_loco_sign_LOG10_P")],
                                         by = "ID")
    indep_hits_full_id_df_stats <- merge(indep_hits_full_id_df_stats,
                                         res_abs[, c("ID", "res_abs_full_sign_LOG10_P")],
                                         by = "ID")
    indep_hits_full_id_df_stats <- merge(indep_hits_full_id_df_stats,
                                         orig[, c("ID", "orig_sign_LOG10_P")],
                                         by = "ID")

    indep_hits_full_id_df_stats <- indep_hits_full_id_df_stats[order(match(indep_hits_full_id_df_stats$ID, minal$ID)),]
    indep_hits_full_id_df_stats <- indep_hits_full_id_df_stats[, c("CHR", "POS", "ID", "rsID", "REF", "ALT", "A1",
                                                                   "snp_pgs_full_sign_LOG10_P",
                                                                   "snp_pgs_loco_sign_LOG10_P",
                                                                   "res_abs_full_sign_LOG10_P",
                                                                   "orig_sign_LOG10_P",
                                                                   "add_full_sign_LOG10_P",
                                                                   "dom_full_sign_LOG10_P",
                                                                   "snp_pgs_wdom_full_sign_LOG10_P")]

    fwrite(indep_hits_full_id_df_stats,
           file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-full-indep-hits-stats.tab"),
           sep = "\t", na = "NA", quote = FALSE)
} else {
    file.create(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-full-indep-hits.tab"))
}
