## Find independent hits from GWAS for each independent hit (LOCO)


library(argparser)
library(data.table)
library(dplyr)
library(pgenlibr)

options(warn=2)  # warnings as errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",     help = "Phenotype code",    nargs = 1)
p <- add_argument(p, "--train_sp", help = "Training samples",         nargs = 1)
p <- add_argument(p, "--test_sp",  help = "Test samples",         nargs = 1)
p <- add_argument(p, "--snp",      help = "Interacting SNP", nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen
train_sp <- argv$train_sp
test_sp  <- argv$test_sp
snp_     <- argv$snp
snp      <- gsub("_", ":", snp_)




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




## function for running stepwise regression
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

stepwise <- function(main_train,
                     main_test,
                     snps_df,
                     snp,
                     snps_prev = NULL,
                     p1) {

    ## check if anchor SNP and others are in same chromosome
    snp_chr  <- as.numeric(strsplit(snp, split = ":")[[1]][1])
    same_chr <- NULL
    ifelse(unique(snps_df$CHR) == snp_chr, same_chr <- TRUE, same_chr <- FALSE)

    ## data.frame to store independent pairwise interactions
    indep_hits <- indep_hits_empty

    ## initial matrix of regressors
    ## check if (absolute) correlation between anchor SNP and its square is <0.999
    incl_sq_1_train <- NULL
    if (abs(cor(main_train[[snp]], main_train[[snp]]^2)) < 0.999) {
        incl_sq_1_train <- TRUE
        X_train <- as.matrix(cbind(main_train[[snp]],
                                   main_train[[snp]]^2))  # anchor SNP + its square
    } else {
        incl_sq_1_train <- FALSE
        X_train <- as.matrix(main_train[[snp]])
    }
    ## (it's ok for training and test to be different in whether or not
    ## squared term is included as test set is only used to compute heritability)
    if (abs(cor(main_test[[snp]], main_test[[snp]]^2)) < 0.999) {
        X_test <- as.matrix(cbind(main_test[[snp]],
                                  main_test[[snp]]^2))
    } else {
        X_test <- as.matrix(main_test[[snp]])
    }

    ## if there are previously identified independent pairwise hits to take into account
    if (length(snps_prev) > 0) {
        for (i in 1:length(snps_prev)) {
            if (abs(cor(main_train[[snps_prev[i]]], main_train[[snps_prev[i]]]^2)) < 0.999) {
                X_train <- cbind(X_train,
                                 main_train[[snps_prev[i]]],
                                 main_train[[snps_prev[i]]]^2,
                                 main_train[[snp]] * main_train[[snps_prev[i]]])
            } else {
                X_train <- cbind(X_train,
                                 main_train[[snps_prev[i]]],
                                 main_train[[snp]] * main_train[[snps_prev[i]]])
            }
            if (abs(cor(main_test[[snps_prev[i]]], main_test[[snps_prev[i]]]^2)) < 0.999) {
                X_test  <- cbind(X_test,
                                 main_test[[snps_prev[i]]],
                                 main_test[[snps_prev[i]]]^2,
                                 main_test[[snp]] * main_test[[snps_prev[i]]])
            } else {
                X_test  <- cbind(X_test,
                                 main_test[[snps_prev[i]]],
                                 main_test[[snp]] * main_test[[snps_prev[i]]])
            }
        }
    }

    snps_main_only <- NULL
    for (i in 1:nrow(snps_df)) {

        ## check that current SNP is not in high LD (r^2 > 0.9) with any of the SNPs included so far
        snps_so_far <- c(snps_prev, indep_hits$ID_2, snps_main_only)
        if (length(snps_so_far) > 0) {
            cor_so_far <- cor(main_train[[snps_df$ID[i]]], main_train[, snps_so_far])
            if (max(cor_so_far^2) > 0.9) { next } 
        }
        
        ## check if correlation between current SNP and its square is <0.999 to see if should include its square
        ### training set
        incl_sq_2_train <- NULL
        if (abs(cor(main_train[[snps_df$ID[i]]], main_train[[snps_df$ID[i]]]^2)) < 0.999) {
            incl_sq_2_train <- TRUE
            X_add_train <- as.matrix(cbind(X_train,
                                           main_train[[snps_df$ID[i]]],
                                           main_train[[snps_df$ID[i]]]^2,  # add next most significant SNP + its square
                                           main_train[[snp]] * main_train[[snps_df$ID[i]]]))
        } else {
            incl_sq_2_train <- FALSE
            X_add_train <- as.matrix(cbind(X_train,
                                           main_train[[snps_df$ID[i]]],
                                           main_train[[snp]] * main_train[[snps_df$ID[i]]]))
        }
        ### test set
        incl_sq_2_test <- NULL
        if (abs(cor(main_test[[snps_df$ID[i]]], main_test[[snps_df$ID[i]]]^2)) < 0.999) {
            incl_sq_2_test <- TRUE
            X_add_test <- as.matrix(cbind(X_test,
                                          main_test[[snps_df$ID[i]]],
                                          main_test[[snps_df$ID[i]]]^2,
                                          main_test[[snp]] * main_test[[snps_df$ID[i]]]))
        } else {
            incl_sq_2_test <- FALSE
            X_add_test <- as.matrix(cbind(X_test,
                                          main_test[[snps_df$ID[i]]],
                                          main_test[[snp]] * main_test[[snps_df$ID[i]]]))
        }


        ## run regression
        reg_train <- lm(main_train[[paste0(phen, ".res.cov")]] ~ X_add_train, singular.ok = FALSE)
        sum_coeff  <- summary(reg_train)$coeff  # just for readability


        ## if interaction is significant
        if (sum_coeff[nrow(sum_coeff), 4] <= p1) {

            remove_snp_test <- TRUE
            remove_sq_test  <- NULL

            ## compute LD with anchor SNP if in same chromosome
            if (same_chr) {
                ld <- cor(main_train[[snp]], main_train[[snps_df$ID[i]]])^2
            } else {
                ld <- NA
            }

            ## check that at least 1 person in test set has both mutations and if so compute h2
            if (sum(X_add_test[, ncol(X_add_test)]) > 0) {

                remove_snp_test <- FALSE
                remove_sq_test  <- FALSE

                ## catch errors
                reg_test <- try(lm(main_test[[paste0(phen, ".res.cov")]] ~ X_add_test, singular.ok = FALSE), silent = TRUE)
                reg_no_int_test <- try(lm(main_test[[paste0(phen, ".res.cov")]] ~ X_add_test[, -ncol(X_add_test)], singular.ok = FALSE), silent = TRUE)

                if (class(reg_test) == "try-error" | class(reg_no_int_test) == "try-error") {

                    ## if there's an error, check if it's due to square of SNP
                    if (incl_sq_2_test) {
                        
                        reg_test_no_dom <- try(lm(main_test[[paste0(phen, ".res.cov")]] ~ X_add_test[, -(ncol(X_add_test) - 1)], singular.ok = FALSE), silent = TRUE)
                        reg_no_int_test_no_dom <- try(lm(main_test[[paste0(phen, ".res.cov")]] ~ X_add_test[, -seq(ncol(X_add_test) - 1, ncol(X_add_test))], singular.ok = FALSE), silent = TRUE)

                        if (class(reg_test_no_dom) == "try-error" | class(reg_no_int_test_no_dom) == "try-error") {
                            print(paste0("Remove ", snps_df$ID[i],
                                         " entirely from test set matrix as it causes a singularity error."))
                            remove_snp_test <- TRUE
                            h2 <- NA
                        } else {
                            print(paste0("Remove the square of ", snps_df$ID[i],
                                         " from test set matrix as it causes a singularity error."))
                            remove_sq_test <- TRUE
                            h2 <- summary(reg_test_no_dom)$r.squared - summary(reg_no_int_test_no_dom)$r.squared
                        }
                    } else {
                        print(paste0("Remove ", snps_df$ID[i],
                                     " entirely from test set matrix as it causes a singularity error."))
                        remove_snp_test <- TRUE  # remove SNP from test set matrix
                        h2 <- NA
                    }
                } else {
                    h2 <- summary(reg_test)$r.squared - summary(reg_no_int_test)$r.squared
                }
            } else {
                h2 <- NA
            }
            
            ## add to list of independent pairwise interactions
            indep_hits_add <- data.frame(
                ## first SNP in pair
                ID_1         = snp,
                A1_1         = snp_a1,
                MAF_1        = snp_maf,
                BETA_1       =                 sum_coeff[2, 1],
                LOG10_P_1    = - ( pnorm( abs( sum_coeff[2, 3] ), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2) ),  # to avoid numerical underflow
                BETA_sq_1    = ifelse(incl_sq_1_train,                 sum_coeff[3, 1],                                                             NA),
                LOG10_P_sq_1 = ifelse(incl_sq_1_train, - ( pnorm( abs( sum_coeff[3, 3] ), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2) ), NA),
                ## 2nd SNP in pair
                ID_2         = snps_df$ID[i],
                A1_2         = snps_df$A1[i],
                MAF_2        = maf$MAF[maf$ID == snps_df$ID[i]],
                BETA_2       = ifelse(incl_sq_2_train,                 sum_coeff[nrow(sum_coeff) - 2, 1],                                                             sum_coeff[nrow(sum_coeff) - 1, 1]),
                LOG10_P_2    = ifelse(incl_sq_2_train, - ( pnorm( abs( sum_coeff[nrow(sum_coeff) - 2, 3] ), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2) ), sum_coeff[nrow(sum_coeff) - 1, 3]),
                BETA_sq_2    = ifelse(incl_sq_2_train,                 sum_coeff[nrow(sum_coeff) - 1, 1],                                                             NA),
                LOG10_P_sq_2 = ifelse(incl_sq_2_train, - ( pnorm( abs( sum_coeff[nrow(sum_coeff) - 1, 3] ), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2) ), NA),
                ## interaction
                BETA_int     =                 sum_coeff[nrow(sum_coeff), 1],
                LOG10_P_int  = - ( pnorm( abs( sum_coeff[nrow(sum_coeff), 3] ), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2) ),
                h2           = h2,
                LD           = ld)
            indep_hits <- rbind(indep_hits, indep_hits_add)

            ## update base matrix
            X_train <- X_add_train

            if (remove_snp_test) {
                next
            } else if (remove_sq_test) {
                X_test <- X_add_test[, -(ncol(X_add_test) - 1)]  # remove square term
            } else {
                X_test <- X_add_test
            }
            
        } else if (sum_coeff[nrow(sum_coeff) - ifelse(incl_sq_2_train, 2, 1), 4] <= p1) {

            ## if the interaction isn't significant but the main effect is, we keep it in the model
            ## (as well as its square if not in too high LD with main variable)
            X_train <- X_add_train[, -ncol(X_add_train)]
            X_test  <- X_add_test[, -ncol(X_add_test)]

            snps_main_only <- c(snps_main_only, snps_df$ID[i])

        }
    }

    return(indep_hits)
}




## load MAF
maf_f <- paste0("../data/imputed-wb_all-stats/minor-alleles/chr", seq(1, 22), ".tab")
maf <- lapply(maf_f, function(x) fread(x, data.table = FALSE))
maf <- do.call(rbind, maf)


## get SNP info
snp_chr <- as.numeric(strsplit(snp, split = ":")[[1]][1])
snp_pos <- as.numeric(strsplit(snp, split = ":")[[1]][2])
## A1
gwas <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/gwas/", phen, ".", train_sp, ".",
                     snp_, ".int.loco.chr", snp_chr, ".", phen, ".res.cov.glm.linear.gz"), data.table = FALSE)
gwas <- gwas[gwas$TEST == "ADD",]
snp_a1  <- gwas$A1[gwas$ID == snp]
## MAF
snp_maf <- maf$MAF[maf$ID == snp]


## load list of independent SNPs in the PGS
clump_f <- paste0("../results/02-pgs/ld-clump/", phen, "/p1_0.05-r2_0.1-kb_500/",
                  train_sp, "-chr", seq(1, 22), "-clumped.tab")
clump <- lapply(clump_f, function(x) fread(x, data.table = FALSE))
clump <- do.call(rbind, clump)
clump <- clump[clump$LOG10_P >= -log10(5e-8),]
p1_pgs <- 0.05 / (nrow(clump) * 2)  # p-value threshold for testing independent PGS loci


## load phenotypes (residuals from basic GWAS)
resid_train <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/", phen, "-", train_sp, "-resid-covars.tab"),
                     data.table = FALSE)
resid_train <- resid_train[!is.na(resid_train[[paste0(phen, ".res.cov")]]),
                           c("IID", paste0(phen, ".res.cov"))]
resid_test  <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/", phen, "-", test_sp, "-resid-covars.tab"),
                     data.table = FALSE)
resid_test <- resid_test[!is.na(resid_test[[paste0(phen, ".res.cov")]]),
                         c("IID", paste0(phen, ".res.cov"))]


## load sample file
sample_f <- fread("../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam", data.table = FALSE)
## load genotype for independent interaction SNP
geno_snp_train <- load_geno(snp, snp_a1, resid_train$IID)
geno_snp_test  <- load_geno(snp, snp_a1, resid_test$IID)




## find independent interaction pairs
indep_hits <- replicate(22, indep_hits_empty, simplify = FALSE)

for (chr in 1:22) {

    ## load interaction GWAS summary stats
    gwas <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/gwas/", phen, ".", train_sp, ".",
                         snp_, ".int.loco.chr", chr, ".", phen, ".res.cov.glm.linear.gz"), data.table = FALSE)
    colnames(gwas)[1] <- "CHR"
    ## keep only interaction summary stats
    gwas <- gwas[gwas$TEST == paste0("ADDx", snp, "_", snp_a1),]

    ## if in same chromosome
    if (chr == snp_chr) {

        ## remove anchor SNP as has NA GWAS results
        gwas <- gwas[gwas$ID != snp,]

        ## remove any SNP with LD (R^2) > 0.1 within 3Mb
        gwas_ld <- gwas[abs(gwas$POS - snp_pos) <= 1.5e6,]
        geno_ld_train <- load_geno(gwas_ld$ID,
                                   gwas_ld$A1,
                                   resid_train$IID)

        cor_df <- data.frame(ID     = gwas_ld$ID,
                             COR_SQ = as.numeric(cor(geno_snp_train[, 2],
                                                     geno_ld_train[, -1])^2))
        
        cor_df_high <- cor_df[cor_df$COR_SQ > 0.1,]
        gwas <- gwas[!(gwas$ID %in% cor_df_high$ID),]

        rm(geno_ld_train)  # memory management
    }

    ## remove SNPs with ERRCODE == "VIF_INFINITE"
    ## (pairs of SNPs where no one has both mutations have NA results and this error code)
    if (sum(gwas$ERRCODE == "VIF_INFINITE") > 0) {
        gwas <- gwas[!(gwas$ERRCODE == "VIF_INFINITE"),]
    }

    ## check that there are no SNPs with NA GWAS results left
    if (sum(is.na(gwas$LOG10_P)) > 0) {
        stop("There are SNPs for which GWAS failed to run.")
    }

    
    ## check if there are GWS hits
    if (sum(gwas$LOG10_P >= -log10(5e-8)) > 0) {  # run stepwise regression
        
        ## filter and sort hits
        gwas_gws <- gwas[gwas$LOG10_P >= -log10(5e-8),]
        gwas_gws <- gwas_gws %>%
            arrange(-LOG10_P)

        ## load genotypes for these SNPs
        geno_gws_train <- load_geno(gwas_gws$ID, gwas_gws$A1, resid_train$IID)
        geno_gws_test  <- load_geno(gwas_gws$ID, gwas_gws$A1, resid_test$IID)

        ## run stepwise regression
        main_train <- resid_train
        main_train <- left_join(main_train, geno_snp_train, by = "IID")  # add anchor SNP
        main_train <- left_join(main_train, geno_gws_train, by = "IID")  # add GWS hits

        main_test <- resid_test
        main_test <- left_join(main_test, geno_snp_test, by = "IID")  # add anchor SNP
        main_test <- left_join(main_test, geno_gws_test, by = "IID")  # add GWS hits

        indep_hits[[chr]] <- stepwise(main_train = main_train,
                                      main_test  = main_test,
                                      snps_df    = gwas_gws,
                                      snp        = snp,
                                      p1         = 5e-8)
    } 


    ## check if there is any additional signal in the independent PGS loci (w/ less stringent p-value threshold)
    gwas_clump <- gwas[gwas$ID %in% clump$ID,]

    if (nrow(indep_hits[[chr]]) > 0) {

        gwas_clump <- gwas_clump[!(gwas_clump$ID %in% indep_hits[[chr]]$ID),]

        if (sum(gwas_clump$LOG10_P >= -log10(p1_pgs)) > 0) {  # run stepwise regression including previous hits

            ## filter and sort hits
            gwas_pgsws <- gwas_clump[gwas_clump$LOG10_P >= -log10(p1_pgs),]
            gwas_pgsws <- gwas_pgsws %>%
                arrange(-LOG10_P)

            ## load genotypes for these + previously identified SNPs
            geno_pgsws_train <- load_geno(c(gwas_pgsws$ID, indep_hits[[chr]]$ID_2),
                                          c(gwas_pgsws$A1, indep_hits[[chr]]$A1_2),
                                          resid_train$IID)
            geno_pgsws_test  <- load_geno(c(gwas_pgsws$ID, indep_hits[[chr]]$ID_2),
                                          c(gwas_pgsws$A1, indep_hits[[chr]]$A1_2),
                                          resid_test$IID)

            ## run stepwise regression
            main_train <- resid_train
            main_train <- left_join(main_train, geno_snp_train, by = "IID")  # add anchor SNP
            main_train <- left_join(main_train, geno_pgsws_train, by = "IID")  # add PGSWS hits

            main_test <- resid_test
            main_test <- left_join(main_test, geno_snp_test, by = "IID")  # add anchor SNP
            main_test <- left_join(main_test, geno_pgsws_test, by = "IID")  # add PGSWS hits

            indep_hits_add <- stepwise(main_train = main_train,
                                       main_test  = main_test,
                                       snps_df    = gwas_pgsws,
                                       snp        = snp,
                                       snps_prev  = indep_hits[[chr]]$ID_2,
                                       p1         = p1_pgs)
            indep_hits[[chr]] <- rbind(indep_hits[[chr]], indep_hits_add)
        } 

    } else if (sum(gwas_clump$LOG10_P >= -log10(p1_pgs)) > 0) {  # run stepwise regression

        ## filter and sort hits
        gwas_pgsws <- gwas_clump[gwas_clump$LOG10_P >= -log10(p1_pgs),]
        gwas_pgsws <- gwas_pgsws %>%
            arrange(-LOG10_P)

        ## load genotypes for these SNPs
        geno_pgsws_train <- load_geno(gwas_pgsws$ID, gwas_pgsws$A1, resid_train$IID)
        geno_pgsws_test  <- load_geno(gwas_pgsws$ID, gwas_pgsws$A1, resid_test$IID)

        ## run stepwise regression
        main_train <- resid_train
        main_train <- left_join(main_train, geno_snp_train, by = "IID")    # add anchor SNP
        main_train <- left_join(main_train, geno_pgsws_train, by = "IID")  # add PGSWS hits

        main_test <- resid_test
        main_test <- left_join(main_test, geno_snp_test, by = "IID")    # add anchor SNP
        main_test <- left_join(main_test, geno_pgsws_test, by = "IID")  # add PGSWS hits

        indep_hits[[chr]] <- stepwise(main_train = main_train,
                                      main_test  = main_test,
                                      snps_df    = gwas_pgsws,
                                      snp        = snp,
                                      p1         = p1_pgs)
    }
}




## bind data.frames for different chromosomes
indep_hits_all <- do.call(rbind, indep_hits)

## export table with independent interacting pairs (or empty file if there are none)
if (nrow(indep_hits_all) > 0) {

    ## remove any interacting pairs with LD > 0.1
    indep_hits_all <- indep_hits_all[is.na(indep_hits_all$LD) |
                                     (!is.na(indep_hits_all$LD) & indep_hits_all$LD <= 0.1),]

    ## round numbers to 8 decimal places
    cols <- c("MAF_1", "BETA_1", "LOG10_P_1", "BETA_sq_1", "LOG10_P_sq_1",
              "MAF_2", "BETA_2", "LOG10_P_2", "BETA_sq_2", "LOG10_P_sq_2",
              "BETA_int", "LOG10_P_int", "h2", "LD")
    indep_hits_all <- as.data.table(indep_hits_all)
    indep_hits_all[,(cols) := round(.SD, 8), .SDcols = cols]

    fwrite(indep_hits_all,
           file = paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", train_sp, "-loco-", snp_, "-pairwise-hits.tab"),
           sep = "\t", na = "NA", quote = FALSE)

} else {

    file.create(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", train_sp, "-loco-", snp_, "-pairwise-hits.tab"))

}
