## desc: Compute PGS for optimal global p-value threshold.
##       Also assess performance on training and test data

library(argparser)
library(biglm)
library(data.table)


p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code",           nargs = 1)
p <- add_argument(p, "--samples", help = "Sample set code", nargs = 1)
p <- add_argument(p, "--r2",      help = "LD clumping r2 parameter", nargs = 1)
p <- add_argument(p, "--kb",      help = "LD clumping kb parameter", nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples
r2      <- argv$r2
kb      <- argv$kb



time <- system.time({
## compute and load score for each chromosome
ft <- TRUE
for (chr in 1:22) {

    ## read coefficient data frame to check if empty
    coeff_chr <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                              "/coeff/p1_opt/chr", chr, "-coeff.tab"),
                       data.table = FALSE)
    
    if (nrow(coeff_chr) > 0) {

        system(paste0("./scripts/02-pgs/02-compute-score.sh ",
                      "../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                      "../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                      "../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam ",
                      "../data/sample-ids/filtered/all-ids.tab ",
                      "../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                      "../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                      "/coeff/p1_opt/chr", chr, "-coeff.tab ",
                      "../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex, "/pgs-all/p1_opt/chr", chr))

        if (isTRUE(ft)) {
            pgs_df <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                                   "/pgs-all/p1_opt/chr", chr, ".sscore"))
            system(paste0("rm ../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                          "/pgs-all/p1_opt/chr", chr, ".sscore"))

            ## rescale score by allele count
            pgs_df[, 5] <- pgs_df[, 5] * pgs_df$ALLELE_CT
            ft <- FALSE
        } else {
            pgs_df_chr <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                                       "/pgs-all/p1_opt/chr", chr, ".sscore"))
            system(paste0("rm ../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                          "/pgs-all/p1_opt/chr", chr, ".sscore"))

            pgs_df_chr[, 5] <- pgs_df_chr[, 5] * pgs_df_chr$ALLELE_CT
            pgs_df[, 5] <- pgs_df[, 5] + pgs_df_chr[, 5]
        }
    }
}

pgs_df <- pgs_df[, c(1, 2, 5)]
colnames(pgs_df)[3] <- "p1_opt_pgs"

## export
pgs_df$p1_opt_pgs <- round(pgs_df$p1_opt_pgs, digits = 8)
fwrite(pgs_df, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                             "/pgs-all/p1_opt/", phen, "-all-pgs0.tab"),
       sep = "\t", na = "NA", quote = FALSE)




## evaluate performance on training and test set and add to table
p1_df <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                      "/pgs-vali/p1_opt/p1_opt-rsq_vali.tab"),
               data.table = FALSE)


## load residual phenotypes: training
phen_train <- fread(paste0("../results/01-gwas/residuals-covars/", phen,
                             "/", phen, "-wb_", sex, "-resid-covars.tab"),
                      data.table = FALSE)
phen_train <- phen_train[, 2:3]
colnames(phen_train)[2] <- "phen_res"

## add PGS
phen_train_pgs <- merge(phen_train, pgs_df[, 2:3], by = "IID")

## regression
reg_train <- biglm(phen_res ~ p1_opt_pgs, data = phen_train_pgs)
p1_df$rsq_train <- round(summary(reg_train)$rsq, digits = 8)


## load residual phenotypes: test
phen_test <- fread(paste0("../results/01-gwas/residuals-covars/", phen,
                             "/", phen, "-other_gb_test_", sex, "-resid-covars.tab"),
                      data.table = FALSE)
phen_test <- phen_test[, 2:3]
colnames(phen_test)[2] <- "phen_res"

## add PGS
phen_test_pgs <- merge(phen_test, pgs_df[, 2:3], by = "IID")

## regression
reg_test <- biglm(phen_res ~ p1_opt_pgs, data = phen_test_pgs)
p1_df$rsq_test <- round(summary(reg_test)$rsq, digits = 8)


## write new p1_opt-rsq_val table
p1_df <- p1_df[, c(1, 3, 2, 4)]
fwrite(p1_df, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                            "/pgs-all/p1_opt/", phen, "-pgs0-rsq.tab"),
       sep = "\t", na = "NA", quote = FALSE)

})
fwrite(list(time[3]), file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/wb_", sex,
                                    "/pgs-all/p1_opt/total-runtime.txt"),
       col.names = FALSE, quote = FALSE)
