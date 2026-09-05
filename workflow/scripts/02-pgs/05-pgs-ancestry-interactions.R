## Run regression with PGS*ancestry component interactions and get residuals

library(argparser)
library(data.table)
library(dplyr)


p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--train_sp", help = "Training samples",         nargs = 1)
p <- add_argument(p, "--r2",      help = "LD clumping r2 parameter", nargs = 1)
p <- add_argument(p, "--kb",      help = "LD clumping kb parameter", nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen
train_sp <- argv$train_sp
r2       <- argv$r2
kb       <- argv$kb



## load residuals after regressing out covars
phen_res_cov_train <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/", phen, "-", train_sp, "-resid-covars.tab"),
                            data.table = FALSE)
phen_res_cov_train <- phen_res_cov_train[!is.na(phen_res_cov_train[, 3]),]
phen_res_cov_train <- phen_res_cov_train[, c(2, 3)]  # remove FID col

covars <- fread("../data/covars/age-sex-batch-centre-ac-all.tab", data.table = FALSE)
covars <- covars[, -1]  # remove FID col

train_sp_df <- fread(paste0("../data/sample-ids/filtered/", train_sp, "-ids.tab"), data.table = FALSE)

## merge
data <- merge(phen_res_cov_train, covars, by = "IID")

## load full PGS (mean-corrected)
pgs_full_mc <- fread(paste0("../results/02-pgs/mean-corrected/", phen,
                            "/r2_", r2, "-kb_", kb, "/", phen, "-", train_sp, "-pgs-full_mc.tab"),
                     data.table = FALSE)
pgs_full_mc <- pgs_full_mc[, 2:3]

## add PGS to covars table
data_full <- merge(data, pgs_full_mc, by = "IID")




## add PGS*AC interactions
ac_st <- which(colnames(data_full) == "G_Anglia")
ac_en <- which(colnames(data_full) == "W_Philippines")

for (col in ac_st:ac_en) {
    data_full[, paste0("pgs_full_mc_x_", colnames(data_full)[col])] <- data_full[[colnames(data_full)[col]]] * data_full$pgs_full_mc
}

## add PGS*assessment centre interactions
centre_st <- which(colnames(data_full) == "centre_11001")
centre_en <- which(colnames(data_full) == "centre_11023")

for (col in centre_st:centre_en) {
    data_full[, paste0("pgs_full_mc_x_", colnames(data_full)[col])] <- data_full[[colnames(data_full)[col]]] * data_full$pgs_full_mc
}




## regression
regress <- lm(as.formula(paste0(phen, ".res.cov ~ .")), data = data_full[, -1])

## save summary stats
sink(paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_", r2, "-kb_", kb, "/",
            phen, "-", train_sp, "-pgs-full_mc-ac-centre-sumstats.tab"))
summary(regress)
sink()

## export residuals
resid_df <- data.frame(IID = data_full$IID,
                       phen_res = resid(regress))

## add IIDs and samples with missing values
resid_df_all <- merge(train_sp_df[, "IID", drop=F], resid_df, by = "IID", all.x = TRUE)
resid_df_all <- resid_df_all[order(match(resid_df_all$IID, train_sp_df$IID)),]  # same order as initially
resid_df_all[, 2] <- round(resid_df_all[, 2], digits = 8)

resid_df_all$FID <- resid_df_all$IID
resid_df_all <- resid_df_all[, c("FID", "IID", "phen_res")]
colnames(resid_df_all) <- c("#FID", "IID", paste0(phen, ".res.cov.pgs_full_mc.ac.ctr"))

## export
fwrite(resid_df_all,
       file = paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_", r2, "-kb_", kb, "/",
                     phen, "-", train_sp, "-resid-covars-pgs-full_mc-ac-centre.tab"),
       sep = "\t", na = "NA", quote = FALSE)


## compute and export absolute residuals
resid_abs_df_all <- resid_df_all
resid_abs_df_all[, 3] <- abs(resid_abs_df_all[, 3])
names(resid_abs_df_all)[3] <- paste0(phen, ".res_abs.cov.pgs_full_mc.ac.ctr")

fwrite(resid_abs_df_all,
       file = paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_", r2, "-kb_", kb, "/",
                     phen, "-", train_sp, "-resid_abs-covars-pgs-full_mc-ac-centre.tab"),
       sep = "\t", na = "NA", quote = FALSE)




## LOCO PGS
for (chr in 1:22) {
    ## load LOCO PGS (mean-corrected)
    pgs_loco_mc <- fread(paste0("../results/02-pgs/mean-corrected/", phen,
                                "/r2_", r2, "-kb_", kb, "/", phen, "-", train_sp, "-pgs-loco", chr, "_mc.tab"),
                         data.table = FALSE)
    pgs_loco_mc <- pgs_loco_mc[, 2:3]

    ## add PGS to covars table
    data_loco <- merge(data, pgs_loco_mc, by = "IID")

 


    ## add PGS*AC interactions
    ac_st <- which(colnames(data_loco) == "G_Anglia")
    ac_en <- which(colnames(data_loco) == "W_Philippines")

    for (col in ac_st:ac_en) {
        data_loco[, paste0("pgs_loco", chr, "_mc_x_", colnames(data_loco)[col])] <- data_loco[[colnames(data_loco)[col]]] * data_loco[[paste0("pgs_loco", chr, "_mc")]]
    }

    ## add PGS*assessment centre interactions
    centre_st <- which(colnames(data_loco) == "centre_11001")
    centre_en <- which(colnames(data_loco) == "centre_11023")

    for (col in centre_st:centre_en) {
        data_loco[, paste0("pgs_loco", chr, "_mc_x_", colnames(data_loco)[col])] <- data_loco[[colnames(data_loco)[col]]] * data_loco[[paste0("pgs_loco", chr, "_mc")]]
    }




    ## regression
    regress <- lm(as.formula(paste0(phen, ".res.cov ~ .")), data = data_loco[, -1])

    ## save summary stats
    sink(paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_", r2, "-kb_", kb, "/",
                phen, "-", train_sp, "-pgs-loco", chr, "_mc-ac-centre-sumstats.tab"))
    print(summary(regress))
    sink()

    ## export residuals
    resid_df <- data.frame(IID = data_loco$IID,
                           phen_res = resid(regress))

    ## add IIDs and samples with missing values
    resid_df_all <- merge(train_sp_df[, "IID", drop=F], resid_df, by = "IID", all.x = TRUE)
    resid_df_all <- resid_df_all[order(match(resid_df_all$IID, train_sp_df$IID)),]  # same order as initially
    resid_df_all[, 2] <- round(resid_df_all[, 2], digits = 8)

    resid_df_all$FID <- resid_df_all$IID
    resid_df_all <- resid_df_all[, c("FID", "IID", "phen_res")]
    colnames(resid_df_all) <- c("#FID", "IID", paste0(phen, ".res.cov.pgs_loco", chr, "_mc.ac.ctr"))

    ## export
    fwrite(resid_df_all,
           file = paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_", r2, "-kb_", kb, "/",
                         phen, "-", train_sp, "-resid-covars-pgs-loco", chr, "_mc-ac-centre.tab"),
           sep = "\t", na = "NA", quote = FALSE)


    ## compute and export absolute residuals
    resid_abs_df_all <- resid_df_all
    resid_abs_df_all[, 3] <- abs(resid_abs_df_all[, 3])
    names(resid_abs_df_all)[3] <- paste0(phen, ".res_abs.cov.pgs_loco", chr, "_mc.ac.ctr")

    fwrite(resid_abs_df_all,
           file = paste0("../results/02-pgs/ancestry-interactions/", phen, "/r2_", r2, "-kb_", kb, "/",
                         phen, "-", train_sp, "-resid_abs-covars-pgs-loco", chr, "_mc-ac-centre.tab"),
           sep = "\t", na = "NA", quote = FALSE)
}
