## desc: Run a linear regression of quantitative phenotype on all the covars
##       in the file data/covars/age-sex-nonlin-ac-batch-all.tab
##
##       Standardise (zero-centred, unit variance) the phenotypes before
##       running the regression.
##
##       Export the residuals.

library(argparser)
library(bestNormalize)
library(data.table)


p <- arg_parser('Argument parser')
p <- add_argument(p, '--phen',    help = 'Phenotype code', nargs = 1)
p <- add_argument(p, '--samples', help = 'Sample code',    nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples




## load data
phen_df <- fread(paste0("../data/phenotypes/clean/", phen, ".tab"), data.table = FALSE)
colnames(phen_df)[1] <- 'IID'
phen_name <- colnames(phen_df)[2]

covars <- fread('../data/covars/age-sex-batch-centre-ac-all.tab', data.table = FALSE)

sample <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)
sample_vec <- sample$IID




## keep only samples in sample file
phen_sp   <- phen_df[phen_df$IID %in% sample_vec,]
covars_sp <- covars[covars$IID %in% sample_vec,]
covars_sp <- covars_sp[, 2:ncol(covars_sp)]  # remove FID column

## merge phen and covars
phen_covars <- merge(phen_sp, covars_sp, by = "IID")

## remove samples for whom the phenotype is NA
phen_covars_nonmiss <- phen_covars[!is.na(phen_covars[,2]),]

## standardise phenotype
phen_covars_nonmiss[, 2] <- as.numeric(scale(phen_covars_nonmiss[, 2]))




## regression
regress <- lm(as.formula(paste(phen, "~ .")), data = phen_covars_nonmiss[, -1])

## extract residuals
resid <- residuals(regress)
resid_df <- data.frame(IID = phen_covars_nonmiss$IID, res = resid)
if (sum(is.na(resid_df$res)) > 0) {
    stop("Some residuals are missing.")
}




## quantile normalisation
resid_df_qn <- resid_df
resid_df_qn$res <- orderNorm(resid_df_qn$res)$x.t

## add IIDs and samples with missing values
resid_df_qn_all <- merge(phen_sp, resid_df_qn, by = 'IID', all.x = TRUE)
resid_df_qn_all <- resid_df_qn_all[order(match(resid_df_qn_all$IID, sample_vec)),]  # same order as initially

resid_df_qn_all <- resid_df_qn_all[, c(1,3)]
resid_df_qn_all[, 2] <- round(resid_df_qn_all[, 2], digits = 8)

## add FIDs
resid_df_qn_all$FID <- resid_df_qn_all$IID
resid_df_qn_all <- resid_df_qn_all[, c('FID', 'IID', 'res')]
colnames(resid_df_qn_all) <- c('#FID', 'IID', paste0(phen_name, '_qn.res.cov'))

## export
fwrite(resid_df_qn_all,
       file = paste0("../results/01-gwas/residuals-covars/", phen, "_qn/", phen, "_qn-", samples, "-resid-covars.tab"),
       sep = '\t', na = 'NA', quote = FALSE)
