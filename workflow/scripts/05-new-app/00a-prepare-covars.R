## desc: Prepare covars table with the following variables:
##       - age (field 21022 - Age at recruitment)
#          >> should be standardised (mean centred w/ unit variance)
##       - sex (field 31)
##       - age * sex, age^2, age^2 * sex
##       - batch dummy variables (field 22000 - Genotype measurement batch)
##       - ancestry components (previously computed by Sile)


library(data.table)
library(dplyr)



## load most recent main UKB dataset
main <- fread('../data/phenotypes/ukb-raw/ukb677925.tab',
              select = c('f.eid', 'f.21022.0.0', 'f.31.0.0', 'f.22000.0.0', 'f.54.0.0'),
              data.table = FALSE)

## load list of withdrawn samples (as of 4 Mar 2024)
withdrawn <- fread('../data/sample-ids/withdrawals/withdraw103076_179_20240301.txt', data.table = FALSE)

## load most recent imputed data sample file (downloaded 4 Mar 2024)
imp_sp <- fread('../data/sample-ids/ukb-103076-imp-auto-s486989.psam', data.table = FALSE)




## keep only samples:
## - who haven't withdrawn
## - for whom imputed genotypes are available
## - who don't have missing fields
covars <- main %>%
    filter(!(f.eid %in% withdrawn$V1),
           f.eid %in% imp_sp$IID,
           !is.na(f.21022.0.0),
           !is.na(f.31.0.0),
           !is.na(f.22000.0.0),
           !is.na(f.54.0.0))




## age
colnames(covars)[2] <- 'age'
covars$age <- scale(covars$age)

## sex
colnames(covars)[3] <- 'sex'

## add non-linear terms: age*sex + age^2 + age^2*sex
covars$age_t_sex    <- covars$age * covars$sex
covars$age_sq       <- covars$age^2
covars$age_sq_t_sex <- covars$age_sq * covars$sex

## re-order columns
covars <- covars %>%
    select(f.eid, age, sex, age_t_sex, age_sq, age_sq_t_sex, f.22000.0.0, f.54.0.0)




## batch dummy variables
batch_ <- factor(covars$f.22000.0.0)
batch_dummy <- as.data.frame(model.matrix(~ batch_ + 0))
colnames(batch_dummy) <- gsub('_-', '_.', colnames(batch_dummy))

## first value dummy (-11 = UKBiLEVEAX_b11) is removed to avoid multicollinearity
covars <- cbind(subset(covars, select=f.eid:age_sq_t_sex),
                batch_dummy[, 2:ncol(batch_dummy)],
                subset(covars, select=f.54.0.0))




## assessment centre
## (we take the first value which corresponds to the centre of the initial assessment visit)
centre_ <- factor(covars$f.54.0.0)
centre_dummy <- as.data.frame(model.matrix(~ centre_ + 0))

## first value dummy (10003 = Stockport (pilot)) is removed to avoid multicollinearity
covars <- cbind(subset(covars, select=f.eid:batch_95),
                centre_dummy[, 2:ncol(centre_dummy)])




## ancestry components
anc <- readRDS('../data/covars/ancestry/v2_487409.rds')

## most common AC (G_Merseyside) is removed to avoid multicollinearity
anc <- anc[, c(1:7, 9:ncol(anc))]

## collate (samples are in same order as imputed sample file)
anc_iid <- cbind(imp_sp[, "IID", drop = FALSE], anc)

## join
covars <- left_join(covars, anc_iid, by = c('f.eid' = 'IID'))




## reorder by original sample file
covars <- covars[match(imp_sp$IID[imp_sp$IID %in% covars$f.eid], covars$f.eid),]
colnames(covars)[1] <- 'IID'
## add FID as first column
covars$'#FID' <- covars$IID
covars <- covars[, c(ncol(covars), 1:(ncol(covars) - 1))]




## export
fwrite(covars, file = '../data/covars/age-sex-batch-centre-ac-all.tab',
       sep = '\t', na = 'NA', quote = FALSE)
