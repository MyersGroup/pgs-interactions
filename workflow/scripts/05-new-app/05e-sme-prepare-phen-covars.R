## desc: Prepare phenotype file and run SME

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

plink2 <- "/path/to/plink2"



## prepare phenotype file with only relevant samples

## load fam file showing samples in order in genetic data files
fam <- fread(paste0("../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_all_chr.fam"),
             data.table = FALSE)
colnames(fam)[2] <- "IID"

## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## load raw (non-QN) phenotype
phen_raw <- fread(paste0("../data/phenotypes/clean/", substring(phen, 1, nchar(phen) - 3), ".tab"), data.table = FALSE)
colnames(phen_raw) <- c("IID", "phen")
phen_raw <- phen_raw[!is.na(phen_raw$phen),]
phen_raw$FID <- phen_raw$IID
phen_raw <- phen_raw[, c(3, 1, 2)]

## keep only WB samples
phen_raw <- phen_raw[phen_raw$IID %in% wb_ids$IID,]

## load covariates from most recent main UKB dataset
covars <- fread("../data/phenotypes/ukb-raw/ukb677925.tab",
                select = c("f.eid", "f.21022.0.0", "f.31.0.0"),
                data.table = FALSE)
colnames(covars) <- c("IID", "age", "sex")

## regress out age and sex
phen_cov <- left_join(phen_raw, covars, by = "IID")
reg <- lm(phen ~ age + sex, data = phen_cov)
resid_df <- data.frame(FID = phen_cov$FID,
                       IID = phen_cov$IID,
                       phen = round(residuals(reg), digits = 8))
                       
## order by fam file
resid_df <- resid_df[order(match(resid_df$IID, fam$IID)),]
fwrite(resid_df,
       file = paste0("../results/05-new-app/sme/output/", phen, "/", phen, "-", snp_, ".pheno"),
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')

## remove samples with missing phenotype values from bed/bim/fam
samples_to_keep <- resid_df[, c("FID", "IID")]
colnames(samples_to_keep)[1] <- "#FID"
fwrite(samples_to_keep,
       file = paste0("../results/05-new-app/sme/output/", phen, "/", phen, "-", snp_, ".samples"),
       row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', na = 'NA')
system(paste0(plink2,
              " --bfile ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_all_chr ",
              "--keep ../results/05-new-app/sme/output/", phen, "/", phen, "-", snp_, ".samples ",
              "--make-bed ",
              "--threads 4 ",
              "--memory 60000 ",
              "--out ../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_all_flt_samples_chr"))
