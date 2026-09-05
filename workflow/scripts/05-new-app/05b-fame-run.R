## desc: Run FAME

library(argparser)
library(data.table)
library(bestNormalize)

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen", help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--snp",  help = "Target SNP",     nargs = 1)
argv <- parse_args(p)

phen <- argv$phen
snp_ <- argv$snp
snp  <- gsub("_", ":", snp_)

samples <- "wb_all"

genie_gxg <- "/path/to/FAME/build/GENIE_GxG"


## make phenotype file with only relevant samples

## load fam file showing samples in order
fam <- fread("../results/05-new-app/fame/genotype-calls/ukb_all_chr.fam",
             data.table = FALSE)
colnames(fam)[2] <- "IID"

## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## load raw (non-QN) phenotype
phen_raw <- fread(paste0("../data/phenotypes/clean/", substring(phen, 1, nchar(phen) - 3), ".tab"), data.table = FALSE)
colnames(phen_raw) <- c("IID", "phen")
phen_raw$FID <- phen_raw$IID
phen_raw <- phen_raw[, c(3, 1, 2)]

## keep only WB samples
phen_raw <- phen_raw[phen_raw$IID %in% wb_ids$IID,]

## quantile normalise without regressing out covariates
phen_qn <- phen_raw
phen_qn$phen_qn <- orderNorm(phen_qn$phen)$x.t

## order by fam file
phen_qn <- phen_qn[order(match(phen_qn$IID, fam$IID)),]

fwrite(phen_qn[, c("FID", "IID", "phen_qn")],
       file = paste0("../results/05-new-app/fame/output/", phen, "/", phen, "-", snp_, ".pheno"),
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')
      


## make covars file with age, sex and top 20 PCs

## load most recent main UKB dataset
covars <- fread("../data/phenotypes/ukb-raw/ukb677925.tab",
                select = c("f.eid", "f.21022.0.0", "f.31.0.0", paste0("f.22009.0.", seq(1, 20))),
                data.table = FALSE)
colnames(covars) <- c("IID", "age", "sex", paste0("pc", seq(1, 20)))
covars <- cbind("FID" = covars$IID, covars)

## keep only WB samples
covars <- covars[covars$IID %in% wb_ids$IID,]

## order by fam file
covars <- covars[order(match(covars$IID, fam$IID)),]

fwrite(covars,
       file = paste0("../results/05-new-app/fame/output/", phen, "/", phen, "-", snp_, ".covar"),
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')



## make annotation file: should include all SNPs in bed file other than those in same chromosome
chr <- as.numeric(sapply(strsplit(snp, split = ":"), `[[`, 1))
bim <- fread("../results/05-new-app/fame/genotype-calls/ukb_all_chr.bim",
             col.names = c("CHR", "rsID", "cM", "POS", "REF", "ALT"),
             data.table = FALSE)
bim$annot <- as.numeric(!(bim$CHR == chr))
fwrite(list(bim$annot),
       file = paste0("../results/05-new-app/fame/output/", phen, "/annot-", snp_, ".tab"),
       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')

## get index of SNP
bim$CHR_POS <- paste0(bim$CHR, "_", bim$POS)
snp_chr_pos <- paste(sapply(strsplit(snp, split = ":"), `[[`, 1), 
                     sapply(strsplit(snp, split = ":"), `[[`, 2), sep = "_")
snp_index <- which(bim$CHR_POS == snp_chr_pos)



## run FAME
system(paste0(genie_gxg,
              " -g ../results/05-new-app/fame/genotype-calls/ukb_all_chr ",
              " -p ../results/05-new-app/fame/output/", phen, "/", phen, "-", snp_, ".pheno ",
              " -c ../results/05-new-app/fame/output/", phen, "/", phen, "-", snp_, ".covar ",
              " -annot ../results/05-new-app/fame/output/", phen, "/annot-", snp_, ".tab ",
              " -gxgbin 0 ",
              " -snp ", snp_index,
              " -k 100 ",
              " -jn 100 ",
              " -o ../results/05-new-app/fame/output/", phen, "/results-", snp_, ".tab ",
              " > ../results/05-new-app/fame/output/", phen, "/results-", snp_, ".log"))
