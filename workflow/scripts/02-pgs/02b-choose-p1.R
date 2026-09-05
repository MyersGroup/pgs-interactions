## desc: Choose optimum global p-value threshold


library(argparser)
library(biglm)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("scripts/misc/fn-ggp1.R")


p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen", help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--train_sp", help = "Training sample code", nargs = 1)
p <- add_argument(p, "--vali_sp", help = "Validation sample code", nargs = 1)
p <- add_argument(p, "--r2", help = "LD clumping r2 parameter", nargs = 1)
p <- add_argument(p, "--kb", help = "LD clumping kb parameter", nargs = 1)
argv <- parse_args(p)

phen         <- argv$phen
train_sp     <- argv$train_sp
vali_sp      <- argv$vali_sp
r2           <- argv$r2
kb           <- argv$kb





time <- system.time({

## load residual phenotypes to check how well PGS predict them
phen_df <- fread(paste0("../results/01-gwas/residuals-covars/", phen,
                        "/", phen, "-", vali_sp, "-resid-covars.tab"),
                 data.table = FALSE)
phen_df <- phen_df[, 2:3]
colnames(phen_df)[2] <- "phen_res"
## remove samples for whom the phenotype is NA
phen_df <- phen_df[!is.na(phen_df[,2]),]

## df to store params and performance for all scores
pgs_comp <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                         "/coeff/p1_grid/pgs-p1-n_snps.tab"),
                  data.table = FALSE)






## load chr1 score components
pgs_df <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                       "/pgs-vali/p1_grid/chr1.sscore"))
ncols <- ncol(pgs_df)
## rescale components by allele count
pgs_df[, 5:ncols] <- pgs_df[, 5:ncols] * pgs_df$ALLELE_CT

## add remaining components
for (chr in 2:22) {
    pgs_df_chr <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                               "/pgs-vali/p1_grid/chr", chr, ".sscore"))
    pgs_df_chr[, 5:ncols] <- pgs_df_chr[, 5:ncols] * pgs_df_chr$ALLELE_CT
    pgs_df[, 5:ncols] <- pgs_df[, 5:ncols] + pgs_df_chr[, 5:ncols]
}






## compute performance and store

## remove individuals with missing phenotype and reorder
pgs_df <- pgs_df[pgs_df$IID %in% phen_df$IID,]
pgs_df <- pgs_df[order(match(pgs_df$IID, phen_df$IID)),]

## confirm same IIDs
if (!isTRUE(all.equal(phen_df$IID, pgs_df$IID))) {
    stop("Different samples in the phenotype and PGS data.")
}


## regression function
reg_col <- function(x) {
    reg_df <- data.frame(phen = phen_df[, 2],
                         pgs  = x)
    regress <- biglm(phen ~ pgs, data = reg_df)
    return(summary(regress)$rsq)
}

pgs_comp$rsq <- apply(pgs_df[, 5:ncols], 2, reg_col)

## export comparison table
pgs_comp_out     <- pgs_comp
pgs_comp_out$rsq <- round(pgs_comp_out$rsq, digits = 8)
fwrite(pgs_comp_out, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                   "/pgs-vali/p1_grid/", vali_sp, "-pgs-comparison.tab"),
       sep = "\t", na = "NA", quote = FALSE)






## plot R2 as a function of score size

## chose optimum p-value threshold
rsq_rg  <- max(pgs_comp$rsq) - pgs_comp$rsq[1]
opt_col <- min(which(pgs_comp$rsq - pgs_comp$rsq[1] >= rsq_rg))
p1_opt  <- pgs_comp$p1[opt_col]
rsq_opt <- pgs_comp$rsq[opt_col]

## get phenotype description
all_desc  <- fread("../data/phenotypes/clean/code-desc-map-all.tab", data.table = FALSE)
phen_desc <- all_desc$description[which(all_desc$field_id == phen)]


## plot
p1_plot <- ggp1(pgs_comp, p1_opt = p1_opt, rsq_opt = rsq_opt,
                title = paste("PGS validation R-squared by p-value threshold:", phen_desc))
ggsave(p1_plot, filename = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                  "/figs/", phen, "-", vali_sp, "-p1_rsq.png"),
       type = "cairo-png", width = 700/120, height = 350/120, units = "in", dpi = 120)






## save optimal p1 and validation R-squared
p1_df <- data.frame(p1_opt = p1_opt,
                    rsq_vali = round(rsq_opt, digits = 8))
fwrite(p1_df, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                            "/pgs-vali/p1_opt/p1_opt-rsq_vali.tab"),
       sep = "\t", na = "NA", quote = FALSE)


## export optimal p-value score coefficients
for (chr in 1:22) {
    coeff_chr <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                              "/coeff/p1_grid/chr", chr, "-coeff.tab"),
                       data.table = FALSE)

    ## keep only optimal pv score
    coeff_chr <- coeff_chr[, c(1:4, 4 + opt_col)]

    ## remove SNPs above p-value threshold
    coeff_chr <- coeff_chr %>%
        filter(LOG10_P >= -log10(p1_opt))

    fwrite(coeff_chr, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                    "/coeff/p1_opt/chr", chr, "-coeff.tab"),
           sep = "\t", na = "NA", quote = FALSE)
}

})
fwrite(list(time[3]), file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                    "/coeff/p1_opt/choose-p1-runtime.txt"),
       col.names = FALSE, quote = FALSE)
