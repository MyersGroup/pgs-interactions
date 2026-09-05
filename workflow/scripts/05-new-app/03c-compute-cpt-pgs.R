library(argparser)
library(data.table)
setDTthreads(2)
library(dplyr)
library(pgenlibr)

options(warn = 2)  # turn warnings into errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",     help = "Phenotype code",   nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen

samples <- "wb_all"

plink2  <- "/path/to/plink2"


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
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## compute standard C+T (R^2 = 0.1) LOCO PGS
pgs_loco_df <- data.frame(IID = wb_ids$IID)
pgs_chrs <- c()
for (chr in 1:22) {

    ## read coefficient data frame to check if empty
    coeff_chr <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_0.1-kb_500/", samples,
                              "/coeff/p1_opt/chr", chr, "-coeff.tab"),
                       data.table = FALSE)
    coeff_chr <- coeff_chr[coeff_chr[, 5] != 0,]
    
    if (nrow(coeff_chr) > 0) {

        pgs_chrs <- c(pgs_chrs, chr)
        system(paste0(plink2,
                      " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                      "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                      "--psam ../data/sample-ids/ukb-103076-imp-auto-s486989.psam ",
                      "--extract ../data/variant-ids/chr", chr, "-var-ids.tab ",
                      "--read-freq ../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                      "--keep ../data/sample-ids/filtered/", samples, "-ids.tab ",
                      "--score ../results/02-pgs/initial/", phen, "/r2_0.1-kb_500/", samples, "/coeff/p1_opt/chr", chr, "-coeff.tab 1 2 5 header-read ",
                      "--threads 4 ",
                      "--memory 64000 ",
                      "--out ../results/05-new-app/initial/", phen, "/pgs-chr", chr))
        pgs_chr <- fread(paste0("../results/05-new-app/initial/", phen, "/pgs-chr", chr, ".sscore"), data.table = FALSE)

        ## rescale score by allele count
        pgs_chr[, 5] <- pgs_chr[, 5] * pgs_chr$ALLELE_CT
        colnames(pgs_chr)[5] <- paste0("chr", chr)

        ## add to data.frame
        pgs_loco_df <- left_join(pgs_loco_df, pgs_chr[, c("IID", paste0("chr", chr))], by = "IID")

        ## remove temporary files
        system(paste0("rm ../results/05-new-app/initial/", phen, "/pgs-chr", chr, ".sscore"))
        system(paste0("rm ../results/05-new-app/initial/", phen, "/pgs-chr", chr, ".log"))

    } else {

        pgs_loco_df[, paste0("chr", chr)] <- 0

    }
}

if (length(pgs_chrs) == 1) {
    for (chr in 1:22) {
        if (chr == pgs_chrs) {
            pgs_loco_df[, paste0("pgs_loco", chr)] <- NA
        } else {
            pgs_loco_df[, paste0("pgs_loco", chr)] <- pgs_loco_df[, paste0("chr", pgs_chrs)]
        }
    }    
} else {
    for (chr in 1:22) {
        pgs_loco_df[, paste0("pgs_loco", chr)] <- rowSums(pgs_loco_df[, paste0("chr", setdiff(pgs_chrs, chr))])
    }
}
pgs_loco_df <- pgs_loco_df[, c("IID", paste0("pgs_loco", 1:22))]
save(pgs_loco_df,
     file = paste0("../results/05-new-app/initial/", phen, "/pgs-loco.RData"))
