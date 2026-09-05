## desc: Function to load genotypes from PGEN files
##
##       Takes as input:
##       - `var_ids`: character vector of SNP IDs (must match PVAR file);
##       - `a1`: character vector of effect alleles for each SNP;
##       - `sample_ids_filtered`: vector of sample IIDs.
##
##       Also requires a PSAM sample file to be loaded.


library(data.table)
library(dplyr)
library(pgenlibr)



## load sample file
sample_f <- fread("../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam", data.table = FALSE)

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
