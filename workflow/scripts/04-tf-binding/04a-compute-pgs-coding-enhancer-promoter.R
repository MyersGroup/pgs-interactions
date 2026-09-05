## Compute PGSs partitioned by coding/H3K4me1/3 annotations

library(argparser)
library(data.table)
setDTthreads(2)
library(dplyr)

options(warn = 2)  # turn warnings into errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples", help = "Sample code",    nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples

source("scripts/misc/fn-load_geno.R")



## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## load original GWAS summary stats
orig_f <- paste0("../results/01-gwas/plink-output/", phen, "/", samples,
                 ".chr", seq(1, 22), ".", phen, ".res.cov.glm.linear.gz")
orig <- lapply(orig_f, function(x) fread(x, header = TRUE, data.table = FALSE))
orig <- do.call("rbind", orig)



## compute PGSs
for (anno in c("coding", "h3k4me1", "h3k4me3")) {

    ## load list of SNPs and weights
    load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/", anno, "-pgs-snps.RData"))
    if (anno == "coding") {
        anno_probs <- coding_probs
    } else if (anno == "h3k4me1") {
        anno_probs <- h3k4me1_probs
    } else if (anno == "h3k4me3") {
        anno_probs <- h3k4me3_probs
    }

    ## check that there are no duplicated tag IDs
    if (sum(duplicated(anno_probs$tag_ID)) > 0) {
        stop("There are duplicated tag SNPs.")
    }

    ## merge with original GWAS results to get BETA coeffs (and A1 allele)
    anno_coeff <- left_join(anno_probs,
                            orig[, c("ID", "A1", "BETA")],
                            by = c("tag_ID" = "ID"))

    ## compute weighted coefficients
    anno_coeff$BETA_w <- anno_coeff$BETA * anno_coeff$post_prob

    ## compute sub-PGS by chromosome
    anno_coeff$CHR <- as.numeric(sapply(strsplit(anno_coeff$tag_ID, split = ":"), "[[", 1))
    anno_coeff_chrs <- sort(unique(anno_coeff$CHR))
    anno_pgs <- data.frame(IID = wb_ids$IID)
    for (chr in anno_coeff_chrs) {

        ## subset coeff df
        anno_coeff_chr <- anno_coeff[anno_coeff$CHR == chr,]

        ## load genotypes
        chr_geno <- load_geno(anno_coeff_chr$tag_ID, anno_coeff_chr$A1, wb_ids$IID)
        
        ## compute score
        anno_pgs[, paste0("chr", chr)] <- as.matrix(chr_geno[, -1]) %*% anno_coeff_chr$BETA_w
    }

    ## compute full PGS
    anno_pgs$pgs_full <- rowSums(anno_pgs[, paste0("chr", anno_coeff_chrs)])

    ## compute LOCO PGSs
    for (chr in 1:22) {

        if (length(anno_coeff_chrs) > 1) {

            if (chr %in% anno_coeff_chrs) {
                all_other_chrs <- anno_coeff_chrs[anno_coeff_chrs != chr]
                anno_pgs[, paste0("pgs_loco", chr)] <- rowSums(anno_pgs[, paste0("chr", all_other_chrs), drop = FALSE])
            } else {
                anno_pgs[, paste0("pgs_loco", chr)] <- anno_pgs$pgs_full
            }

        } else {

            if (chr == anno_coeff_chrs) {
                ## if all SNPs are in the same chromosome, we can't compute a LOCO PGS for that chromosome 
                anno_pgs[, paste0("pgs_loco", chr)] <- NA
            } else {
                anno_pgs[, paste0("pgs_loco", chr)] <- anno_pgs$pgs_full
            }

        }
    }
    
    anno_pgs <- anno_pgs[, c("IID", "pgs_full", paste0("pgs_loco", 1:22))]
    save(anno_pgs,
         file = paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/", anno, "-pgs.RData"))
}
