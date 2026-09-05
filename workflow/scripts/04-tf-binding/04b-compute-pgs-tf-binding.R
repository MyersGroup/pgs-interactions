## Compute PGSs partitioned by TF binding sites

library(argparser)
library(data.table)
setDTthreads(2)
library(dplyr)

options(warn = 2)  # turn warnings into errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",     help = "Phenotype code",   nargs = 1)
p <- add_argument(p, "--samples",  help = "Sample code",      nargs = 1)
p <- add_argument(p, "--tf_batch", help = "TF batch indices", nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen
samples  <- argv$samples
tf_batch <- argv$tf_batch

source("scripts/misc/fn-load_geno.R")



## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## load original GWAS summary stats
orig_f <- paste0("../results/01-gwas/plink-output/", phen, "/", samples,
                 ".chr", seq(1, 22), ".", phen, ".res.cov.glm.linear.gz")
orig <- lapply(orig_f, function(x) fread(x, header = TRUE, data.table = FALSE))
orig <- do.call("rbind", orig)

## load list of independent LOCO interaction hits
indep_hits <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                           phen, "-", samples, "-loco-indep-hits.tab"), data.table = FALSE)
indep_hits_chrs <- sort(unique(indep_hits$CHR))

## load HOCOMOCO data
hocomoco_raw <- scan("../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt", what = "char")
hocomoco <- vector("list", length(grep(">", hocomoco_raw)))
starts <- grep(">", hocomoco_raw)
starts <- c(starts, length(hocomoco_raw) + 1)
for (i in 1:length(hocomoco)) {
    start <- starts[i] + 1
    end   <- starts[i + 1] - 1
    hocomoco[[i]]$name  <- hocomoco_raw[starts[i]]
    hocomoco[[i]]$nameh <- gsub(">", "", hocomoco_raw[starts[i]])
    hocomoco[[i]]$motif <- t(matrix(as.double(hocomoco_raw[start:end]), nrow = 4, byrow = TRUE))
}
namesho <- unlist(lapply(1:length(hocomoco), function(i) hocomoco[[i]]$nameh))
names(hocomoco) <- namesho

## load list with weights for each motif
load(paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-snps.RData"))

## extract batch indices
tf_batch_vec <- strsplit(tf_batch, split = "_")[[1]]
tf_batch_ind <- as.numeric(tf_batch_vec)
tf_batch_names <- namesho[seq(tf_batch_ind[1], tf_batch_ind[2])]


## compute PGSs
pgs_ls <- list()
for (motif in tf_batch_names) {

    ## load list of SNPs and weights
    motif_probs <- tf_probs_ls[[motif]]

    ## check that a set of SNPs is available for this motif
    if (class(motif_probs) == "character") {
        if (motif_probs == "Unavailable") {
            pgs_ls[[motif]] <- "Unavailable"
            next
        }
    }

    ## check that there are no duplicated tag IDs
    if (sum(duplicated(motif_probs$tag_ID)) > 0) {
        stop("There are duplicated tag SNPs.")
    }

    ## merge with original GWAS results to get BETA coeffs (and A1 allele)
    motif_coeff <- left_join(motif_probs,
                             orig[, c("ID", "A1", "BETA")],
                             by = c("tag_ID" = "ID"))

    ## compute weighted coefficient
    motif_coeff$BETA_w <- motif_coeff$BETA * motif_coeff$post_prob

    ## compute sub-PGS by chromosome
    motif_coeff$CHR <- as.numeric(sapply(strsplit(motif_coeff$tag_ID, split = ":"), "[[", 1))
    motif_coeff_chrs <- sort(unique(motif_coeff$CHR))
    motif_pgs <- data.frame(IID = wb_ids$IID)
    for (chr in motif_coeff_chrs) {

        ## subset coeff df
        motif_coeff_chr <- motif_coeff[motif_coeff$CHR == chr,]

        ## load genotypes
        chr_geno <- load_geno(motif_coeff_chr$tag_ID, motif_coeff_chr$A1, wb_ids$IID)
        
        ## compute score
        motif_pgs[, paste0("chr", chr)] <- as.matrix(chr_geno[, -1, drop = FALSE]) %*% motif_coeff_chr$BETA_w
    }

    ## compute full PGS
    motif_pgs$pgs_full <- rowSums(motif_pgs[, paste0("chr", motif_coeff_chrs)])

    ## compute LOCO PGSs
    for (chr in 1:22) {

        if (length(motif_coeff_chrs) > 1) {

            if (chr %in% motif_coeff_chrs) {
                all_other_chrs <- motif_coeff_chrs[motif_coeff_chrs != chr]
                motif_pgs[, paste0("pgs_loco", chr)] <- rowSums(motif_pgs[, paste0("chr", all_other_chrs), drop = FALSE])
            } else {
                motif_pgs[, paste0("pgs_loco", chr)] <- motif_pgs$pgs_full
            }

        } else {

            if (chr == motif_coeff_chrs) {
                ## if all SNPs are in the same chromosome, we can't compute a LOCO PGS for that chromosome 
                motif_pgs[, paste0("pgs_loco", chr)] <- NA
            } else {
                motif_pgs[, paste0("pgs_loco", chr)] <- motif_pgs$pgs_full
            }
        }
    }

    ## store in list, keeping only LOCO PGSs for chrs for which we have an LHS hit to save disk space
    pgs_ls[[motif]] <- motif_pgs[, c("IID", "pgs_full", paste0("pgs_loco", indep_hits_chrs))]
}

save(pgs_ls,
     file = paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-", tf_batch, ".RData"))
