## desc: Test whether SNP*TF-PGS interactions merely reflect pairwise
##       interactions with a SNP in the PGS

library(data.table)
setDTthreads(8)
library(plyr)
library(dplyr)
library(pgenlibr)

options(warn=2)  # turn warnings into errors

samples <- "wb_all"

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

## load list of independent TF interactions
load("../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.RData")

## prepare table to store new results
indep_hits_agg_retest <- indep_hits_agg[, c(1:6, 12)]
indep_hits_agg_retest <- cbind(indep_hits_agg_retest,
                               data.frame(matrix(NA, nrow = nrow(indep_hits_agg_retest), ncol = 13)))
colnames(indep_hits_agg_retest)[8:20] <- c("n_tag_snps",
                                           "snp_BETA", "snp_LOG10_P",
                                           "snp_sq_BETA", "snp_sq_LOG10_P",
                                           "pgs_loco_BETA", "pgs_loco_LOG10_P",
                                           "int_pgs_loco_BETA", "int_pgs_loco_LOG10_P",
                                           "pgs_tf_BETA", "pgs_tf_LOG10_P",
                                           "int_pgs_tf_BETA", "int_pgs_tf_LOG10_P")

## prepare list to store removed driver SNPs
removed_snps <- vector(mode = "list", length = nrow(indep_hits_agg))

for (i in 1:nrow(indep_hits_agg)) {

    phen   <- indep_hits_agg$field_id[i]

    hit_id  <- indep_hits_agg$ID[i]
    hit_a1  <- indep_hits_agg$A1[i]
    hit_chr <- as.numeric(sapply(strsplit(indep_hits_agg$ID[i], split = ":"), "[[", 1))

    motif  <- indep_hits_agg$motif[i]


    ## load list of TF PGS SNPs and weights for this phenotype
    load(paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-snps.RData"))

    ## get SNPs and weights for this motif
    motif_probs <- tf_probs_ls[[motif]]

    ## remove SNPs from hit SNP's chromosome
    motif_probs$CHR <- as.numeric(sapply(strsplit(motif_probs$hit_ID, split = ":"), `[[`, 1))
    motif_probs <- motif_probs[motif_probs$CHR != hit_chr,]


    ## load pairwise GWAS results
    pairwise_f <- paste0("../results/03-interaction-gwas/indep-hits/", phen,
                             "/snp-pgs/gwas/", phen, ".", samples, ".", gsub(":", "_", hit_id), ".",
                             "int.loco.chr", setdiff(1:22, hit_chr), ".", phen, ".res.cov.glm.linear.gz")
    pairwise <- lapply(pairwise_f, function(x) fread(x, header = TRUE, data.table = FALSE))
    pairwise <- do.call("rbind", pairwise)

    ## keep only SNP*SNP term
    pairwise <- pairwise[pairwise$TEST == paste0("ADDx", hit_id, "_", hit_a1),]

    ## add pairwise interaction p-values to TF PGS SNPs/coeffs df
    motif_probs <- left_join(motif_probs, pairwise[, c("ID", "LOG10_P")], by = c("tag_ID" = "ID"))

    ## save data.frame of TF PGS SNPs/coeffs df to then check if any are removed
    motif_probs_backup <- motif_probs

    ## remove any tag SNP with log10(pv) >= 5 and all other tags of the corresponding hit SNP
    for (hit_snp in unique(motif_probs_backup$hit_ID)) {

        ## retrieve highest p-value across all its tags
        max_pv <- max(motif_probs$LOG10_P[motif_probs$hit_ID == hit_snp])

        if (max_pv >= 5) {
            motif_probs <- motif_probs[motif_probs$hit_ID != hit_snp,]
            removed_snps[[i]] <- c(removed_snps[[i]], hit_snp)
        } else {
            next
        }
    }


    ## if any SNP was removed, repeat TF-PGS interaction test
    if (nrow(motif_probs) < nrow(motif_probs_backup)) {

        ## load original GWAS summary stats
        orig_f <- paste0("../results/01-gwas/plink-output/", phen, "/", samples,
                         ".chr", seq(1, 22), ".", phen, ".res.cov.glm.linear.gz")
        orig <- lapply(orig_f, function(x) fread(x, header = TRUE, data.table = FALSE))
        orig <- do.call("rbind", orig)

        ## merge with original GWAS results to get BETA coeffs (and A1 allele)
        motif_coeff <- left_join(motif_probs,
                                 orig[, c("ID", "A1", "BETA")],
                                 by = c("tag_ID" = "ID"))

        ## compute weighted coefficient
        motif_coeff$BETA_w <- motif_coeff$BETA * motif_coeff$post_prob

        ## count number of tag SNPs not in SNP's chromosome
        indep_hits_agg_retest$n_tag_snps[i] <- nrow(motif_coeff)


        ## compute sub-PGS by chromosome
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

        ## compute LOCO PGS for LHS hit's chromosome
        ## (since we had removed all SNPs in same chromosome as LHS hit, we just need to sum all other chromosomes)
        motif_pgs[, paste0("pgs_loco", hit_chr)] <- rowSums(motif_pgs[, paste0("chr", motif_coeff_chrs)])


        ## load simple INT phenotype
        resid_int <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-", samples, "-resid-covars.tab"),
                           data.table = FALSE)

        ## add genotype of LHS SNP
        geno <- load_geno(hit_id, hit_a1, wb_ids$IID)
        main <- left_join(resid_int[, 2:3], geno, by = "IID")
        colnames(main)[2:3] <- c("phen", "geno")

        ## check whether correlation between genotype and its square is < 0.999
        sq_cor_low <- (cor(geno[, 2], geno[, 2]^2) < 0.999)
        if (sq_cor_low) main$geno_sq <- main$geno^2

        ## add standard C+T (R^2 = 0.1) LOCO PGS
        load(paste0("../results/05-new-app/initial/", phen, "/pgs-loco.RData"))
        main <- left_join(main, pgs_loco_df[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
        colnames(main)[length(colnames(main))] <- "pgs_loco"

        ## add TF PGS
        main <- left_join(main, motif_pgs[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
        colnames(main)[ncol(main)] <- "pgs_tf"


        ## run regression
        if (sq_cor_low) {
            reg <- lm(phen ~ geno + geno_sq + geno*pgs_loco + geno*pgs_tf, data = main)
        } else {
            reg <- lm(phen ~ geno + geno*pgs_loco + geno*pgs_tf, data = main)
        }

        ## add sumstats to data.frame
        coeffs <- summary(reg)$coeff
        indep_hits_agg_retest$snp_BETA[i]             <- coeffs["geno", 1]
        indep_hits_agg_retest$snp_LOG10_P[i]          <- - (pt(abs(coeffs["geno", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))  # to avoid numerical underflow
        indep_hits_agg_retest$snp_sq_BETA[i]          <- ifelse(sq_cor_low, coeffs["geno_sq", 1], NA)
        indep_hits_agg_retest$snp_sq_LOG10_P[i]       <- ifelse(sq_cor_low,
                                                                - (pt(abs(coeffs["geno_sq", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                                                NA)
        indep_hits_agg_retest$pgs_loco_BETA[i]        <- coeffs["pgs_loco", 1]
        indep_hits_agg_retest$pgs_loco_LOG10_P[i]     <- - (pt(abs(coeffs["pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
        indep_hits_agg_retest$int_pgs_loco_BETA[i]    <- coeffs["geno:pgs_loco", 1]
        indep_hits_agg_retest$int_pgs_loco_LOG10_P[i] <- - (pt(abs(coeffs["geno:pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
        indep_hits_agg_retest$pgs_tf_BETA[i]          <- coeffs["pgs_tf", 1]
        indep_hits_agg_retest$pgs_tf_LOG10_P[i]       <- - (pt(abs(coeffs["pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))
        indep_hits_agg_retest$int_pgs_tf_BETA[i]      <- coeffs["geno:pgs_tf", 1]
        indep_hits_agg_retest$int_pgs_tf_LOG10_P[i]   <- - (pt(abs(coeffs["geno:pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2))

    } else {
        next
    }
}

## export
indep_hits_agg_retest_out <- indep_hits_agg_retest
cols <- colnames(indep_hits_agg_retest_out)[9:20]

indep_hits_agg_retest_out <- as.data.table(indep_hits_agg_retest_out)
indep_hits_agg_retest_out[,(cols) := round(.SD, 8), .SDcols = cols]
indep_hits_agg_retest_out <- as.data.frame(indep_hits_agg_retest_out)

fwrite(indep_hits_agg_retest_out,
       file = "../results/05-new-app/tf-gwas/aggregate/int-gwas-indep-hits-aggregate-retest-pairwise.tab",
       sep = "\t", na = "NA", quote = FALSE)
save(indep_hits_agg_retest, removed_snps,
     file = "../results/05-new-app/tf-gwas/aggregate/int-gwas-indep-hits-aggregate-retest-pairwise.RData")
