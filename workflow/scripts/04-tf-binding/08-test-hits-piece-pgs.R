## desc: Test independent hits again now with redefined TF PGS which is
##       comprised only of SNPs that are part of the LOCO PGS (i.e., hit SNPs)

library(data.table)
library(dplyr)
library(jsonlite)

options(warn = 2)  # turn warnings into errors

samples <- "wb_all"

source("scripts/misc/fn-get_post.R")
source("scripts/misc/fn-load_geno.R")
source("scripts/misc/fn-test_enrich_mat.R")


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

## load table with all independent hits
load("../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.RData")

## prepare table to store new results
indep_hits_agg_retest <- indep_hits_agg[, c(1:12, 15:27)]
indep_hits_agg_retest[, 13:ncol(indep_hits_agg_retest)] <- NA
colnames(indep_hits_agg_retest)[c(13, 18:25)] <-
    c("n_hit_snps",
      "pgs_loco_diff_BETA",
      "pgs_loco_diff_LOG10_P",
      "int_pgs_loco_diff_BETA",
      "int_pgs_loco_diff_LOG10_P",
      "pgs_tf_hits_BETA",
      "pgs_tf_hits_LOG10_P",
      "int_pgs_tf_hits_BETA",
      "int_pgs_tf_hits_LOG10_P")


## loop over table rows grouped by phenotype
for (phen in unique(indep_hits_agg$field_id)) {

    ## load basic phenotype (residuals after regressing out covariates)
    phen_df <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/",
                            phen, "-", samples, "-resid-covars.tab"), data.table = FALSE)
    phen_df <- phen_df[!is.na(phen_df[, 3]),]

    ## load standard C+T (R^2 = 0.1) LOCO PGS
    load(paste0("../results/04-tf-binding/gwas/", phen, "/pgs-loco.RData"))

    ## load list of SNPs in C+T PGS so we can remove GWS hit SNPs whose p-value is > p1_opt
    pgs_loco_coeff <- data.frame()
    for (chr in 1:22) {

        ## read coefficient data frame to check if empty
        coeff_chr <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_0.1-kb_500/", samples,
                                  "/coeff/p1_opt/chr", chr, "-coeff.tab"),
                           data.table = FALSE)
        coeff_chr <- coeff_chr[coeff_chr[, 5] != 0,]
    
        if (nrow(coeff_chr) > 0) {
            pgs_loco_coeff <- rbind(pgs_loco_coeff, coeff_chr)
        }
    }
    ## read optimal p-value threshold p1_optimal
    p1_df <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_0.1-kb_500/wb_all",
                          "/pgs-all/p1_opt/", phen, "-pgs0-rsq.tab"), data.table = FALSE)
    p1_opt <- p1_df$p1_opt

    ## load list of TF PGS SNPs and weights for this phenotype
    load(paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-snps.RData"))

    ## load original GWAS summary stats
    orig_f <- paste0("../results/01-gwas/plink-output/", phen, "/", samples,
                     ".chr", seq(1, 22), ".", phen, ".res.cov.glm.linear.gz")
    orig <- lapply(orig_f, function(x) fread(x, header = TRUE, data.table = FALSE))
    orig <- do.call("rbind", orig)

    for (i in which(indep_hits_agg$field_id == phen)) {

        lhs_snp <- indep_hits_agg$ID[i]
        lhs_chr <- as.numeric(strsplit(lhs_snp, ":")[[1]][1])
        lhs_a1  <- indep_hits_agg$A1[i]

        motif <- indep_hits_agg$motif[i]

        ## get SNPs and weights for this motif
        motif_probs <- tf_probs_ls[[motif]]

        ## remove hit SNPs that are not present in the optimal PGS
        ## (this will only happen if the optimal p-value threshold is lower than the
        ## GWS level used to define the initial set of hits for which to find tags)
        if (p1_opt < 5e-8) {
            motif_probs <- motif_probs[motif_probs$hit_ID %in% pgs_loco_coeff$ID,]
        }

        ## keep only hit SNPs and assign them the sum of weights of their tags
        hit_ID_order <- motif_probs$hit_ID
        motif_probs <- motif_probs %>%
            group_by(hit_ID) %>%
            summarise(post_prob = sum(post_prob)) %>%
            as.data.frame()
        motif_probs <- motif_probs[order(match(motif_probs$hit_ID, hit_ID_order)),]
        n_hit_snps <- nrow(motif_probs)

        ## merge with original GWAS results to get BETA coeffs (and A1 allele)
        motif_coeff <- left_join(motif_probs,
                                 orig[, c("ID", "A1", "BETA")],
                                 by = c("hit_ID" = "ID"))

        ## compute weighted coefficient
        motif_coeff$BETA_w <- motif_coeff$BETA * motif_coeff$post_prob

        ## compute sub-PGS by chromosome
        motif_coeff$CHR <- as.numeric(sapply(strsplit(motif_coeff$hit_ID, split = ":"), "[[", 1))
        motif_coeff_chrs <- sort(unique(motif_coeff$CHR))
        motif_pgs <- data.frame(IID = wb_ids$IID)
        for (chr in motif_coeff_chrs) {

            ## subset coeff df
            motif_coeff_chr <- motif_coeff[motif_coeff$CHR == chr,]

            ## load genotypes
            chr_geno <- load_geno(motif_coeff_chr$hit_ID, motif_coeff_chr$A1, wb_ids$IID)
            
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

        ## keep only chromosome of LHS hit
        motif_pgs <- motif_pgs[, c("IID", "pgs_full", paste0("pgs_loco", lhs_chr))]

        ## load genotypes
        geno <- load_geno(lhs_snp, lhs_a1, wb_ids$IID)
        ## check whether correlation between genotype and its square is < 0.999
        sq_cor_low <- (cor(geno[, 2], geno[, 2]^2) < 0.999)

        ## prepare data.frame with phenotype, SNP and standard LOCO PGS
        main <- left_join(phen_df[, 2:3], geno, by = "IID")
        colnames(main)[2:3] <- c("phen", "geno")
        if (sq_cor_low) main$geno_sq <- main$geno^2
        main <- left_join(main, pgs_loco_df[, c("IID", paste0("pgs_loco", lhs_chr))], by = "IID")
        colnames(main)[length(colnames(main))] <- "pgs_loco"

        ## add motif PGS
        main_motif <- left_join(main, motif_pgs[, c("IID", paste0("pgs_loco", lhs_chr))], by = "IID")
        colnames(main_motif)[ncol(main_motif)] <- "pgs_tf"
        ## subtract it from LOCO PGS
        main_motif$pgs_loco_diff <- main_motif$pgs_loco - main_motif$pgs_tf

        ## run regression
        if (sq_cor_low) {
            reg <- lm(phen ~ geno + geno_sq + geno*pgs_loco_diff + geno*pgs_tf, data = main_motif)
        } else {
            reg <- lm(phen ~ geno + geno*pgs_loco_diff + geno*pgs_tf, data = main_motif)
        }

        ## add sumstats to results data.frame
        coeffs <- summary(reg)$coeff
        sumstats_vec <- c(n_hit_snps = n_hit_snps,
                          snp_BETA                  = coeffs["geno", 1],
                          snp_LOG10_P               = - (pt(abs(coeffs["geno", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),  # to avoid numerical underflow
                          snp_sq_BETA               = ifelse(sq_cor_low, coeffs["geno_sq", 1], NA),
                          snp_sq_LOG10_P            = ifelse(sq_cor_low,
                                                             - (pt(abs(coeffs["geno_sq", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                                             NA),
                          pgs_loco_diff_BETA        = coeffs["pgs_loco_diff", 1],
                          pgs_loco_diff_LOG10_P     = - (pt(abs(coeffs["pgs_loco_diff", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                          int_pgs_loco_diff_BETA    = coeffs["geno:pgs_loco_diff", 1],
                          int_pgs_loco_diff_LOG10_P = - (pt(abs(coeffs["geno:pgs_loco_diff", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                          pgs_tf_hits_BETA          = coeffs["pgs_tf", 1],
                          pgs_tf_hits_LOG10_P       = - (pt(abs(coeffs["pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                          int_pgs_tf_hits_BETA      = coeffs["geno:pgs_tf", 1],
                          int_pgs_tf_hits_LOG10_P   = - (pt(abs(coeffs["geno:pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)))
        indep_hits_agg_retest[i, 13:25] <- sumstats_vec
    }
}


## add motif cluster information
cluster_list <- fread("../data/tf-binding/hocomoco/v13/cluster_list.tsv", data.table = FALSE)
indep_hits_agg_retest$representative_motif <- NA
indep_hits_agg_retest$primary_family <- NA

for (i in 1:nrow(indep_hits_agg_retest)) {

    motif <- indep_hits_agg_retest$motif[i]

    ## find motif in clusters list
    motif_ind <- which(grepl(motif, cluster_list$Clustered_Motifs))
    if (length(motif_ind) > 1) { stop("There are multiple matches in clusters list.") }

    ## add Representative Motif and Primary Family (from TFClass)
    indep_hits_agg_retest$representative_motif[i] <- cluster_list$Representative_Motif[motif_ind]
    indep_hits_agg_retest$primary_family[i] <- cluster_list$Primary_Family[motif_ind]
}

## reorder columns
indep_hits_agg_retest <- cbind(indep_hits_agg_retest[, 1:12],
                               indep_hits_agg_retest[, 26:27],
                               indep_hits_agg_retest[, 13:25])


## export
save(indep_hits_agg_retest,
     file = "../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate-retest.RData")

## round numbers before exporting as text file
cols <- c("MAF",
          "snp_BETA",               "snp_LOG10_P",
          "snp_sq_BETA",            "snp_sq_LOG10_P",
          "pgs_loco_diff_BETA",     "pgs_loco_diff_LOG10_P",
          "int_pgs_loco_diff_BETA", "int_pgs_loco_diff_LOG10_P",
          "pgs_tf_hits_BETA",       "pgs_tf_hits_LOG10_P",
          "int_pgs_tf_hits_BETA",   "int_pgs_tf_hits_LOG10_P")
indep_hits_agg_retest_out <- as.data.table(indep_hits_agg_retest)
indep_hits_agg_retest_out[, (cols) := round(.SD, 8), .SDcols = cols]
fwrite(indep_hits_agg_retest_out,
       file = "../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate-retest.tab",
       sep = "\t", na = "NA", quote = FALSE)
