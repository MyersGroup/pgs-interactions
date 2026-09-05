## Run GWAS of interaction with partitioned PGSs

library(argparser)
library(data.table)
setDTthreads(4)
library(dplyr)

options(warn = 2)  # turn warnings into errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples", help = "Sample code",    nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples

source("scripts/misc/fn-load_geno.R")

plink2 <- "/path/to/plink2"



## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## load basic phenotype (residuals after regressing out covariates)
phen_df <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/",
                        phen, "-", samples, "-resid-covars.tab"), data.table = FALSE)
phen_df <- phen_df[!is.na(phen_df[, 3]),]

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
                      "--psam ../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam ",
                      "--extract ../data/variant-ids/chr", chr, "-var-ids.tab ",
                      "--read-freq ../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                      "--keep ../data/sample-ids/filtered/", samples, "-ids.tab ",
                      "--score ../results/02-pgs/initial/", phen, "/r2_0.1-kb_500/", samples, "/coeff/p1_opt/chr", chr, "-coeff.tab 1 2 5 header-read ",
                      "--threads 4 ",
                      "--memory 64000 ",
                      "--out ../results/04-tf-binding/gwas/", phen, "/pgs-chr", chr))
        pgs_chr <- fread(paste0("../results/04-tf-binding/gwas/", phen, "/pgs-chr", chr, ".sscore"), data.table = FALSE)

        ## rescale score by allele count
        pgs_chr[, 5] <- pgs_chr[, 5] * pgs_chr$ALLELE_CT
        colnames(pgs_chr)[5] <- paste0("chr", chr)

        ## add to data.frame
        pgs_loco_df <- left_join(pgs_loco_df, pgs_chr[, c("IID", paste0("chr", chr))], by = "IID")

        ## remove temporary files
        system(paste0("rm ../results/04-tf-binding/gwas/", phen, "/pgs-chr", chr, ".sscore"))
        system(paste0("rm ../results/04-tf-binding/gwas/", phen, "/pgs-chr", chr, ".log"))

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
     file = paste0("../results/04-tf-binding/gwas/", phen, "/pgs-loco.RData"))


## load list of independent LOCO interaction hits
indep_hits <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                           phen, "-", samples, "-loco-indep-hits.tab"), data.table = FALSE)


## run interaction GWAS
gwas <- list()

for (hit in indep_hits$ID) {

    hit_chr <- indep_hits$CHR[indep_hits$ID == hit]

    ## load genotypes
    geno <- load_geno(hit, indep_hits$A1[indep_hits$ID == hit], wb_ids$IID)
    ## check whether correlation between genotype and its square is < 0.999
    sq_cor_low <- (cor(geno[, 2], geno[, 2]^2) < 0.999)

    ## prepare results data.frame
    gwas[[hit]] <- data.frame(matrix(ncol = 13, nrow = 0))
    colnames(gwas[[hit]]) <- c("motif",
                               "snp_BETA",          "snp_LOG10_P",
                               "snp_sq_BETA",       "snp_sq_LOG10_P",
                               "pgs_loco_BETA",     "pgs_loco_LOG10_P",
                               "int_pgs_loco_BETA", "int_pgs_loco_LOG10_P",
                               "pgs_tf_BETA",       "pgs_tf_LOG10_P",
                               "int_pgs_tf_BETA",   "int_pgs_tf_LOG10_P")

    ## prepare data.frame with phenotype, SNP and standard LOCO PGS
    main <- left_join(phen_df[, 2:3], geno, by = "IID")
    colnames(main)[2:3] <- c("phen", "geno")
    if (sq_cor_low) main$geno_sq <- main$geno^2
    main <- left_join(main, pgs_loco_df[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
    colnames(main)[length(colnames(main))] <- "pgs_loco"

    ## coding and H3K4me1/3 PGSs
    for (anno in c("coding", "h3k4me1", "h3k4me3")) {

        ## load PGS
        load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/", anno, "-pgs.RData"))

        ## check that a LOCO PGS is available for this hit's chromosome
        ## (it won't be if all PGS SNPs fall on this chromosome)
        if (sum(!is.na(anno_pgs[paste0("pgs_loco", hit_chr)])) == 0) {
            add_row <- data.frame(anno, matrix(rep(NA, 12), ncol = 12, nrow = 1))
            colnames(add_row) <- colnames(gwas[[hit]])
            gwas[[hit]] <- rbind(gwas[[hit]], add_row)
            next
        }

        ## add annotation PGS to main data.frame
        main_anno <- left_join(main, anno_pgs[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
        colnames(main_anno)[ncol(main_anno)] <- "pgs_tf"

        ## run regression
        if (sq_cor_low) {
            reg <- lm(phen ~ geno + geno_sq + geno*pgs_loco + geno*pgs_tf, data = main_anno)
        } else {
            reg <- lm(phen ~ geno + geno*pgs_loco + geno*pgs_tf, data = main_anno)
        }

        ## add sumstats to results data.frame
        coeffs <- summary(reg)$coeff
        add_df <- data.frame(motif                = anno,
                             snp_BETA             = coeffs["geno", 1],
                             snp_LOG10_P          = - (pt(abs(coeffs["geno", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),  # to avoid numerical underflow
                             snp_sq_BETA          = ifelse(sq_cor_low, coeffs["geno_sq", 1], NA),
                             snp_sq_LOG10_P       = ifelse(sq_cor_low,
                                                           - (pt(abs(coeffs["geno_sq", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                                           NA),
                             pgs_loco_BETA        = coeffs["pgs_loco", 1],
                             pgs_loco_LOG10_P     = - (pt(abs(coeffs["pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                             int_pgs_loco_BETA    = coeffs["geno:pgs_loco", 1],
                             int_pgs_loco_LOG10_P = - (pt(abs(coeffs["geno:pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                             pgs_tf_BETA          = coeffs["pgs_tf", 1],
                             pgs_tf_LOG10_P       = - (pt(abs(coeffs["pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                             int_pgs_tf_BETA      = coeffs["geno:pgs_tf", 1],
                             int_pgs_tf_LOG10_P   = - (pt(abs(coeffs["geno:pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)))
        gwas[[hit]] <- rbind(gwas[[hit]], add_df)
    }


    ## TFs
    tf_ind_df <- data.frame(st  = seq(1, 1601, 50),
                            en  = c(seq(50, 1600, 50), 1611))
    tf_ind_df$ind <- paste(sprintf("%04d", tf_ind_df$st), sprintf("%04d", tf_ind_df$en), sep = "_")
    
    for (i in 1:nrow(tf_ind_df)) {

        ## load PGS batch
        tf_batch <- tf_ind_df$ind[i]
        load(paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-", tf_batch, ".RData"))

        for (motif in names(pgs_ls)) {

            ## get PGS data.frame
            motif_pgs <- pgs_ls[[motif]]

            ## check that a PGS is available for this motif
            if (class(motif_pgs) == "character") {
                if (motif_pgs == "Unavailable") {
                    add_row <- data.frame(motif, matrix(rep(NA, 12), ncol = 12, nrow = 1))
                    colnames(add_row) <- colnames(gwas[[hit]])
                    gwas[[hit]] <- rbind(gwas[[hit]], add_row)
                    next
                }
            }

            ## check that a LOCO PGS is available for this hit's chromosome
            ## (it won't be if all PGS SNPs fall on this chromosome)
            if (sum(!is.na(motif_pgs[paste0("pgs_loco", hit_chr)])) == 0) {
                add_row <- data.frame(motif, matrix(rep(NA, 12), ncol = 12, nrow = 1))
                colnames(add_row) <- colnames(gwas[[hit]])
                gwas[[hit]] <- rbind(gwas[[hit]], add_row)
                next
            }

            ## add annotation PGS to main data.frame
            main_motif <- left_join(main, motif_pgs[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
            colnames(main_motif)[ncol(main_motif)] <- "pgs_tf"

            ## run regression
            if (sq_cor_low) {
                reg <- lm(phen ~ geno + geno_sq + geno*pgs_loco + geno*pgs_tf, data = main_motif)
            } else {
                reg <- lm(phen ~ geno + geno*pgs_loco + geno*pgs_tf, data = main_motif)
            }

            ## add sumstats to results data.frame
            coeffs <- summary(reg)$coeff
            add_df <- data.frame(motif                = motif,
                                 snp_BETA             = coeffs["geno", 1],
                                 snp_LOG10_P          = - (pt(abs(coeffs["geno", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),  # to avoid numerical underflow
                                 snp_sq_BETA          = ifelse(sq_cor_low, coeffs["geno_sq", 1], NA),
                                 snp_sq_LOG10_P       = ifelse(sq_cor_low,
                                                               - (pt(abs(coeffs["geno_sq", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                                               NA),
                                 pgs_loco_BETA        = coeffs["pgs_loco", 1],
                                 pgs_loco_LOG10_P     = - (pt(abs(coeffs["pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                 int_pgs_loco_BETA    = coeffs["geno:pgs_loco", 1],
                                 int_pgs_loco_LOG10_P = - (pt(abs(coeffs["geno:pgs_loco", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                 pgs_tf_BETA          = coeffs["pgs_tf", 1],
                                 pgs_tf_LOG10_P       = - (pt(abs(coeffs["pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                 int_pgs_tf_BETA      = coeffs["geno:pgs_tf", 1],
                                 int_pgs_tf_LOG10_P   = - (pt(abs(coeffs["geno:pgs_tf", 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)))
            gwas[[hit]] <- rbind(gwas[[hit]], add_df)
        }
    }
}

save(gwas,
     file = paste0("../results/04-tf-binding/gwas/", phen, "/int-gwas-sumstats.RData"))
