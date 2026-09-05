## Run stepwise regression to identify independent hits

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
p1      <- 1e-5  # significance threshold for stepwise regression

source("scripts/misc/fn-load_geno.R")



## load list of WB sample IIDs
wb_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"), data.table = FALSE)

## load basic phenotype (residuals after regressing out covariates)
phen_df <- fread(paste0("../results/01-gwas/residuals-covars/", phen, "/",
                        phen, "-", samples, "-resid-covars.tab"), data.table = FALSE)
phen_df <- phen_df[!is.na(phen_df[, 3]),]

## load standard C+T (R^2 = 0.1) LOCO PGS
load(paste0("../results/04-tf-binding/gwas/", phen, "/pgs-loco.RData"))

## load list of independent LOCO interaction hits
indep_hits_loco <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/",
                                phen, "-", samples, "-loco-indep-hits.tab"), data.table = FALSE)

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
## process 50 TFs at a time
tf_ind_df <- data.frame(st  = seq(1, 1601, 50),
                        en  = c(seq(50, 1600, 50), 1611))
tf_ind_df$ind <- paste(sprintf("%04d", tf_ind_df$st), sprintf("%04d", tf_ind_df$en), sep = "_")

## load list with weights for each motif (used to count SNPs in each TF PGS)
load(paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-snps.RData"))

## load interaction GWAS results
load(paste0("../results/04-tf-binding/gwas/", phen, "/int-gwas-sumstats.RData"))



## run stepwise regression
indep_hits <- list()
indep_hits_stepwise_sumstats <- list()

for (hit in indep_hits_loco$ID) {

    hit_chr <- indep_hits_loco$CHR[indep_hits_loco$ID == hit]

    ## get GWAS results for this LHS hit
    gwas_hit <- gwas[[hit]]
    
    ## keep only GWS hits
    gwas_hit <- gwas_hit[!is.na(gwas_hit$int_pgs_tf_LOG10_P),]
    gwas_hit <- gwas_hit[gwas_hit$int_pgs_tf_LOG10_P >= -log10(p1),]

    if (nrow(gwas_hit) == 0) {  # if no GWS hits, move to next LHS SNP
        next
    } else if (nrow(gwas_hit) == 1) {  # if exactly one GWS hit, add to list and move to next LHS SNP
        indep_hits[[hit]] <- c(indep_hits[[hit]], gwas_hit$motif[1])
        ## add sumstats from stepwise regression
        indep_hits_stepwise_sumstats[[hit]] <- gwas_hit[1,]
        next
    }

    ## load genotypes
    geno <- load_geno(hit, indep_hits_loco$A1[indep_hits_loco$ID == hit], wb_ids$IID)
    ## check whether correlation between genotype and its square is < 0.999
    sq_cor_low <- (cor(geno[, 2], geno[, 2]^2) < 0.999)

    ## prepare data.frame with phenotype, SNP and standard LOCO PGS
    main <- left_join(phen_df[, 2:3], geno, by = "IID")
    colnames(main)[2:3] <- c("phen", "geno")
    if (sq_cor_low) main$geno_sq <- main$geno^2
    main <- left_join(main, pgs_loco_df[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
    colnames(main)[length(colnames(main))] <- "pgs_loco"
    main$geno_x_pgs_loco <- main$geno * main$pgs_loco

    ## order hits in decreasing order of significance
    gwas_hit <- gwas_hit %>% arrange(-int_pgs_tf_LOG10_P)

    ## add top hit to list
    motif_top <- gwas_hit$motif[1]
    indep_hits[[hit]] <- c(indep_hits[[hit]], motif_top)
    ## add sumstats from stepwise regression
    indep_hits_stepwise_sumstats[[hit]] <- gwas_hit[1,]
    
    ## load PGS for most significant motif
    if (motif_top %in% c("coding", "h3k4me1", "h3k4me3")) {

        load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/", motif_top, "-pgs.RData"))
        motif_top_pgs <- anno_pgs

    } else {

        motif_top_ind <- ceiling(which(namesho == motif_top) / 50)
        load(paste0("../results/04-tf-binding/tf-binding/", phen,
                    "/tf-pgs-", tf_ind_df$ind[motif_top_ind], ".RData"))
        motif_top_pgs <- pgs_ls[[motif_top]]
    }

    ## join to main data.frame
    main_motif <- left_join(main, motif_top_pgs[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
    colnames(main_motif)[ncol(main_motif)] <- paste0("pgs_tf_", motif_top)
    main_motif[, paste0("geno_x_pgs_tf_", motif_top)] <- main_motif$geno * main_motif[, paste0("pgs_tf_", motif_top)]

    ## run stepwise regression    
    for (motif in gwas_hit$motif[-1]) {

        ## load PGS for this motif
        if (motif %in% c("coding", "h3k4me1", "h3k4me3")) {

            load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/", motif, "-pgs.RData"))
            motif_pgs <- anno_pgs

        } else {

            motif_ind <- ceiling(which(namesho == motif) / 50)
            load(paste0("../results/04-tf-binding/tf-binding/", phen,
                        "/tf-pgs-", tf_ind_df$ind[motif_ind], ".RData"))
            motif_pgs <- pgs_ls[[motif]]
        }

        ## add to temporary matrix
        main_motif_add <- left_join(main_motif, motif_pgs[, c("IID", paste0("pgs_loco", hit_chr))], by = "IID")
        colnames(main_motif_add)[ncol(main_motif_add)] <- paste0("pgs_tf_", motif)
        main_motif_add[, paste0("geno_x_pgs_tf_", motif)] <- main_motif_add$geno * main_motif_add[, paste0("pgs_tf_", motif)]
       
        ## regression
        reg <- lm(main_motif_add$phen ~ as.matrix(main_motif_add[, 3:ncol(main_motif_add)]))
        coeffs <- summary(reg)$coeff

        if (coeffs[nrow(coeffs), 4] <= p1) {
            main_motif <- main_motif_add
            ## add to list of independent hits
            indep_hits[[hit]] <- c(indep_hits[[hit]], motif)
            ## add sumstats from stepwise regression
            add_df <- data.frame(motif                = motif,
                                 snp_BETA             = coeffs[2, 1],
                                 snp_LOG10_P          = - (pt(abs(coeffs[2, 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),  # to avoid numerical underflow
                                 snp_sq_BETA          = ifelse(sq_cor_low, coeffs[3, 1], NA),
                                 snp_sq_LOG10_P       = ifelse(sq_cor_low,
                                                               - (pt(abs(coeffs[3, 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                                               NA),
                                 pgs_loco_BETA        = coeffs[3 + ifelse(sq_cor_low, 1, 0), 1],
                                 pgs_loco_LOG10_P     = - (pt(abs(coeffs[3 + ifelse(sq_cor_low, 1, 0), 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                 int_pgs_loco_BETA    = coeffs[4 + ifelse(sq_cor_low, 1, 0), 1],
                                 int_pgs_loco_LOG10_P = - (pt(abs(coeffs[4 + ifelse(sq_cor_low, 1, 0), 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                 pgs_tf_BETA          = coeffs[nrow(coeffs) - 1, 1],
                                 pgs_tf_LOG10_P       = - (pt(abs(coeffs[nrow(coeffs) - 1, 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)),
                                 int_pgs_tf_BETA      = coeffs[nrow(coeffs), 1],
                                 int_pgs_tf_LOG10_P   = - (pt(abs(coeffs[nrow(coeffs), 3]), df = df.residual(reg), lower.tail = FALSE, log.p = TRUE) / log(10) + log10(2)))
            indep_hits_stepwise_sumstats[[hit]] <- rbind(indep_hits_stepwise_sumstats[[hit]], add_df)
        }
    }
}

save(indep_hits, indep_hits_stepwise_sumstats,
     file = paste0("../results/04-tf-binding/gwas/", phen, "/int-gwas-indep-hits.RData"))
