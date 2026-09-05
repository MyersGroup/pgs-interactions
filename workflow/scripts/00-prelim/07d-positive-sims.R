## make simulated data with interactions

library(data.table)
library(dplyr)
library(ggplot2)
library(pgenlibr)




## load coefficients from null simulation with 1000 SNPs and no clustering/transformation
null_sim_f <- paste0("../data/phenotypes/simulations-prep/coeff/sim-coeff-1k_rg-chr", seq(1, 22), ".tab")
null_sim   <- lapply(null_sim_f, function(x) fread(x, data.table = FALSE))
null_sim   <- do.call("rbind", null_sim)

## add CHR and POS variables
null_sim$CHR <- as.numeric(
    sapply(strsplit(null_sim$ID, split = ":"),
           "[[", 1)
)
null_sim$POS <- as.numeric(
    sapply(strsplit(null_sim$ID, split = ":"),
           "[[", 2)
)




## load list of all SNPs and compute distance to closest SNP with direct effect

## function to compute distance to closest SNP in same chr
closest_snp <- function(x) {
    pos <- x[2]

    dist <- min(
        abs(pos - null_sim_chr$POS)
    )

    return(dist)
}

all_vars <- list()
a1 <- list()
for (chr in 1:22) {

    ## load data
    all_vars[[chr]] <- fread(paste0("../data/variant-ids/chr", chr, "-var-ids.tab"),
                             header = FALSE,
                             col.names = "ID",
                             data.table = FALSE)

    ## add CHR and POS variables
    all_vars[[chr]]$CHR <- as.numeric(
        sapply(strsplit(all_vars[[chr]]$ID, split = ":"),
               "[[", 1)
    )
    all_vars[[chr]]$POS <- as.numeric(
        sapply(strsplit(all_vars[[chr]]$ID, split = ":"),
               "[[", 2)
    )

    ## add A1 variable
    a1[[chr]] <- fread(paste0("../data/imputed-wb_all-stats/minor-alleles/chr", chr, ".tab"),
                       data.table = FALSE)
    all_vars[[chr]] <- left_join(all_vars[[chr]], a1[[chr]][, c("ID", "A1", "MAF")], by = "ID")

    ## compute distance to closest SNP with direct effect
    null_sim_chr <- null_sim[null_sim$CHR == chr,]
    dist_vec <- apply(all_vars[[chr]][, 2:3], 1, closest_snp)
    all_vars[[chr]]$dist_close <- dist_vec

}
all_vars <- do.call("rbind", all_vars)
a1 <- do.call("rbind", a1)




## choose left side of (SNP, SNP) interaction pairs

## tag SNPs with MAF >= 5%
null_sim <- left_join(null_sim, a1[, c("ID", "MAF")], by = "ID")
null_sim$maf5 <- (null_sim$MAF >= 0.05)
null_top <- null_sim[null_sim$maf5,]
## remove indels (removes 48 SNPs)
null_top$REF <- sapply(strsplit(null_top$ID, split = ":"), "[[", 3)
null_top$ALT <- sapply(strsplit(null_top$ID, split = ":"), "[[", 4)
null_top <- null_top[(nchar(null_top$REF) == 1 &
                      nchar(null_top$ALT) == 1),]

set.seed(92783)

## randomly choose 8 SNPs with direct effect
## these must be at least 500kb away from one another
inter_dir <- null_top[sample(1:nrow(null_top), 1), c("ID", "A1", "MAF", "BETA")]
## remaining independent hits for next sampling round
lt_inter_chr <- as.numeric(strsplit(inter_dir$ID, split = ":")[[1]][1])
lt_inter_pos <- as.numeric(strsplit(inter_dir$ID, split = ":")[[1]][2])
null_top_ind <- null_top[(null_top$CHR != lt_inter_chr) |
                         (abs(null_top$POS - lt_inter_pos) > 500000),]

for (i in 1:7) {
    inter_dir <- rbind(inter_dir,
                       null_top_ind[sample(1:nrow(null_top_ind), 1), c("ID", "A1", "MAF", "BETA")])

    ## remaining independent hits for next sampling round
    lt_inter_chr <- as.numeric(strsplit(last(inter_dir$ID), split = ":")[[1]][1])
    lt_inter_pos <- as.numeric(strsplit(last(inter_dir$ID), split = ":")[[1]][2])
    null_top_ind <- null_top_ind[(null_top_ind$CHR != lt_inter_chr) |
                                 (abs(null_top_ind$POS - lt_inter_pos) > 500000),]
}

## randomly choose 8 SNPs that do not have a direct effect on the trait
## these must be at least 500kb away from any SNP with a direct effect and from one another
## and have MAF >= 5%
all_vars$maf5 <- (all_vars$MAF >= 0.05)
all_vars_no <- all_vars[(all_vars$dist_close > 500000) &
                       all_vars$maf5,]

inter_no <- all_vars_no[sample(1:nrow(all_vars_no), 1), c("ID", "A1", "MAF")]
lt_inter_chr <- as.numeric(strsplit(inter_no$ID, split = ":")[[1]][1])
lt_inter_pos <- as.numeric(strsplit(inter_no$ID, split = ":")[[1]][2])
all_vars_no_ind <- all_vars_no[(all_vars_no$CHR != lt_inter_chr) |
                               (abs(all_vars_no$POS - lt_inter_pos) > 500000),]

for (i in 1:7) {
    inter_no <- rbind(inter_no,
                      all_vars_no_ind[sample(1:nrow(all_vars_no_ind), 1), c("ID", "A1", "MAF")])

    lt_inter_chr <- as.numeric(strsplit(last(inter_no$ID), split = ":")[[1]][1])
    lt_inter_pos <- as.numeric(strsplit(last(inter_no$ID), split = ":")[[1]][2])
    all_vars_no_ind <- all_vars_no_ind[(all_vars_no_ind$CHR != lt_inter_chr) |
                                       (abs(all_vars_no_ind$POS - lt_inter_pos) > 500000),]
}
inter_no$BETA <- 0

## randomly choose 6 SNPs that do not have a direct effect on the trait but which will have a biased interaction
## these must be at least 500kb away from any of the previous SNPs and from one another
inter_bia <- data.frame()
for (i in 1:6) {
    inter_bia <- rbind(inter_bia,
                       all_vars_no_ind[sample(1:nrow(all_vars_no_ind), 1), c("ID", "A1", "MAF")])

    lt_inter_chr <- as.numeric(strsplit(last(inter_bia$ID), split = ":")[[1]][1])
    lt_inter_pos <- as.numeric(strsplit(last(inter_bia$ID), split = ":")[[1]][2])
    all_vars_no_ind <- all_vars_no_ind[(all_vars_no_ind$CHR != lt_inter_chr) |
                                       (abs(all_vars_no_ind$POS - lt_inter_pos) > 500000),]
}
inter_bia$BETA <- 0


## choose whether interactions increase or decrease original effect
inter_dir$BETA.int.sign <- sample(c(-1, 1), nrow(inter_dir), replace = TRUE)
inter_no$BETA.int.sign  <- sample(c(-1, 1), nrow(inter_no),  replace = TRUE)
inter_bia$BETA.int.sign <- sample(c(-1, 1), nrow(inter_bia), replace = TRUE)




## choose right side of (SNP, SNP) interaction pairs

inter <- data.frame(ID.1            = as.character(),
                    A1.1            = as.character(),
                    MAF.1           = as.numeric(),
                    BETA.1          = as.numeric(),
                    BETA.1.int.sign = as.numeric(),
                    type            = as.character(),
                    affects         = as.character(),
                    ID.2            = as.character(),
                    A1.2            = as.character(),
                    MAF.2           = as.numeric(),
                    BETA.2          = as.numeric())

colnames(inter_dir) <- c("ID.1", "A1.1", "MAF.1", "BETA.1", "BETA.1.int.sign")
inter_dir$type <- "Direct effect"
inter_dir$affects <- NA
inter_dir$ID.2 <- NA
inter_dir$A1.2 <- NA
inter_dir$MAF.2 <- NA
inter_dir$BETA.2 <- NA

colnames(inter_no) <- c("ID.1", "A1.1", "MAF.1", "BETA.1", "BETA.1.int.sign")
inter_no$type <- "No direct effect"
inter_no$affects <- NA
inter_no$ID.2 <- NA
inter_no$A1.2 <- NA
inter_no$MAF.2 <- NA
inter_no$BETA.2 <- NA

colnames(inter_bia) <- c("ID.1", "A1.1", "MAF.1", "BETA.1", "BETA.1.int.sign")
inter_bia$type <- "Biased inter. sign"
inter_bia$affects <- NA
inter_bia$ID.2 <- NA
inter_bia$A1.2 <- NA
inter_bia$MAF.2 <- NA
inter_bia$BETA.2 <- NA


## 1 SNP affected
## choose 2 SNPs with direct effect
inter_add <- inter_dir[sample(1:nrow(inter_dir), 2),]
inter_dir <- inter_dir[!(inter_dir$ID.1 %in% inter_add$ID.1),]
## choose 2 SNPs without direct effect
inter_add <- rbind(inter_add,
                   inter_no[sample(1:nrow(inter_no), 2),])
inter_no  <- inter_no[!(inter_no$ID.1 %in% inter_add$ID.1),]
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "1 SNP"
## convert each row to separate data.frame in list
## note that split orders list elements alphabetically
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNP is affected
for (i in 1:length(inter_add)) {
    null_sim_exc <- null_sim[!(null_sim$ID %in% inter_add[[i]]$ID.1),]
    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- null_sim_exc[sample(1:nrow(null_sim_exc), 1),
                                                                           c("ID", "A1", "MAF", "BETA")]
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))

## choose 2 SNPs without direct effect but with biased interaction sign
inter_add <- inter_bia[sample(1:nrow(inter_bia), 2),]
inter_bia  <- inter_bia[!(inter_bia$ID.1 %in% inter_add$ID.1),]
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "1 SNP"
## convert each row to separate data.frame in list
## note that split orders list elements alphabetically
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNP is affected
for (i in 1:length(inter_add)) {
    null_sim_pos <- null_sim[null_sim$BETA > 0,]  # chosen interacting SNP must have positive additive effect size
    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- null_sim_pos[sample(1:nrow(null_sim_pos), 1),
                                                                           c("ID", "A1", "MAF", "BETA")]
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))




## 2% of SNPs affected
## choose 2 SNPs with direct effect
inter_add <- inter_dir[sample(1:nrow(inter_dir), 2),]
inter_dir <- inter_dir[!(inter_dir$ID.1 %in% inter_add$ID.1),]
## choose 2 SNPs without direct effect
inter_add <- rbind(inter_add,
                   inter_no[sample(1:nrow(inter_no), 2),])
inter_no  <- inter_no[!(inter_no$ID.1 %in% inter_add$ID.1),]
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "2% of SNPs"
## convert each row to separate data.frame in list
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNPs are affected
for (i in 1:length(inter_add)) {
    ## repeat row 2% * 1000 times
    inter_add[[i]] <- do.call("rbind",
                              replicate(0.02 * 1000, inter_add[[i]], simplify = FALSE))
    ## draw interacting SNPs
    null_sim_exc <- null_sim[!(null_sim$ID %in% inter_add[[i]]$ID.1),]
    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- null_sim_exc[sample(1:nrow(null_sim_exc), 0.02 * 1000),
                                                                           c("ID", "A1", "MAF", "BETA")]
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))

## choose 2 SNPs without direct effect but with biased interaction sign
inter_add <- inter_bia[sample(1:nrow(inter_bia), 2),]
inter_bia  <- inter_bia[!(inter_bia$ID.1 %in% inter_add$ID.1),]
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "2% of SNPs"
## convert each row to separate data.frame in list
## note that split orders list elements alphabetically
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNPs are affected
for (i in 1:length(inter_add)) {
    ## repeat row 2% * 1000 times
    inter_add[[i]] <- do.call("rbind",
                              replicate(0.02 * 1000, inter_add[[i]], simplify = FALSE))
    ## draw interacting SNPs
    null_sim_pos <- null_sim[null_sim$BETA > 0,]
    null_sim_neg <- null_sim[null_sim$BETA < 0,]

    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- rbind(null_sim_pos[sample(1:nrow(null_sim_pos),
                                                                                        0.8 * 0.02 * 1000),
                                                                                 c("ID", "A1", "MAF", "BETA")],
                                                                    null_sim_neg[sample(1:nrow(null_sim_neg),
                                                                                        0.2 * 0.02 * 1000),
                                                                                 c("ID", "A1", "MAF", "BETA")])
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))




## 10% of SNPs affected
## choose 2 SNPs with direct effect
inter_add <- inter_dir[sample(1:nrow(inter_dir), 2),]
inter_dir <- inter_dir[!(inter_dir$ID.1 %in% inter_add$ID.1),]
## choose 2 SNPs without direct effect
inter_add <- rbind(inter_add,
                   inter_no[sample(1:nrow(inter_no), 2),])
inter_no  <- inter_no[!(inter_no$ID.1 %in% inter_add$ID.1),]
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "10% of SNPs"
## convert each row to separate data.frame in list
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNPs are affected
for (i in 1:length(inter_add)) {
    ## repeat row 2% * 1000 times
    inter_add[[i]] <- do.call("rbind",
                              replicate(0.1 * 1000, inter_add[[i]], simplify = FALSE))
    ## draw interacting SNPs
    null_sim_exc <- null_sim[!(null_sim$ID %in% inter_add[[i]]$ID.1),]
    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- null_sim_exc[sample(1:nrow(null_sim_exc), 0.1 * 1000),
                                                                           c("ID", "A1", "MAF", "BETA")]
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))

## choose 2 SNPs without direct effect but with biased interaction sign
inter_add <- inter_bia[sample(1:nrow(inter_bia), 2),]
inter_bia  <- inter_bia[!(inter_bia$ID.1 %in% inter_add$ID.1),]
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "10% of SNPs"
## convert each row to separate data.frame in list
## note that split orders list elements alphabetically
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNPs are affected
for (i in 1:length(inter_add)) {
    ## repeat row 2% * 1000 times
    inter_add[[i]] <- do.call("rbind",
                              replicate(0.1 * 1000, inter_add[[i]], simplify = FALSE))
    ## draw interacting SNPs
    null_sim_pos <- null_sim[null_sim$BETA > 0,]
    null_sim_neg <- null_sim[null_sim$BETA < 0,]

    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- rbind(null_sim_pos[sample(1:nrow(null_sim_pos),
                                                                                        0.8 * 0.1 * 1000),
                                                                                 c("ID", "A1", "MAF", "BETA")],
                                                                    null_sim_neg[sample(1:nrow(null_sim_neg),
                                                                                        0.2 * 0.1 * 1000),
                                                                                 c("ID", "A1", "MAF", "BETA")])
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))




## 100% of SNPs affected
## take remaining 2 SNPs with direct effect
inter_add <- inter_dir
## take remaining 2 SNPs without direct effect
inter_add <- rbind(inter_add,
                   inter_no)
## add info on how many SNPs each LHS SNP affects
inter_add$affects <- "100% of SNPs"
## convert each row to separate data.frame in list
inter_add <- split(inter_add, inter_add$ID.1)
## choose which SNPs is affected
for (i in 1:length(inter_add)) {
    ## repeat row 100% * 1000 times
    inter_add[[i]] <- do.call("rbind",
                              replicate(1000, inter_add[[i]], simplify = FALSE))
    ## 'draw' interacting SNPs (interaction with all SNPs with direct effect)
    inter_add[[i]][, c("ID.2", "A1.2", "MAF.2", "BETA.2")] <- null_sim[, c("ID", "A1", "MAF", "BETA")]
}
## add to main table
inter <- rbind(inter,
               do.call("rbind", inter_add))




## add CHR and POS to table of interactions
inter$CHR.1 <- as.numeric(sapply(strsplit(inter$ID.1, split = ":"), "[[", 1))
inter$CHR.2 <- as.numeric(sapply(strsplit(inter$ID.2, split = ":"), "[[", 1))
inter$POS.1 <- as.numeric(sapply(strsplit(inter$ID.1, split = ":"), "[[", 2))
inter$POS.2 <- as.numeric(sapply(strsplit(inter$ID.2, split = ":"), "[[", 2))

## rearrange cols and order rows by CHR + POS
inter <- inter[, c("CHR.1", "POS.1", "ID.1", "A1.1", "MAF.1", "BETA.1", "BETA.1.int.sign",
                   "CHR.2", "POS.2", "ID.2", "A1.2", "MAF.2", "BETA.2", "type", "affects")]
inter <- inter %>%
    arrange(CHR.1, POS.1, CHR.2, POS.2)




## compute interaction coefficients
## the effect size of the interaction term will be x% * (direct effect size of the right-hand-side SNP)
inter <- inter %>%
    mutate(BETA.int.001 = BETA.2 * BETA.1.int.sign * 0.01,
           BETA.int.010 = BETA.2 * BETA.1.int.sign * 0.1,
           BETA.int.050 = BETA.2 * BETA.1.int.sign * 0.5,
           BETA.int.100 = BETA.2 * BETA.1.int.sign * 1)

## reduce effect size by 60% for SNPs that interact with 100% of SNPs in PGS /100% effect size only/
inter$BETA.int.100[inter$affects == "100% of SNPs"] <- 0.4 * inter$BETA.int.100[inter$affects == "100% of SNPs"]

## make version of coefficients that's scaled by 1/(2p(1-p))^0.25
## this will give a var. of beta^2*(2p(1-p))^0.5 for ech SNP
inter <- inter %>%
    mutate(BETA.int.001.a5 = BETA.int.001 / ((2 * MAF.2 * (1 - MAF.2))^0.25),
           BETA.int.010.a5 = BETA.int.010 / ((2 * MAF.2 * (1 - MAF.2))^0.25),
           BETA.int.050.a5 = BETA.int.050 / ((2 * MAF.2 * (1 - MAF.2))^0.25),
           BETA.int.100.a5 = BETA.int.100 / ((2 * MAF.2 * (1 - MAF.2))^0.25))

## export
fwrite(inter,
       file = "../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab",
       sep = "\t", na = NA, quote = FALSE)




## load genotypes
## make list of all SNPs to load
snps <- data.frame(ID = unique(c(null_sim$ID,
                                 inter$ID.1,
                                 inter$ID.2)),
                   stringsAsFactors = FALSE)

snps$CHR <- as.numeric(sapply(strsplit(snps$ID, split = ":"), "[[", 1))
snps$POS <- as.numeric(sapply(strsplit(snps$ID, split = ":"), "[[", 2))
snps$REF <- as.character(sapply(strsplit(snps$ID, split = ":"), "[[", 3))
snps$ALT <- as.character(sapply(strsplit(snps$ID, split = ":"), "[[", 4))
snps <- left_join(snps, a1[, c("ID", "A1")], by = "ID")

snps <- snps %>%
    select(CHR, POS, ID, REF, ALT, A1) %>%
    arrange(CHR, POS)

## load sample file
sample_f <- fread("../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam", data.table = FALSE)

## load list of all samples considered in this analysis
all_ids <- fread("../data/sample-ids/filtered/all-ids.tab", data.table = FALSE)

geno <- list()
for (chr in unique(snps$CHR)) {

    f.pgen <- paste0("../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen")
    f.pvar <- paste0("../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar")

    pvar <- pgenlibr::NewPvar(f.pvar)
    pgen <- pgenlibr::NewPgen(f.pgen, pvar=pvar)

    ## load dosages from PGEN
    var_ids <- snps$ID[snps$CHR == chr]
    var_num <- rep(NA_real_, length(var_ids))
    for (i in 1:length(var_ids)) {
        var_num[i] <- pgenlibr::GetVariantsById(pvar = pvar, id = var_ids[i])
    }

    geno[[chr]] <- as.matrix(ReadList(pgen, var_num))

    ## add column names
    colnames(geno[[chr]]) <- paste0("g.", var_ids)

    ## keep only samples from all_ids file
    geno[[chr]] <- geno[[chr]][(sample_f$IID %in% all_ids$IID),]

    ## flip allele where needed to match A1
    ref_snps <- (snps$A1[snps$CHR == chr] == snps$REF[snps$CHR == chr])
    ref_snps_cols <- seq(1, ncol(geno[[chr]]))[ref_snps]
    geno[[chr]][, ref_snps_cols] <- 2 - geno[[chr]][, ref_snps_cols]

}
geno <- do.call("cbind", geno)




## mean-centre all SNPs
wb_samples <- fread("../data/sample-ids/filtered/wb_all-ids.tab", data.table = FALSE)

for (i in 1:ncol(geno)) {

    ## compute mean of SNP in White British samples
    snp.geno.wb.mean <- mean(geno[all_ids$IID %in% wb_samples$IID, i])

    ## mean-centre SNP
    geno[, i]  <- geno[, i] - snp.geno.wb.mean

}




## compute basic additive score (without interactions)

if (!all.equal(paste0("g.", null_sim$ID), colnames(geno)[snps$ID %in% null_sim$ID])) {
    stop("Missing SNPs or different order between genotype matrix and `snps` data.frame.")
}

sim_pgs <- all_ids
add_cols <- seq(1, ncol(geno))[snps$ID %in% null_sim$ID]
sim_pgs$phen_add <- geno[, add_cols] %*% null_sim$BETA
sim_pgs$phen_add_a5 <- geno[, add_cols] %*% null_sim$BETA.a5  # scaled

## ## load PGS to double check data
## ## (note that this works only before SNPs are mean-centred above)
## ## load chr1 score
## pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
##                        "1k", "_", "rg", "-chr1.sscore"),
##                 data.table = FALSE)
## ## rescale components by allele count
## pgs_df[, 5] <- pgs_df[, 5] * pgs_df$ALLELE_CT
## ## scaled
## pgs_df[, 6] <- pgs_df[, 6] * pgs_df$ALLELE_CT
## ## add remaining components
## for (chr in 2:22) {
##   pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
##                              "1k", "_", "rg", "-chr", chr, ".sscore"),
##                       data.table = FALSE)
##   pgs_df_chr[, 5] <- pgs_df_chr[, 5] * pgs_df_chr$ALLELE_CT
##   pgs_df[, 5] <- pgs_df[, 5] + pgs_df_chr[, 5]
##   ## scaled
##   pgs_df_chr[, 6] <- pgs_df_chr[, 6] * pgs_df_chr$ALLELE_CT
##   pgs_df[, 6] <- pgs_df[, 6] + pgs_df_chr[, 6]
## }
## ## > cor(pgs_df$SCORE1_AVG, sim_pgs$phen_add)
## ##      [,1]
## ## [1,]    1
## ## > cor(pgs_df$SCORE2_AVG, sim_pgs$phen_add_a5)
## ##      [,1]
## ## [1,]    1




## compute pairwise variables corresponding to interaction terms SNP1*SNP2

geno_inter <- matrix(data = NA_real_,
                     nrow = nrow(geno),
                     ncol = nrow(inter))

for (i in 1:nrow(inter)) {

    snp1 <- inter$ID.1[i]
    snp2 <- inter$ID.2[i]

    ## get index of cols in geno matrix
    snp1.col_ind <- which(colnames(geno) == paste0("g.", snp1))
    snp2.col_ind <- which(colnames(geno) == paste0("g.", snp2))

    ## multiply
    geno_inter[, i] <- geno[, snp1.col_ind] * geno[, snp2.col_ind]

}

colnames(geno_inter) <- paste0("g.", inter$ID.1, "_", inter$ID.2)




## compute interactions
effect_sz <- c("001", "010", "050", "100")
alphas <- c("0", "-0.5")
inter_set <- c("all", "no100")

for (e in effect_sz) {

    for (i in inter_set) {  # make version of the trait without the SNPs that interact with the full PGS

        for (a in alphas) {

            if (i == "all") {

                if (a == "0") {
                    sim_pgs$phen_inter_all <- geno_inter %*% inter[[paste0("BETA.int.", e)]]
                    colnames(sim_pgs)[ncol(sim_pgs)] <- paste0("phen_inter_", e, "_all")
                } else {
                    sim_pgs$phen_inter_all <- geno_inter %*% inter[[paste0("BETA.int.", e, ".a5")]]
                    colnames(sim_pgs)[ncol(sim_pgs)] <- paste0("phen_inter_", e, "_all_a5")
                }

            } else {

                if (a == "0") {
                    no100_ind <- seq(1, nrow(inter))[inter$affects != "100% of SNPs"]
                    sim_pgs$phen_inter_no100 <- geno_inter[, no100_ind] %*% inter[no100_ind, paste0("BETA.int.", e)]
                    colnames(sim_pgs)[ncol(sim_pgs)] <- paste0("phen_inter_", e, "_no100")
                } else {
                    no100_ind <- seq(1, nrow(inter))[inter$affects != "100% of SNPs"]
                    sim_pgs$phen_inter_no100 <- geno_inter[, no100_ind] %*% inter[no100_ind, paste0("BETA.int.", e, ".a5")]
                    colnames(sim_pgs)[ncol(sim_pgs)] <- paste0("phen_inter_", e, "_no100_a5")
                }
            }
        }
    }
}


## sum additive and interaction component
sim_phen <- sim_pgs

for (e in effect_sz) {

    for (i in inter_set) {

        for (a in alphas) {

            if (a == "0") {

                sim_phen[[paste0("phen_all_", e, "_", i)]] <- sim_phen$phen_add + sim_phen[[paste0("phen_inter_", e, "_", i)]]

            } else {

                sim_phen[[paste0("phen_all_", e, "_", i, "_a5")]] <- sim_phen$phen_add_a5 + sim_phen[[paste0("phen_inter_", e, "_", i, "_a5")]]

            }
        }
    }
}


## export additive and interaction components separately and without noise
fwrite(sim_phen,
       file = "../data/phenotypes/simulations-prep/f.sim_int_1k_add.tab",
       sep = '\t', na = NA, quote = FALSE)


## add noise and export
h2 <- 0.6

for (e in effect_sz) {
    for (i in inter_set) {
        for (a in alphas) {

            if (a == "0") {

                cat("# Effect size:", e, "--", i, "interactions -- alpha =", a, "\n")
                cat("Maximum additive heritablity achievable: ")
                cat(cor(sim_phen[["phen_add"]], sim_phen[[paste0("phen_all_", e, "_", i)]])^2, "\n")

                ## add noise to achieve 60% heritability
                var_add <- var(sim_phen[["phen_add"]])
                var_int <- var(sim_phen[[paste0("phen_inter_", e, "_", i)]])
                var_err <- (1 - h2) / h2 * var_add - var_int
                sim_phen[[paste0("phen_all_", e, "_", i)]] <- sim_phen[[paste0("phen_all_", e, "_", i)]] +
                    rnorm(nrow(sim_phen), mean = 0, sd = sqrt(var_err))

                cat("Additive heritability achieved: ")
                cat(cor(sim_phen[["phen_add"]], sim_phen[[paste0("phen_all_", e, "_", i)]])^2, "\n\n")

                ## export
                sim_phen_out <- sim_phen[, c("IID", paste0("phen_all_", e, "_", i))]
                colnames(sim_phen_out) <- c("#IID", paste0("f.sim_int_1k_e", e, "_", i, "_", h2))
                sim_phen_out[, 2] <- round(sim_phen_out[, 2], digits = 8)

                fwrite(sim_phen_out,
                       file = paste0("../data/phenotypes/clean/f.sim_int_1k_e", e, "_", i, "_", h2, ".tab"),
                       sep = '\t', na = NA, quote = FALSE)

            } else {

                cat("# Effect size:", e, "--", i, "interactions -- alpha =", a, "\n")
                cat("Maximum additive heritablity achievable: ")
                cat(cor(sim_phen[["phen_add_a5"]], sim_phen[[paste0("phen_all_", e, "_", i, "_a5")]])^2, "\n")

                ## add noise to achieve 60% heritability
                var_add <- var(sim_phen[["phen_add_a5"]])
                var_int <- var(sim_phen[[paste0("phen_inter_", e, "_", i, "_a5")]])
                var_err <- (1 - h2) / h2 * var_add - var_int
                sim_phen[[paste0("phen_all_", e, "_", i, "_a5")]] <- sim_phen[[paste0("phen_all_", e, "_", i, "_a5")]] +
                    rnorm(nrow(sim_phen), mean = 0, sd = sqrt(var_err))

                cat("Additive heritability achieved: ")
                cat(cor(sim_phen[["phen_add_a5"]], sim_phen[[paste0("phen_all_", e, "_", i, "_a5")]])^2, "\n\n")

                ## export
                sim_phen_out <- sim_phen[, c("IID", paste0("phen_all_", e, "_", i, "_a5"))]
                colnames(sim_phen_out) <- c("#IID", paste0("f.sim_int_1k_e", e, "_", i, "_a5_", h2))
                sim_phen_out[, 2] <- round(sim_phen_out[, 2], digits = 8)

                fwrite(sim_phen_out,
                       file = paste0("../data/phenotypes/clean/f.sim_int_1k_e", e, "_", i, "_a5_", h2, ".tab"),
                       sep = '\t', na = NA, quote = FALSE)

            }
        }
    }
    cat("-----\n\n")
}
