## make simulated data with interactions

library(data.table)
setDTthreads(8)
library(dplyr)
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


## load list of simulated interactions
inter <- fread("../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab", data.table = FALSE)


## load minor alleles
a1_f <- paste0("../data/imputed-wb_all-stats/minor-alleles/chr", seq(1, 22), ".tab")
a1 <- lapply(a1_f, function(x) fread(x, data.table = FALSE))
a1 <- do.call("rbind", a1)




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
sample_f <- fread("../data/sample-ids/ukb-103076-imp-auto-s486989.psam", data.table = FALSE)

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

set.seed(389331)

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
