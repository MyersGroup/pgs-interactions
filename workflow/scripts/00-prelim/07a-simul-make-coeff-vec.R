## generate sets of SNPs and coefficients for simulated phenotypes

library(data.table)
library(dplyr)
library(tidyr)
library(hash)

set.seed(9248203)


## load list of variants
f_vars <- paste0("../data/variant-ids/chr", seq(1:22), "-var-ids.tab")
vars   <- lapply(f_vars, function(x) fread(x, col.names = "ID", header = FALSE, data.table = FALSE))
vars   <- do.call("rbind", vars)

## get chromosome and position from var ID strings
vars <- separate(vars,
                 col = ID,
                 into = c("CHR", "POS", NA, NA),
                 sep = ":",
                 remove = FALSE,
                 convert = TRUE)

## get effect allele (A1) from minor allele stats
f_maf <- paste0("../data/imputed-wb_all-stats/minor-alleles/chr", seq(1:22), ".tab")
maf   <- lapply(f_maf, function(x) fread(x, data.table = FALSE))
maf   <- do.call("rbind", maf)
vars  <- left_join(vars, maf[, c("ID", "A1", "MAF")], by = "ID")  # merge preserving order



## simulate coefficients
n_causal_snps <- c(100, 1000, 10000)
n_causal_snps_hash <- hash(n_causal_snps, c("100", "1k", "10k"))
architecture  <- c("rg", "cl")

for (ncs in n_causal_snps) {
    for (arch in architecture) {

        ## sample causal SNPs
        if (arch == "rg") {

            causal <- vars[sample(1:nrow(vars), size = ncs),]
            causal <- causal[order(match(causal$ID, vars$ID)),]

        } else if (arch == "cl") {

            ## 1/2 of SNPs are clustered
            nrg = ncs / 2 * 0.2  # number of regions
            regions <- vars[sample(1:nrow(vars), size = nrg),]
            regions <- regions[order(match(regions$ID, vars$ID)),]
            regions <- apply(regions, 1,
                             function(x) vars[vars$CHR == as.numeric(x["CHR"]) &
                                              abs(vars$POS - as.numeric(x["POS"])) <= 500000,])

            ncs_per_region <- rmultinom(n = 1, size = ncs / 2, prob = rep(1 / nrg, nrg))
            causal <- list()
            for (i in seq_along(ncs_per_region)) {
                if (i == 1) {
                    region_not_causal <- regions[[i]]
                } else {
                    causal_so_far     <- do.call("rbind", causal[seq(1, i-1)])
                    region_not_causal <- regions[[i]][!(regions[[i]]$ID %in% causal_so_far$ID),]
                }

                causal[[i]] <- region_not_causal[sample(1:nrow(region_not_causal), size = ncs_per_region[i]),] 
            }
            causal <- do.call("rbind", causal)

            ## sample 2nd half
            not_in_causal <- vars[!(vars$ID %in% causal$ID),]
            causal <- rbind(causal,
                            not_in_causal[sample(1:nrow(not_in_causal), size = ncs / 2),])
            causal <- causal[order(match(causal$ID, vars$ID)),]

        }

        if (sum(duplicated(causal$ID)) > 0) {
            stop("There are duplicate SNPs.")
        }


        ## sample effect sizes
        causal$BETA <- round(rnorm(nrow(causal)), digits = 6)

        ## make alternative version of effect sizes where these are scaled by 1/(2p(1-p))^0.25
        causal <- causal %>%
            mutate(BETA.a5 = round(BETA / ((2 * MAF * (1 - MAF))^0.25), digits = 6))


        ## export
        for (chr in 1:22) {
            if (chr %in% unique(causal$CHR)) {
                causal_out <- causal[causal$CHR == chr, c("ID", "A1", "BETA", "BETA.a5")]
                fwrite(causal_out,
                       paste0("../data/phenotypes/simulations-prep/coeff/sim-coeff-",
                              values(n_causal_snps_hash, toString(ncs)), "_",
                              arch, "-chr", chr, ".tab"),
                       sep = "\t", na = "NA", quote = FALSE)
            } else {  # make coeff file with a single SNP with coeff 0 if no causal SNP in this chromosome
                causal_out <- vars[vars$CHR == chr, c("ID", "A1")][1,]
                causal_out$BETA <- 0
                causal_out$BETA.a5 <- 0
                fwrite(causal_out,
                       paste0("../data/phenotypes/simulations-prep/coeff/sim-coeff-",
                              values(n_causal_snps_hash, toString(ncs)), "_",
                              arch, "-chr", chr, ".tab"),
                       sep = "\t", na = "NA", quote = FALSE)
            }
        }
    } 
}
