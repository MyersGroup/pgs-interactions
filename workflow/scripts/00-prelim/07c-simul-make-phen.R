## sum PGS for different chr and add transf. + nois

library(data.table)
library(hash)

set.seed(283619)


n_causal_snps_str  <- c("100", "1k", "10k")
architecture       <- c("rg", "cl")
architecture_hash  <- hash(architecture, c("regular", "clustered"))

transf_type        <- c("nt", "sg")  # no transformation / sigmoid
transf_timing      <- c("bn", "an")  # before/after noise
transf_timing_hash <- hash(transf_timing, c("before add. noise", "after add. noise"))
heritability       <- c(0.3, 0.6)


for (ncs in n_causal_snps_str) {
    for (arch in architecture) {

        ## load chr1 score
        pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                               ncs, "_", arch, "-chr1.sscore"),
                        data.table = FALSE)

        ## rescale components by allele count
        pgs_df[, 5] <- pgs_df[, 5] * pgs_df$ALLELE_CT

        ## add remaining components
        for (chr in 2:22) {
            pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                                       ncs, "_", arch, "-chr", chr, ".sscore"),
                               data.table = FALSE)
            pgs_df_chr[, 5] <- pgs_df_chr[, 5] * pgs_df_chr$ALLELE_CT
            pgs_df[, 5] <- pgs_df[, 5] + pgs_df_chr[, 5]
        }
        pgs_df <- pgs_df[, c(1, 2, 5)]
        

        ## add transformation + noise
        for (h2 in heritability) {
            for (tr in transf_type) {
                if (tr == "sg") {
                    for (tr_t in transf_timing) {

                        phen_df <- pgs_df

                        if (tr_t == "bn") {

                            ## transform
                            phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3]))

                            ## add environmental noise
                            colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                            phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                                rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                        } else {

                            ## add environmental noise
                            colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                            phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                                rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                            ## transform
                            phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3]))

                        }

                        ## export
                        phen_df_out <- phen_df[, 2:3]
                        colnames(phen_df_out)[1] <- "#IID"
                        phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                        fwrite(phen_df_out,
                               file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2, ".tab"),
                               sep = '\t', na = 'NA', quote = FALSE)
                    }        
                } else {
                    
                    phen_df <- pgs_df

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, "_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                    ## export
                    phen_df_out <- phen_df[, 2:3]
                    colnames(phen_df_out)[1] <- "#IID"
                    phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                    fwrite(phen_df_out,
                           file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, "_", h2, ".tab"),
                           sep = '\t', na = 'NA', quote = FALSE)
                }
            }
        }
    }
}




## scaled sigmoids
set.seed(948293)

for (ncs in n_causal_snps_str) {
    for (arch in architecture) {

        ## load chr1 score
        pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                               ncs, "_", arch, "-chr1.sscore"),
                        data.table = FALSE)

        ## rescale components by allele count
        pgs_df[, 5] <- pgs_df[, 5] * pgs_df$ALLELE_CT

        ## add remaining components
        for (chr in 2:22) {
            pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                                       ncs, "_", arch, "-chr", chr, ".sscore"),
                                data.table = FALSE)
            pgs_df_chr[, 5] <- pgs_df_chr[, 5] * pgs_df_chr$ALLELE_CT
            pgs_df[, 5] <- pgs_df[, 5] + pgs_df_chr[, 5]
        }
        pgs_df <- pgs_df[, c(1, 2, 5)]
        

        ## add transformation + noise
        for (h2 in heritability) {
            tr <- "sgsc"

            for (tr_t in transf_timing) {

                phen_df <- pgs_df

                if (tr_t == "bn") {

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                } else {

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))
                    noise <- phen_df[, 3]

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                }

                ## export
                phen_df_out <- phen_df[, 2:3]
                colnames(phen_df_out)[1] <- "#IID"
                phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                fwrite(phen_df_out,
                       file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, tr_t, "_", h2, ".tab"),
                       sep = '\t', na = 'NA', quote = FALSE)
            }        
        }
    }
}




## make versions with scaled betas
set.seed(583920)

for (ncs in n_causal_snps_str) {
    for (arch in architecture) {

        ## load chr1 score
        pgs_df <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                               ncs, "_", arch, "-chr1.sscore"),
                        data.table = FALSE)

        ## rescale components by allele count
        pgs_df[, 6] <- pgs_df[, 6] * pgs_df$ALLELE_CT

        ## add remaining components
        for (chr in 2:22) {
            pgs_df_chr <- fread(paste0("../data/phenotypes/simulations-prep/pgs/sim-pgs-",
                                       ncs, "_", arch, "-chr", chr, ".sscore"),
                               data.table = FALSE)
            pgs_df_chr[, 6] <- pgs_df_chr[, 6] * pgs_df_chr$ALLELE_CT
            pgs_df[, 6] <- pgs_df[, 6] + pgs_df_chr[, 6]
        }
        pgs_df <- pgs_df[, c(1, 2, 6)]
        

        ## add transformation + noise
        for (h2 in heritability) {

            ## no transformation
            tr <- "nt"
                   
            phen_df <- pgs_df

            ## add environmental noise
            colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, "_a5_", h2)
            phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

            ## export
            phen_df_out <- phen_df[, 2:3]
            colnames(phen_df_out)[1] <- "#IID"
            phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

            fwrite(phen_df_out,
                   file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, "_a5_", h2, ".tab"),
                   sep = '\t', na = 'NA', quote = FALSE)


            ## scaled sigmoids
            tr <- "sgsc"

            for (tr_t in transf_timing) {

                phen_df <- pgs_df

                if (tr_t == "bn") {

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                    ## add environmental noise
                    transform <- as.numeric(scale(phen_df[, 3]))  # store transformation before noise but after rescaling for plotting
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_a5_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))

                } else {

                    ## add environmental noise
                    colnames(phen_df)[3] <- paste0("f.sim_", ncs, "_", arch, "_", tr, tr_t, "_a5_", h2)
                    phen_df[, 3] <- as.numeric(scale(phen_df[, 3])) +  # standardise prior to adding noise
                        rnorm(nrow(phen_df), mean = 0, sd = sqrt((1 - h2) / h2))
                    noise <- phen_df[, 3]

                    ## transform
                    stdv <- sd(phen_df[, 3])
                    phen_df[, 3] <- 1 / (1 + exp(-phen_df[, 3] / stdv))

                }

                ## export
                phen_df_out <- phen_df[, 2:3]
                colnames(phen_df_out)[1] <- "#IID"
                phen_df_out[, 2] <- round(phen_df_out[, 2], digits = 8)

                fwrite(phen_df_out,
                       file = paste0("../data/phenotypes/clean/f.sim_", ncs, "_", arch, "_", tr, tr_t, "_a5_", h2, ".tab"),
                       sep = '\t', na = 'NA', quote = FALSE)
            }        
        }
    }
}
