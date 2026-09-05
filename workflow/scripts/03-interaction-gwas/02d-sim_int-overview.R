## Summarise SNP*PGS interaction results for positive simulations

library(data.table)
setDTthreads(4)
library(dplyr)


## load true interaction coefficients
inter_coeff <- fread("../data/phenotypes/simulations-prep/coeff/sim-coeff-inter.tab", data.table = FALSE)

## keep only LHS SNPs
inter_coeff <- inter_coeff %>%
    distinct(ID.1, type, affects, .keep_all = TRUE) %>%
    mutate(type = factor(type, levels = c("Direct effect", "No direct effect", "Biased inter. sign")),
           affects = factor(affects, levels = c("1 SNP", "2% of SNPs", "10% of SNPs", "100% of SNPs"))) %>%
    select("CHR.1", "POS.1", "ID.1", "A1.1", "MAF.1", "BETA.1", "BETA.1.int.sign", "type", "affects") %>%
    arrange(type, affects, CHR.1, POS.1)
colnames(inter_coeff)[5] <- "MAF"



## load SNP*PGS interaction results and compute signed log p-values
effect_sz <- c("001", "010", "050", "100")
alphas <- c("", "_a5")
inter_set <- c("all", "no100")

for (a in alphas) {

    main <- list()

    for (is in inter_set) {
        for (e in effect_sz) {
            
            ## full (QN)
            full_qn <- lapply(paste0("../results/03-interaction-gwas/plink-output/f.sim_int_1k_e", e, "_", is, a, "_0.6_qn/snp-pgs/",
                                     "f.sim_int_1k_e", e, "_", is, a, "_0.6_qn.wb_all.full.chr", seq(1, 22),
                                     ".f.sim_int_1k_e", e, "_", is, a, "_0.6_qn.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz"),
                              function(x) fread(x, data.table = FALSE))
            for (chr in 1:22) {
                full_qn[[chr]] <- full_qn[[chr]][full_qn[[chr]]$TEST == "ADDxpgs_full_mc",]
            }
            full_qn <- do.call(rbind, full_qn)
            ## filter
            full_qn <- full_qn %>%
                filter(ID %in% inter_coeff$ID.1) %>%
                mutate("e{e}_int_full_qn_sign_LOG10_P" := sign(BETA) * LOG10_P) %>%
                select(ID, paste0("e", e, "_int_full_qn_sign_LOG10_P"))

            ## LOCO (QN)
            loco_qn <- lapply(paste0("../results/03-interaction-gwas/plink-output/f.sim_int_1k_e", e, "_", is, a, "_0.6_qn/snp-pgs/",
                                     "f.sim_int_1k_e", e, "_", is, a, "_0.6_qn.wb_all.loco.chr", seq(1, 22),
                                     ".f.sim_int_1k_e", e, "_", is, a, "_0.6_qn.res.cov.pgs_full_mc.ac.ctr.glm.linear.gz"),
                              function(x) fread(x, data.table = FALSE))
            for (chr in 1:22) {
                loco_qn[[chr]] <- loco_qn[[chr]][loco_qn[[chr]]$TEST == paste0("ADDxpgs_loco", chr, "_mc"),]
            }
            loco_qn <- do.call(rbind, loco_qn)
            ## filter
            loco_qn <- loco_qn %>%
                filter(ID %in% inter_coeff$ID.1) %>%
                mutate("e{e}_int_loco_qn_sign_LOG10_P" := sign(BETA) * LOG10_P) %>%
                select(ID, paste0("e", e, "_int_loco_qn_sign_LOG10_P"))


            ## merge
            all <- left_join(full_qn, loco_qn, by = "ID")

            ## add to main data.frame
            if (e == "001") {
                main[[is]] <- all
            } else {
                main[[is]] <- left_join(main[[is]], all, by = "ID")
            }
        }

        ## add variant information and bind into single data.frame
        main[[is]]$inter_set <- is

        ## merge with table of LHS SNPs
        main[[is]] <- left_join(inter_coeff, main[[is]], by = c("ID.1" = "ID"))


        ## add closest independent hits
        ## function to find closest hit
        get_closest <- function(x, df) {
            chr <- as.numeric(strsplit(x, split = ":")[[1]][1])
            pos <- as.numeric(strsplit(x, split = ":")[[1]][2])

            df_filt <- df[df$CHR == chr,]
            if (nrow(df_filt) == 0) {
                return(NA)
            } else {
                min_dist <- min(abs(pos - df_filt$POS))
                return(min_dist)
            }
        }

        for (e in c("001", "010", "050", "100")) {

            ## load independent hits
            ind_full_qn <- data.frame()
            f <- paste0("../results/03-interaction-gwas/indep-hits/f.sim_int_1k_e", e, "_", is, a, "_0.6_qn/snp-pgs/hits/",
                        "f.sim_int_1k_e", e, "_", is, a, "_0.6_qn-wb_all-full-indep-hits.tab")
            if (file.size(f) > 0) {
                ind_full_qn <- fread(f, data.table = FALSE)
            }

            ind_loco_qn <- data.frame()
            f <- paste0("../results/03-interaction-gwas/indep-hits/f.sim_int_1k_e", e, "_", is, a, "_0.6_qn/snp-pgs/hits/",
                        "f.sim_int_1k_e", e, "_", is, a, "_0.6_qn-wb_all-loco-indep-hits.tab")
            if (file.size(f) > 0) {
                ind_loco_qn <- fread(f, data.table = FALSE)
            }
            ind_full_qn <- rbind(ind_loco_qn, ind_full_qn)


            ## add closest hit
            main[[is]][, paste0("e", e, "_int_full_qn_closest_hit")] <- NA
            if (nrow(ind_full_qn) > 0) {
                for (j in 1:nrow(main[[is]])) {
                    main[[is]][j, paste0("e", e, "_int_full_qn_closest_hit")] <- get_closest(main[[is]]$ID.1[j],
                                                                                             ind_full_qn)
                }
            }
            main[[is]][, paste0("e", e, "_int_loco_qn_closest_hit")] <- NA
            if (nrow(ind_loco_qn) > 0) {
                for (j in 1:nrow(main[[is]])) {
                    main[[is]][j, paste0("e", e, "_int_loco_qn_closest_hit")] <- get_closest(main[[is]]$ID.1[j],
                                                                                             ind_loco_qn)
                }
            }
        }
    }

    inter_coeff_res <- do.call(rbind, main)



    ## reorder cols and export
    inter_coeff_res <- inter_coeff_res[, c("inter_set", "ID.1", "A1.1", "MAF", "BETA.1",
                                           "BETA.1.int.sign", "type", "affects",
                                           "e001_int_full_qn_sign_LOG10_P", "e001_int_full_qn_closest_hit",
                                           "e001_int_loco_qn_sign_LOG10_P", "e001_int_loco_qn_closest_hit",
                                           "e010_int_full_qn_sign_LOG10_P", "e010_int_full_qn_closest_hit",
                                           "e010_int_loco_qn_sign_LOG10_P", "e010_int_loco_qn_closest_hit",
                                           "e050_int_full_qn_sign_LOG10_P", "e050_int_full_qn_closest_hit",
                                           "e050_int_loco_qn_sign_LOG10_P", "e050_int_loco_qn_closest_hit",
                                           "e100_int_full_qn_sign_LOG10_P", "e100_int_full_qn_closest_hit",
                                           "e100_int_loco_qn_sign_LOG10_P", "e100_int_loco_qn_closest_hit")]

    fwrite(inter_coeff_res,
           file = paste0("../results/03-interaction-gwas/indep-hits/sim_int-overview/sim_int", a, "-lhs-sumstats.tab"),
           sep = "\t", na = "NA", quote = FALSE)
}
