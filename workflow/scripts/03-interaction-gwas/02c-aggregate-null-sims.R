## Make tables aggregating results of SNP*PGS interaction testing on null sims


library(data.table)
setDTthreads(12)
library(dplyr)


samples <- "wb_all"
p1 <- 5e-8
n_snps <- 12690793


sims_df <- data.frame(field_id     = character(),
                      n_causal     = character(),
                      arch         = character(),
                      transf       = character(),
                      alpha        = numeric(),
                      h2           = numeric(),
                      n_gws_full   = numeric(),
                      fpr_full     = numeric(),
                      n_indep_full = numeric(),
                      n_gws_loco   = numeric(),
                      fpr_loco     = numeric(),
                      n_indep_loco = numeric())


for (a in c("", "_a5")) {
    for (h2 in c(0.3, 0.6)) {

        print(paste(a, h2))

        sims_df_out <- sims_df

        ## these combinations of the remaining variables will be present in the same table
        for (ncs in c("100", "1k", "10k")) {
            for (arch in c("rg", "cl")) {
                for (transf in c("nt", "sgscbn", "sgscan")) {

                    ## get phenotype code
                    phen <- paste0("f.sim_", ncs, "_", arch, "_", transf, a, "_", h2, "_qn")

                    ## load interaction GWAS results
                    gwas_full_f <- paste0("../results/03-interaction-gwas/plink-output/", phen, "/snp-pgs/",
                                          phen, ".", samples, ".full.chr", seq(1, 22), ".", phen, ".res.cov.pgs_full_mc.ac.ctr.glm.linear.gz")
                    gwas_full <- lapply(gwas_full_f, function(x) fread(x, data.table = FALSE))
                    gwas_full <- do.call("rbind", gwas_full)
                    gwas_full <- gwas_full[gwas_full$TEST == "ADDxpgs_full_mc",]

                    gwas_loco_f <- paste0("../results/03-interaction-gwas/plink-output/", phen, "/snp-pgs/",
                                          phen, ".", samples, ".loco.chr", seq(1, 22), ".", phen, ".res.cov.pgs_full_mc.ac.ctr.glm.linear.gz")
                    gwas_loco <- lapply(gwas_loco_f, function(x) fread(x, data.table = FALSE))
                    gwas_loco <- do.call("rbind", gwas_loco)
                    gwas_loco <- gwas_loco[grepl("ADDxpgs_loco", gwas_loco$TEST),]

                    ## count number of GWS hits
                    gwas_full_gws <- sum(gwas_full$LOG10_P >= -log10(p1))
                    gwas_loco_gws <- sum(gwas_loco$LOG10_P >= -log10(p1))

                    if (file.size(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-full-indep-hits.tab")) > 0) {

                        ## load number of independent hits
                        indep_full <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-full-indep-hits-stats.tab"),
                                            data.table = FALSE)
                        gwas_full_indep <- nrow(indep_full)

                    } else {
                        gwas_full_indep <- 0
                    }

                    if (file.size(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-loco-indep-hits.tab")) > 0) {

                        ## load number of independent hits
                        indep_loco <- fread(paste0("../results/03-interaction-gwas/indep-hits/", phen, "/snp-pgs/hits/", phen, "-", samples, "-loco-indep-hits-stats.tab"),
                                            data.table = FALSE)
                        gwas_loco_indep <- nrow(indep_loco)

                    } else {
                        gwas_loco_indep <- 0
                    }

                    ## add to table
                    sims_df_add <- data.frame(field_id     = phen,
                                              n_causal     = ncs,
                                              arch         = arch,
                                              transf       = transf,
                                              alpha        = ifelse(a == "", 0, -0.5),
                                              h2           = h2,
                                              n_gws_full   = gwas_full_gws,
                                              fpr_full     = gwas_full_gws / n_snps,
                                              n_indep_full = gwas_full_indep,
                                              n_gws_loco   = gwas_loco_gws,
                                              fpr_loco     = gwas_loco_gws / n_snps,
                                              n_indep_loco = gwas_loco_indep)
                    sims_df_out <- rbind(sims_df_out, sims_df_add)

                }
            }
        }

        ## add total FPR
        sims_df_total_fpr <- data.frame(field_id     = NA,
                                        n_causal     = NA,
                                        arch         = NA,
                                        transf       = NA,
                                        alpha        = NA,
                                        h2           = NA,
                                        n_gws_full   = NA,
                                        fpr_full     = sum(sims_df_out$n_gws_full) / (nrow(sims_df_out) * n_snps),
                                        n_indep_full = NA,
                                        n_gws_loco   = NA,
                                        fpr_loco     = sum(sims_df_out$n_gws_loco) / (nrow(sims_df_out) * n_snps),
                                        n_indep_loco = NA)
        sims_df_out <- rbind(sims_df_out, sims_df_total_fpr)

        ## export
        fwrite(sims_df_out,
               paste0("../results/03-interaction-gwas/indep-hits/sim_null-overview/sim_null-fp-qn", ifelse(a == "", "", "-a5"), "-", h2, ".tab"),
               sep = "\t", na = "NA", quote = FALSE)
    }
}
