## Aggregate all independent SNP*TF-PGS hits

library(data.table)
setDTthreads(2)
library(dplyr)
library(tidyr)

options(warn = 2)  # turn warnings into errors



## load list of all LOCO QN SNP*PGS hits
indep_hits_loco <- fread("../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab",
                         data.table = FALSE)

## prepare output data.frame
main <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(main) <- c("field_id",
                    "ID",
                    "motif",
                    "n_hit_snps",
                    "snp_BETA",          "snp_LOG10_P",
                    "snp_sq_BETA",       "snp_sq_LOG10_P",
                    "pgs_loco_BETA",     "pgs_loco_LOG10_P",
                    "int_pgs_loco_BETA", "int_pgs_loco_LOG10_P",
                    "pgs_tf_BETA",       "pgs_tf_LOG10_P",
                    "int_pgs_tf_BETA",   "int_pgs_tf_LOG10_P")


## loop over trait-SNP LHS pairs and concatenate independent hits
for (phen in unique(indep_hits_loco$field_id)) {

    ## load independent hits
    load(paste0("../results/04-tf-binding/gwas/", phen, "/int-gwas-indep-hits.RData"))

    if (length(indep_hits) == 0) { next }

    ## retrieve GWAS sumstats
    load(paste0("../results/04-tf-binding/gwas/", phen, "/int-gwas-sumstats.RData"))

    for (hit in names(indep_hits)) {
        gwas_hit <- gwas[[hit]]
        gwas_hit <- gwas_hit[gwas_hit$motif %in% indep_hits[[hit]],]

        ## add to main data.frame
        main <- rbind(main,
                      data.frame(field_id = phen,
                                 ID = hit,
                                 gwas_hit))
    }
    rm(indep_hits, indep_hits_stepwise_sumstats)
}


## add number of SNPs in each TF PGS
main$n_tag_snps <- NA
for (phen in unique(main$field_id)) {

    for (motif in main$motif[main$field_id == phen]){

        ## load PGS coefficients to count number of SNPs
        if (motif %in% c("coding", "h3k4me1", "h3k4me3")) {
            load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/", motif, "-pgs-snps.RData"))
            if (motif == "coding") {
                motif_probs <- coding_probs
            } else if (motif == "h3k4me1") {
                motif_probs <- h3k4me1_probs
            } else if (motif == "h3k4me3") {
                motif_probs <- h3k4me3_probs
            }
        } else {
            load(paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-snps.RData"))
            motif_probs <- tf_probs_ls[[motif]]
        }
        ## count number of tags in motif PGS
        n_tag_snps <- length(unique(motif_probs$tag_ID))

        ## add to table
        main$n_tag_snps[main$field_id == phen & main$motif == motif] <- n_tag_snps
        rm(n_tag_snps)
    }
}


## add rsIDs
map_posid_rsid_f <- paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", seq(1, 22), ".tab")
map_posid_rsid <- lapply(map_posid_rsid_f, function(x) fread(x, data.table = FALSE))
map_posid_rsid <- do.call("rbind", map_posid_rsid)
main <- left_join(main, map_posid_rsid, by = c("ID" = "POSID"))


## add A1 allele and MAF
maf <- lapply(paste0("../data/imputed-wb_all-stats/minor-alleles/chr", seq(1, 22), ".tab"),
              function(x) fread(x, data.table = FALSE))
maf <- do.call("rbind", maf)
main <- left_join(main, maf[, c("ID", "A1", "MAF")], by = "ID")


## add external annotations
## VEP: (most severe) variant consequence
vep_f <- paste0("../data/annotations/vep/anno/vep-chr", seq(1, 22), ".most_severe.tab.gz")
vep <- lapply(vep_f, function(x) fread(x, skip = "#Uploaded_variation", data.table = FALSE))
vep <- do.call("rbind", vep)
vep <- vep[, c("#Uploaded_variation", "Consequence")]
colnames(vep) <- c("ID", "VEP_Consequence")
main <- left_join(main, vep[, c("ID", "VEP_Consequence")], by = "ID")

## Annovar: affected/closest gene(s) from Gencode and RefSeq
## Gencode
annovar_var_fn_f <- paste0("../data/annotations/annovar/gene-anno/annovar-chr", seq(1, 22), ".ensGene.variant_function.gz")
annovar_var_fn <- lapply(annovar_var_fn_f,
                         function(x) fread(x,
                                           header = FALSE,
                                           sep = "\t",
                                           col.names = c("ensGene_Func", "ensGene_Gene", "Comments"),
                                           data.table = FALSE))
annovar_var_fn <- do.call("rbind", annovar_var_fn)
annovar_var_fn <- annovar_var_fn %>%
    separate(Comments, c("POS", "ID"), "comments: ") %>%
    select(-POS) %>%
    separate(ID, c("ID", "rsID"), ",") %>%
    select(-rsID)
main <- left_join(main, annovar_var_fn, by = "ID")
rm(annovar_var_fn)

## RefSeq
annovar_var_fn_f <- paste0("../data/annotations/annovar/gene-anno/annovar-chr", seq(1, 22), ".refGene.variant_function.gz")
annovar_var_fn <- lapply(annovar_var_fn_f,
                         function(x) fread(x,
                                           header = FALSE,
                                           sep = "\t",
                                           col.names = c("refGene_Func", "refGene_Gene", "Comments"),
                                           data.table = FALSE))
annovar_var_fn <- do.call("rbind", annovar_var_fn)
annovar_var_fn <- annovar_var_fn %>%
    separate(Comments, c("POS", "ID"), "comments: ") %>%
    select(-POS) %>%
    separate(ID, c("ID", "rsID"), ",") %>%
    select(-rsID)
main <- left_join(main, annovar_var_fn, by = "ID")


## add trait descriptions
names_clean_qn <- fread("../data/phenotypes/clean/code-desc-map-real.tab", data.table = FALSE)
colnames(names_clean_qn)[4] <- "desc_clean"
main <- left_join(main, names_clean_qn[, c("field_id", "desc_clean")], by = "field_id")


## add motif cluster information
cluster_list <- fread("../data/tf-binding/hocomoco/v13/cluster_list.tsv", data.table = FALSE)
main$representative_motif <- NA
main$primary_family <- NA

for (i in 1:nrow(main)) {

    motif <- main$motif[i]

    ## clusters are only available for TFs
    if (motif %in% c("coding", "h3k4me1", "h3k4me3")) { next }

    ## find motif in clusters list
    motif_ind <- which(grepl(motif, cluster_list$Clustered_Motifs))
    if (length(motif_ind) > 1) { stop("There are multiple matches in clusters list.") }

    ## add Representative Motif and Primary Family (from TFClass)
    main$representative_motif[i] <- cluster_list$Representative_Motif[motif_ind]
    main$primary_family[i] <- cluster_list$Primary_Family[motif_ind]
}


## export
indep_hits_agg <- main[, c("field_id", "desc_clean",
                           "ID", "rsID", "A1", "MAF",
                           "VEP_Consequence", "ensGene_Func", "ensGene_Gene", "refGene_Func", "refGene_Gene",
                           "motif", "representative_motif", "primary_family",
                           "n_tag_snps",
                           "snp_BETA",          "snp_LOG10_P",
                           "snp_sq_BETA",       "snp_sq_LOG10_P",
                           "pgs_loco_BETA",     "pgs_loco_LOG10_P",
                           "int_pgs_loco_BETA", "int_pgs_loco_LOG10_P",
                           "pgs_tf_BETA",       "pgs_tf_LOG10_P",
                           "int_pgs_tf_BETA",   "int_pgs_tf_LOG10_P")]
save(indep_hits_agg,
     file = "../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.RData")

## round numbers before exporting as text file
cols <- c("MAF",
          "snp_BETA",          "snp_LOG10_P",
          "snp_sq_BETA",       "snp_sq_LOG10_P",
          "pgs_loco_BETA",     "pgs_loco_LOG10_P",
          "int_pgs_loco_BETA", "int_pgs_loco_LOG10_P",
          "pgs_tf_BETA",       "pgs_tf_LOG10_P",
          "int_pgs_tf_BETA",   "int_pgs_tf_LOG10_P")
indep_hits_agg_out <- as.data.table(indep_hits_agg)
indep_hits_agg_out[,(cols) := round(.SD, 8), .SDcols = cols]
fwrite(indep_hits_agg_out,
       file = "../results/04-tf-binding/gwas/aggregate/int-gwas-indep-hits-aggregate.tab",
       sep = "\t", na = "NA", quote = FALSE)
