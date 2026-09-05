## Find tags of all LHS/RHS (LOCO QN) hits for fine-mapping

library(data.table)
setDTthreads(4)
library(dplyr)
library(stringr)

options(warn = 2)

plink2 <- "/path/to/plink2"


## load LHS hits
lhs <- fread("../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab", data.table = FALSE)

## load pairwise hits
indep_hits <- readRDS("../results/03-interaction-gwas/indep-hits/aggregate/pairwise-hits-all.RData")
## apply "high confidence"/"moderate confidence" filters
indep_hits <- indep_hits[indep_hits$LOG10_P_int >= -log10(5e-8) &
                         str_detect(indep_hits$field_id, "_qn") &
                         indep_hits$PGS == "LOCO" &
                         (indep_hits$CHR_1 != indep_hits$CHR_2) &
                         (indep_hits$freq_int >= 0.001),]

## make list of unique hits
hits <- unique(c(lhs$ID, indep_hits$ID_1, indep_hits$ID_2))

## order by position
map_posid_rsid_f <- paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", seq(1, 22), ".tab")
map_posid_rsid <- lapply(map_posid_rsid_f, function(x) fread(x, data.table = FALSE))
map_posid_rsid <- do.call("rbind", map_posid_rsid)
hits <- hits[order(match(hits, map_posid_rsid$POSID))]

## add chromosome
hits <- data.frame(ID = hits)
hits$CHR <- sapply(strsplit(hits$ID, split = ":"), "[[", 1)

## find tags
for (chr in 1:22) {

    ## export list to be read by Plink
    write.table(hits$ID[hits$CHR == chr],
                file = paste0("../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-chr", chr, ".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    ## find tags
    system(paste0(plink2,
                  " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                  "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                  "--psam ../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam ",
                  "--extract ../data/variant-ids/chr", chr, "-var-ids.tab ",
                  "--read-freq ../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                  "--keep ../data/sample-ids/filtered/wb_all-ids.tab ",
                  "--r2-unphased ",
                  "--ld-window-kb 500 ",
                  "--ld-window-r2 0.8 ",
                  "--ld-snp-list ../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-chr", chr, ".txt ",
                  "--threads 4 ",
                  "--memory 64000 ",
                  "--out ../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-chr", chr, "-tags"))
}

## load tags
tags <- fread("../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-chr1-tags.vcor",
              data.table = FALSE)
for (chr in 2:22) {
    tags_chr <- fread(paste0("../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-chr", chr, "-tags.vcor"),
                      data.table = FALSE)
    tags <- rbind(tags, tags_chr)
}

## make table symmetric: if B is a tag of A, A should be a tag of B
tags_ab <- tags
tags_ba <- tags[, c("CHROM_B", "POS_B", "ID_B",
                    "#CHROM_A", "POS_A", "ID_A",
                    "UNPHASED_R2")]
colnames(tags_ba) <- colnames(tags_ab)
tags <- rbind(tags_ab, tags_ba)
tags <- tags[tags$ID_A %in% hits$ID,]


## add MAF
maf <- lapply(paste0("../data/imputed-wb_all-stats/minor-alleles/chr", seq(1, 22), ".tab"),
              function(x) fread(x, data.table = FALSE))
maf <- do.call("rbind", maf)
tags <- left_join(tags, maf[, c("ID", "MAF")], by = c("ID_B" = "ID"))

## add functional information
## VEP: (most severe) variant consequence
vep_f <- paste0("../data/annotations/vep/anno/vep-chr", seq(1, 22), ".most_severe.tab.gz")
vep <- lapply(vep_f, function(x) fread(x, skip = "#Uploaded_variation", data.table = FALSE))
vep <- do.call("rbind", vep)
vep <- vep[, c("#Uploaded_variation", "Consequence")]
colnames(vep) <- c("ID", "VEP_Consequence")
## LHS
tags <- left_join(tags, vep, by = c("ID_B" = "ID"))

## Annovar: affected/closest gene(s) from Gencode and refSeq
## Gencode
annovar_var_fn_f <- paste0("../data/annotations/annovar/gene-anno/annovar-chr", seq(1, 22), ".ensGene.variant_function.gz")
annovar_var_fn <- lapply(annovar_var_fn_f,
                         function(x) fread(x,
                                           header = FALSE,
                                           sep = "\t",
                                           col.names = c("ensGene_Func", "ensGene_Gene", "Comments"),
                                           data.table = FALSE))
annovar_var_fn <- do.call("rbind", annovar_var_fn)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$Comments, split = "comments: "), `[[`, 2)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$ID, split = ","), `[[`, 1)
## LHS
tags <- left_join(tags, annovar_var_fn[, c("ID", "ensGene_Func", "ensGene_Gene")], by = c("ID_B" = "ID"))
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
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$Comments, split = "comments: "), `[[`, 2)
annovar_var_fn$ID <- sapply(strsplit(annovar_var_fn$ID, split = ","), `[[`, 1)
## LHS
tags <- left_join(tags, annovar_var_fn[, c("ID", "refGene_Func", "refGene_Gene")], by = c("ID_B" = "ID"))
rm(annovar_var_fn)

## add VEP rank information and make "exonic" indicator
vep_consq <- fread("../data/annotations/vep/vep-var-consq-rel113_202410.csv", data.table = FALSE)
## add explicit ranking
vep_consq$VEP_Rank <- seq(1, nrow(vep_consq))
## add rank to main table
tags <- left_join(tags, vep_consq[, c("SO term", "VEP_Rank")], by = c("VEP_Consequence" = "SO term"))
## order the tags of each hit first by rank (more severe first) then by R^2 (descending) if there are ties
tags <- tags %>%
    arrange(ordered(ID_A, unique(ID_A)), VEP_Rank, desc(UNPHASED_R2))

## add rsIDs
tags <- left_join(tags, map_posid_rsid, by = c("ID_A" = "POSID"))
colnames(tags)[ncol(tags)] <- "rsID_A"
tags <- left_join(tags, map_posid_rsid, by = c("ID_B" = "POSID"))
colnames(tags)[ncol(tags)] <- "rsID_B"

## export
tags_out <- tags
colnames(tags_out)[c(1, 4, 8)] <- c("CHR", "CHR_B", "MAF_B")
tags_out <- tags_out %>%
    select(CHR,
           POS_A, ID_A, rsID_A,
           UNPHASED_R2,
           POS_B, ID_B, rsID_B, MAF_B,
           VEP_Consequence, VEP_Rank,
           ensGene_Func, ensGene_Gene,
           refGene_Func, refGene_Gene)
fwrite(tags_out,
       "../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-loco-qn-tags-0.8.tab",
       sep = "\t", na = "NA", quote = FALSE)


## keep only most severe tag for each hit
tags_ind <- tags_out %>%
    group_by(ID_A) %>%
    slice_head(n = 1) %>%
    as.data.frame()
## restore original order
tags_ind <- tags_ind[order(match(tags_ind$ID_A, map_posid_rsid$POSID)),]
## make indicator of whether each hit tags an "exonic" variant
tags_ind$VEP_Exonic <- (tags_ind$VEP_Rank <= 26)
## export
tags_ind <- tags_ind %>%
    select(CHR,
           POS_A, ID_A, rsID_A,
           UNPHASED_R2,
           POS_B, ID_B, rsID_B, MAF_B,
           VEP_Consequence, VEP_Exonic, VEP_Rank,
           ensGene_Func, ensGene_Gene,
           refGene_Func, refGene_Gene)
fwrite(tags_ind,
       "../results/03-interaction-gwas/indep-hits/aggregate/tags/hits-loco-qn-tags-0.8-most-severe.tab",
       sep = "\t", na = "NA", quote = FALSE)
