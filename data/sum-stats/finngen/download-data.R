library(data.table)
setDTthreads(1)
library(dplyr)


system("wget https://storage.googleapis.com/finngen-public-data-r12/lab_values/analysis_documentation/Kanta_labs_GWAS_results_v2_summary.txt")
system("chmod 440 Kanta_labs_GWAS_results_v2_summary.txt")


## load list of phenotypes + short descriptions
names <- fread("../../phenotypes/clean/code-desc-map-real-qn-paper.tab", data.table = FALSE)

## load table of SNP*PGS hits
hits_loco <- fread("../../../results/03-interaction-gwas/indep-hits/aggregate/hits-loco-qn.tab", data.table = FALSE)

## find traits in "Blood biochemistry", "Blood count" or "Urine assays"
## categories for which we have found SNP*PGS interactions
lab_traits <- names[names$category %in% c("Blood biochemistry", "Blood count", "Urine assays"),]
lab_traits$int <- FALSE
for (i in 1:nrow(lab_traits)) {
    phen <- lab_traits$field_id[i]
    if (phen %in% hits_loco$field_id) {
        lab_traits$int[i] <- TRUE
    }
}
lab_traits <- lab_traits[lab_traits$int,]
lab_traits$int <- NULL

## match to FinnGen phenotypes
fg_ukb <- fread("ukb-finngen-lab-values-correspondence.csv", data.table = FALSE)
lab_traits <- left_join(lab_traits, fg_ukb[, c("field_id", "ukb_units", "OMOPID", "phenostring", "fg_units")], by = "field_id")
lab_traits <- lab_traits[!is.na(lab_traits$OMOPID),]

## export
fwrite(lab_traits,
       file = "ukb-finngen-lab-values-match.tab",
       sep = "\t", na = "NA", quote = FALSE)


## download summary stats
for (omopid in lab_traits$OMOPID) {
    system(paste0("wget https://storage.googleapis.com/finngen-public-data-r12/lab_values/summary_stats/finngen_R12_", omopid, ".gz"))
    system(paste0("chmod 440 finngen_R12_", omopid, ".gz"))
}
