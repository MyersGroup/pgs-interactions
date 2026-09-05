## Prepare Relate results for downstream use

library(data.table)
setDTthreads(4)

options(warn = 2)  # turn warnings into errors


## load ancestral/derived allele info from Relate and save into a single file for all chromosomes
allsnprelate <- matrix(nrow = 0, ncol = 6)
colnames(allsnprelate) <- c("rs.id", "is_not_mapping", "is_flipped",
                            "age_begin", "age_end", "ancestral_allele.alternative_allele")

for(chr in 1:22){

    relate_chr <- fread(paste0("../data/1000G/relate/mut/1000GP_Phase3_mask_prene_chr", chr, ".mut.gz"),
                        col.names = c("snp", "pos_of_snp", "dist", "rs.id", "tree_index",
                                      "branch_indices", "is_not_mapping", "is_flipped", "age_begin",
                                      "age_end", "ancestral_allele.alternative_allele", "V12"),
                        data.table = FALSE)
    newmat <- relate_chr[, c("rs.id", "is_not_mapping", "is_flipped",
                             "age_begin", "age_end", "ancestral_allele.alternative_allele")]
    allsnprelate <- rbind(allsnprelate, newmat) 

}

save(allsnprelate,
     file = "../data/1000G/relate/mut/all-snps-relate.RData")
