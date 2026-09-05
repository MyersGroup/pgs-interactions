## Make alternative PVAR files with variants IDs of the form CHR:POS:REF:ALT
## rather than rsIDs (e.g. rs1234567) as is the default, since the latter are
## not unique.
## Also add INFO column with imputation score.


library(data.table)
library(dplyr)


for (chr in 1:22) {

    ## load "original" .pvar file resulting from BGEN>PGEN conversion
    pvar <- fread(paste0("../ukb_imp_chr", chr, "_v3.pvar"), data.table = FALSE)


    ## load MFI file provided with the imputed data
    mfi <- fread(paste0("PATH_TO_MFI/ukb_mfi_chr", chr, "_v3.txt"), data.table = FALSE)
    colnames(mfi) <- c("Alternate_id", "RS_id", "POS", "REF", "ALT", "MAF", "Minor_Allele", "INFO")


    ## add INFO column to pvar df
    pvar_mfi <- left_join(subset(pvar, select = -ID),  # remove rsID column
                          mfi[, c("POS", "REF", "ALT", "INFO")],  # keep only pos/ref/alt and INFO
                          by = c("POS", "REF", "ALT"))


    ## make new variant ID: chrom:pos:ref:alt
    pvar_mfi$ID <- paste0(chr, ":",
                          pvar_mfi$POS, ":",
                          pvar_mfi$REF, ":",
                          pvar_mfi$ALT)


    ## reorder cols
    pvar_mfi <- pvar_mfi[, c("#CHROM", "POS", "ID", "REF", "ALT", "INFO")]


    ## recode INFO score to VCF 4.3 format (https://samtools.github.io/hts-specs/VCFv4.3.pdf)
    pvar_mfi$INFO <- as.character(pvar_mfi$INFO)
    pvar_mfi$INFO[is.na(pvar_mfi$INFO)] <- "."
    pvar_mfi$INFO <- paste0("IIS=", pvar_mfi$INFO)


    ## check that final order of POS, REF and ALT is the same as in the original pvar files
    if(!all.equal(pvar_mfi$POS, pvar$POS) |
       !all.equal(pvar_mfi$REF, pvar$REF) |
       !all.equal(pvar_mfi$ALT, pvar$ALT)) {
        stop(paste("chr", chr, "variants are not in their original order."))
    }

    
    ## export
    fwrite(pvar_mfi, file = paste0("ukb_imp_POSID_INFO_chr", chr, "_v3.pvar"),
           quote = FALSE, sep = "\t", na = ".")


    ## add first line with details on INFO column (required by VCF 4.3 format)
    system(paste0("echo '##INFO=<ID=IIS,Number=1,Type=Float,Description=\"Imputation quality score\">' | cat - ",
                  "ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                  "> temp && mv temp ",
                  "ukb_imp_POSID_INFO_chr", chr, "_v3.pvar"))

}
