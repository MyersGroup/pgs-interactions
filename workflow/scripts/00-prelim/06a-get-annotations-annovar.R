## Extract annotations from Annovar

library(data.table)
library(dplyr)
library(tidyr)

options(scipen=999)  # important for a few positions


## make input files
for (chr in 1:22) {
    
    var_chr <- fread(paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", chr, ".tab"),
                     header = TRUE, data.table = FALSE)

    var_chr <- var_chr %>%
        separate(POSID,
                 into = c("CHR", "START", "REF", "ALT"),
                 sep  = ":",
                 remove = FALSE,
                 convert = TRUE)
    var_chr$END <- var_chr$START

    ## variants with multiple bases in REF must span that number of bases from START to END
    ## otherwise they're not compatible with Annovar
    var_chr$REF_LN <- nchar(var_chr$REF)
    var_chr$ALT_LN <- nchar(var_chr$ALT)
    var_chr$END[var_chr$REF_LN > 1] <- var_chr$END[var_chr$REF_LN > 1] + var_chr$REF_LN[var_chr$REF_LN > 1] - 1
    
    var_chr$COMMENTS <- paste0("comments: ", var_chr$POSID, ",", var_chr$rsID)
    var_chr <- var_chr[, c("CHR", "START", "END", "REF", "ALT", "COMMENTS")]
    
    fwrite(var_chr,
           file = paste0("../data/annotations/annovar/inputs/chr", chr, ".avinput"),
           sep = " ", col.names = FALSE, quote = FALSE)
}


## retrieve annotations from local databases downloaded from
##   https://annovar.openbioinformatics.org/en/latest/user-guide/download/
##
## - refGene (date: 20211019; last update was 2020-08-17 at UCSC)
## - GENCODE (date: 20241008; last update was 2024-05-13 at UCSC)

for (chr in 1:22) {

    system(paste0("/path/to/annovar/table_annovar.pl ",
                  "../data/annotations/annovar/inputs/chr", chr, ".avinput ",
                  "/path/to/annovar/humandb/ ",
                  "-build hg19 ",
                  "-protocol refGene,ensGene ",
                  "-operation g,g ",
                  "-out ../data/annotations/annovar/gene-anno/annovar-chr", chr))
    
    system(paste0("gzip ../data/annotations/annovar/gene-anno/annovar-chr", chr, ".refGene.variant_function"))
    system(paste0("gzip ../data/annotations/annovar/gene-anno/annovar-chr", chr, ".ensGene.variant_function"))
    system(paste0("gzip ../data/annotations/annovar/gene-anno/annovar-chr", chr, ".hg19_multianno.txt"))
}
