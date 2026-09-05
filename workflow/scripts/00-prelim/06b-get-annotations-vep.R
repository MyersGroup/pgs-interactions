## Extract annotations from Ensembl VEP

library(argparser)
library(data.table)
library(dplyr)
library(tidyr)

p <- arg_parser('Argument parser')
p <- add_argument(p, '--chr', help = 'Chromosome', nargs = 1)
argv <- parse_args(p)

chr <- argv$chr

options(scipen=999)  # important for a few positions



## make input file
var_chr <- fread(paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", chr, ".tab"),
                 header = TRUE, data.table = FALSE)

var_chr <- var_chr %>%
    separate(POSID,
             into = c("CHR", "START", "REF", "ALT"),
             sep  = ":",
             remove = FALSE,
             convert = TRUE)
var_chr$END <- var_chr$START
var_order <- var_chr$POSID  # save IDs in order for sorting later

## deal with indels
var_chr$REF_LN <- nchar(var_chr$REF)
var_chr$ALT_LN <- nchar(var_chr$ALT)
var_chr$MIN_LN <- pmin(var_chr$REF_LN, var_chr$ALT_LN)

var_chr_indel <- var_chr %>% filter(REF_LN > 1 | ALT_LN > 1)
var_chr <- var_chr %>% filter(REF_LN == 1 & ALT_LN == 1)

for (i in 1:nrow(var_chr_indel)) {
    ref_ln <- var_chr_indel$REF_LN[i]
    alt_ln <- var_chr_indel$ALT_LN[i]
    min_ln <- var_chr_indel$MIN_LN[i]

    if (alt_ln > ref_ln) {  # insert: VEP requires start = end + 1
        var_chr_indel$ALT[i] <- substr(var_chr_indel$ALT[i],
                                       start = 1 + min_ln,
                                       stop  = alt_ln)
        var_chr_indel$REF[i] <- "-"
        var_chr_indel$END[i] <- var_chr_indel$END[i] + min_ln
        var_chr_indel$START[i] <- var_chr_indel$END[i] + 1

    } else if (alt_ln < ref_ln) {  # del
        var_chr_indel$REF[i] <- substr(var_chr_indel$REF[i],
                                       start = 1 + min_ln,
                                       stop  = nchar(var_chr_indel$REF[i]))
        var_chr_indel$ALT[i] <- "-"
        var_chr_indel$START[i] <- var_chr_indel$START[i] + min_ln
        var_chr_indel$END[i]   <- var_chr_indel$START[i] + nchar(var_chr_indel$REF[i]) - 1
    }
}

var_chr <- rbind(var_chr, var_chr_indel)
var_chr <- var_chr[match(var_order, var_chr$POSID),]

var_chr$ALLELE <- paste0(var_chr$REF, "/", var_chr$ALT)
var_chr$STRAND <- "+"
var_chr <- var_chr[, c("CHR", "START", "END", "ALLELE", "STRAND", "POSID")]

fwrite(var_chr,
       file = paste0("../data/annotations/vep/inputs/chr", chr, ".tab.gz"),
       sep = "\t", col.names = FALSE, quote = FALSE)



## download annotations (Ensembl database version: 113)
## most severe consequence for each variant
system(paste0("vep --verbose ",
              "--input_file ../data/annotations/vep/inputs/chr", chr, ".tab.gz ",
              "--format ensembl ",
              "--cache --dir_cache ../data/annotations/vep/cache ",
              "--assembly GRCh37 ",
              "--most_severe ",
              "--output_file ../data/annotations/vep/anno/vep-chr", chr, ".most_severe.tab.gz ",
              "--tab ",
              "--force_overwrite ",
              "--compress_output gzip"))

## 'everything'
system(paste0("vep --verbose ",
              "--input_file ../data/annotations/vep/inputs/chr", chr, ".tab.gz ",
              "--format ensembl ",
              "--cache --dir_cache ../data/annotations/vep/cache ",
              "--assembly GRCh37 ",
              "--everything ",
              "--output_file ../data/annotations/vep/anno/vep-chr", chr, ".everything.tab.gz ",
              "--tab ",
              "--force_overwrite ",
              "--compress_output gzip"))
