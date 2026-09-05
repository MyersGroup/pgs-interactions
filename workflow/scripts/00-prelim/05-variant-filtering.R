## desc: Filter variants for further analysis by:
##       - MAF >= 0.001
##       - MAC >= 25
##       - genotype missing rate <= 0.05
##       - HWE p-value >= 1e-10
##       - INFO imputation score >= 0.8
##       Make POSID -> rsID mapping file.
##       Make MAF reference file.


library(data.table)
library(dplyr)
library(stringr)



## initialise vector of filtered variants across all chrs
all_chr_vars <- character()

## initialise empty POSID -> rsID mapping df
id_map <- data.frame(POSID = character(),
                     rsID  = character())

## initialise empty MAF reference file
maf_ref <- data.frame(ID = character(),
                      MAF   = numeric())


for (chr in 1:22) {
    freq    <- fread(paste0('../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr', chr, '_v3.afreq'), data.table = FALSE)
    counts  <- fread(paste0('../data/imputed-wb_all-stats/allele-counts/ukb_imp_chr', chr, '_v3.acount'), data.table = FALSE)
    missing <- fread(paste0('../data/imputed-wb_all-stats/missing-rates/ukb_imp_chr', chr, '_v3.vmiss'), data.table = FALSE)
    hardy   <- fread(paste0('../data/imputed-wb_all-stats/hardy-eq/ukb_imp_chr', chr, '_v3.hardy'), data.table = FALSE)
    info    <- fread(paste0('../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr', chr, '_v3.pvar'), data.table = FALSE)


    ## convert INFO field from 'IIS=0.1234' to numeric field
    info$score   <- substr(info$INFO, 5, nchar(info$INFO))
    info$score[info$score == '.'] <- NA
    info$score <- as.numeric(info$score)


    ## MAF >= 0.001  <=>  ALT_AF >= 0.001 & ALT_AF <= 0.999
    vars_freq    <- freq$ID[!is.na(freq$ALT_FREQS) &
                            freq$ALT_FREQS >= 0.001 &
                            freq$ALT_FREQS <= 0.999]
    ## MAC >= 25     <=>  ALT_AC >= 25 & ALT_AC <= OBS_CT - 25
    vars_counts  <- counts$ID[!is.na(counts$ALT_CTS) &
                              counts$ALT_CTS >= 25 &
                              counts$ALT_CTS <= counts$OBS_CT - 25]
    ## genotype missing rate <= 0.05
    vars_missing <- missing$ID[!is.na(missing$F_MISS) &
                               missing$F_MISS <= 0.05]
    ## HWE p-value >= 1e-10
    vars_hardy   <- hardy$ID[!is.na(hardy$MIDP) &
                             hardy$MIDP >= 1e-10]
    ## INFO imputation score >= 0.8
    vars_info    <- info$ID[!is.na(info$score) &
                            info$score >= 0.8]


    ## intersect
    vars_filter <- Reduce(intersect, list(vars_freq, vars_counts, vars_missing, vars_hardy, vars_info))
    ## export
    fwrite(list(vars_filter), file = paste0('../data/variant-ids/chr', chr, '-var-ids.tab'),
           row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', na = 'NA')



    ## make POSID -> rsID mapping df
    pvar <- fread(paste0('../data/imputed-genotypes/ukb_imp_chr', chr, '_v3.pvar'), data.table = FALSE)
    colnames(pvar)[c(1, 3)] <- c('CHROM', 'rsID')
    id_map_chr <- pvar %>%
        mutate(POSID = paste0(CHROM, ":", POS, ":", REF, ":", ALT)) %>%
        select(POSID, rsID) %>%
        filter(POSID %in% vars_filter)
    fwrite(id_map_chr, file = paste0('../data/variant-ids/map-posid-rsid/posid-rsid-chr', chr, '.tab'),
           row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')



    ## make file with minor allele/A1 and MAF
    maf_ref_chr <- freq %>%
        mutate(POS = as.numeric(word(ID, 2, sep = ':')),
               A1 = case_when(ALT_FREQS <= 0.5 ~ ALT,
                              ALT_FREQS >  0.5 ~ REF),
               MAF = pmin(ALT_FREQS, 1 - ALT_FREQS)) %>%
        select(ID, POS, REF, ALT, A1, ALT_FREQS, MAF) %>%
        filter(ID %in% vars_filter)
    fwrite(maf_ref_chr, file = paste0('../data/imputed-wb_all-stats/minor-alleles/chr', chr, '.tab'),
           row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t', na = 'NA')
}
