library(data.table)
setDTthreads(4)
library(dplyr)
library(tidyr)


## load list of traits to analyse
lab_traits <- fread("../data/sum-stats/finngen/ukb-finngen-lab-values-match.tab", data.table = FALSE)

## make list of all SNPs across all traits to convert from hg38 to hg19
gwas <- fread(paste0("../data/sum-stats/finngen/finngen_R12_", lab_traits$OMOPID[1], ".gz"), data.table = FALSE)
main <- gwas[, 1:4]
main$ID <- paste(gwas[,1], gwas[,2], gwas[,3], gwas[,4], sep = ":")

for (omopid in lab_traits$OMOPID[-1]) {
    print(omopid)

    ## load sumstats
    gwas <- fread(paste0("../data/sum-stats/finngen/finngen_R12_", omopid, ".gz"), data.table = FALSE)
    gwas <- gwas[, 1:4]
        
    ## keep only new variants
    gwas$ID <- paste(gwas[,1], gwas[,2], gwas[,3], gwas[,4], sep = ":")
    gwas <- gwas[!(gwas$ID %in% main$ID),]

    ## add to main
    main <- rbind(main, gwas)
}

## sort by physical position
main <- main %>%
    arrange(`#chrom`, pos)

## export
fwrite(main,
       file = "../data/sum-stats/finngen/variant-ids-hg38.tab.gz",
       sep = "\t", na = "NA", quote = FALSE)

## write VCF file
vcf <- main[, c("#chrom", "pos", "ID", "ref", "alt")]
colnames(vcf) <- c("#CHROM", "POS", "ID", "REF", "ALT")

## remove `chr23`
vcf <- vcf[vcf$`#CHROM` %in% 1:22,]

## change chromosome notation from 1, 2, etc. to chr1, chr2, etc. to match reference sequence fasta files
vcf$`#CHROM` <- paste0("chr", vcf$`#CHROM`)

## add remaining mandatory VCF fields
vcf$QUAL <- "."
vcf$FILTER <- "."
vcf$INFO <- "."

## export
fwrite(vcf,
       file = "../data/sum-stats/finngen/variant-ids-hg38.vcf",
       sep = "\t", na = "NA", quote = FALSE)

## add header line
system("sed -i '1i##fileformat=VCFv4.3' ../data/sum-stats/finngen/variant-ids-hg38.vcf")

## must compress with bgzip and build index before converting
system("bgzip ../data/sum-stats/finngen/variant-ids-hg38.vcf")
system("tabix -p vcf ../data/sum-stats/finngen/variant-ids-hg38.vcf.gz")



## convert using the liftover plugin for bcftools
system(paste0("bcftools plugin ../../../software/score/liftover.so ",
              "-Ou ../data/sum-stats/finngen/variant-ids-hg38.vcf.gz ",
              "-- -s ../data/liftOver/hg38.fa ",
              "-f ../data/liftOver/hg19.fa ",
              "-c ../data/liftOver/hg38ToHg19.over.chain ",
              "| bcftools sort -Oz -o ../data/sum-stats/finngen/variant-ids-hg19.vcf.gz -W=tbi"))
## $ zcat ../data/sum-stats/finngen/variant-ids-hg19.vcf.gz | awk '!(/^#/)' | wc -l
## 19176052
## > nrow(vcf) - 19176052
## [1] 81854
## only 81854 out of 19257906 variants missing from converted file (0.425%)



## make list with only variants present in UKB imputed data

## load list of UKB imputed variants
pvar_f <- paste0("../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", seq(1, 22), "_v3.pvar")
pvar <- lapply(pvar_f, function(x)  fread(x, data.table = FALSE))
pvar <- do.call("rbind", pvar)

## load converted VCF file
vcf_hg19 <- fread("../data/sum-stats/finngen/variant-ids-hg19.vcf.gz",
                  skip = "#CHROM", data.table = FALSE)

## filter
vcf_hg19$ID_hg19 <- paste(gsub("chr", "", vcf_hg19$`#CHROM`), vcf_hg19$POS, vcf_hg19$REF, vcf_hg19$ALT, sep = ":")
vcf_hg19_flt <- vcf_hg19[vcf_hg19$ID_hg19 %in% pvar$ID,]
## > nrow(vcf_hg19_flt) / nrow(vcf_hg19)
## [1] 0.8257379
## 82.57% of FinnGen variants are present in the UKB imputed dataset

## remove non-unique matches (two hg38 variants assigned to the same hg19 variant)
hg38_hg19_ukb <- vcf_hg19_flt[, c("ID", "ID_hg19")]
colnames(hg38_hg19_ukb)[1] <- "ID_hg38"
hg38_hg19_ukb <- hg38_hg19_ukb[!duplicated(hg38_hg19_ukb$ID_hg19) & !duplicated(hg38_hg19_ukb$ID_hg19, fromLast = TRUE),]

## remove variants mapped to different chromosome in two assemblies
hg38_hg19_ukb$CHR_hg38 <- sapply(strsplit(hg38_hg19_ukb$ID_hg38, split = ":"), `[[`, 1)
hg38_hg19_ukb$CHR_hg19 <- sapply(strsplit(hg38_hg19_ukb$ID_hg19, split = ":"), `[[`, 1)
hg38_hg19_ukb <- hg38_hg19_ukb[hg38_hg19_ukb$CHR_hg38 == hg38_hg19_ukb$CHR_hg19,]
hg38_hg19_ukb$CHR_hg38 <- NULL
hg38_hg19_ukb$CHR_hg19 <- NULL

## export
fwrite(hg38_hg19_ukb,
       file = "../data/sum-stats/finngen/variant-ids-hg38-hg19-ukb.tab.gz",
       sep = "\t", na = "NA", quote = FALSE)



## make BIM file with variants in 1KGP EUR reference panel that overlap with intersection of FinnGen and UKB
## load info on SNPs in LD reference panel for PRS-CS
snpinfo_1kg <- fread("../../../software/PRScs/ldblk_1kg_eur/snpinfo_1kg_hm3", data.table = FALSE)

## add ID of form CHR:POS:REF:ALT
snpinfo_1kg$ID_hg19 <- paste(snpinfo_1kg$CHR, snpinfo_1kg$BP, snpinfo_1kg$A2, snpinfo_1kg$A1, sep = ":")

## intersect 1KGP with (intersection of) UKB and FinnGen variants
snpinfo_1kg_fg_ukb <- snpinfo_1kg[snpinfo_1kg$ID_hg19 %in% hg38_hg19_ukb$ID_hg19,]

## make BIM file
snpinfo_1kg_fg_ukb_bim <- snpinfo_1kg_fg_ukb
snpinfo_1kg_fg_ukb_bim$cM <- 0
snpinfo_1kg_fg_ukb_bim <- snpinfo_1kg_fg_ukb_bim[, c("CHR", "SNP", "cM", "BP", "A1", "A2")]
fwrite(snpinfo_1kg_fg_ukb_bim,
       file = "../results/05-new-app/finngen/bim/fg-ukb-1kg.bim",
       sep = "\t", na = "NA", quote = FALSE, col.names = FALSE)
