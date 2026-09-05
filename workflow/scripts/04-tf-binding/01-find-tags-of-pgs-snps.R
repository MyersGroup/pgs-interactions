## Make df with PGS SNPs and their tags with Bayes Factors

library(argparser)
library(data.table)
setDTthreads(4)
library(dplyr)
library(rtracklayer)

options(warn = 2)  # turn warnings into errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples", help = "Sample code",    nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples

plink2 <- "/path/to/plink2"



## find all tags (R2 >= 0.75 in 1Mb window) of independent PGS SNPs
pgs_chrs <- c()
system(paste0("mkdir -p ../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps"))
system(paste0("mkdir -p ../results/04-tf-binding/pgs-snps-tags/", phen, "/tags"))
for (chr in 1:22) {

    ## load independent hits from original GWAS
    clump <- fread(paste0("../results/02-pgs/ld-clump/", phen, "/p1_0.05-r2_0.1-kb_500/",
                          samples, "-chr", chr, "-clumped.tab"), data.table = FALSE)

    ## keep only GWS hits
    clump <- clump[clump$LOG10_P >= -log10(5e-8),]

    if (nrow(clump) == 0) {
        next
    } else {
        pgs_chrs <- c(pgs_chrs, chr)
    }
    
    ## export list to be read by Plink
    write.table(unique(clump$ID),
                file = paste0("../results/04-tf-binding/pgs-snps-tags/", phen, "/",
                              "pgs-snps/pgs-snps-chr", chr, ".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    ## find tags
    system(paste0(plink2,
                  " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                  "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                  "--psam ../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam ",
                  "--extract ../data/variant-ids/chr", chr, "-var-ids.tab ",
                  "--read-freq ../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                  "--keep ../data/sample-ids/filtered/", samples, "-ids.tab ",
                  "--r2-unphased ",
                  "--ld-window-kb 500 ",
                  "--ld-window-r2 0.75 ",
                  "--ld-snp-list ../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps/pgs-snps-chr", chr, ".txt ",
                  "--threads 4 ",
                  "--memory 64000 ",
                  "--out ../results/04-tf-binding/pgs-snps-tags/", phen, "/tags/pgs-snps-tags-chr", chr))   
}


## load tags
tags <- fread(paste0("../results/04-tf-binding/pgs-snps-tags/", phen, "/tags/pgs-snps-tags-chr", pgs_chrs[1], ".vcor"), data.table = FALSE)
if (length(pgs_chrs) > 1) {
    for (chr in pgs_chrs[-1]) {
        tags_chr <- fread(paste0("../results/04-tf-binding/pgs-snps-tags/", phen, "/tags/pgs-snps-tags-chr", chr, ".vcor"), data.table = FALSE)
        tags <- rbind(tags, tags_chr)
    }
}
## remove temporary files
for (chr in pgs_chrs) {
    system(paste0("rm ../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps/pgs-snps-chr", chr, ".txt"))
    system(paste0("rm ../results/04-tf-binding/pgs-snps-tags/", phen, "/tags/pgs-snps-tags-chr", chr, ".vcor"))
}
system(paste0("rmdir ../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps"))

##  keep at most 100 tags for each hit SNP (ordered by strength of LD)
tags_top100 <- tags %>%
    group_by(ID_A) %>%
    slice_max(order_by = UNPHASED_R2, n = 100, with_ties = FALSE)  # this reorders rows
tags <- tags[paste(tags$ID_A, tags$ID_B, sep = "_") %in% paste(tags_top100$ID_A, tags_top100$ID_B, sep = "_"),]


## prepare output data.frame
pgs_snps_tags <- tags[, c("ID_A", "ID_B", "CHROM_B")]
colnames(pgs_snps_tags) <- c("hit_ID", "tag_ID", "tag_CHR")

## add each hit SNP as its first own tag
new <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(new) <- colnames(pgs_snps_tags)
for (hit in unique(pgs_snps_tags$hit_ID)) {
    new <- rbind(new,
                 data.frame(hit_ID = hit,
                            tag_ID = hit,
                            tag_CHR = sapply(strsplit(hit, split = ":"), "[[", 1)))
    new <- rbind(new,
                 pgs_snps_tags[pgs_snps_tags$hit_ID == hit,])
}
pgs_snps_tags <- new
    

## load CHR:POS:REF:ALT -> rsID map
posid_rsid_f <- paste0("../data/variant-ids/map-posid-rsid/posid-rsid-chr", seq(1, 22), ".tab")
posid_rsid <- lapply(posid_rsid_f, function(x) fread(x, header = TRUE, data.table = FALSE))
posid_rsid <- do.call("rbind", posid_rsid)
colnames(posid_rsid)[1] <- "ID"

## add rsIDs
pgs_snps_tags <- left_join(pgs_snps_tags, posid_rsid, by = c("hit_ID" = "ID"))
colnames(pgs_snps_tags)[4] <- "hit_rsID"
pgs_snps_tags <- left_join(pgs_snps_tags, posid_rsid, by = c("tag_ID" = "ID"))
colnames(pgs_snps_tags)[5] <- "tag_rsID"


## add position of tags in *GRCh38* to match with epigenetic data and for scoring motifs

## extract POS, REF, ALT information for each tag
pgs_snps_tags$tag_POS <- as.numeric(sapply(strsplit(pgs_snps_tags$tag_ID, split = ":"), "[[", 2))
pgs_snps_tags$tag_REF <- sapply(strsplit(pgs_snps_tags$tag_ID, split = ":"), "[[", 3)
pgs_snps_tags$tag_ALT <- sapply(strsplit(pgs_snps_tags$tag_ID, split = ":"), "[[", 4)

## make GenomicRanges object with necessary info
tags_pos <- data.frame(chrom  = paste0("chr", pgs_snps_tags$tag_CHR),
                       start  = pgs_snps_tags$tag_POS,
                       end    = pgs_snps_tags$tag_POS,
                       strand = rep("+", nrow(pgs_snps_tags)),
                       tag_ID = pgs_snps_tags$tag_ID)
tags_pos <- unique(tags_pos)
GR <- GenomicRanges::makeGRangesFromDataFrame(tags_pos, keep.extra.columns = TRUE)

## convert from hg19 to hg38
chain_hg19ToHg38 <- rtracklayer::import.chain("../data/liftOver/hg19ToHg38.over.chain")
GRhg38 <- rtracklayer::liftOver(GR, chain_hg19ToHg38)
GRhg38 <- as.data.frame(unlist(GRhg38))
colnames(GRhg38)[2] <- "tag_POS_hg38"

## join (note that this removes a small number of SNPs that couldn't be mapped across assemblies)
pgs_snps_tags <- inner_join(pgs_snps_tags, GRhg38[, c("tag_ID", "tag_POS_hg38")], by = "tag_ID")


## add ancestral/derived allele info from Relate
load("../data/1000G/relate/mut/all-snps-relate.RData")
allsnprelate$tag_anc <- sapply(strsplit(allsnprelate$ancestral_allele.alternative_allele, split = "/"), "[[", 1)
allsnprelate$tag_der <- sapply(strsplit(allsnprelate$ancestral_allele.alternative_allele, split = "/"), "[[", 2)
pgs_snps_tags <- left_join(pgs_snps_tags, allsnprelate[, c("rs.id", "tag_anc", "tag_der")], by = c("tag_rsID" = "rs.id"))


## compute Bayes Factors

## load original GWAS summary stats
orig_f <- paste0("../results/01-gwas/plink-output/", phen, "/", samples, ".chr", seq(1, 22), ".", phen, ".res.cov.glm.linear.gz")
orig <- lapply(orig_f, function(x) fread(x, header = TRUE, data.table = FALSE))
orig <- do.call("rbind", orig)

## add p-value to main data.frame
pgs_snps_tags <- left_join(pgs_snps_tags, orig[, c("ID", "LOG10_P")], by = c("tag_ID" = "ID"))

## compute z-scores and from these BFs (in favour of the alternative hypothesis)
K <- 1
pgs_snps_tags$z2 <- qchisq(10^(-pgs_snps_tags$LOG10_P), df = 1, lower.tail = FALSE)
pgs_snps_tags$bf <- 1 / sqrt(1 + K) * exp((pgs_snps_tags$z2 / 2) * K / (1 + K))

## large values become nan or infinite, so make nan/infinite equal to maximum value
maxval <- max(pgs_snps_tags$bf[!is.na(pgs_snps_tags$bf) & !is.infinite(pgs_snps_tags$bf)], na.rm = TRUE)
pgs_snps_tags$bf[is.nan(pgs_snps_tags$bf) | is.infinite(pgs_snps_tags$bf)] <- maxval

## compute sum of BFs for each hit SNP
pgs_snps_tags <- pgs_snps_tags %>%
    group_by(hit_ID) %>%
    mutate(bftots = sum(bf)) %>%
    as.data.frame()

## posterior probability is ratio of BF to sum of BFs for the corresponding hit SNP
pgs_snps_tags$post_prob <- pgs_snps_tags$bf / pgs_snps_tags$bftots

## reorder cols and export
pgs_snps_tags <- pgs_snps_tags[, c("hit_ID", "hit_rsID", "tag_ID", "tag_rsID",
                                   "tag_CHR", "tag_POS", "tag_REF", "tag_ALT", "tag_POS_hg38",
                                   "tag_anc", "tag_der", "post_prob")]
save(pgs_snps_tags,
     file = paste0("../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps-tags.RData"))
