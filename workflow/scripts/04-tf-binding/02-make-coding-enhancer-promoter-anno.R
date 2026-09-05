## Prepare annotations on overlap with coding, enhancer and promoter sequences 

library(argparser)
library(data.table)
setDTthreads(4)
library(dplyr)

options(warn = 2)  # turn warnings into errors

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples", help = "Sample code",    nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples

source("scripts/misc/fn-get_overlap.R")
source("scripts/misc/fn-test_enrich_mat.R")
source("scripts/misc/fn-get_post.R")



## load data.frame with PGS hits and their tags and BFs
load(paste0("../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps-tags.RData"))


## read in H3K4me1 peaks and make annotation vectors
h3k4me1_datasets <- fread("../data/epigenome/list-blueprint-h3k4me1.txt",
                          col.names = "dataset", header = FALSE, data.table = FALSE)
    
h3k4me1_annot <- list()
for (ds in h3k4me1_datasets$dataset) {

    ## read peaks
    h3k4me1_peaks <- fread(paste0("../data/epigenome/bed/", ds, ".bed"), data.table = FALSE)

    ## extract chr, start, end
    target <- h3k4me1_peaks[, c(1, 2, 3)]
    target <- target[target[, 1] %in% paste("chr", 1:22, sep = ""),]  # keep autosomes only

    ## prepare TRUE/FALSE vector of whether each tag is in an epigenetic interval
    h3k4me1_logic <- rep(FALSE, nrow(pgs_snps_tags))

    for (chr in 1:22) {

        ## extract positions of tags for this chromome from pgs_snps_tags
        tags_pos_chr <- pgs_snps_tags[pgs_snps_tags$tag_CHR == chr, "tag_POS_hg38"]
        chr_indices <- which(pgs_snps_tags$tag_CHR == chr)

        ## peaks: make matrix with start/end positions from above
        peaks_chr <- as.matrix(target[target[,1] == paste0("chr", chr), 2:3])

        ## check overlap of each tag with the enhancer peaks
        overlap <- get_overlap_simple(peaks_chr, tags_pos_chr)

        ## final output: indices of rows in target matrix for all tag SNPs in current chromosome
        h3k4me1_logic[chr_indices] <- overlap
    }

    ## make annotation vector
    h3k4me1_annot[[ds]] <- rep(0, nrow(pgs_snps_tags))
    h3k4me1_annot[[ds]][h3k4me1_logic] <- 1
}


## read in H3K4me3 peaks and make annotation vectors
h3k4me3_datasets <- fread("../data/epigenome/list-blueprint-h3k4me3.txt",
                          col.names = "dataset", header = FALSE, data.table = FALSE)

h3k4me3_annot <- list()
for (ds in h3k4me3_datasets$dataset) {

    ## read peaks
    h3k4me3_peaks <- fread(paste0("../data/epigenome/bed/", ds, ".bed"), data.table = FALSE)

    ## extract chr, start, end
    target <- h3k4me3_peaks[, c(1, 2, 3)]
    target <- target[target[, 1] %in% paste("chr", 1:22, sep = ""),]  # keep autosomes only

    ## prepare TRUE/FALSE vector of whether each tag is in an epigenetic interval
    h3k4me3_logic <- rep(FALSE, nrow(pgs_snps_tags))

    for (chr in 1:22) {

        ## extract positions of tags for this chromome from pgs_snps_tags
        tags_pos_chr <- pgs_snps_tags[pgs_snps_tags$tag_CHR == chr, "tag_POS_hg38"]
        chr_indices <- which(pgs_snps_tags$tag_CHR == chr)

        ## peaks: make matrix with start/end positions from above
        peaks_chr <- as.matrix(target[target[,1] == paste0("chr", chr), 2:3])

        ## check overlap of each tag with the enhancer peaks
        overlap <- get_overlap_simple(peaks_chr, tags_pos_chr)

        ## final output: indices of rows in target matrix for all tag SNPs in current chromosome
        h3k4me3_logic[chr_indices] <- overlap
    }

    ## make annotation vector
    h3k4me3_annot[[ds]] <- rep(0, nrow(pgs_snps_tags))
    h3k4me3_annot[[ds]][h3k4me3_logic] <- 1
}


## H3K4me1: test enrichment of annotation in each dataset
qqannot_h3k4me1 <- list()
for (ds in h3k4me1_datasets$dataset) {
    qqannot_h3k4me1[[ds]] <- test_enrich_mat(pgs_snps_tags$post_prob,
                                             as.matrix(h3k4me1_annot[[ds]]),
                                             pgs_snps_tags[, c("hit_ID", "tag_ID")])
}

## find the dataset for which the likelihood is maximised
qqannot_h3k4me1_lhood_max <- data.frame(dataset = h3k4me1_datasets$dataset,
                                        lhood_max = rep(NA, length(h3k4me1_datasets$dataset)))
for (ds in h3k4me1_datasets$dataset) {
    vv2b <- apply(qqannot_h3k4me1[[ds]]$lhoodsu, 1, which.max)
    qqannot_h3k4me1_lhood_max$lhood_max[qqannot_h3k4me1_lhood_max$dataset == ds] <-
        qqannot_h3k4me1[[ds]]$lhoodsu[vv2b]
}
h3k4me1_ds_best <- qqannot_h3k4me1_lhood_max$dataset[which.max(qqannot_h3k4me1_lhood_max$lhood_max)]

## H3K4me3: test enrichment of annotation in each dataset
qqannot_h3k4me3 <- list()
for (ds in h3k4me3_datasets$dataset) {
    qqannot_h3k4me3[[ds]] <- test_enrich_mat(pgs_snps_tags$post_prob,
                                             as.matrix(h3k4me3_annot[[ds]]),
                                             pgs_snps_tags[, c("hit_ID", "tag_ID")])
}

## find the dataset for which the likelihood is maximised
qqannot_h3k4me3_lhood_max <- data.frame(dataset = h3k4me3_datasets$dataset,
                                        lhood_max = rep(NA, length(h3k4me3_datasets$dataset)))
for (ds in h3k4me3_datasets$dataset) {
    vv2b <- apply(qqannot_h3k4me3[[ds]]$lhoodsu, 1, which.max)
    qqannot_h3k4me3_lhood_max$lhood_max[qqannot_h3k4me3_lhood_max$dataset == ds] <-
        qqannot_h3k4me3[[ds]]$lhoodsu[vv2b]
}
h3k4me3_ds_best <- qqannot_h3k4me3_lhood_max$dataset[which.max(qqannot_h3k4me3_lhood_max$lhood_max)]

## export
save(h3k4me1_datasets, h3k4me1_annot, qqannot_h3k4me1, h3k4me1_ds_best,
     h3k4me3_datasets, h3k4me3_annot, qqannot_h3k4me3, h3k4me3_ds_best,
     file = paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/",
                   "annot/h3k4me1_3-annot.RData"))




## make combined annotation from VEP

## load (most severe) variant consequence
vep_f <- paste0("../data/annotations/vep/anno/vep-chr", seq(1, 22), ".most_severe.tab.gz")
vep <- lapply(vep_f, function(x) fread(x, skip = "#Uploaded_variation", data.table = FALSE))
vep <- do.call("rbind", vep)
vep <- vep[, c("#Uploaded_variation", "Consequence")]
colnames(vep) <- c("ID", "VEP_Consequence")

## check if each tag SNP's most severe consequence is in list of top 24
top_consq <- c("transcript_ablation",
               "splice_acceptor_variant",
               "splice_donor_variant",
               "stop_gained",
               "frameshift_variant",
               "stop_lost",
               "start_lost",
               "transcript_amplification",
               "feature_elongation",
               "feature_truncation",
               "inframe_insertion",
               "inframe_deletion",
               "missense_variant",
               "protein_altering_variant",
               "splice_donor_5th_base_variant",
               "splice_region_variant",
               "splice_donor_region_variant",
               "splice_polypyrimidine_tract_variant",
               "incomplete_terminal_codon_variant",
               "start_retained_variant",
               "stop_retained_variant",
               "synonymous_variant",
               "coding_sequence_variant",
               "mature_miRNA_variant")
vep$top_consq <- (vep$VEP_Consequence %in% top_consq)

## make annotation vector
vep_pgs_snps_tags <- left_join(pgs_snps_tags[, c("tag_ID"), drop = FALSE],
                               vep[, c("ID", "top_consq")],
                               by = c("tag_ID" = "ID"))
functional_annot <- rep(0, nrow(pgs_snps_tags))
functional_annot[vep_pgs_snps_tags$top_consq] <- 1

## export
save(functional_annot,
     file = paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/",
                   "annot/functional-annot.RData"))




## test for enrichment of SNPs within each 'confidence set' (e.g., tags associated
## with a driver SNP) in these three annotations separately
qqannot <- test_enrich_mat(pgs_snps_tags$post_prob,
                           cbind(functional_annot, h3k4me1_annot[[h3k4me1_ds_best]], h3k4me3_annot[[h3k4me3_ds_best]]),
                           pgs_snps_tags)

## take the alpha values that maximise the predictive power of each annotation and
## build single combined annotation based on the three types
annot_comb <- cbind(1, (1 + functional_annot * (qqannot$oru[1] - 1)) *
                       (1 + h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1)) *
                       (1 + h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1)) - 1)

## get posterior probability using new combined annotation
## (in using alphavalues = c(1, 2), we are using the factormat supplied as input, but with each tag having a unique driver)
newpost <- get_post(pgs_snps_tags$post_prob,
                    annot_comb,
                    pgs_snps_tags[, c("hit_ID", "tag_ID")],
                    alphavalues = c(1, 2))


## coding annotations
## keep only posterior probabilities for coding SNPs
coding_probs <- data.frame(hit_ID = pgs_snps_tags$hit_ID,
                           tag_ID = pgs_snps_tags$tag_ID,
                           functional_annot = functional_annot,
                           post_prob = newpost$post_prob)
coding_probs <- coding_probs[coding_probs$functional_annot == 1 &
                             !is.na(coding_probs$post_prob),]  # also discard repeated tag SNPs
## keep highest PP SNPs accounting for 80% of posterior density in total
temp <- order(coding_probs$post_prob, decreasing = TRUE)
frac <- cumsum(coding_probs$post_prob[temp]) / sum(coding_probs$post_prob, na.rm = TRUE)
cutoff <- coding_probs$post_prob[temp[which(frac > 0.8)[1]]]
coding_probs <- coding_probs[coding_probs$post_prob >= cutoff, c("hit_ID", "tag_ID", "post_prob")]
## export
save(coding_probs,
     file = paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/coding-pgs-snps.RData"))


## H3K4me1 annotations
## keep only posterior probabilities for SNPs in H3K4me1 marks
h3k4me1_probs <- data.frame(hit_ID = pgs_snps_tags$hit_ID,
                            tag_ID = pgs_snps_tags$tag_ID,
                            h3k4me1_annot = h3k4me1_annot[[h3k4me1_ds_best]],
                            post_prob = newpost$post_prob)
h3k4me1_probs <- h3k4me1_probs[h3k4me1_probs$h3k4me1_annot == 1 &
                               !is.na(h3k4me1_probs$post_prob),]  # also discard repeated tag SNPs
## keep highest PP SNPs accounting for 80% of posterior density in total
temp <- order(h3k4me1_probs$post_prob, decreasing = TRUE)
frac <- cumsum(h3k4me1_probs$post_prob[temp]) / sum(h3k4me1_probs$post_prob, na.rm = TRUE)
cutoff <- h3k4me1_probs$post_prob[temp[which(frac > 0.8)[1]]]
h3k4me1_probs <- h3k4me1_probs[h3k4me1_probs$post_prob >= cutoff, c("hit_ID", "tag_ID", "post_prob")]
## export
save(h3k4me1_probs,
     file = paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/h3k4me1-pgs-snps.RData"))


## H3K4me3 annotations
## keep only posterior probabilities for SNPs in H3K4me3 marks
h3k4me3_probs <- data.frame(hit_ID = pgs_snps_tags$hit_ID,
                            tag_ID = pgs_snps_tags$tag_ID,
                            h3k4me3_annot = h3k4me3_annot[[h3k4me3_ds_best]],
                            post_prob = newpost$post_prob)
h3k4me3_probs <- h3k4me3_probs[h3k4me3_probs$h3k4me3_annot == 1 &
                               !is.na(h3k4me3_probs$post_prob),]  # also discard repeated tag SNPs
## keep highest PP SNPs accounting for 80% of posterior density in total
temp <- order(h3k4me3_probs$post_prob, decreasing = TRUE)
frac <- cumsum(h3k4me3_probs$post_prob[temp]) / sum(h3k4me3_probs$post_prob, na.rm = TRUE)
cutoff <- h3k4me3_probs$post_prob[temp[which(frac > 0.8)[1]]]
h3k4me3_probs <- h3k4me3_probs[h3k4me3_probs$post_prob >= cutoff, c("hit_ID", "tag_ID", "post_prob")]
## export
save(h3k4me3_probs,
     file = paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/h3k4me3-pgs-snps.RData"))
