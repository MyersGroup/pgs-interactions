## Prepare annotations on overlap with TF binding motifs

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

source("scripts/misc/fn-test_enrich_mat.R")
source("scripts/misc/fn-get_post.R")
source("scripts/misc/fn-score.R")
source("scripts/misc/fn-score_all.R")
source("scripts/misc/fn-sample_motif.R")



## load data.frame with PGS hits and their tags and BFs
load(paste0("../results/04-tf-binding/pgs-snps-tags/", phen, "/pgs-snps-tags.RData"))

## load processed enhancer/promoter annotations
load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/",
            "annot/h3k4me1_3-annot.RData"))

## load VEP annotations
load(paste0("../results/04-tf-binding/coding-epigenetics/", phen, "/",
            "annot/functional-annot.RData"))


## get best alpha values for coding/epigenetic annotations as before
qqannot <- test_enrich_mat(pgs_snps_tags$post_prob,
                           cbind(functional_annot, h3k4me1_annot[[h3k4me1_ds_best]], h3k4me3_annot[[h3k4me3_ds_best]]),
                           pgs_snps_tags)


## load HOCOMOCO data
hocomoco_raw <- scan("../data/tf-binding/hocomoco/v13/H13CORE_jaspar_format.txt", what = "char")
hocomoco <- vector("list", length(grep(">", hocomoco_raw)))
starts <- grep(">", hocomoco_raw)
starts <- c(starts, length(hocomoco_raw) + 1)
for (i in 1:length(hocomoco)) {
    start <- starts[i] + 1
    end   <- starts[i + 1] - 1
    hocomoco[[i]]$name  <- hocomoco_raw[starts[i]]
    hocomoco[[i]]$nameh <- gsub(">", "", hocomoco_raw[starts[i]])
    hocomoco[[i]]$motif <- t(matrix(as.double(hocomoco_raw[start:end]), nrow = 4, byrow = TRUE))
}
namesho <- unlist(lapply(1:length(hocomoco), function(i) hocomoco[[i]]$nameh))
names(hocomoco) <- namesho




## for each Hocomoco motif, take sample of bases from each position 10,000 times
## according to PWM and make score, then compute quantiles (100) of these scores
set.seed(58201)
motsquantile <- matrix(nrow = length(hocomoco), ncol = 101)
for (i in 1:length(hocomoco)) {
    a <- sample_motif(namesho[i], n = 10000)
    motsquantile[i,] <- quantile(a, seq(0, 1, .01))
}
motsquantile[, 1] <- 0




## STEP 1: obtain sample set of tag SNPs falling in H3K4me1 marks with the same
## frequency spectrum and total number of tags as the overall set of tags

## load Relate results for GBR population
load("../data/1000G/relate/allele_ages_EUR/allele_ages_GBR.RData")
allele_ages_gbr <- as.data.frame(allele_ages); rm(allele_ages)
## extract rsIDs
allele_ages_gbr$REF <- sapply(strsplit(allele_ages_gbr$ID, split = ":"), "[[", 3)
allele_ages_gbr$ALT <- sapply(strsplit(allele_ages_gbr$ID, split = ":"), "[[", 4)
allele_ages_gbr$rsID <- sapply(strsplit(allele_ages_gbr$ID, split = ":"), "[[", 1)
## make ID of form CHR:POS_REF_ALT when no rsID is available
cols_no_rsid <- (!startsWith(allele_ages_gbr$rsID, "rs"))
allele_ages_gbr$rsID[cols_no_rsid] <- paste0(
    allele_ages_gbr$CHR[cols_no_rsid], ":",
    allele_ages_gbr$BP[cols_no_rsid],  "_",
    allele_ages_gbr$REF[cols_no_rsid], "_",
    allele_ages_gbr$ALT[cols_no_rsid])
## compute derived frequency spectrum but using only SNPs found in our list of tags
allele_ages_gbr_in_pgs <- allele_ages_gbr[allele_ages_gbr$rsID %in% pgs_snps_tags$tag_rsID,]
target <- table(allele_ages_gbr_in_pgs[, "DAF"])

## load best H3K4me1 dataset
h3k4me1_bed <- fread(paste0("../data/epigenome/bed/", h3k4me1_ds_best, ".bed"), data.table = FALSE)

## find overlaps between Relate and those in H3K4me1 dataset
## prepare Relate results
df <- data.frame(chrom  = paste0("chr", allele_ages_gbr$CHR),
                 start  = allele_ages_gbr$BP,
                 end    = allele_ages_gbr$BP,
                 strand = rep("+", nrow(allele_ages_gbr)))
HR <- GenomicRanges::makeGRangesFromDataFrame(df); rm(df)
## prepare H3K4me1 data
df <- data.frame(chrom  = h3k4me1_bed[, 1],
                 start  = h3k4me1_bed[, 2],
                 end    = h3k4me1_bed[, 3],
                 strand = rep("+", nrow(h3k4me1_bed)))
GR <- GenomicRanges::makeGRangesFromDataFrame(df); rm(df)
chain <- rtracklayer::import.chain("../data/liftOver/hg38ToHg19.over.chain")  # convert from hg38 to hg19
GRhg19 <- rtracklayer::liftOver(GR, chain)
GRhg19 <- unlist(GRhg19)
## find overlaps
vv <- GenomicRanges::findOverlaps(GRhg19, HR)
checkmat <- cbind(from(vv), to(vv))
h3k4me1overlap <- rep(0, nrow(allele_ages_gbr))
h3k4me1overlap[checkmat[, 2]] <- checkmat[, 1]  # row corresponds to Relate df, value is row in H3K4me1 data

## sample from SNPs in Relate results overlapping H3K4me1 marks to obtain set of SNPs
## with approximately the same derived allele frequency spectrum as the set of tag SNP
set.seed(48201)
tosamplefrom <- allele_ages_gbr[h3k4me1overlap > 0,]  # keep only Relate SNPs that fall in an H3K4me1 interval
cols <- match(tosamplefrom[, "DAF"], names(target))
ourrows <- vector(length = 0)
for (i in 1:length(target)) {
    ss <- which(cols == i)
    m <- target[i]
    ## if there are not enough SNPs to sample from, take all that are available
    if (m > length(ss)) {
        m <- length(ss)
    }
    ourrows <- c(ourrows, sample(which(cols == i), m))
}
newsamples <- tosamplefrom[sort(ourrows),]
h3k4me1samples <- newsamples

## find overlaps between SNPs in Relate (1000GP) results and those in H3K4me3 dataset
h3k4me3_bed <- fread(paste0("../data/epigenome/bed/", h3k4me3_ds_best, ".bed"), data.table = FALSE)
## prepare H3K4me3 data
df <- data.frame(chrom  = h3k4me3_bed[, 1],
                 start  = h3k4me3_bed[, 2],
                 end    = h3k4me3_bed[, 3],
                 strand = rep("+", nrow(h3k4me3_bed)))
GR <- GenomicRanges::makeGRangesFromDataFrame(df); rm(df)
GRhg19 <- rtracklayer::liftOver(GR, chain)
GRhg19 <- unlist(GRhg19)
## find overlaps
vv <- GenomicRanges::findOverlaps(GRhg19, HR)
checkmat <- cbind(from(vv), to(vv))
h3k4me3overlap <- rep(0, nrow(allele_ages_gbr))
h3k4me3overlap[checkmat[, 2]] <- checkmat[, 1]

## sample from SNPs in Relate results *not* overlapping either methylation mark to obtain set
## of SNPs with approximately the same derived allele frequency spectrum as the set of tag SNP
set.seed(732028)
tosamplefrom <- allele_ages_gbr[h3k4me1overlap == 0 & h3k4me3overlap == 0 ,]
cols <- match(tosamplefrom[, "DAF"], names(target))
ourrows <- vector(length = 0)
for (i in 1:length(target)) {
    ss <- which(cols == i)
    m <- target[i]
    if (m > length(ss)) {
        m <- length(ss)
    }
    ourrows <- c(ourrows, sample(which(cols == i), m))
}
newsamples <- tosamplefrom[sort(ourrows),]
neithersamples <- newsamples




## STEP 2: for each SNP in set of 'neither SNPs', obtain 150bp surrounding sequence
## and score it for how closely it is bound by each TF

library(BSgenome.Hsapiens.UCSC.hg19)

## SNPs not overlapping either methylation mark (neithersamples)
## get sequences +-150bp around each SNP
curmat <- as.data.frame(neithersamples)
mm <- vector(length = 0)
for(chrom in 1:22){
    starts <- as.double(curmat[as.double(curmat[, "CHR"]) == chrom, "BP"]) - 75
    ends   <- as.double(curmat[as.double(curmat[, "CHR"]) == chrom, "BP"]) + 75
    r1 <- IRanges::IRanges(start = starts, end = ends)

    chrcur <- Biostrings::unmasked(BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chrom]])
    myseqs <- IRanges::Views(chrcur, start = start(r1), end = end(r1))

    tempset <- as.character(myseqs)
    mm <- c(mm, tempset)
}

## transform vector of strings into matrix with one base per column
refseqs <- mm
refsurrounds <- matrix(nrow = length(mm), ncol = 151)
for (i in 1:nchar(mm[1])) {
    refsurrounds[, i] <- substring(refseqs, i, i)
}

## convert bases to integers and score
sequencemat <- refsurrounds
sequencemat[sequencemat == "A"] <- 1
sequencemat[sequencemat == "C"] <- 2
sequencemat[sequencemat == "G"] <- 3
sequencemat[sequencemat == "T"] <- 4
sequencemat[sequencemat == "N"] <- 5
sequencemat <- matrix(as.double(sequencemat), nrow = nrow(sequencemat))
maxscoremat <- score_all(sequencemat = sequencemat, hocomoco = hocomoco, namesho = namesho)

## for each motif, compute percentile where each sequence falls based
## on distribution of scores for that motif obtained by sampling
## 10,000 random sequences above
qrefscoremath <- maxscoremat
y <- seq(0, 1, .01)
for (i in 1:ncol(qrefscoremath)) {
	qrefscoremath[, i] <- approx(x = motsquantile[i,], y = y, xout = qrefscoremath[, i], ties = mean)$y	
}
qrefscoremath[is.na(qrefscoremath)] <- 1
colnames(qrefscoremath) <- namesho

maxscoremat_neither <- maxscoremat; rm(maxscoremat)
qrefscoremath_neither <- qrefscoremath; rm(qrefscoremath)


## repeat for sample of SNPs overlapping H3K4me1 marks
curmat <- h3k4me1samples
mm <- vector(length = 0)
colnames(curmat) <- c("CHR","BP")
for (chrom in 1:22) {
    starts <- as.double(curmat[as.double(curmat[, "CHR"]) == chrom, "BP"]) - 75
    ends   <- as.double(curmat[as.double(curmat[, "CHR"]) == chrom, "BP"]) + 75
    r1 <- IRanges::IRanges(start = starts, end = ends)

    chrcur <- Biostrings::unmasked(Hsapiens[[chrom]])
    myseqs <- IRanges::Views(chrcur, start = start(r1), end = end(r1))

    tempset <- as.character(myseqs)
    mm <- c(mm,tempset)
}

refseqs <- mm
refsurrounds <- matrix(nrow = length(mm), ncol = 151)
for (i in 1:nchar(mm[1])) {
    refsurrounds[, i] <- substring(refseqs, i, i)
}

sequencemat <- refsurrounds
sequencemat[sequencemat=="A"] <- 1
sequencemat[sequencemat=="C"] <- 2
sequencemat[sequencemat=="G"] <- 3
sequencemat[sequencemat=="T"] <- 4
sequencemat[sequencemat=="N"] <- 5
sequencemat <- matrix(as.double(sequencemat), nrow = nrow(sequencemat))
maxscoremat <- score_all(sequencemat = sequencemat, hocomoco = hocomoco, namesho = namesho)

qrefscoremath <- maxscoremat
y <- seq(0, 1, .01)
for (i in 1:ncol(qrefscoremath)) {
    qrefscoremath[, i] <- approx(x = motsquantile[i,], y = y, xout = qrefscoremath[, i], ties = mean)$y	
}
qrefscoremath[is.na(qrefscoremath)] <- 1
colnames(qrefscoremath) <- namesho

maxscoremat_h3k4me1 <- maxscoremat; rm(maxscoremat)
qrefscoremath_h3k4me1 <- qrefscoremath; rm(qrefscoremath)


## for each motif, compute odds ratios of having SNPs in 9 top percentiles of
## score for 'neither' and H3K4me1 sets separately;
## then compute ratio of odds ratios
nmat <- qrefscoremath_neither
hmat <- qrefscoremath_h3k4me1

overalltempres <- matrix(nrow = 0, ncol = length(namesho))
qmaxmat <- qrefscoremath_h3k4me1

fracset <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25)

for (frac in fracset[1:length(fracset)]) {
    qmaxmat2 <- qmaxmat * 0
    for (i in 1:length(namesho)) {
        thresh <- quantile(qmaxmat[, i], 1 - frac, na.rm = TRUE)
        qmaxmat2[qmaxmat[, i] > thresh, i] <- 1
    }

    v <- matrix(0L, nrow = 4, ncol = ncol(nmat))
    for (i in 1:length(namesho)) {
        thresh <- quantile(nmat[, i], 1 - frac, na.rm = TRUE)
        v[2, i] <- sum(nmat[, i] > thresh)
        v[4, i] <- sum(hmat[, i] > thresh)  # notice that the same threshold is used
    }
    v[1,] <- nrow(nmat) - v[2,]
    v[3,] <- nrow(hmat) - v[4,]

    oddsratio <- v[1,] * v[4,] / v[2,] / v[3,]
    overalltempres <- rbind(overalltempres, oddsratio)
}
overallmath3k4me1random_150 <- overalltempres

## replace Inf ratios with NaN so they are not considered for the interpolation below
overallmath3k4me1random_150[is.infinite(overallmath3k4me1random_150)] <- NaN

## make matrix of 9 (high) quantiles (fracset) for each motif
nmat <- qrefscoremath_neither
quantilemat <- matrix(nrow = length(fracset), ncol = length(namesho))
for (i in 1:length(namesho)) {
    quantilemat[, i] <- quantile(nmat[, i], 1 - fracset, na.rm = TRUE)
}
## for each motif, compute average of scores above threshold
meanquantilemat <- quantilemat * 0
temp <- t(nmat)
for (i in 1:nrow(quantilemat)) {
    temp2 <- temp * as.double(temp > quantilemat[i,])  # keep only temp elements > corresponding quantile, the others become 0
    meanquantilemat[i,] <- rowSums(temp2) / rowSums(temp2 != 0)
}
## retrieve odds ratios only from overallmath3k4me1random_150 above
oddsmat <- rbind(overallmath3k4me1random_150, 1)  # add bottom rows of 1s
colnames(oddsmat) <- namesho
meanquantilemat <- rbind(meanquantilemat, colMeans(nmat))  # add bottom rows with average score across all tag SNPs for each motif
quantilemat <- rbind(quantilemat, 0)

save(oddsmat, meanquantilemat, quantilemat,
     file = paste0("../results/04-tf-binding/tf-binding/", phen, "/neithermaxscorequantilematq150.RData"))




## STEP 3: score the tag SNPs for overlap with each TF
detach("package:BSgenome.Hsapiens.UCSC.hg19", unload = TRUE)
library(BSgenome.Hsapiens.UCSC.hg38)

## replace REF/ALT alleles with ancestral/derived from Relate if available
ref <- pgs_snps_tags[, "tag_REF"]
alt <- pgs_snps_tags[, "tag_ALT"]
ref[!is.na(pgs_snps_tags[, "tag_anc"])] <- pgs_snps_tags[!is.na(pgs_snps_tags[, "tag_anc"]), "tag_anc"]
alt[!is.na(pgs_snps_tags[, "tag_anc"])] <- pgs_snps_tags[!is.na(pgs_snps_tags[, "tag_anc"]), "tag_der"]

## extract surround sequences (30bp on either side)
mm <- vector(length = 0)
for(chr in 1:22){
    ## get start/end positions 30bp from each tag (extended on RHS if deletion, i.e. REF has more than one base)
    r1 <- IRanges(start = as.double(pgs_snps_tags[as.double(pgs_snps_tags[, "tag_CHR"]) == chr, "tag_POS_hg38"]) - 30,
                  end   = as.double(pgs_snps_tags[as.double(pgs_snps_tags[, "tag_CHR"]) == chr, "tag_POS_hg38"]) + 29 + nchar(ref[as.double(pgs_snps_tags[, "tag_CHR"]) == chr]))

    ## get sequences from reference genome
    chrcur <- Biostrings::unmasked(BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr]])
    myseqs <- IRanges::Views(chrcur, start = start(r1), end = end(r1))

    ## convert to character and append to previous chromosome
    tempset <- as.character(myseqs)
    mm <- c(mm, tempset)
}

## replace REF/ALT sequences from reference genome with original REF/ALT
allseqs <- mm
refseqs <- allseqs
refseqs <- paste0(substring(refseqs, 1, 30), ref, substring(refseqs, 31 + nchar(ref), nchar(refseqs)))
altseqs <- allseqs
altseqs <- paste0(substring(altseqs, 1, 30), alt, substring(altseqs, 31 + nchar(ref), nchar(altseqs)))

## keep only first 61bp of sequences
refseqs <- substring(refseqs, 1, 61)
altseqs <- substring(altseqs, 1, 61)

## make matrices with one column per base
sequencemat <- matrix(nrow = length(mm), ncol = 61)
for (i in 1:61) {
    sequencemat[, i] <- substring(refseqs, i, i)
}
refsurrounds <- sequencemat

sequencemat <- matrix(nrow = length(mm), ncol = 61)
for (i in 1:61) {
    sequencemat[, i] <- substring(altseqs, i, i)
}
altsurrounds <- sequencemat

## score REF allele
sequencemat <- refsurrounds
sequencemat[sequencemat == "A"] <- 1
sequencemat[sequencemat == "C"] <- 2
sequencemat[sequencemat == "G"] <- 3
sequencemat[sequencemat == "T"] <- 4
sequencemat[sequencemat == "N"] <- 5
sequencemat <- matrix(as.double(sequencemat), nrow = nrow(sequencemat))
transformedseq <- sequencemat
transformedseq[is.na(transformedseq)] <- 5

## score ALT allele
sequencemat <- altsurrounds
sequencemat[sequencemat == "A"] <- 1
sequencemat[sequencemat == "C"] <- 2
sequencemat[sequencemat == "G"] <- 3
sequencemat[sequencemat == "T"] <- 4
sequencemat[sequencemat == "N"] <- 5
sequencemat <- matrix(as.double(sequencemat), nrow = nrow(sequencemat))
transformedseqalt <- sequencemat
transformedseqalt[is.na(transformedseqalt)] <- 5

## score motifs
refscoremath <- matrix(nrow = nrow(transformedseq), ncol = length(hocomoco))
altscoremath <- refscoremath

for (mot in 1:length(hocomoco)) {
   
    ## score REF and ALT, standard and complement sequences
    scoremat <- score(transformedseq, namesho[mot], transformed = TRUE)
    scorematalt <- score(transformedseqalt, namesho[mot], transformed = TRUE)
    scorematcomp <- score(transformedseq, namesho[mot], transformed = TRUE, comp = TRUE)
    scorematcompalt <- score(transformedseqalt, namesho[mot], transformed = TRUE, comp = TRUE)

    ## get highest score between standard and complement seqs. for REF and ALT separately
    bestscore <- scoremat[, 1]
    bestalt <- scorematalt[, 1]
    for (i in 2:ncol(scoremat)) {
	bestscore[scoremat[, i] > bestscore] <- scoremat[scoremat[, i] > bestscore, i]
	bestalt[scorematalt[, i] > bestalt]  <- scorematalt[scorematalt[, i] > bestalt, i]
    }
    for (i in 1:ncol(scoremat)) {
	bestscore[scorematcomp[, i] > bestscore] <- scorematcomp[scorematcomp[, i] > bestscore, i]
	bestalt[scorematcompalt[, i] > bestalt]  <- scorematcompalt[scorematcompalt[, i] > bestalt, i]
    }
    refscoremath[, mot] <- bestscore
    altscoremath[, mot] <- bestalt
}
rownames(refscoremath) <- pgs_snps_tags$tag_ID
colnames(refscoremath) <- namesho
rownames(altscoremath) <- pgs_snps_tags$tag_ID
colnames(altscoremath) <- namesho

qrefscoremath <- refscoremath
qaltscoremath <- altscoremath
y <- seq(0, 1, .01)
for (i in 1:ncol(qrefscoremath)) {
	qrefscoremath[, i] <- approx(x = motsquantile[i,], y = y, xout = qrefscoremath[, i], ties = mean)$y	
	qaltscoremath[, i] <- approx(x = motsquantile[i,], y = y, xout = qaltscoremath[, i], ties = mean)$y
}
qrefscoremath[is.na(qrefscoremath)] <- 1
qaltscoremath[is.na(qaltscoremath)] <- 1
rownames(qrefscoremath) <- pgs_snps_tags$tag_ID
colnames(qrefscoremath) <- namesho
rownames(qaltscoremath) <- pgs_snps_tags$tag_ID
colnames(qaltscoremath) <- namesho

save(qrefscoremath, qaltscoremath, ref, alt,
     file = paste0("../results/04-tf-binding/tf-binding/", phen, "/motif-scores.RData"))

## test for motifs enriched in higher weight relative to lower weight cases
oddsmat <- rbind(oddsmat, 1)
meanquantilemat <- rbind(1, meanquantilemat)
meanquantilemat <- rbind(meanquantilemat, 0)
oddsmat <- rbind(oddsmat[1,], oddsmat)

qmaxscoremat <- qrefscoremath
qmaxscoremat[qrefscoremath < qaltscoremath] <- qaltscoremath[qrefscoremath < qaltscoremath]

hitor <- qmaxscoremat
meanquantilemat[is.na(meanquantilemat)] <- 1
for(i in 1:length(namesho)){
    hitor[, i] <- approx(x = meanquantilemat[, i], xout = qmaxscoremat[, i], y = oddsmat[, i], ties = mean)$y
}




## STEP 4: fit model by MLE to obtain alphas
qq  <- test_enrich_mat(pgs_snps_tags$post_prob, hitor, pgs_snps_tags[, c("hit_ID", "tag_ID")],
                       mode = "R", alpharange = seq(0.1, 4, 0.05))
qq2 <- test_enrich_mat(pgs_snps_tags$post_prob, hitor, pgs_snps_tags[, c("hit_ID", "tag_ID")],
                       mode = "R", alpharange = c(1, seq(4.25, 6, 0.25)))

qqnew <- qq
qqnew$lhoodsu <- cbind(qq$lhoodsu, qq2$lhoodsu[, 2:ncol(qq2$lhoodsu)])
qqnew$alphas <- c(qq$alphas, qq2$alphas[2:ncol(qq2$lhoodsu)])
cond <- (qq2$psu < qq$psu)
cond[is.na(cond)] <- FALSE
qqnew$psu[cond] <- qq2$psu[cond]
qqnew$oru[cond] <- qq2$oru[cond]

tf_probs_ls <- list()
for (k in 1:length(namesho)){

    ## if the inferred alpha for this motif is <= 1, skip it
    if (qqnew$oru[k] <= 1) {
        tf_probs_ls[[namesho[k]]] <- "Unavailable"
        next
    }

    ## STEP 5: build combined annotation for each motif
    extra <- (hitor[, k] - 1) * (qqnew$oru[k] - 1)
    annot_comb <- cbind(1, functional_annot * (qqannot$oru[1] - 1) +
                           h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1) +
                           h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1) +
                           functional_annot * (qqannot$oru[1] - 1) * h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1) +
                           functional_annot * (qqannot$oru[1] - 1) * h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1) +
                           h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1) * h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1) +
                           extra +
                           h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1) * extra +
                           h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1) * extra +
                           functional_annot * (qqannot$oru[1] - 1) * h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1) * h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1) +
                           h3k4me1_annot[[h3k4me1_ds_best]] * (qqannot$oru[2] - 1) * h3k4me3_annot[[h3k4me3_ds_best]] * (qqannot$oru[3] - 1) * extra)
    

    ## STEP 6: find posterior assuming we combine all the alpha values partially
    ## multiplicatively as above to weight SNPs
    curpost <- get_post(pgs_snps_tags$post_prob,
                        annot_comb,
                        pgs_snps_tags[, c("hit_ID", "tag_ID")],
                        alphavalues = c(1, 2))


    ## STEP 7: compute final weights: weight by quality of motif match
    ## (best one) and posterior probability of being causal
    tf_probs <- data.frame(hit_ID = pgs_snps_tags$hit_ID,
                           tag_ID = pgs_snps_tags$tag_ID,
                           extra_annot = extra,
                           post_prob = curpost$post_prob)
    tf_probs$post_prob <- (tf_probs$extra / (1 + tf_probs$extra)) * tf_probs$post_prob
    ## remove cases where probability is negative
    tf_probs$post_prob[tf_probs$post_prob < 0] <- 0
    ## also discard repeated tag SNPs
    tf_probs <- tf_probs[!is.na(tf_probs$post_prob),]

    ## if there are fewer than 2 SNPs with positive probabilities, skip this motif
    if (sum(tf_probs$post_prob > 0) < 2) {
        tf_probs_ls[[namesho[k]]] <- "Unavailable"
        next
    }
    
    ## keep highest PP SNPs accounting for 80% of posterior density in total
    temp <- order(tf_probs$post_prob, decreasing = TRUE)
    frac <- cumsum(tf_probs$post_prob[temp]) / sum(tf_probs$post_prob, na.rm = TRUE)
    cutoff <- tf_probs$post_prob[temp[which(frac > 0.8)[1]]]
    tf_out <- tf_probs[tf_probs$post_prob >= cutoff, c("hit_ID", "tag_ID", "post_prob")]

    ## add to list
    tf_probs_ls[[namesho[k]]] <- tf_out
}

save(tf_probs_ls,
     file = paste0("../results/04-tf-binding/tf-binding/", phen, "/tf-pgs-snps.RData"))
