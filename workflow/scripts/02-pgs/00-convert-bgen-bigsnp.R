## Convert BGEN files to bigSNP format

library(argparser)
library(data.table)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)


p <- arg_parser("Argument parser")
p <- add_argument(p, "--chr", help = "Chromosome", nargs = 1)
p <- add_argument(p, '--threads', help = 'Number of CPU cores', nargs = 1)
argv <- parse_args(p)

chr  <- argv$chr
threads  <- as.numeric(argv$threads)



time <- system.time({
    
## load list of filtered SNPs
snps <- fread(paste0("../data/variant-ids/chr", chr, "-var-ids.tab"),
              col.names = "ID", header = FALSE, data.table = FALSE)

## load file with minor allele/A1 and MAF
maf_ref_chr <- fread(paste0("../data/imputed-wb_all-stats/minor-alleles/chr", chr, ".tab"), data.table = FALSE)

## ensure original order of SNPs is maintained
## (ordering by position is not enough because of multiallelic SNPs)
snps$order <- rows_along(snps)
snps <- merge(snps, maf_ref_chr[, c("ID", "POS", "REF", "ALT", "A1")], by = "ID")
snps <- snps[order(snps$order),]
snps <- snps[, -2]


## load sample lists
ukb_sample <- fread("../data/sample-ids/ukb22828_c1_b0_v3_s487256.sample", data.table = FALSE)
ukb_sample <- ukb_sample[-1, 1:2]
colnames(ukb_sample) <- c("#FID", "IID")



## load genotypes for all samples
bgen = paste0("/path/to/data/imputed-genotypes/bgen/ukb_imp_chr", chr, "_v3.bgen")
variant_ids = gsub(":", "_", snps$ID)

geno_rds_chr <- snp_readBGEN(bgen,
                             backingfile = paste0("../data/imputed-genotypes/bigsnp/chr", chr),
                             list_snp_id = list(variant_ids),
                             read_as = "dosage",
                             ncores = threads - 1)
geno <- snp_attach(geno_rds_chr)

## swap ref/alt counts where needed to match GWAS
ref_snps <- (snps$A1 == snps$REF)
ref_snps_cols <- seq(1, ncol(geno$genotypes))[ref_snps]
ref_snps_n <- length(ref_snps_cols)

if (ref_snps_n > 5000) {  # split into batches of 5000 to avoid crashes
    if ((ref_snps_n %% 5000) != 0) {
        ## Seq of intervals into which to slice data for faster processing
        st <- seq(1, ref_snps_n, 5000)
        en <- c(seq(5000, ref_snps_n, 5000), ref_snps_n)
    } else {
        st <- seq(1, ref_snps_n, 5000)
        en <- seq(5000, ref_snps_n, 5000)
    }

    geno$genotypes[, ref_snps_cols[st[1]:en[1]]] <- 207 - round(100 * geno$genotypes[, ref_snps_cols[st[1]:en[1]]])
    if (length(st) > 1) {
        for (t in 2:length(st)) {
            geno$genotypes[, ref_snps_cols[st[t]:en[t]]] <- 207 - round(100 * geno$genotypes[, ref_snps_cols[st[t]:en[t]]])
        }
    }
} else {
    geno$genotypes[, ref_snps_cols] <- 207 - round(100 * geno$genotypes[, ref_snps_cols])
}



## add fam table and transform map table from tibble (tbl_df) to data.frame
## (these changes seems to be required by snp_subset function)
geno$fam <- ukb_sample
geno$map <- as.data.frame(geno$map)
## add the usual variant IDs
geno$map$ID <- paste0(chr, ":", geno$map$physical.pos, ":",
                      geno$map$allele1, ":", geno$map$allele2)
## add A1 column
geno$map$A1 <- snps$A1

## check that order of SNPs remains correct
if (sum(geno$map$ID != snps$ID) > 0) {
    stop("SNPs in bigSNP object are in different order than list of filtered SNPs.")
}

## overwrite RDS file to save changes to fam and map dfs
saveRDS(geno, paste0("../data/imputed-genotypes/bigsnp/chr", chr, ".rds"))


})
fwrite(list(time[3]), file = paste0("../data/imputed-genotypes/bigsnp/logs/chr", chr, "-runtime.txt"),
       col.names = FALSE, quote = FALSE)
