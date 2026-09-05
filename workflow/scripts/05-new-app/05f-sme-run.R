## desc: Prepare phenotype file and run SME

library(argparser)
library(data.table)
setDTthreads(4)
library(dplyr)
library(smer)

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen", help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--snp",  help = "Target SNP",     nargs = 1)
argv <- parse_args(p)

phen <- argv$phen
snp_ <- argv$snp
snp  <- gsub("_", ":", snp_)

samples <- "wb_all"

plink2 <- "/path/to/plink2"



## make masking file with variants in LOCO PGS

## load bim file
bim <- fread(paste0("../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_all_flt_samples_chr.bim"),
             col.names = c("CHR", "ID", "cM", "POS", "REF", "ALT"), data.table = FALSE)

## load coefficients of trait's iterative PGS
r2 <- 0.9; kb <- 500
coeff <- list()
for (chr in 1:22) {
    coeff_f <- paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", samples,
                      "/coeff/final/", phen, "-coeff-chr", chr, ".tab")
    if (file.exists(coeff_f)) {
        coeff[[chr]] <- fread(coeff_f, data.table = FALSE)
    }
}
coeff <- do.call("rbind", coeff)

## code below is adapted from: https://cran.r-project.org/web/packages/smer/vignettes/tutorial-create-mask-file.html
## Group names
gxg_h5_group <- "gxg"
ld_h5_group <- "ld"

## Focal SNP (still in 1-based R indexing)
focal_snp <- which(bim$ID == snp)

## Data (still in 1-based R indexing)
include_gxg_snps <- which(bim$ID %in% coeff$ID)
snp_chr <- as.numeric(sapply(strsplit(snp, split = ":"), "[[", 1))
exclude_ld_snps <- which(bim$ID %in% coeff$ID & bim$CHR == snp_chr & bim$ID != snp)

## Dataset names
dataset_name_pattern <- "%s/%s"
## 0-based index!
gxg_dataset <- sprintf(dataset_name_pattern, gxg_h5_group, focal_snp - 1)
ld_dataset <- sprintf(dataset_name_pattern, ld_h5_group, focal_snp - 1)

## Create an empty HDF5 file
hdf5_file <- paste0("../results/05-new-app/sme/output/", phen, "/", phen, "-", snp_, ".mask")
create_hdf5_file(hdf5_file)

## Write LD data
write_hdf5_dataset(hdf5_file, ld_dataset, exclude_ld_snps - 1)  # 0-based index!

## Write GXG data
write_hdf5_dataset(hdf5_file, gxg_dataset, include_gxg_snps - 1)





## run SME
sme_result <- sme(
    plink_file  = paste0("../results/05-new-app/sme/genotype-calls/", phen, "-", snp_, "-ukb_all_flt_samples_chr"),
    pheno_file  = paste0("../results/05-new-app/sme/output/", phen, "/", phen, "-", snp_, ".pheno"),
    mask_file   = hdf5_file,
    gxg_indices = focal_snp,
    chunk_size  = 10,
    n_randvecs  = 10,
    n_blocks    = 100,
    n_threads   = 4,
    rand_seed   = 982102
)

save(sme_result,
     file = paste0("../results/05-new-app/sme/output/", phen, "/", phen, "-", snp_, "-sme-result.RData"))
