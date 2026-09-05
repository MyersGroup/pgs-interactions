## desc: Prepare text files with coefficients of LD clump index variants
##       to be fed into PLINK for PGS computation

library(argparser)
library(data.table)
library(dplyr)
library(tidyr)



p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",    help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples", help = "Sample set code", nargs = 1)
p <- add_argument(p, "--p1_max",  help = "LD clumping maximum p1 parameter value", nargs = 1)
p <- add_argument(p, "--r2",      help = "LD clumping r2 parameter value", nargs = 1)
p <- add_argument(p, "--kb",      help = "LD clumping kb parameter value", nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen
samples <- argv$samples
p1_max  <- argv$p1_max
r2      <- argv$r2
kb      <- argv$kb




time <- system.time({

## p1 value grid: 1000 values equally spaced on the log scale between 1e-9 and 0.05
p1s <- exp(seq(log(1e-9), log(0.05), length.out = 1000))


## load clumping data and add BETA and P-VALUE
clump <- fread(paste0("../results/02-pgs/ld-clump/", phen,
                      "/p1_", p1_max, "-r2_", r2, "-kb_", kb,
                      "/", samples, "-chr1-clumped.tab"))
clump$CHR <- 1

for (chr in 2:22) {
    clump_chr <- fread(paste0("../results/02-pgs/ld-clump/", phen,
                              "/p1_", p1_max, "-r2_", r2, "-kb_", kb,
                              "/", samples, "-chr", chr, "-clumped.tab"))
    clump_chr$CHR <- chr
    clump <- rbind(clump, clump_chr)
}


## make score vector table
coeff <- clump
coeff <- coeff[order(-coeff$LOG10_P),]

cols <- rep(0, length(p1s))
for (i in 1:length(p1s)) {
    cols[i] <- max(which(coeff$LOG10_P >= -log10(p1s[i])))
}

n <- nrow(coeff)
for (r in 1:length(cols)) {
    new_col_nm <- paste0("S_", format(p1s[r], scientific = FALSE))
    coeff[, (new_col_nm) := c(BETA[1:cols[r]], rep(0, n - cols[r]))]
}


## make table with p-value threshold and #SNPs for each score
pgs_nsnp <- data.frame(r2 = r2,
                       kb = kb,
                       p1 = format(p1s, scientific = FALSE),
                       n_snps = cols)
fwrite(pgs_nsnp, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", samples,
                               "/coeff/p1_grid/pgs-p1-n_snps.tab"),
       sep = "\t", na = "NA", quote = FALSE)


## export by chromosome
for (chr in 1:22) {
    coeff_chr <- coeff[coeff$CHR == chr, -5]
    coeff_chr <- coeff_chr[order(match(coeff_chr$ID, clump$ID)),]

    fwrite(coeff_chr, file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", samples,
                                    "/coeff/p1_grid/chr", chr, "-coeff.tab"),
           sep = "\t", na = "NA", quote = FALSE)
}

})
fwrite(list(time[3]), file = paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", samples,
                                    "/coeff/p1_grid/make-coeff-vec-runtime.txt"),
       col.names = FALSE, quote = FALSE)
