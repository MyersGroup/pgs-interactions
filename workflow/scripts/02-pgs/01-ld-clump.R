## desc: Run LD clumping on PLINK2 --glm log10 output


library(argparser)
library(data.table)
library(bigsnpr)
## disable parallel BLAS if enabled by default
## see `?assert_cores`
options(bigstatsr.check.parallel.blas = FALSE)


p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",     help = "Phenotype code",           nargs = 1)
p <- add_argument(p, "--chr",      help = "Chromosome",               nargs = 1)
p <- add_argument(p, "--train_sp", help = "Training samples",         nargs = 1)
p <- add_argument(p, '--p1',       help = 'LD clumping p1 parameter', nargs = 1)
p <- add_argument(p, '--r2',       help = 'LD clumping r2 parameter', nargs = 1)
p <- add_argument(p, '--kb',       help = 'LD clumping kb parameter', nargs = 1)
p <- add_argument(p, '--threads',  help = 'Number of CPU cores',      nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen
chr      <- argv$chr
train_sp <- argv$train_sp
p1       <- as.numeric(argv$p1)
r2       <- as.numeric(argv$r2)
kb       <- as.numeric(argv$kb)
threads  <- as.numeric(argv$threads)




time <- system.time({

## load GWAS summary stats
gwas <- fread(paste0("../results/01-gwas/plink-output/", phen, "/", train_sp, ".chr", chr, ".", phen, ".res.cov.glm.linear.gz"),
              data.table = FALSE)

## load imputed data
geno <- snp_attach(paste0("../data/imputed-genotypes/bigsnp/chr", chr, ".rds"))

## confirm SNPs are in same order
if (sum(gwas$ID != geno$map$ID) > 0) {
    stop("SNPs in bigSNP object are in different order than in GWAS results.")
}


## load UKB sample file
ukb_sample <- fread("../data/sample-ids/ukb22828_c1_b0_v3_s487256.sample", data.table = FALSE)
ukb_sample <- ukb_sample[-1, 1:2]
colnames(ukb_sample) <- c("#FID", "IID")

## get training IDs for whom phenotype is not missing
phen_train <- fread(paste0('../results/01-gwas/residuals-covars/', phen,
                           '/', phen, '-', train_sp, '-resid-covars.tab'),
                    data.table = FALSE)
colnames(phen_train)[3] <- 'phen_res'
phen_train <- phen_train[!is.na(phen_train$phen_res),]
train_nonmiss_ids <- phen_train[, c('#FID', 'IID')]
train_sp_ind <- seq(1, nrow(ukb_sample))[ukb_sample$IID %in% train_nonmiss_ids$IID]


## clump
exc_pv <- seq(1, nrow(gwas))[gwas$LOG10 < -log10(p1)]
ind_keep <- snp_clumping(geno$genotypes,
                         ind.row = train_sp_ind,
                         S = gwas$LOG10_P,
                         thr.r2 = r2,
                         size = kb,
                         exclude = exc_pv,
                         ncores = threads - 1,
                         infos.chr = as.numeric(geno$map$chromosome),
                         infos.pos = geno$map$physical.pos)


clump_out <- gwas[ind_keep, c("ID", "A1", "BETA", "LOG10_P")]
fwrite(clump_out, file = paste0("../results/02-pgs/ld-clump/", phen, "/p1_", p1, "-r2_", r2, "-kb_", kb, "/",
                                train_sp, "-chr", chr, "-clumped.tab"),
       sep = '\t', na = 'NA', quote = FALSE)

})
fwrite(list(time[3]), file = paste0("../results/02-pgs/ld-clump/", phen, "/p1_", p1, "-r2_", r2, "-kb_", kb, "/",
                                    train_sp, "-chr", chr, "-runtime.txt"),
       col.names = FALSE, quote = FALSE)
