library(argparser)
library(data.table)

options(warn = 2)


p <- arg_parser('Argument parser')
p <- add_argument(p, '--phen',    help = 'Phenotype code', nargs = 1)
argv <- parse_args(p)

phen    <- argv$phen

plink2 <- "/path/to/plink2"



## load list of all samples
all_sp_ids <- fread("../data/sample-ids/filtered/all-ids.tab", data.table = FALSE)

## initialise data.frame with score for each chromosome
pgs_bychr <- all_sp_ids

## compute score for each chromosome
r2 <- 0.9; kb <- 500; train_sp <- "wb_all"

for (chr in 1:22) {

    ## check if chromosome is not empty
    if (file.exists(paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                           "/coeff/final/", phen, "-coeff-chr", chr, ".tab"))) {

        ## compute score with Plink
        system(paste0(plink2,
                      " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                      "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                      "--psam ../data/sample-ids/ukb-103076-imp-auto-s486989.psam ",
                      "--keep ../data/sample-ids/filtered/all-ids.tab ",
                      "--read-freq ../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                      "--score ../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/coeff/final/", phen, "-coeff-chr", chr, ".tab header-read ",
                      "--threads 2 ",
                      "--memory 30000 ",
                      "--out ../results/05-new-app/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/pgs/pgs-chr", chr))

        ## load score and delete temporary files
        pgs_chr <- fread(paste0("../results/05-new-app/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/pgs/pgs-chr", chr, ".sscore"), data.table = FALSE)
        system(paste0("rm ../results/05-new-app/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/pgs/pgs-chr", chr, ".*"))

        ## rescale score by allele count
        pgs_chr[, 5] <- pgs_chr[, 5] * pgs_chr$ALLELE_CT
        
        ## add to main df
        pgs_bychr[[paste0("pgs_chr", chr)]] <- pgs_chr[, 5]
        
    } else {

        ## add to main df
        pgs_bychr[[paste0("pgs_chr", chr)]] <- 0

    }
}

pgs_bychr_out <- as.data.table(pgs_bychr)
cols <- colnames(pgs_bychr_out)[3:ncol(pgs_bychr_out)]
pgs_bychr_out[, (cols) := round(.SD, digits = 8), .SDcols = cols]
fwrite(pgs_bychr_out,
       file = paste0("../results/05-new-app/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/pgs/", phen, "-all-pgs-bychr.tab"),
       sep = "\t", na = "NA", quote = FALSE)



## compute full PGS
pgs_full <- all_sp_ids
pgs_full$pgs_full <- rowSums(pgs_bychr[, paste0("pgs_chr", 1:22)])
pgs_full_out <- pgs_full
pgs_full_out$pgs_full <- round(pgs_full_out$pgs_full, digits = 8)
fwrite(pgs_full_out,
       file = paste0("../results/05-new-app/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/pgs/", phen, "-all-pgs-full.tab"),
       sep = "\t", na = "NA", quote = FALSE)



## compute leave-one-chromosome-out PGS
pgs_loco <- all_sp_ids
for (chr in 1:22) {
    if (isTRUE(all.equal(pgs_bychr[[paste0("pgs_chr", chr)]], rep(0, nrow(pgs_bychr))))) {
        pgs_loco[[paste0("pgs_loco", chr)]] <- pgs_full$pgs_full
    } else {
        pgs_loco[[paste0("pgs_loco", chr)]] <- rowSums(pgs_bychr[, paste0("pgs_chr", setdiff(1:22, chr))])
    }
}
pgs_loco_out <- as.data.table(pgs_loco)
cols <- colnames(pgs_loco_out)[3:ncol(pgs_loco_out)]
pgs_loco_out[, (cols) := round(.SD, digits = 8), .SDcols = cols]
fwrite(pgs_loco_out,
       file = paste0("../results/05-new-app/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/pgs/", phen, "-all-pgs-loco.tab"),
       sep = "\t", na = "NA", quote = FALSE)
