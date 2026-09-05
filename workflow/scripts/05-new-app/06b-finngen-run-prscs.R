library(argparser)
library(data.table)
setDTthreads(4)
library(dplyr)

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen", help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--omopid", help = "FinnGen phenotype code", nargs = 1)
argv <- parse_args(p)

phen <- argv$phen
omopid <- argv$omopid

samples <- "wb_all"
ref_dir <- "/path/to/PRScs/ldblk_1kg_eur"
plink2  <- "/path/to/plink2"


## find sample size
lab_docs <- fread("../data/sum-stats/finngen/Kanta_labs_GWAS_results_v2_summary.txt", data.table = FALSE)
sample_sz <- lab_docs$N_total[lab_docs$OMOPID == omopid]

## prepare summary statistics file
gwas <- fread(paste0("../data/sum-stats/finngen/finngen_R12_", omopid, ".gz"), data.table = FALSE)

## add hg38 ID and subset to variants in Finn Gen & UKB
hg38_hg19_ukb <- fread("../data/sum-stats/finngen/variant-ids-hg38-hg19-ukb.tab.gz", data.table = FALSE)
gwas$ID_hg38 <- paste(gwas[,1], gwas[,2], gwas[,3], gwas[,4], sep = ":")
gwas <- inner_join(gwas, hg38_hg19_ukb, by = "ID_hg38")

## add rsID from 1KG dataset
snpinfo_1kg_fg_ukb <- fread("../results/05-new-app/finngen/bim/fg-ukb-1kg.bim",
                            col.names = c("CHR", "SNP", "cM", "BP", "A1", "A2"), data.table = FALSE)
snpinfo_1kg_fg_ukb$ID_hg19 <- paste(snpinfo_1kg_fg_ukb$CHR, snpinfo_1kg_fg_ukb$BP, snpinfo_1kg_fg_ukb$A2, snpinfo_1kg_fg_ukb$A1, sep = ":")
gwas <- inner_join(gwas, snpinfo_1kg_fg_ukb[, c("SNP", "ID_hg19")], by = "ID_hg19")

## export
gwas_out <- gwas[, c("SNP", "alt", "ref", "beta", "sebeta")]
colnames(gwas_out) <- c("SNP", "A1", "A2", "BETA", "SE")
fwrite(gwas_out,
       file = paste0("../results/05-new-app/finngen/sumstats/", phen, "-", omopid, "-sumstats.tab"),
       sep = "\t", na = "NA", quote = FALSE)



## run PRS-CS
system("export MKL_NUM_THREADS=4")
system("export export NUMEXPR_NUM_THREADS=4")
system("export OMP_NUM_THREADS=4")
system(paste0("python /path/to/PRScs/PRScs.py ",
              "--ref_dir=", ref_dir, " ",
              "--bim_prefix=/path/to/results/05-new-app/finngen/bim/fg-ukb-1kg ",
              "--sst_file=/path/to/results/05-new-app/finngen/sumstats/", phen, "-", omopid, "-sumstats.tab ",
              "--n_gwas=", sample_sz, " ",
              "--out_dir=../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, " ",
              "--seed 674839"))



## compute PGS

## load list of all samples
all_ids <- fread("../data/sample-ids/filtered/all-ids.tab", data.table = FALSE)

## initialise data.frame with score for each chromosome
pgs_bychr <- data.frame(IID = all_ids$IID)

for (chr in 1:22) {

    ## check if chromosome is not empty
    if (file.exists(paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt"))) {

        ## load coeff
        coeff <- fread(paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt"),
                       col.names = c("CHR", "rsID", "POS", "ALT", "REF", "BETA"), data.table = FALSE)

        ## prepare file for scoring with Plink
        coeff$ID <- paste(coeff$CHR, coeff$POS, coeff$REF, coeff$ALT, sep = ":")
        coeff$A1 <- coeff$ALT
        fwrite(coeff[, c("ID", "A1", "BETA")],
               file = paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "_pst_eff_a1_b0.5_phiauto_chr", chr, "-coeff.txt"),
               sep = "\t", na = "NA", quote = FALSE)
              
        ## compute score with Plink
        system(paste0(plink2,
                      " --pgen ../data/imputed-genotypes/ukb_imp_chr", chr, "_v3.pgen ",
                      "--pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr", chr, "_v3.pvar ",
                      "--psam ../data/sample-ids/ukb-103076-imp-auto-s486989.psam ",
                      "--keep ../data/sample-ids/filtered/all-ids.tab ",
                      ## "--read-freq ../data/imputed-wb_all-stats/allele-freq/ukb_imp_chr", chr, "_v3.afreq ",
                      "--score ../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "_pst_eff_a1_b0.5_phiauto_chr", chr, "-coeff.txt header-read ",
                      "--threads 4 ",
                      "--memory 60000 ",
                      "--out ../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-pgs-chr", chr))

        ## load score and delete temporary files
        pgs_chr <- fread(paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-pgs-chr", chr, ".sscore"),
                         data.table = FALSE)
        system(paste0("rm ../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-pgs-chr", chr, ".*"))
        system(paste0("rm ../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "_pst_eff_a1_b0.5_phiauto_chr", chr, "-coeff.txt"))

        ## rescale score by allele count
        pgs_chr[, 5] <- pgs_chr[, 5] * pgs_chr$ALLELE_CT
        
        ## add to main df
        pgs_bychr[[paste0("pgs_chr", chr)]] <- pgs_chr[, 5]
        
    } else {

        ## add to main df
        pgs_bychr[[paste0("pgs_chr", chr)]] <- 0

    }
}

## compute full PGS
pgs_full <- all_ids
pgs_full$pgs_full <- rowSums(pgs_bychr[, paste0("pgs_chr", 1:22)])
pgs_full_out <- pgs_full
pgs_full_out$pgs_full <- round(pgs_full_out$pgs_full, digits = 8)
fwrite(pgs_full_out,
       file = paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-all-pgs-full.tab"),
       sep = "\t", na = "NA", quote = FALSE)

## compute leave-one-chromosome-out PGS
pgs_loco <- all_ids
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
       file = paste0("../results/05-new-app/finngen/PRScs/", phen, "/", phen, "-", omopid, "-all-pgs-loco.tab"),
       sep = "\t", na = "NA", quote = FALSE)
