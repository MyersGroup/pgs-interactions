# desc: Improve initial PGS with stepwise regression

library(argparser)
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tictoc)

source("scripts/misc/fn-ggmanh.R")
source("scripts/misc/fn-ggqq.R")
source("scripts/misc/fn-ggp1.R")


p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",     help = "Phenotype code",           nargs = 1)
p <- add_argument(p, "--train_sp", help = "Training samples",         nargs = 1)
p <- add_argument(p, "--vali_sp",  help = "Validation sample code",   nargs = 1)
p <- add_argument(p, "--test_sp",  help = "Test sample code",         nargs = 1)
p <- add_argument(p, "--r2",       help = "LD clumping r2 parameter", nargs = 1)
p <- add_argument(p, "--kb",       help = "LD clumping kb parameter", nargs = 1)
p <- add_argument(p, '--threads',  help = 'Number of CPU cores',      nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen
train_sp <- argv$train_sp
vali_sp  <- argv$vali_sp
test_sp  <- argv$test_sp
r2       <- argv$r2
kb       <- argv$kb
threads  <- as.numeric(argv$threads)



runtime_df <- data.frame(task = character(),
                         time = numeric())

## STEP 0 #######################################################################

tic()

## regress out PGS and get residuals --------------------------------------------

## phenotype training data
phen_train <- fread(paste0("../results/01-gwas/residuals-covars/", phen,
                           "/", phen, "-", train_sp, "-resid-covars.tab"),
                    data.table = FALSE)
colnames(phen_train)[3] <- "phen_res"
phen_train <- phen_train[!is.na(phen_train$phen_res),]
train_nonmiss_ids <- phen_train[, c("#FID", "IID")]

## phenotype validation data
phen_vali <- fread(paste0("../results/01-gwas/residuals-covars/", phen,
                          "/", phen, "-", vali_sp, "-resid-covars.tab"),
                   data.table = FALSE)
colnames(phen_vali)[3] <- "phen_res"
phen_vali <- phen_vali[!is.na(phen_vali$phen_res),]
vali_nonmiss_ids <- phen_vali[, c("#FID", "IID")]

## PGS
pgs_cum <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                        "/pgs-all/p1_opt/", phen, "-all-pgs0.tab"),
                 data.table = FALSE)
colnames(pgs_cum)[3] <- "pgs0"

## merge
phen_pgs_train <- merge(phen_train, pgs_cum[, c("IID", "pgs0")], by = "IID")
phen_pgs_vali  <- merge(phen_vali, pgs_cum[, c("IID", "pgs0")], by = "IID")

## regress
reg_train <- lm(phen_res ~ pgs0, data = phen_pgs_train)
reg_vali  <- lm(phen_res ~ pgs0, data = phen_pgs_vali)

## save PGS coefficients
pgs_scaling <- data.frame(step = numeric(),
                          coeff_pgs_old_train = numeric(),
                          coeff_pgs_cpn_train = numeric())

resid_train <- cbind(phen_pgs_train[, c("#FID", "IID")], residuals(reg_train))
colnames(resid_train)[3] <- "pgs_res0"
resid_train <- resid_train[order(match(resid_train$IID, train_nonmiss_ids$IID)),]

resid_vali <- cbind(phen_pgs_vali[, c("#FID", "IID")], residuals(reg_vali))
colnames(resid_vali)[3] <- "pgs_res"
resid_vali <- resid_vali[order(match(resid_vali$IID, vali_nonmiss_ids$IID)),]




## get set of original PGS SNPs -------------------------------------------------

## load optimal p1 threshold
p1_opt_0 <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                         "/pgs-all/p1_opt/", phen, "-pgs0-rsq.tab"),
                  data.table = FALSE)
p1_opt <- p1_opt_0$p1_opt


## load list of initial GWAS SNPs and coefficients
coeff_cum    <- list()
coeff        <- list()
chr_pgs_orig <- numeric()
for (chr in 1:22) {
    coeff_cum[[chr]] <- fread(paste0("../results/02-pgs/initial/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                     "/coeff/p1_opt/chr", chr, "-coeff.tab"),
                           data.table = FALSE)
    colnames(coeff_cum[[chr]])[5] <- "step0"
    coeff_cum[[chr]] <- coeff_cum[[chr]][, c("ID", "A1", "step0")]

    ## make df with chrom, pos, ref/alt alleles for use below (CHR + POS needed for Manhattan plots)
    coeff[[chr]] <- coeff_cum[[chr]][, c(1:2)]
    coeff[[chr]] <- coeff[[chr]] %>%
        mutate(CHR = chr,
               POS = as.numeric(word(ID, 2, sep = ":")),
               REF = word(ID, 3, sep = ":")) %>%
        select(CHR, ID, POS, REF, A1)

    if (nrow(coeff[[chr]]) > 0) {
        chr_pgs_orig <- c(chr_pgs_orig, chr)
    }
}




## load genotypes for all samples -----------------------------------------------

geno       <- list()
geno_train <- list()
geno_vali  <- list()
geno_rest  <- list()

## UKB sample file
ukb_sample <- fread("../data/sample-ids/ukb22828_c1_b0_v3_s487256.sample", data.table = FALSE)
ukb_sample <- ukb_sample[-1, 1:2]
colnames(ukb_sample) <- c("#FID", "IID")

all_sp_ids <- fread("../data/sample-ids/filtered/all-ids.tab", data.table = FALSE)
all_sp_ind <- rows_along(ukb_sample)[ukb_sample$IID %in% all_sp_ids$IID]

train_sp_ids <- ukb_sample[ukb_sample$IID %in% train_nonmiss_ids$IID,]
train_sp_ind <- rows_along(ukb_sample)[ukb_sample$IID %in% train_sp_ids$IID]

vali_sp_ids <- ukb_sample[ukb_sample$IID %in% vali_nonmiss_ids$IID,]
vali_sp_ind <- rows_along(ukb_sample)[ukb_sample$IID %in% vali_sp_ids$IID]

rest_sp_ids <- all_sp_ids[!(all_sp_ids$IID %in% c(train_nonmiss_ids$IID, vali_nonmiss_ids$IID)),]
rest_sp_ind <- seq(1, nrow(ukb_sample))[ukb_sample$IID %in% rest_sp_ids$IID]


for (chr in chr_pgs_orig) {
    print(paste0("Loading genotypes: chr ", chr))
    print("Memory used so far:")
    print(object.size(x=lapply(ls(), get)), units="Mb")
    
    bgen = paste0("/path/to/data/imputed-genotypes/bgen/ukb_imp_chr", chr, "_v3.bgen")
    variant_ids = gsub(":", "_", coeff[[chr]]$ID)

    geno_rds_chr <- snp_readBGEN(bgen,
                                 backingfile = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/geno/chr", chr),
                                 list_snp_id = list(variant_ids),
                                 read_as = "dosage",
                                 ncores = threads - 1)
    geno[[chr]] <- snp_attach(geno_rds_chr)

    ## swap ref/alt counts where needed to match GWAS
    ref_snps <- (coeff[[chr]]$A1 == coeff[[chr]]$REF)
    ref_snps_cols <- seq(1, ncol(geno[[chr]]$genotypes))[ref_snps]
    geno[[chr]]$genotypes[, ref_snps_cols] <- 207 - round(100 * geno[[chr]]$genotypes[, ref_snps_cols])

    ## add fam table transform map table from tibble (tbl_df) to data.frame
    ## (these changes seems to be required for snp_subset function below)
    geno[[chr]]$fam <- ukb_sample
    geno[[chr]]$map <- as.data.frame(geno[[chr]]$map)
    ## add the usual variant IDs
    geno[[chr]]$map$ID <- paste0(chr, ":", geno[[chr]]$map$physical.pos, ":",
                                 geno[[chr]]$map$allele1, ":", geno[[chr]]$map$allele2)
    ## add A1 column
    geno[[chr]]$map$A1 <- coeff[[chr]]$A1

    ## confirm SNPs are in same order in initial PGS coefficient table as in bigSNP object
    if (sum(geno[[chr]]$map$ID != coeff[[chr]]$ID) > 0) {
        stop(paste0("Chromosome ", chr, ": SNPs in bigSNP object are in different order than in table of PGS coeffs."))
    }

    ## overwrite RDS file to save changes to fam and map dfs
    saveRDS(geno[[chr]], paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/geno/chr", chr, ".rds"))


    ## subset train
    geno_train_rds_chr <- snp_subset(geno[[chr]],
                                     ind.row = train_sp_ind,
                                     backingfile = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/geno/chr", chr, "-train"))
    geno_train[[chr]] <- snp_attach(geno_train_rds_chr)


    ## read validation data into memory as data.frame
    geno_vali[[chr]] <- list()  # geno_vali becomes list of lists temporarily
    geno_vali[[chr]]$genotypes <- geno[[chr]]$genotypes[vali_sp_ind, ]
    geno_vali[[chr]]$map <- geno[[chr]]$map

    cols_var_ids <- geno_vali[[chr]]$map$ID  # save SNP IDs before overwriting geno_vali[[chr]] (becomes df)
    geno_vali[[chr]] <- cbind(vali_sp_ids, geno_vali[[chr]]$genotypes[])
    colnames(geno_vali[[chr]])[-(1:2)] <- cols_var_ids


    ## subset rest
    geno_rest_rds_chr <- snp_subset(geno[[chr]],
                                    ind.row = rest_sp_ind,
                                    backingfile = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/geno/chr", chr, "-rest"))
    geno_rest[[chr]] <- snp_attach(geno_rest_rds_chr)


    ## delete original bk and RDS files
    system(paste0("rm ../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/geno/chr", chr, ".*"))
}

print("Finished loading genotypes")
print("Memory used so far:")
print(object.size(x=lapply(ls(), get)), units="Mb")

runtime <- toc()
runtime_df <- rbind(runtime_df,
                    data.frame(task = "Read data", time = as.numeric(runtime$toc - runtime$tic)))




## ITERATION ####################################################################

## df to store thresholds and validation performance
p1_opt_comp <- data.frame(step = 0,
                          p1_opt = round(p1_opt, digits = 8),
                          rsq_train = p1_opt_0$rsq_train,
                          rsq_vali  = p1_opt_0$rsq_vali,
                          rsq_test  = p1_opt_0$rsq_test)
                    
## get phenotype descriptions for plotting
all_desc  <- fread("../data/phenotypes/clean/code-desc-map-all.tab", data.table = FALSE)
phen_desc <- all_desc$desc_full[which(all_desc$field_id == phen)]

s <- 0
while (TRUE) {
    tic()
    s <- s + 1
    print(paste("Step", s))

    if (s == 1){
        p1s <- exp(seq(log(5e-8), log(1), length.out = 100))
    } else {
        p1s <- exp(seq(log(5e-8), log(1), length.out = 100))
        p1_opt_pos <- which(p1s == p1_opt)
        p1s <- p1s[1:p1_opt_pos]
    }

    gwas <- list()
    clump <- list()
    geno_score_vali <- list()
    geno_score_train <- list()

    for (chr in chr_pgs_orig) {
        ## run GWAS -------------------------------------------------------------
        linreg <- big_univLinReg(geno_train[[chr]]$genotypes,
                                 resid_train[, ncol(resid_train)],
                                 ncores = threads - 1)

        ## rerun regressions that failed
        if (sum(is.nan(linreg$estim)) > 0) {
            nan_snp_ind <- rows_along(linreg)[is.nan(linreg$estim)]
            for (i in nan_snp_ind) {
                reg <- lm(resid_train[, ncol(resid_train)] ~
                              geno_train[[chr]]$genotypes[, i])
                linreg[i, 1:2] <- summary(reg)$coeff[2, 1:2]
                linreg[i, 3]   <- linreg[i, 1] / linreg[i, 2]
            }
        }

        gwas[[chr]] <- coeff[[chr]]
        if (sum(geno_train[[chr]]$map$ID != gwas[[chr]]$ID) > 0) {
            stop(paste0("Chromosome ", chr, ": SNPs in training data bigSNP object are in different order than in table of GWAS coeffs."))
        }

        gwas[[chr]]$BETA <- linreg$estim
        gwas[[chr]]$LOG10_P <- -(log10(2) + pnorm(-abs(linreg$score), log = TRUE) / log(10))


        ## LD clumping ----------------------------------------------------------
        ind_keep <- snp_clumping(geno_train[[chr]]$genotypes,
                                 S = gwas[[chr]]$LOG10_P,
                                 thr.r2 = 0.1,
                                 size = 500,
                                 ncores = threads - 1,
                                 infos.chr = as.numeric(geno_train[[chr]]$map$chromosome),
                                 infos.pos = geno_train[[chr]]$map$physical.pos)
        clump[[chr]] <- gwas[[chr]][ind_keep,]


        ## compute PGS for different p1 thresholds ------------------------------
        coeff_chr_cpn <- clump[[chr]] %>%
            arrange(-LOG10_P)

        cols <- rep(0, length(p1s))
        for (i in 1:length(p1s)) {
            if (sum(coeff_chr_cpn$LOG10_P >= -log10(p1s[i])) > 0) {
                cols[i] <- max(which(coeff_chr_cpn$LOG10_P >= -log10(p1s[i])))
            }
        }
        
        geno_score_vali[[chr]] <- geno_vali[[chr]][, c("#FID", "IID")]
        geno_score_train[[chr]] <- train_sp_ids

        if (cols[1] == 0) {
            geno_score_vali[[chr]]$new <- 0
            colnames(geno_score_vali[[chr]])[ncol(geno_score_vali[[chr]])] <- paste0("S_", format(p1s[1], scientific = FALSE))

            geno_score_train[[chr]]$new <- 0
            colnames(geno_score_train[[chr]])[ncol(geno_score_train[[chr]])] <- paste0("S_", format(p1s[1], scientific = FALSE))
        } else {
            new_cols <- seq(1, cols[1])
            new_snps <- coeff_chr_cpn[new_cols, c("ID", "BETA", "POS")] %>%
                arrange(POS)
            
            keep_cols <- c("#FID", "IID", new_snps$ID)

            ## validation
            geno_score_vali_add <- geno_vali[[chr]][, keep_cols]
            snp_cols <- new_snps[, "ID"]
            for (k in 1:length(snp_cols)) {
                geno_score_vali_add[[snp_cols[k]]] <- geno_score_vali_add[[snp_cols[k]]] * new_snps[k, "BETA"]
            }
            geno_score_vali_add$tot <- rowSums(geno_score_vali_add[, 3:ncol(geno_score_vali_add), drop = FALSE])

            geno_score_vali[[chr]]$new <- geno_score_vali_add$tot
            colnames(geno_score_vali[[chr]])[ncol(geno_score_vali[[chr]])] <- paste0("S_", format(p1s[1], scientific = FALSE))

            ## training
            score <- big_prodVec(geno_train[[chr]]$genotypes,
                                 new_snps$BETA,
                                 ind.col = rows_along(geno_train[[chr]]$map)[geno_train[[chr]]$map$ID %in% new_snps$ID],
                                 ncores = threads - 1)
            geno_score_train[[chr]]$new <- score
            colnames(geno_score_train[[chr]])[ncol(geno_score_train[[chr]])] <- paste0("S_", format(p1s[1], scientific = FALSE))
        }
       
        if (length(p1s) > 1) {
            for (i in 2:length(p1s)) {
                if (cols[i] == cols[[i-1]]) {
                    last_col <- ncol(geno_score_vali[[chr]])
                    geno_score_vali[[chr]]$new <- geno_score_vali[[chr]][, last_col]
                    colnames(geno_score_vali[[chr]])[ncol(geno_score_vali[[chr]])] <- paste0("S_", format(p1s[i], scientific = FALSE))

                    last_col <- ncol(geno_score_train[[chr]])
                    geno_score_train[[chr]]$new <- geno_score_train[[chr]][, last_col]
                    colnames(geno_score_train[[chr]])[ncol(geno_score_train[[chr]])] <- paste0("S_", format(p1s[i], scientific = FALSE))
                } else {
                    new_cols <- seq(cols[i-1] + 1, cols[i])
                    new_snps <- coeff_chr_cpn[new_cols, c("ID", "BETA", "POS")] %>%
                        arrange(POS)
                    keep_cols <- c("#FID", "IID", new_snps$ID)

                    ## validation
                    geno_score_vali_add <- geno_vali[[chr]][, keep_cols]

                    snp_cols <- new_snps[, "ID"]
                    for (k in 1:length(snp_cols)) {
                        geno_score_vali_add[[snp_cols[k]]] <- geno_score_vali_add[[snp_cols[k]]] * new_snps[k, "BETA"]
                    }
                    geno_score_vali_add$tot <- rowSums(geno_score_vali_add[, 3:ncol(geno_score_vali_add), drop = FALSE])

                    last_col <- ncol(geno_score_vali[[chr]])
                    geno_score_vali[[chr]]$new <- geno_score_vali_add$tot + geno_score_vali[[chr]][, last_col]
                    colnames(geno_score_vali[[chr]])[ncol(geno_score_vali[[chr]])] <- paste0("S_", format(p1s[i], scientific = FALSE))

                    ## training
                    score_add <- big_prodVec(geno_train[[chr]]$genotypes,
                                             new_snps$BETA,
                                             ind.col = rows_along(geno_train[[chr]]$map)[geno_train[[chr]]$map$ID %in% new_snps$ID],
                                             ncores = threads - 1)

                    last_col <- ncol(geno_score_train[[chr]])
                    geno_score_train[[chr]]$new <- score_add + geno_score_train[[chr]][, last_col]
                    colnames(geno_score_train[[chr]])[ncol(geno_score_train[[chr]])] <- paste0("S_", format(p1s[i], scientific = FALSE))
                }
            }
        }
    }

    ## make Manhattan and Q-Q plots
    gwas_allchr <- bind_rows(gwas)
    manhattan <- ggmanh(gwas_allchr, title = paste0("PGS iteration ", s, ": ", phen_desc))
    ggsave(manhattan, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                        "/figs/", phen, "-step", s, "-manh.png"),
           type = "cairo-png", width = 800/120, height = 400/120, units = "in", dpi = 120)

    if (s == 1) {
        fwrite(gwas_allchr, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                          "/figs/", phen, "-step", s, "-gwas.tab"),
               sep = "\t", na = "NA", quote = FALSE)
    }

    qq <- ggqq(gwas_allchr, title = paste0("PGS iteration ", s, ": ", phen))
    ggsave(qq, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                 "/figs/", phen, "-step", s, "-qq.png"),
           type = "cairo-png", width = 180/120, height = 200/120, units = "in", dpi = 120)


    ## exit while loop if no SNP is genome-wide significant
    if (sum(gwas_allchr$LOG10_P >= -log10(5e-8)) == 0) {
        runtime <- toc()
        runtime_df <- rbind(runtime_df,
                            data.frame(task = "Partial iteration", time = as.numeric(runtime$toc - runtime$tic)))
        break
    }


    ## choose optimal p1 threshold ----------------------------------------------
    ## sum PGS components across chromosomes
    pgs_cpn_step_vali <- geno_score_vali[[chr_pgs_orig[1]]]
    for (chr in chr_pgs_orig[-1]) {
        pgs_cpn_step_vali[, 3:ncol(pgs_cpn_step_vali)] <- pgs_cpn_step_vali[, 3:ncol(pgs_cpn_step_vali)] +
            geno_score_vali[[chr]][, 3:ncol(geno_score_vali[[chr]])]
    } 

    pgs_cpn_step_train <- geno_score_train[[chr_pgs_orig[1]]]
    for (chr in chr_pgs_orig[-1]) {
        pgs_cpn_step_train[, 3:ncol(pgs_cpn_step_train)] <- pgs_cpn_step_train[, 3:ncol(pgs_cpn_step_train)] +
            geno_score_train[[chr]][, 3:ncol(geno_score_train[[chr]])]
    } 

    ## merge original phenotype, latest PGS and new scores for validation set
    phen_pgs_vali_cpn <- merge(phen_vali, pgs_cum[, c(1, 2, ncol(pgs_cum))], by = c("#FID", "IID"))
    phen_pgs_vali_cpn <- merge(phen_pgs_vali_cpn, pgs_cpn_step_vali, by = c("#FID", "IID"))
    
    phen_pgs_train_cpn <- merge(phen_train, pgs_cum[, c(1, 2, ncol(pgs_cum))], by = c("#FID", "IID"))
    phen_pgs_train_cpn <- merge(phen_pgs_train_cpn, pgs_cpn_step_train, by = c("#FID", "IID"))

    ## compute R-squared for each threshold
    pgs_cpn_step_comp <- data.frame(p1 = p1s,
                                    rsq_vali = NA,
                                    coeff_pgs_old = NA,
                                    coeff_pgs_cpn = NA)
        
    for (i in 1:length(p1s)) {
        phen_pgs_train_cpn_reg <- phen_pgs_train_cpn[, c(3, 4, 4 + i)]
        reg_train_cpn <- lm(phen_pgs_train_cpn_reg$phen_res ~
                               phen_pgs_train_cpn_reg[, 2] +
                               phen_pgs_train_cpn_reg[, 3])
        pgs_cpn_step_comp$coeff_pgs_old[i] <- summary(reg_train_cpn)$coeff[2, 1]
        pgs_cpn_step_comp$coeff_pgs_cpn[i] <- summary(reg_train_cpn)$coeff[3, 1]
 
        phen_pgs_vali_cpn_reg <- phen_pgs_vali_cpn[, c(3, 4, 4 + i)]
        phen_pgs_vali_cpn_reg$old_new <- pgs_cpn_step_comp$coeff_pgs_old[i] * phen_pgs_vali_cpn_reg[, 2] +
            pgs_cpn_step_comp$coeff_pgs_cpn[i] * phen_pgs_vali_cpn_reg[, 3]

        reg_vali_cpn <- lm(phen_pgs_vali_cpn_reg$phen_res ~ 
                               phen_pgs_vali_cpn_reg$old_new)
        ## reg_vali_cpn <- lm(phen_pgs_vali_cpn_reg$phen_res ~ 
        ##                        phen_pgs_vali_cpn_reg$pgs0 + phen_pgs_vali_cpn_reg[, 3])
 
        pgs_cpn_step_comp$rsq_vali[i] <- summary(reg_vali_cpn)$r.squared
    }
   
    ## optimal p-value threshold
    opt_col <- min(which(pgs_cpn_step_comp$rsq_vali == max(pgs_cpn_step_comp$rsq_vali)))
    p1_opt  <- pgs_cpn_step_comp$p1[opt_col]  # overwrites previous value
    rsq_opt <- pgs_cpn_step_comp$rsq_vali[opt_col]

    pgs_scaling_add <- data.frame(step = s,
                                  coeff_pgs_old_train = pgs_cpn_step_comp$coeff_pgs_old[opt_col],
                                  coeff_pgs_cpn_train = pgs_cpn_step_comp$coeff_pgs_cpn[opt_col])
    pgs_scaling <- rbind(pgs_scaling, pgs_scaling_add)

    ## plot R2 as a function of score size
    p1_plot <- ggp1_iter(pgs_cpn_step_comp, p1_opt = p1_opt, rsq_opt = rsq_opt,
                         title = paste0("PGS iteration ", s, " valid. perf. by p-value: ", phen_desc))
    ggsave(p1_plot, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                      "/figs/", phen, "-step", s, "-", vali_sp, "-p1_rsq.png"),
           type = "cairo-png", width = 700/120, height = 350/120, units = "in", dpi = 120)

 
    ## compute optimal PGS for training + validation samples --------------------
    chr_pgs_iter <- numeric()
    geno_score_train_opt <- list()
    geno_score_vali_opt  <- list()
    for (chr in chr_pgs_orig) {
        ## make score vectors
        coeff_chr_cpn <- clump[[chr]] %>%
            filter(LOG10_P >= -log10(p1_opt))
        
        if (nrow(coeff_chr_cpn) == 0) next
        chr_pgs_iter <- c(chr_pgs_iter, chr)

        ## training
        keep_cols_train <- c(1:2, 2 + opt_col)
        geno_score_train_opt[[chr]] <- geno_score_train[[chr]][, keep_cols_train]
        colnames(geno_score_train_opt[[chr]])[ncol(geno_score_train_opt[[chr]])] <- "tot"

        ## validation
        keep_cols_vali <- c(1:2, 2 + opt_col)
        geno_score_vali_opt[[chr]] <- geno_score_vali[[chr]][, keep_cols_vali]
        colnames(geno_score_vali_opt[[chr]])[ncol(geno_score_vali_opt[[chr]])] <- "tot"
    }
   
    geno_score_train_allchr <- cbind(geno_score_train_opt[[chr_pgs_iter[1]]][, 1:2],
                                     rowSums(sapply(geno_score_train_opt[chr_pgs_iter], `[[`, 3), na.rm = TRUE))
    colnames(geno_score_train_allchr)[3] <- "score"
    
    geno_score_vali_allchr <- cbind(geno_score_vali_opt[[chr_pgs_iter[1]]][, 1:2],
                                    rowSums(sapply(geno_score_vali_opt[chr_pgs_iter], `[[`, 3), na.rm = TRUE))
    colnames(geno_score_vali_allchr)[3] <- "score"

    ## add optimal PGS component to component df
    pgs_cum <- merge(pgs_cum, rbind(geno_score_train_allchr, geno_score_vali_allchr),
                     by = c("#FID", "IID"))
    colnames(pgs_cum)[ncol(pgs_cum)] <- paste0("pgs_cum", s)

    ## get old/new coeffs from training set
    phen_pgs_train_cum <- merge(phen_train, pgs_cum, by = c("#FID", "IID"))
    reg_train_cum      <- lm(phen_pgs_train_cum$phen_res ~
                                 phen_pgs_train_cum[, ncol(phen_pgs_train_cum) - 1] +
                                 phen_pgs_train_cum[, ncol(phen_pgs_train_cum)])

    ## compute new (cumulative) PGS by multiplying old PGS by *training* weight
    pgs_cum[, ncol(pgs_cum)] <- pgs_scaling$coeff_pgs_old_train[s] * pgs_cum[, ncol(pgs_cum) - 1] +
        pgs_scaling$coeff_pgs_cpn_train[s] * pgs_cum[, ncol(pgs_cum)]

    ## measure accumulated + weighted performance of all components so far on validation set
    phen_pgs_vali_cum     <- merge(phen_vali, pgs_cum, by = c("#FID", "IID"))
    reg_vali_cum          <- lm(phen_pgs_vali_cum$phen_res ~
                                    phen_pgs_vali_cum[, ncol(phen_pgs_vali_cum)])

    ## store residuals
    resid_train <- merge(resid_train, cbind(phen_pgs_train_cum[, c("#FID", "IID")], residuals(reg_train_cum)),
                         by = c("#FID", "IID"))
    colnames(resid_train)[ncol(resid_train)] <- paste0("pgs_res", s)
    resid_train <- resid_train[order(match(resid_train$IID, train_nonmiss_ids$IID)),]

    resid_vali <- merge(resid_vali, cbind(phen_pgs_vali_cum[, c("#FID", "IID")], residuals(reg_vali_cum)),
                         by = c("#FID", "IID"))
    colnames(resid_vali)[ncol(resid_vali)] <- paste0("pgs_res", s)
    resid_vali <- resid_vali[order(match(resid_vali$IID, vali_nonmiss_ids$IID)),]
    
    ## record performance
    p1_opt_comp <- rbind(p1_opt_comp,
                         c(s,
                           round(p1_opt, digits = 8),
                           round(summary(reg_train_cum)$r.squared, digits = 8),
                           round(summary(reg_vali_cum)$r.squared, digits = 8),
                           NA))
    fwrite(p1_opt_comp, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                      "/", phen, "-pgs-performance.tab"),
           sep = "\t", na = "NA", quote = FALSE)
 

    ## add component and accumulated coefficients to tables ---------------------
    for (chr in chr_pgs_orig) {
        ## remove SNPs above p-value threshold from LD clumping results
        coeff_chr_cum_step <- clump[[chr]] %>%
            filter(LOG10_P >= -log10(p1_opt)) %>%
            select(ID, BETA)
        colnames(coeff_chr_cum_step)[2] <- "coeff"

        if (nrow(coeff_chr_cum_step) == 0) {
            ## add to component table
            coeff_cum[[chr]]$new <- coeff_cum[[chr]][, ncol(coeff_cum[[chr]])] * pgs_scaling$coeff_pgs_old_train[s]
            colnames(coeff_cum[[chr]])[ncol(coeff_cum[[chr]])] <- paste0("step", s)
            next
        }

        ## add to component table
        coeff_cum[[chr]] <- merge(coeff_cum[[chr]], coeff_chr_cum_step, by = "ID", all = TRUE)
        colnames(coeff_cum[[chr]])[ncol(coeff_cum[[chr]])] <- paste0("step", s)
        coeff_cum[[chr]][[paste0("step", s)]][is.na(coeff_cum[[chr]][[paste0("step", s)]])] <- 0 
        coeff_cum[[chr]][[paste0("step", s)]] <- coeff_cum[[chr]][[paste0("step", s - 1)]] * pgs_scaling$coeff_pgs_old_train[s] +
            coeff_cum[[chr]][[paste0("step", s)]] * pgs_scaling$coeff_pgs_cpn_train[s]
        
        ## make sure original order is kept
        coeff_cum[[chr]] <- coeff_cum[[chr]][order(match(coeff_cum[[chr]]$ID, geno_train[[chr]]$map$ID)),]
    }

    runtime <- toc()
    runtime_df <- rbind(runtime_df,
                        data.frame(task = paste("Iteration", s), time = as.numeric(runtime$toc - runtime$tic)))
}




## export coefficients ----------------------------------------------------------

## export coefficient components table
for (chr in chr_pgs_orig) {
    coeff_cum_chr_out <- as.data.table(coeff_cum[[chr]])
    cols <- names(coeff_cum_chr_out)[3:ncol(coeff_cum_chr_out)]
    coeff_cum_chr_out[, (cols) := round(.SD, digits = 8), .SDcols = cols]

    fwrite(coeff_cum_chr_out, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                           "/coeff/cpn/", phen, "-coeff-cpn-chr", chr, ".tab"),
           sep = "\t", na = "NA", quote = FALSE)
}


## compute and export final coefficients
coeff_final <- list()
for (chr in chr_pgs_orig) {
    coeff_final[[chr]] <- coeff_cum[[chr]][, c(1, 2, ncol(coeff_cum[[chr]]))]
    colnames(coeff_final[[chr]])[3] <- "BETA"

    coeff_final_chr_out <- coeff_final[[chr]]
    coeff_final_chr_out$BETA <- round(coeff_final_chr_out$BETA, digits = 8)
    fwrite(coeff_final_chr_out, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                           "/coeff/final/", phen, "-coeff-chr", chr, ".tab"),
           sep = "\t", na = "NA", quote = FALSE)
}


## export table of old/new PGS coefficients
pgs_scaling_out <- as.data.table(pgs_scaling)
cols <- names(pgs_scaling_out)[2:ncol(pgs_scaling_out)]
pgs_scaling_out[, (cols) := round(.SD, digits = 8), .SDcols = cols]
fwrite(pgs_scaling_out, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                      "/coeff/", phen, "-pgs-scaling.tab"),
       sep = "\t", na = "NA", quote = FALSE)




## export scores ----------------------------------------------------------------

tic()

## compute PGS for each chromosome for all samples
geno_score_bychr_train <- train_sp_ids
for (chr in 1:22) {
    if (chr %in% chr_pgs_orig) {
        score <- big_prodVec(geno_train[[chr]]$genotypes,
                             coeff_final[[chr]]$BETA,
                             ncores = threads - 1)
    } else {
        score <- 0
    }
    geno_score_bychr_train$new <- score
    colnames(geno_score_bychr_train)[ncol(geno_score_bychr_train)] <- paste0("pgs_chr", chr)
}

geno_score_bychr_vali <- vali_sp_ids
for (chr in 1:22) {
    if (chr %in% chr_pgs_orig) {
        geno_score_bychr_vali$new <- as.matrix(geno_vali[[chr]][, c(-1, -2)]) %*% coeff_final[[chr]]$BETA
    } else {
        geno_score_bychr_vali$new <- 0
    }
    colnames(geno_score_bychr_vali)[ncol(geno_score_bychr_vali)] <- paste0("pgs_chr", chr)
}

geno_score_bychr_rest <- rest_sp_ids
for (chr in 1:22) {
    if (chr %in% chr_pgs_orig) {
        score <- big_prodVec(geno_rest[[chr]]$genotypes,
                             coeff_final[[chr]]$BETA,
                             ncores = threads - 1)
    } else {
        score <- 0
    }
    geno_score_bychr_rest$new <- score
    colnames(geno_score_bychr_rest)[ncol(geno_score_bychr_rest)] <- paste0("pgs_chr", chr)
}

geno_score_bychr <- rbind(geno_score_bychr_train,
                          geno_score_bychr_vali,
                          geno_score_bychr_rest)
geno_score_bychr <- geno_score_bychr[order(match(geno_score_bychr$IID, all_sp_ids$IID)),]

geno_score_bychr_out <- as.data.table(geno_score_bychr)
cols <- names(geno_score_bychr_out)[3:ncol(geno_score_bychr_out)]
geno_score_bychr_out[, (cols) := round(.SD, digits = 8), .SDcols = cols]
fwrite(geno_score_bychr_out, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                               "/pgs/", phen, "-all-pgs-bychr.tab"),
       sep = "\t", na = "NA", quote = FALSE)


## compute full PGS
geno_score_full <- all_sp_ids
geno_score_full$pgs_full <- rowSums(geno_score_bychr[, 3:ncol(geno_score_bychr)])

geno_score_full_out <- geno_score_full
geno_score_full_out$pgs_full <- round(geno_score_full_out$pgs_full, digits = 8)
fwrite(geno_score_full_out, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                          "/pgs/", phen, "-all-pgs-full.tab"),
       sep = "\t", na = "NA", quote = FALSE)


## compute leave-one-chromosome-out PGS
geno_score_loco <- all_sp_ids
for (chr in 1:22) {
    if (chr %in% chr_pgs_orig) {
        geno_score_loco$new <- rowSums(geno_score_bychr[, setdiff(1:22, chr) + 2])
    } else {
        geno_score_loco$new <- geno_score_full$pgs_full
    }
    colnames(geno_score_loco)[ncol(geno_score_loco)] <- paste0("pgs_loco", chr)
}

geno_score_loco_out <- as.data.table(geno_score_loco)
cols <- names(geno_score_loco_out)[3:ncol(geno_score_loco_out)]
geno_score_loco_out[, (cols) := round(.SD, digits = 8), .SDcols = cols]
fwrite(geno_score_loco_out, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                          "/pgs/", phen, "-all-pgs-loco.tab"),
       sep = "\t", na = "NA", quote = FALSE)


## delete folder with bigSNP files
system(paste0("rm -rf ../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/geno"))




## evaluate performance on test set ---------------------------------------------

phen_test <- fread(paste0("../results/01-gwas/residuals-covars/", phen,
                          "/", phen, "-", test_sp, "-resid-covars.tab"),
                   data.table = FALSE)
colnames(phen_test)[3] <- "phen_res"
phen_test <- phen_test[!is.na(phen_test$phen_res),]

## evaluate on test set
phen_pgs_test_full <- merge(phen_test, geno_score_full, by = c("#FID", "IID"))
reg_test_full <- lm(phen_res ~ pgs_full, data = phen_pgs_test_full)
p1_opt_comp$rsq_test[nrow(p1_opt_comp)] <- round(summary(reg_test_full)$r.squared, digits = 8)
fwrite(p1_opt_comp, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                  "/", phen, "-pgs-performance.tab"),
       sep = "\t", na = "NA", quote = FALSE)





## make plot showing evolution of rsq with iteration steps ----------------------

colours <- c("Training" = "#0073AD", "Validation" = "#D35F27", "Test" = "#04A472")
p1_rsq_plt <- ggplot(data = p1_opt_comp, aes(x = step)) +
    ## train
    geom_point(aes(y = rsq_train, colour="Training"), size = 0.5) +
    geom_line(aes(y = rsq_train, colour="Training"), size = 0.5) +
    ## test
    geom_point(data = p1_opt_comp[c(1, nrow(p1_opt_comp)),],
               aes(y = rsq_test, colour="Test"), size = 1) +
    geom_segment(data = p1_opt_comp[c(1, nrow(p1_opt_comp)),],
                 aes(xend=step, y = 0, yend = rsq_test, colour = "Test"), size = 0.5) +
    geom_label(data = p1_opt_comp[c(1, nrow(p1_opt_comp)),],
               aes(label=round(rsq_test, digits = 4), y=rsq_test),
               size=2.5, nudge_y = -0.2 * max(p1_opt_comp$rsq_test, na.rm = TRUE)) +
    ## vali
    geom_point(aes(y = rsq_vali, colour="Validation"), size = 0.5) +
    geom_line(aes(y = rsq_vali, colour="Validation"), size = 0.5) +
    ylim(0, NA) +
    ggtitle(paste0("PGS performance evolution: ", phen_desc)) +
    labs(x = "Iteration",
         y = "R-squared",
         colour = "Legend") +
    scale_colour_manual(values = colours) +
    theme_bw() +
    theme(plot.title = element_text(size = 9),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7),
          legend.position = c(0.5, 0.08),
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.background = element_rect(fill = "white", color = "black", size = 0.1),
          legend.margin=margin(c(1, 10, 1, 1)))
ggsave(p1_rsq_plt, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                     "/figs/", phen, "-iter_step-rsq.png"),
       type = "cairo-png", width = 800/120, height = 400/120, units = "in", dpi = 120)
ggsave(p1_rsq_plt, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                     "/figs/", phen, "-iter_step-rsq.pdf"),
       device = cairo_pdf, width = 14, height = 7, units = "cm")



runtime <- toc()
runtime_df <- rbind(runtime_df,
                    data.frame(task = "Export results", time = as.numeric(runtime$toc - runtime$tic)))

## plot runtime per iteration
colours <- c("Training" = "#0073AD", "Validation" = "#D35F27", "Test" = "#04A472")
runtime_df_plt <- runtime_df[2:(nrow(runtime_df) - 2),]
runtime_df_plt$iter <- seq(1:nrow(runtime_df_plt))
runtime_plt <- ggplot(data = runtime_df_plt, aes(x = iter, y = time)) +
    ## train
    geom_point(colour="#0073AD", size = 0.5) +
    geom_line(colour="#0073AD", size = 0.5) +
    ylim(0, NA) +
    ggtitle(paste0("Runtime per iteration: ", phen_desc)) +
    labs(x = "Iteration",
         y = "Time in seconds") +
    theme_bw() +
    theme(plot.title = element_text(size = 9),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7))
ggsave(runtime_plt, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                      "/figs/", phen, "-iter_step-runtime.png"),
       type = "cairo-png", width = 800/120, height = 400/120, units = "in", dpi = 120)
ggsave(runtime_plt, filename = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp,
                                      "/figs/", phen, "-iter_step-runtime.pdf"),
       device = cairo_pdf, width = 14, height = 7, units = "cm")


runtime_df$time <- round(runtime_df$time, digits = 3)
fwrite(runtime_df, file = paste0("../results/02-pgs/iterative/", phen, "/r2_", r2, "-kb_", kb, "/", train_sp, "/runtime.txt"),
       sep = "\t", na = "NA", quote = FALSE)
