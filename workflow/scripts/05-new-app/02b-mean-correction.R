## Adjust PGS to mean phenotype value within PGS quantiles

library(argparser)
library(data.table)
library(dplyr)
library(Hmisc)

p <- arg_parser("Argument parser")
p <- add_argument(p, "--phen",     help = "Phenotype code", nargs = 1)
p <- add_argument(p, "--samples" , help = "Samples",        nargs = 1)
argv <- parse_args(p)

phen     <- argv$phen
samples  <- argv$samples

r2       <- 0.9
kb       <- 500


## load sample IIDs
sample_ids <- fread(paste0("../data/sample-ids/filtered/", samples, "-ids.tab"),
                    data.table = FALSE)

## load residuals after regressing out covars
phen_res_cov <- fread(paste0("../results/05-new-app/residuals-covars/", phen, "/", phen, "-", samples, "-resid-covars.tab"),
                      data.table = FALSE)
colnames(phen_res_cov)[3] <- "phen_res_cov"
phen_res_cov <- phen_res_cov[!is.na(phen_res_cov$phen_res_cov),]

## load full PGS
pgs_full <- fread(paste0("../results/05-new-app/iterative/", phen, "/r2_0.9-kb_500/wb_all/pgs/", phen, "-all-pgs-full.tab"),
                  data.table = FALSE)
## load LOCO PGS
pgs_loco <- fread(paste0("../results/05-new-app/iterative/", phen, "/r2_0.9-kb_500/wb_all/pgs/", phen, "-all-pgs-loco.tab"),
                  data.table = FALSE)

## merge
phen_pgs <- left_join(phen_res_cov, pgs_full, by = c("#FID", "IID"))
phen_pgs <- left_join(phen_pgs, pgs_loco, by = c("#FID", "IID"))



## add mean correction

## full PGS

## slice
pgs_full_mc <- phen_pgs %>%
    select("#FID", IID, phen_res_cov, pgs_full) %>%
    arrange(pgs_full)

## grouping
if (samples == "wb_all") {
    group_vec <- c(rep(1:20, each = 50),
                   20 + rep(1:(ceiling(nrow(phen_pgs)/1000) - 1), each = 1000))
    group_vec <- group_vec[1:(nrow(phen_pgs) - 1000)]
    last <- group_vec[length(group_vec)]
    group_vec <- c(group_vec,
                   last + rep(1:20, each = 50))
} else {
    group_vec <- c(rep(1, each = floor(nrow(phen_pgs) %% 100 / 2)),
                   rep(2:((nrow(phen_pgs) - nrow(phen_pgs) %% 100) / 100 + 1), each = 100),
                   rep(((nrow(phen_pgs) - nrow(phen_pgs) %% 100) / 100 + 2), each = ceiling(nrow(phen_pgs) %% 100 / 2)))
}

## mean y / median x per bin 
pgs_full_mc$group <- group_vec
pgs_full_mc <- pgs_full_mc %>%
    group_by(group) %>%
    mutate(phen_res_cov_mean = mean(phen_res_cov)) %>%
    as.data.frame()

## linear interpolation
x_int <- pgs_full_mc %>%
    group_by(group) %>%
    dplyr::summarise(pgs_full_median = median(pgs_full))
x_int <- as.numeric(x_int$pgs_full_median)
y_int <- pgs_full_mc %>%
    group_by(group) %>%
    dplyr::summarise(phen_res_cov_mean = mean(phen_res_cov))
y_int <- as.numeric(y_int$phen_res_cov_mean)
pgs_full_mc$phen_res_cov_int <- approxExtrap(x = x_int,
                                             y = y_int,
                                             xout = pgs_full_mc$pgs_full)$y

## export
pgs_full_mc_out <- pgs_full_mc[, c("#FID", "IID", "phen_res_cov_int")]
colnames(pgs_full_mc_out)[3] <- "pgs_full_mc"
pgs_full_mc_out$pgs_full_mc <- round(pgs_full_mc_out$pgs_full_mc, digits = 8)

## add IIDs and samples with missing values
pgs_full_mc_out <- left_join(sample_ids, pgs_full_mc_out, by = c("#FID", "IID"))
fwrite(pgs_full_mc_out,
       file = paste0("../results/05-new-app/mean-corrected/", phen,
                     "/r2_", r2, "-kb_", kb, "/", phen, "-", samples, "-pgs-full_mc.tab"),
       sep = "\t", na = "NA", quote = FALSE)



## LOCO PGS
for (chr in 1:22) {
    ## slice
    pgs_loco_mc <- phen_pgs[, c("#FID", "IID", "phen_res_cov", paste0("pgs_loco", chr))]
    colnames(pgs_loco_mc)[4] <- "pgs_loco"
    pgs_loco_mc <- pgs_loco_mc[with(pgs_loco_mc, order(pgs_loco)),]

    ## grouping
    if (samples == "wb_all") {
        group_vec <- c(rep(1:20, each = 50),
                       20 + rep(1:(ceiling(nrow(phen_pgs)/1000) - 1), each = 1000))
        group_vec <- group_vec[1:(nrow(phen_pgs) - 1000)]
        last <- group_vec[length(group_vec)]
        group_vec <- c(group_vec,
                       last + rep(1:20, each = 50))
    } else {
        group_vec <- c(rep(1, floor(nrow(phen_pgs) %% 100 / 2)),
                       rep(2:((nrow(phen_pgs) - nrow(phen_pgs) %% 100) / 100 + 1), each = 100),
                       rep(((nrow(phen_pgs) - nrow(phen_pgs) %% 100) / 100 + 2), ceiling(nrow(phen_pgs) %% 100 / 2)))
    }

    ## mean y / median x per bin 
    pgs_loco_mc$group <- group_vec
    pgs_loco_mc <- pgs_loco_mc %>%
        group_by(group) %>%
        mutate(phen_res_cov_mean = mean(phen_res_cov)) %>%
        as.data.frame()

    ## linear interpolation
    x_int <- pgs_loco_mc %>%
        group_by(group) %>%
        dplyr::summarise(pgs_loco_median = median(pgs_loco))
    x_int <- as.numeric(x_int$pgs_loco_median)
    y_int <- pgs_loco_mc %>%
        group_by(group) %>%
        dplyr::summarise(phen_res_cov_mean = mean(phen_res_cov))
    y_int <- as.numeric(y_int$phen_res_cov_mean)
    pgs_loco_mc$phen_res_cov_int <- approxExtrap(x = x_int,
                                                 y = y_int,
                                                 xout = pgs_loco_mc$pgs_loco)$y

    ## export
    pgs_loco_mc_out <- pgs_loco_mc[, c("#FID", "IID", "phen_res_cov_int")]
    pgs_loco_mc_out$phen_res_cov_int <- round(pgs_loco_mc_out$phen_res_cov_int, digits = 8)
    colnames(pgs_loco_mc_out)[3] <- paste0("pgs_loco", chr, "_mc")

    ## add IIDs and samples with missing values
    pgs_loco_mc_out <- left_join(sample_ids, pgs_loco_mc_out, by = c("#FID", "IID"))

    ## add full PGS to df
    pgs_loco_mc_out <- left_join(pgs_loco_mc_out, pgs_full_mc_out, by = c("#FID", "IID"))
    
    fwrite(pgs_loco_mc_out,
           file = paste0("../results/05-new-app/mean-corrected/", phen,
                         "/r2_", r2, "-kb_", kb, "/", phen, "-", samples, "-pgs-loco", chr, "_mc.tab"),
           sep = "\t", na = "NA", quote = FALSE)
}
