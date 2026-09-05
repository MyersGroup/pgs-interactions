## desc: Prepare files with sample IDs

library(data.table)


## load covariates file
covars <- fread('../data/covars/age-sex-batch-centre-ac-all.tab', data.table = FALSE)

## create list of all samples
all_ids_df <- covars[, c('#FID', 'IID')]
fwrite(all_ids_df, file = '../data/sample-ids/filtered/all-ids.tab',
       sep = '\t', na = 'NA', quote = FALSE)



## make list of white British IIDs
wb_ind_mar <- readRDS('../data/covars/ancestry/wb_ids.rds')
imp_sp <- fread('../data/sample-ids/ukb-27960-imp-auto-s487256-20210614.psam', data.table = FALSE)
wb_ids_mar <- imp_sp$IID[wb_ind_mar]

## keep only those in covariates file
wb_ids_df <- covars[covars$IID %in% wb_ids_mar, c('#FID', 'IID')]
fwrite(wb_ids_df, file = '../data/sample-ids/filtered/wb_all-ids.tab',
       sep = '\t', na = 'NA', quote = FALSE)



## make validation and test sets of White British individuals not included in J Marchini's 340k
## compute share of British ancestry (excluding Ireland)
covars$british <- rowSums(subset(covars, select=G_Anglia:G_Noi))

## keep only White British samples who are not in the set of WB IIDs above (32287 samples)
wb_ids_tv_df <- covars[covars$british >= 0.9 & !(covars$IID %in% wb_ids_df$IID), c('#FID', 'IID')]

## load relatedness table and keep only 2nd degree relationships or closer
ukb_rel <- fread('../data/relatedness/ukb_rel_a27960_s488224.dat', data.table = FALSE)
ukb_rel_2 <- ukb_rel[ukb_rel$Kinship >= 0.0884,]

## remove samples related (2nd degree or closer) with anyone in training set
wb_ids_tv_related_w_train <- numeric()
for (i in 1:nrow(wb_ids_tv_df)) {
    if (wb_ids_tv_df$IID[i] %in% ukb_rel_2$ID1) {  # is test/vali IID in relatedness tbl (col 1)?
        rel_i <- which(ukb_rel_2$ID1 == wb_ids_tv_df$IID[i])  # get row index in relatedness tbl
        for (j in rel_i) {
            if (ukb_rel_2$ID2[j] %in% wb_ids_df$IID) {  # is corresponding ID2 in training set?
                wb_ids_tv_related_w_train <- c(wb_ids_tv_related_w_train, wb_ids_tv_df$IID[i])  # add to exclusion list
                next
            }
        }
    } else if (wb_ids_tv_df$IID[i] %in% ukb_rel_2$ID2) {
        rel_i <- which(ukb_rel_2$ID2 == wb_ids_tv_df$IID[i])
        for (j in rel_i) {
            if (ukb_rel_2$ID1[j] %in% wb_ids_df$IID) {
                wb_ids_tv_related_w_train <- c(wb_ids_tv_related_w_train, wb_ids_tv_df$IID[i])
                next
            }
        } 
    }
}

wb_ids_tv_unrel_w_train <- wb_ids_tv_df[!(wb_ids_tv_df$IID %in% wb_ids_tv_related_w_train),]
## 19694 samples remaining

## relatedness matrix of these samples (among themselves)
ukb_rel_2_tv <- ukb_rel_2[ukb_rel_2$ID1 %in% wb_ids_tv_unrel_w_train$IID &
                          ukb_rel_2$ID2 %in% wb_ids_tv_unrel_w_train$IID,]

## remove test/vali samples related (2nd degree or closer) among themselves
freqs <- as.data.frame(table(c(ukb_rel_2_tv$ID1, ukb_rel_2_tv$ID2)))
colnames(freqs) <- c("IID", "freq")
freqs$IID <- as.numeric(levels(freqs$IID))[freqs$IID]

freqs <- freqs[freqs$freq > 1,]
freqs <- freqs[order(-freqs$freq, freqs$IID),]

ukb_rel_2_tv_trim <- ukb_rel_2_tv
freqs_trim <- freqs
i <- 1
while (nrow(freqs_trim) > 0) {
    ukb_rel_2_tv_trim <- ukb_rel_2_tv_trim[ukb_rel_2_tv_trim$ID1 != freqs_trim$IID[1] &
                                           ukb_rel_2_tv_trim$ID2 != freqs_trim$IID[1],]

    freqs_trim <- as.data.frame(table(c(ukb_rel_2_tv_trim$ID1, ukb_rel_2_tv_trim$ID2)))
    colnames(freqs_trim) <- c("IID", "freq")
    freqs_trim$IID <- as.numeric(levels(freqs_trim$IID))[freqs_trim$IID]

    freqs_trim <- freqs_trim[freqs_trim$freq > 1,]
    freqs_trim <- freqs_trim[order(-freqs_trim$freq, freqs_trim$IID),]
    i <- i + 1
}

## then once everyone appears only once we remove half at random
wb_ids_tv_unrel_all <- wb_ids_tv_unrel_w_train[!(wb_ids_tv_unrel_w_train$IID %in% ukb_rel_2_tv_trim$ID2),]
## 19071 samples remaining

## sample 10k to make validation set
set.seed(394857)
wb_ids_vali_cols <- sort(sample(1:nrow(wb_ids_tv_unrel_all), 10000))
wb_ids_vali_df   <- wb_ids_tv_unrel_all[wb_ids_vali_cols,]
## match original order
wb_ids_vali_df   <- wb_ids_vali_df[order(match(wb_ids_vali_df$IID, covars$IID)),]
fwrite(wb_ids_vali_df, file = '../data/sample-ids/filtered/other_gb_vali_all-ids.tab',
       sep = '\t', na = 'NA', quote = FALSE)

## remaining 9071 are test set
wb_ids_test_df <- wb_ids_tv_unrel_all[setdiff(1:nrow(wb_ids_tv_unrel_all), wb_ids_vali_cols),]
## match original order
wb_ids_test_df   <- wb_ids_test_df[order(match(wb_ids_test_df$IID, covars$IID)),]
fwrite(wb_ids_test_df, file = '../data/sample-ids/filtered/other_gb_test_all-ids.tab',
       sep = '\t', na = 'NA', quote = FALSE)
