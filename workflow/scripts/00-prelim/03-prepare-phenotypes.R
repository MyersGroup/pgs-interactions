## desc: Prepare text files with phenotypes for each sample
##       Include traits chosen based on heritability and
##       sample size as well as 5 additional traits:
##       - f.leglength: leg length (height - sitting height)
##       - f.leglength_rel: relative leg length (leg length / height)
##       - f.20015_rel: relative sitting height (sitting height / height)
##       - f.waist_hip: waist-to-hip ratio (waist circumference / hip circumference)
##       - f.fev1_fvc: FEV1/FVC ratio



library(data.table)
library(dplyr)

options(scipen=999)  # disable scientific notation




## load list of samples to consider
all_sp <- fread('../data/sample-ids/filtered/all-ids.tab', data.table = FALSE)
all_sp_ids <- all_sp$IID


## load raw UKB data files, keep only relevant samples

## ID 45188, Basket ID 2010858, 4 Feb 2021 
## load column names only due to file size
phen_45188_cols <- fread('../data/phenotypes/ukb-raw/ukb45188.tab',
                         header = FALSE, nrows = 1, data.table = FALSE)
phen_45188_cols <- unname(unlist(phen_45188_cols))  # get vector

## ID 26618, Basket ID 2002132, 12 Feb 2019
phen_26618 <- fread('../data/phenotypes/ukb-raw/ukb26618.tab', data.table = FALSE)
phen_26618 <- phen_26618[phen_26618$f.eid %in% all_sp_ids,]
## order by all-ids.tab file (matches order of original sample file)
phen_26618 <- phen_26618[match(all_sp_ids, phen_26618$f.eid),]
phen_26618_cols <- colnames(phen_26618)


## load table of traits
trait_tbl <- fread("../data/phenotypes/clean/code-desc-map-real.tab", data.table = FALSE)
traits <- trait_tbl$ukb_id[!is.na(trait_tbl$ukb_id)]






## process traits in dataset 45188
## this dataset is more recent and therefore should be prioritised
## only traits which aren't available in this dataset are taken from ID 26618

## list of traits available in phenotype table
phen_45188_ids_list <- strsplit(phen_45188_cols, split = '.', fixed = TRUE)
## get 2nd and 4th item from each element in list
phen_45188_ids <- sapply(phen_45188_ids_list, '[', 2)
phen_45188_array <- sapply(phen_45188_ids_list, '[', 4)

## initialise empty vector of traits which are not available in phenotype table
na_traits_45188 <- vector()
## initialise empty vector of traits with more than one array index
array_traits_45188 <- vector()

## function to get right-most element which isn't NA
get_right <- function(x) {
    if (!is.numeric(x)) { stop('input not a numeric vector') }
    right_index <- which(!is.na(x))  # index of right-most non-NA/NaN element
    if (length(right_index) != 0) {
        return(x[[max(right_index)]])
    } else {  # if all elements are NA return NA
        return(NA)
    }
}

## process traits with only one array index
for (t in traits) {
    print(t)

    if (t %in% phen_45188_ids) { 
        cols <- which(phen_45188_ids == t)  # get column indices

        ## load data for these columns only
        phen_45188_t <- fread('../data/phenotypes/ukb-raw/ukb45188.tab',
                              select = c(1, cols), data.table = FALSE)
        ## keep only relevant IIDs
        phen_45188_t <- phen_45188_t[phen_45188_t$f.eid %in% all_sp_ids,]
        ## order by all-ids.tab file (matches order of original sample file)
        phen_45188_t <- phen_45188_t[match(all_sp_ids, phen_45188_t$f.eid),]


        ## check that there are no negative values indicating no response, etc
        ## field 78 (Heel bone mineral density (BMD) T-score, automated) is excluded
        ## because it can have meaningful negative values
        if (t != 78 & any(phen_45188_t[, 2:ncol(phen_45188_t)] < 0, na.rm = TRUE)) { 
            stop(paste('phenotype', t, 'has negative values.')) 
        }


        ## check that there is only one array index per instance (visit)
        if (identical(as.numeric(phen_45188_array[cols]), rep(0, length(cols)))) {

            ## get right-most value if there is more than one
            if (length(cols) > 1) {
                latest <- apply(phen_45188_t[, 2:ncol(phen_45188_t)], 1, get_right)
            } else {
                latest <- phen_45188_t[, 2]
            }
                
            ## add eid column
            latest_df <- data.frame(phen_45188_t$f.eid, latest)
            colnames(latest_df) <- c('#IID', paste0('f.', t))
            ## export
            fwrite(latest_df, file = paste0('../data/phenotypes/clean/f.', t, '.tab'),
                   sep = '\t', na = 'NA', quote = FALSE)

        } else {
            array_traits_45188 <- c(array_traits_45188, t)
        }
    } else {
        na_traits_45188 <- c(na_traits_45188, t)
    }
}


## 3 traits with multiple array indices:
##   102 - Pulse rate, automated reading
##   4079 - Diastolic blood pressure, automated reading
##   4080 - Systolic blood pressure, automated reading
##
## Blood pressure is measure twice at each visit and pulse rate is
## measured at the same time.  We take the *average* of the two
## measurements from the most recent visit. If only one measurement is
## available at the most recent visit, we take that

for (t in array_traits_45188) {
    cols <- which(phen_45188_ids == t)  # get column indices

    ## load data for these columns only
    phen_45188_t <- fread('../data/phenotypes/ukb-raw/ukb45188.tab',
                          select = c(1, cols), data.table = FALSE)
    ## keep only relevant IIDs
    phen_45188_t <- phen_45188_t[phen_45188_t$f.eid %in% all_sp_ids,]
    ## order by all-ids.tab file (matches order of original sample file)
    phen_45188_t <- phen_45188_t[match(all_sp_ids, phen_45188_t$f.eid),]


    ## average of reads from each visit/instance
    df <- data.frame(inst.1 = rowMeans(phen_45188_t[, 2:3], na.rm=TRUE),
                     inst.2 = rowMeans(phen_45188_t[, 4:5], na.rm=TRUE),
                     inst.3 = rowMeans(phen_45188_t[, 6:7], na.rm=TRUE),
                     inst.4 = rowMeans(phen_45188_t[, 8:9], na.rm=TRUE))
    ## keep most recent non-NA value
    latest <- apply(df, 1, get_right)
    ## add eid column
    latest_df <- data.frame(phen_45188_t$f.eid, latest)
    colnames(latest_df) <- c('#IID', paste0('f.', t))
    ## export
    fwrite(latest_df, file = paste0('../data/phenotypes/clean/f.', t, '.tab'),
           sep = '\t', na = 'NA', quote = FALSE)
}


## add leg length = standing height (50) - sitting height (20015)
df_50 <- fread('../data/phenotypes/clean/f.50.tab')
df_20015 <- fread('../data/phenotypes/clean/f.20015.tab')
df_leglength <- merge(df_50, df_20015, by='#IID')
df_leglength <- df_leglength[match(all_sp_ids, df_leglength$'#IID'),]  ## reorder IIDs
## compute leg length
df_leglength$f.leglength <- df_leglength$f.50 - df_leglength$f.20015

## add relative leg length = leg length / standing height (50)
df_leglength$f.leglength_rel <- df_leglength$f.leglength / df_leglength$f.50
## add relative sitting height = sitting height / standing height (50)
df_leglength$f.20015_rel <- df_leglength$f.20015 / df_leglength$f.50

## export
df_leglength_out <- df_leglength[, c('#IID', 'f.leglength')]
fwrite(df_leglength_out, file = '../data/phenotypes/clean/f.leglength.tab',
       sep = '\t', na = 'NA', quote = FALSE)
df_leglength_rel_out <- df_leglength[, c('#IID', 'f.leglength_rel')]
df_leglength_rel_out$f.leglength_rel <- round(df_leglength_rel_out$f.leglength_rel, digits = 8)
fwrite(df_leglength_rel_out, file = '../data/phenotypes/clean/f.leglength_rel.tab',
       sep = '\t', na = 'NA', quote = FALSE)
df_20015_rel_out <- df_leglength[, c('#IID', 'f.20015_rel')]
df_20015_rel_out$f.20015_rel <- round(df_20015_rel_out$f.20015_rel, digits = 8)
fwrite(df_20015_rel_out, file = '../data/phenotypes/clean/f.20015_rel.tab',
       sep = '\t', na = 'NA', quote = FALSE)


## add waist-to-hip ratio
df_48 <- fread('../data/phenotypes/clean/f.48.tab', data.table = FALSE)
df_49 <- fread('../data/phenotypes/clean/f.49.tab', data.table = FALSE)
df_waist_hip <- merge(df_48, df_49, by='#IID')
df_waist_hip <- df_waist_hip[match(all_sp_ids, df_waist_hip$'#IID'),]  ## reorder IIDs
df_waist_hip$f.waist_hip <- df_waist_hip$f.48 / df_waist_hip$f.49

## export
df_waist_hip_out <- df_waist_hip[, c('#IID', 'f.waist_hip')]
df_waist_hip_out$f.waist_hip <- round(df_waist_hip_out$f.waist_hip, digits = 8)
fwrite(df_waist_hip_out, file = '../data/phenotypes/clean/f.waist_hip.tab',
       sep = '\t', na = 'NA', quote = FALSE)


## add FEV1/FVC ratio
df_20150 <- fread('../data/phenotypes/clean/f.20150.tab', data.table = FALSE)
df_20151 <- fread('../data/phenotypes/clean/f.20151.tab', data.table = FALSE)
df_fev1_fvc <- merge(df_20150, df_20151, by='#IID')
df_fev1_fvc <- df_fev1_fvc[match(all_sp_ids, df_fev1_fvc$'#IID'),]  ## reorder IIDs
df_fev1_fvc$f.fev1_fvc <- df_fev1_fvc$f.20150 / df_fev1_fvc$f.20151

## export
df_fev1_fvc_out <- df_fev1_fvc[, c('#IID', 'f.fev1_fvc')]
df_fev1_fvc_out$f.fev1_fvc <- round(df_fev1_fvc_out$f.fev1_fvc, digits = 8)
fwrite(df_fev1_fvc_out, file = '../data/phenotypes/clean/f.fev1_fvc.tab',
       sep = '\t', na = 'NA', quote = FALSE)






## process traits in dataset with ID 26618

## only four traits are not available in dataset 45188:
## - 23105: Basal metabolic rate
## - 20016: Fluid intelligence score
## - 399:   Number of incorrect matches in round
## - 20023: Mean time to correctly identify matches

## list of traits available in phenotype table
phen_26618_ids_list <- strsplit(phen_26618_cols, split = '.', fixed = TRUE)
## get 2nd and 4th item from each element in list
phen_26618_ids <- sapply(phen_26618_ids_list, '[', 2)
phen_26618_array <- sapply(phen_26618_ids_list, '[', 4)

## initialise empty vector of traits which are not available in phenotype table
na_traits_26618 <- vector()
## initialise empty vector of traits with more than one array index
array_traits_26618 <- vector()

## process traits with only one array index
for (t in na_traits_45188) {
    print(t)

    if (t %in% phen_26618_ids) { 
        cols <- which(phen_26618_ids == t)  # get column indices

        ## check that there are no negative values indicating no response, etc
        if (any(phen_26618[, cols] < 0, na.rm = TRUE)) { 
            stop(paste('phenotype', t, 'has negative values.')) 
        }

        ## check that there is only one array index per instance (visit)
        if (identical(as.numeric(phen_26618_array[cols]), rep(0, length(cols)))) {
            latest <- apply(phen_26618[, cols], 1, get_right)
            ## add eid column
            latest_df <- data.frame(phen_26618$f.eid, latest)
            colnames(latest_df) <- c('#IID', paste0('f.', t))
            ## export
            fwrite(latest_df, file = paste0('../data/phenotypes/clean/f.', t, '.tab'),
                   sep = '\t', na = 'NA', quote = FALSE)
        } else {
            array_traits_26618 <- c(array_traits_26618, t)
        }
    } else {
        na_traits_26618 <- c(na_traits_26618, t)
    }
}


## 1 trait with multiple array indices:
##   399 - Number of incorrect matches in round
## Two rounds of pair matching are conducted in each visit and the
## number of incorrect matches is recorded for each. The docs for this
## category at
## https://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=100030 say that
## a third round was conducted in the pilot, which explains why the
## array has 3 indices. We *sum* the number of incorrect matches in the
## first two rounds from the most recent visit
cols_399   <- which(phen_26618_ids == 399)  # get column indices
## check that there are no negative values indicating no response, etc
if (any(phen_26618[, cols_399] < 0, na.rm = TRUE)) {
    stop('phenotype 399 has negative values.')
}
cols_399_1 <- cols_399[1:2]
cols_399_2 <- cols_399[4:5]
cols_399_3 <- cols_399[7:8]
## make df with sum of pairs of cols (include NAs)
df_399 <- data.frame(f.eid  = phen_26618$f.eid,
                     inst.1 = rowSums(phen_26618[, cols_399_1], na.rm=FALSE),
                     inst.2 = rowSums(phen_26618[, cols_399_2], na.rm=FALSE),
                     inst.3 = rowSums(phen_26618[, cols_399_3], na.rm=FALSE))
## keep most recent non-NA value
latest_399 <- apply(df_399[, 2:4], 1, get_right)
## add eid column
latest_df_399 <- data.frame(phen_26618$f.eid, latest_399)
colnames(latest_df_399) <- c('#IID', 'f.399')
## export
fwrite(latest_df_399, file = '../data/phenotypes/clean/f.399.tab',
       sep = '\t', na = 'NA', quote = FALSE)
