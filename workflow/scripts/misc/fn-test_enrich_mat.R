## desc: Function to test enrichment of a set of SNPs for set of annotations
##
##       Takes as input:
##       - likelihood_ratios: numeric vector of posterior probabilities for
##       set of SNPs of interest (tag SNPs in snpinfo below)
##       - factormat: 0/1 matrix with one row per SNP and one column per source
##       of annotation;
##       - snpinfo: data.frame with info on driver and tag SNPs (must have
##       "hit_ID" and "tag_ID" cols)
##       - mode: whether factormat is a matrix of relative (R) or
##       unconditional (U) weights
##       - alpharange: sequence of alpha values to compute model likelihoods for


library(dplyr)


test_enrich_mat <- function(likelihood_ratios,
                            factormat,
                            snpinfo,
                            mode = "U",
                            alpharange = NULL) {

    ## set alpharange to default if not provided
    if(is.null(alpharange)) alpharange <- seq(0.1, 4, 0.05)

    ## if factormat is matrix of relative weights
    if(mode == "R") factormat <- factormat - 1

    ## extract drivers and tags from snpinfo
    drivers <- snpinfo[, "hit_ID"]
    tags    <- snpinfo[, "tag_ID"]

    ## remove drivers/tags with missing LRs or annotations
    cond <- (!is.na(likelihood_ratios) & !is.na(rowSums(factormat)))
    tags <- tags[cond]
    drivers <- drivers[cond]
    likelihood_ratios <- likelihood_ratios[cond]
    factormat <- factormat[cond, , drop = FALSE]


    ## for each driver SNP, sum its LRs (posterior probabilities)
    aa <- table(drivers)  # this orders drivers alphabetically rather than by position
    indices <- match(drivers, names(aa))  # this allows restoring positional order
    tots <- rep(0, length(aa))
    names(tots) <- names(aa)
    for (i in 1:length(indices)) {  # loop through driver SNPs in standard (positional) order
        if (!is.na(likelihood_ratios[i])) {
            tots[indices[i]] <- tots[indices[i]] + likelihood_ratios[i]
        }
    }
    totsold <- tots  # sums to 1 (LRs are ratio of BF for each tag to sum of BFs for corresponding driver)

    
    ## associate each element of tags to a unique row (driver): choose driver with fewest number of other tags
    counts <- aa[indices]  # number of tags for each driver but with repeated rows (in positional order)
    main <- data.frame(tag    = tags,
                       driver = drivers,
                       count  = as.numeric(counts),
                       factormat,
                       lr = likelihood_ratios)
    main <- main %>% arrange(count)  # sort tags by number of other tags for that driver SNP (increasing)

    untags <- unique(tags)
    tagrows <- rep(0, length(untags))
    tagrows <- match(untags, main$tag)  # find first (and therefore smallest) match (i.e., implicitly find driver SNP with fewest number of other tags)
    tagrows <- sort(tagrows)  # sort so we have positional ordering again

    newtags      <- main[tagrows, 1]
    newdrivers   <- main[tagrows, 2]
    newfactormat <- as.matrix(main[tagrows, 4:(4 + ncol(factormat) - 1)], rownames = FALSE)
    newlikelihood_ratios <- main[tagrows, 4 + ncol(factormat)]

    
    ## sum likelihoods for each driver as before
    ## LRs may now not sum to 1 for drivers that 'lost' a tag to another driver
    aa <- table(newdrivers)
    indices <- match(newdrivers, names(aa))
    tots <- rep(0, length(aa))
    names(tots) <- names(aa)
    for (i in 1:length(indices)) {
        if (!is.na(newlikelihood_ratios[i])) {
            tots[indices[i]] <- tots[indices[i]] + newlikelihood_ratios[i]
        }
    }
    
    ## keep those that are convincingly a single region:
    ## if one of a driver's tags accounted for more than 5% of its total
    ## likelihood and was moved to a different driver, this suggests
    ## the two drivers' regions overlap
    keepset <- which(tots > 0.95 * totsold[names(tots)])


    ## for each values of alpha, compute model likelihood
    lhoodsu <- matrix(ncol = 0, nrow = ncol(factormat))
    alphas  <- vector(length = 0)
    totsnew <- matrix(0, ncol = length(aa), nrow = ncol(newfactormat))
    denom   <- matrix(0, ncol = length(aa), nrow = ncol(newfactormat))
    newfactormat <- t(newfactormat)

    for (alpha in alpharange) {

        ## alpha represents whether annotations should increase (alpha > 1) or decrease (alpha < 1) likelihood of each tag
        v1 <- 1 + (alpha - 1) * newfactormat
        v1[v1 < 0] <- 0.001

        ## multiply weights v1 by posterior probability for each SNP (separately for each type of annotation)
        newbg <- t(t(v1) * newlikelihood_ratios)

        ## totals for each driver
        totsnew <- totsnew * 0
        denom <- denom * 0
        colnames(totsnew) <- aa

        ## calculate numerator likelihood of data summed over weighted SNPs (sum PPs for each driver)
        for (i in 1:length(indices)) totsnew[, indices[i]] <- totsnew[, indices[i]] + newbg[, i]

        ## calculate normalising constant as conditioning on one of SNPs being causal (sum annotation weights for each driver)
        for (i in 1:length(indices)) denom[, indices[i]] <- denom[, indices[i]] + v1[, i]
        ## so denominator is region-specific conditional on one of these sites being causal

        ## for driver SNPs that convincingly form a single region from above
        ## use denominator summed over all sites as probability a site is causal
        ## so prior is we have all SNPs concatenated then choose one proportional to its weight
        ## likelihood is sum over all SNPs for region then product across regions
        temp <- totsnew[, keepset, drop = FALSE] / rowSums(denom[, keepset, drop = FALSE])
        l2 <- rowSums(log(temp))
        
        lhoodsu <- cbind(lhoodsu, l2)
        alphas <- c(alphas, alpha)
    }

    ## calculate some basic statistics
    ## unconditional model: causal snps are drawn from overall pool of SNPs independently with probabilities proportional to fitted weights
    vv2b <- apply(lhoodsu, 1, which.max)
    psu <- pchisq(2 * (lhoodsu[cbind(1:nrow(lhoodsu), vv2b)] - lhoodsu[, which(alphas == 1)]), lower.tail = FALSE, df = 1)  # LR test
    oru <- alphas[vv2b]  # alpha values that maximise the likelihood for each kind of annotation

    ## return all alphas, likelihoods, p-values, and alpha values that maximise likelihood for each kind of annotation
    return(list(alphas  = alphas,
                lhoodsu = lhoodsu,
                psu = psu,
                oru = oru))
}
