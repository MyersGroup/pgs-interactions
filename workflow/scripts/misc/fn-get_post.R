## desc: Get updated posterior probabilities using a single annotation
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
##       - alpharange: sequence of alpha values that now give weights on each
##       column of factormat


get_post <- function(likelihood_ratios,
                     factormat,
                     snpinfo,
                     mode = "U",
                     alphavalues) {

    ## if factormat is matrix of relative weights
    if (mode == "R") factormat <- factormat - 1

    ## extract drivers and tags from snpinfo
    drivers <- snpinfo[, "hit_ID"]
    tags <- snpinfo[, "tag_ID"]

    ## remove drivers/tags with missing LRs or annotations
    cond <- (!is.na(likelihood_ratios) & !is.na(rowSums(factormat)))
    tags <- tags[cond]
    drivers <- drivers[cond]
    likelihood_ratios <- likelihood_ratios[cond]
    factormat <- factormat[cond, , drop = FALSE]

    ## for each driver SNP, sum its LRs
    aa <- table(drivers)  # this orders drivers alphabetically rather than by position
    indices <- match(drivers, names(aa))  # this allows restoring positional order
    tots <- rep(0, length(aa))
    names(tots) <- names(aa)
    for (i in 1:length(indices)) {  # loop through driver SNPs in standard (positional) order
        if (!is.na(likelihood_ratios[i])) {
            tots[indices[i]] <- tots[indices[i]] + likelihood_ratios[i]
        }
    }
    totsold <- tots


    ## associate each element of tags to a unique row (driver): choose driver with fewest number of other tags
    counts <- aa[indices]  # number of tags for each driver but with repeated rows (in positional order)
    main <- data.frame(tag    = tags,
                       driver = drivers,
                       count  = as.numeric(counts),
                       factormat,
                       lr = likelihood_ratios)
    main <- main %>% arrange(count)  # sort tags by number of other tags for that driver snp (increasing)

    untags <- unique(tags)
    tagrows <- rep(0, length(untags))
    tagrows <- match(untags, main$tag)  # find first (and therefore smallest) match (i.e., implicitly find driver snp with fewest number of other tags)
    tagrows <- sort(tagrows)  # sort so we have positional ordering again

    newtags      <- main[tagrows, 1]
    newdrivers   <- main[tagrows, 2]
    newfactormat <- as.matrix(main[tagrows, 4:(4 + ncol(factormat) - 1)], rownames = FALSE)
    newlikelihood_ratios <- main[tagrows, 4 + ncol(factormat)]
    
   
    ## sum likelihoods for each driver as before
    ## lrs may now not sum to 1 for drivers that 'lost' a tag to another driver
    aa <- table(newdrivers)
    indices <- match(newdrivers, names(aa))
    tots <- rep(0, length(aa))
    names(tots) <- names(aa)
    for (i in 1:length(indices)) {
        if (!is.na(newlikelihood_ratios[i])) {
            tots[indices[i]] <- tots[indices[i]] + newlikelihood_ratios[i]
        }
    }

    ## 'old' posterior terms (computed based directly on input, removing only some tags so each has a unique driver)
    oldpostprobs <- newlikelihood_ratios / tots[indices]


    ## re-weight posterior probabilities by weights in factormat
    terms <- t(1 + t(newfactormat) * (alphavalues - 1))
    overallterms <- apply(terms, 1, prod)
    overallterms <- overallterms * newlikelihood_ratios
    totsnew <- rep(0, length(aa))
    names(totsnew) <- names(aa)
    for (i in 1:length(indices)) {
        totsnew[indices[i]] <- totsnew[indices[i]] + overallterms[i]
    }

    postprobs <- overallterms / totsnew[indices]


    ## return new posteriors for unique set of tags
    postprobs_df <- data.frame(hit_ID = newdrivers,
                               tag_ID = newtags,
                               post_prob = postprobs)
    out <- left_join(snpinfo[ ,c("hit_ID", "tag_ID")], postprobs_df, by = c("hit_ID", "tag_ID"))
 
    return(out)
}
