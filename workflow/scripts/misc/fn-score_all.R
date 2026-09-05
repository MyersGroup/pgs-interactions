## desc: Function to score a vector of sequences for disruption of a set of TF binding motifs

score_all <- function(sequencemat,
                      hocomoco,
                      namesho) {

    ## make skeleton output matrix
    maxscoremat <- matrix(nrow = nrow(sequencemat), ncol = length(namesho))
    colnames(maxscoremat) <- namesho
    maxscorematcomp <- maxscoremat

    minr <- 10000
    maxr <- 1

    for (i in 1:length(namesho)) {

        ## get Position Count Matrix for this motif
        motifname <- namesho[i]
        ourmot <- hocomoco[[motifname]]$motif

        ## find minimum number of rows in PCM across all motifs
        if (nrow(ourmot) < minr) {
            minr <- nrow(ourmot)
        }

        ## find maximum number of rows in PCM across all motifs
        if (nrow(ourmot) > maxr) {
            maxr <- nrow(ourmot)
        }
    }

    ## now set up sequences
    left <- 1
    right <- maxr
    ## get maximum number of "windows" to test (based on shortest sequence motif)
    leftmax <- ncol(sequencemat) - (minr - 1)

    ## make matrix with:
    ## - one group of rows per window
    ## - within each window, one row per SNP
    ## - one column per base, where number of cols is length of longest motif
    ## - contents are chunks of sequencemat, in a moving window left-to-right way
    ##   (width of window is length of longest motif; once it no longer fits it starts shrinking down to length of shortest motif)
    ##   (contents always start from the left; once length of sequences is shorter than ncols, have NAs)
    tempmat <- matrix(nrow = nrow(sequencemat) * leftmax, ncol = maxr)
    for (i in 1:leftmax) {
        ll <- left + i - 1
        rr <- right + i - 1
        if (rr > ncol(sequencemat)) {
            rr <- ncol(sequencemat)
        }
        tempmat[(1 + nrow(sequencemat) * (i - 1)):(nrow(sequencemat) * i), 1:(rr - ll + 1)] <- sequencemat[, ll:rr]
    }

    ## score all motifs in hocomoco
    for(c in 1:length(namesho)) {

        motifname <- namesho[c]
        ourmot <- hocomoco[[motifname]]$motif
        ourmot <- cbind(ourmot, 1)  # add fifth column of 1s
        ourmot[ourmot == 0] <- 1
        ourmot <- ourmot / rowSums(ourmot)  # compute row-wise frequencies (each row sums to 1)

        ## start with first window, from first position in sequences with width = length of motif
        scores <- rep(1, length = nrow(tempmat))
        for (k in 1:nrow(ourmot)) {  # loop over length of current motif
            vv <- ourmot[k,]  # get frequencies
            ## update scores vector by multiplying current score by frequency of each base at current position
            ## if base at this position has high frequency in our motif, score will be higher
            scores <- scores * vv[tempmat[, k]]
        }
        scores <- scores[!is.na(scores)]
        maxscores <- scores[1:nrow(sequencemat)]  # keep only first window / group of observations from tempmat

        ## move window along length of sequencemat
        left <- 1
        right <- nrow(ourmot)
        leftmax <- ncol(sequencemat) - (right - left)  # right-most starting point for moving window
        ## find maximum frequency starting from each position in sequence
        for (i in 2:leftmax) {
            temp <- scores[((i - 1) * nrow(sequencemat) + 1):(i * nrow(sequencemat))]
            ## if new position scores higher than current, keep that
            maxscores[temp > maxscores] <- temp[temp > maxscores]
        }

        maxscoremat[, c] <- maxscores
    }

    ## score reverse complement of each sequence
    for (c in 1:length(namesho)) {

        motifname <- namesho[c]
        ourmot <- hocomoco[[motifname]]$motif
        ourmot <- cbind(ourmot, 1)
        ourmot[ourmot == 0] <- 1
        ourmot <- ourmot / rowSums(ourmot)

        ourmot <- ourmot[nrow(ourmot):1,]  # reverse order of rows
        ourmot <- ourmot[, c(4:1, 5)]  # reverse order of columns (but keep last column of 1s)

        scores <- rep(1, length = nrow(tempmat))
        for (k in 1:nrow(ourmot)) {
            vv <- ourmot[k,]
            scores <- scores * vv[tempmat[, k]]
        }
        scores <- scores[!is.na(scores)]
        maxscores <- scores[1:nrow(sequencemat)]

        left <- 1
        right <- nrow(ourmot)
        leftmax <- ncol(sequencemat) - (right - left)
        for (i in 2:leftmax) {
            temp <- scores[((i - 1) * nrow(sequencemat) + 1):(i * nrow(sequencemat))]
            maxscores[temp > maxscores] <- temp[temp > maxscores]
        }

        maxscorematcomp[, c] <- maxscores
    }

    ## keep highest of regular / reverse complement scores
    maxscoremat[maxscorematcomp > maxscoremat] <- maxscorematcomp[maxscorematcomp > maxscoremat]

    return(maxscoremat)
}
