## function to build score
score <- function(sequencemat,
                  motifname,
                  transformed = FALSE,
                  comp = FALSE) {

    if (transformed == FALSE) {
        sequencemat[sequencemat == "A"] <- 1
        sequencemat[sequencemat == "C"] <- 2
        sequencemat[sequencemat == "G"] <- 3
        sequencemat[sequencemat == "T"] <- 4
        sequencemat <- matrix(as.double(sequencemat), nrow = nrow(sequencemat))
    }

    ourmot <- hocomoco[[motifname]]$motif

    ourmot <- cbind(ourmot, 1)
    ourmot[ourmot == 0] <- 1
    ourmot <- ourmot / rowSums(ourmot)

    if (comp == TRUE) {
        ourmot <- ourmot[nrow(ourmot):1,]
        ourmot <- ourmot[, c(4:1, 5)]
    }

    res <- matrix(1, nrow = nrow(sequencemat), ncol = nrow(ourmot))
    starts <- (31 - nrow(ourmot) + 1):31  # window with width=length of motif, must overlap central position
    ends <- 31:(31 + nrow(ourmot) - 1)
    curscores <- matrix(nrow = nrow(sequencemat), ncol = nrow(ourmot))
    index <- matrix(nrow = nrow(sequencemat), ncol = 2)
    index[, 1] <- 1
    for (curpos in 1:ncol(res)) {
        res[, curpos] <- 1
        temp <- sequencemat[, starts[curpos]:ends[curpos]]  # extract sequence with same length as motif
        for (curcol in 1:ncol(curscores)) {
            index[, 1] <- curcol
            index[, 2] <- temp[, curcol]  # matrix with first col=position in motif; second col=nucleotide at that position for each sequence
            res[, curpos] <- res[, curpos] * ourmot[index]  # multiply current score by weight from that nucleotide
        }
    }
    return(res)
}
