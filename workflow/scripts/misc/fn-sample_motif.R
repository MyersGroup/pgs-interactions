## desc: Function to sample base of each position of a motif n times

sample_motif <- function(motifname,
                         n = 1000,
                         prob = NULL) {

    ## get Position Count Matrix for this motif
    ourmot <- hocomoco[[motifname]]$motif

    ## compute row-wise frequencies (each row sums to 1)
    ourmot[ourmot == 0] <- 1
    ourmot <- ourmot / rowSums(ourmot)

    ## for each row/position in length of motif:
    ## - sample n bases and multiply their probability
    scores <- rep(1, n)
    for (i in 1:nrow(ourmot)) {
        if (is.null(prob)) {
            scores <- scores * sample(ourmot[i,], prob = ourmot[i,], size = n, replace = TRUE)
        } else {
            scores <- scores * sample(ourmot[i,], prob = prob, size = n, replace = TRUE)
        }
    }
    return(scores)
}
