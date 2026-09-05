## desc: Function to check overlap between a vector of positions and a set of intervals
##
##       Takes as input a matrix object (df) with two columns for:
##       - start position of interval
##       - end position of interval
##       and a vector (v) with the single physical positions whose overlap with
##       intervals of interest is to be checked.
##
##       Returns vector of same length as v where the elements are the row index
##       in df of the interval containing the corresponding element of v, and 0
##       if no match is found.


get_overlap <- function(df, v) {

    ## make output vector with length of matrix whose intervals are to be checked
    ## for overlap with second matrix representing annotations
    out <- vector("numeric", length = length(v))

    ## proceed row by row
    for (i in 1:length(v)) {

        ## check if there are any intervals in df that contain the *starting* position for this row
        if (length(which(v[i] >= df[, 1] & v[i] <= df[, 2])) > 0) {

            ## output indices of intervals in df inside wich the current position in df fits
            out_indices <- which(v[i] >= df[, 1] & v[i] <= df[, 2])

            ## check that there is overlap with *only one* interval
            if (length(out_indices) > 1) {
                stop(paste("Position in row", i, "overlaps with more than one interval."))
            } else {
                out[i] <- out_indices
            }
        }
    }

    return(out)
}



## return TRUE/FALSE vector of whether each element in v overlaps any interval in df
get_overlap_simple <- function(df, v) {

    ## make output vector with length of matrix whose intervals are to be checked
    ## for overlap with second matrix representing annotations
    out <- rep(FALSE, length = length(v))

    ## proceed row by row
    for (i in 1:length(v)) {

        ## check if there are any intervals in df that contain the *starting* position for this row
        if (length(which(v[i] >= df[, 1] & v[i] <= df[, 2])) > 0) {

            out[i] <- TRUE

        }
    }

    return(out)
}
