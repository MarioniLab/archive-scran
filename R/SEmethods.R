.getUsedMatrix <- function(x, assay="counts", get.spikes=FALSE) {
    cur.mat <- assayDataElement(x, assay)
    if (!get.spikes) {
        nokeep <- is.spike(x)
        if (!is.null(nokeep) && any(nokeep)) { 
            cur.mat <- cur.mat[!nokeep,,drop=FALSE]
        }
    }
    return(cur.mat)
}

.spikeSubset <- function(x, get.spikes) {
    if (!get.spikes) {
        nokeep <- is.spike(x)
        if (!is.null(nokeep) && any(nokeep)) {
            return(!nokeep)
        }
    } 
    return(NULL)
}

.subset_to_index <- function(subset, names) {
    if (is.logical(subset)) { 
        subset <- which(subset)
    } else if (is.character(subset)) {
        subset <- match(subset, names)
        if (any(is.na(subset))) { 
            stop("missing names in subset vector")
        }
    }
    as.integer(subset)
}
