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
