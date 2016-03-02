setGeneric("computeSpikeFactors", function(x, ...) { standardGeneric("computeSpikeFactors") })

setMethod("computeSpikeFactors", "SCESet", function(x) 
# Uses the total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
# last modified 29 February 2016
{
    out <- colSums(spikes(x))
    if (any(out < 1e-8)) { 
        warning("zero spike-in counts during spike-in normalization")
    } 
    out <- log(out)
    sizeFactors(x) <- exp(out - mean(out, na.rm=TRUE))
    x
})

setGeneric("spikes", function(x, ...) standardGeneric("spikes"))

setMethod("spikes", "SCESet", function(x, type=c("counts", "exprs")) {
    type <- match.arg(type)
    cur.assay <- assayDataElement(x, type)[is.spike(x),,drop=FALSE]
    return(cur.assay)
})

setGeneric("isSpike", function(x) standardGeneric("isSpike"))

is.spike <- function(x) { fData(x)$is_feature_spike }

setMethod("isSpike", "SCESet", function(x) {
    keep <- is.spike(x)
    if (is.null(keep)) { stop("set 'isSpike(x)' to identify spike-in rows") }
    return(keep)
}

setGeneric("isSpike<-", function(x, value) standardGeneric("isSpike<-"))
setReplaceMethod("isSpike", "SCESet", function(x, value) {
    if (!is.logical(value)) { stop("'isSpike(x)' must be a logical vector") }
    fData(x)$is_feature_spike <- value
    return(x) 
})

