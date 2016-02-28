setGeneric("normalizeBySpikes", function(x, ...) { standardGeneric("normalizeBySpikes") })

setMethod("normalizeBySpikes", "SCESet", function(x) 
# Uses the total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
{
    out <- colSums(spikes(x))
    if (any(out < 1e-8)) { 
        warning("zero spike-in counts during spike-in normalization")
    } 
    out <- log(out)
    exp(out - mean(out, na.rm=TRUE))
})

setGeneric("spikes", function(x, ...) standardGeneric("spikes"))

setMethod("spikes", "SCESet", function(x, type=c("counts", "exprs")) {
    type <- match.arg(type)
    cur.assay <- assayDataElement(x, type)[is.spike(x),,drop=FALSE]
    return(cur.assay)
})

setGeneric("isSpike", function(x) standardGeneric("isSpike"))

setMethod("isSpike", "SCESet", is.spike)
is.spike <- function(x) { 
    keep <- fData(x)$is_feature_spike 
    if (is.null(keep)) { stop("set 'isSpike(x)' to identify spike-in rows") }
    return(keep)
}

setGeneric("isSpike<-", function(x, value) standardGeneric("isSpike<-"))
setReplaceMethod("isSpike", "SCESet", function(x, value) {
    fData(x)$is_feature_spike <- value
    return(x) 
})

