setGeneric("normalizeBySpikes", function(x, ...) { standardGeneric("normalizeBySpikes") })

setMethod("normalizeBySpikes", "SCESet", function(x) 
# Uses the total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
{
    out <- log(colSums(spikes(x)))
    if (any(is.na(out))) { 
        warning("zero spike-in counts during spike-in normalization")
    }
    exp(out - mean(out, na.rm=TRUE))
})

setGeneric("spikes", function(x, ...) standardGeneric("spikes"))

setMethod("spikes", "SCESet", function(x, type=c("counts", "norm_exprs")) {
    type <- match.arg(type)
    cur.assay <- assayDataElement(x, type)[is.spike(x),,drop=FALSE]
    return(cur.assay)
})

is.spike <- function(x) { 
    keep <- fData(x)$is_feature_spike 
    if (is.null(keep)) { stop("'is_feature_spike' not set to identify spike-in rows") }
    return(keep)
}
