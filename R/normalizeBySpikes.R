setGeneric("normalizeBySpikes", function(x, ...) { standardGeneric("normalizeBySpikes") })

setMethod("normalizeBySpikes", "SummarizedExperiment0", function(x) 
# Uses the total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
{
    if (is.null(x$spikes)) { stop("no 'spikes' are present in 'colData'") }
    out <- log(sapply(x$spikes, FUN=sum))
    exp(out - mean(out))
})

setGeneric("spikes", function(x, ...) { standardGeneric("spikes") })

setMethod("spikes", "SummarizedExperiment0", function(x, type=c("counts", "exprs")) {
    type <- match.arg(type)
    if (type=="counts") {
        if (is.null(x$spikes)) { stop("no 'spikes' are present in 'colData'") }
        cur.assay <- do.call(cbind, x$spikes)
    } else if (type=="exprs") {
        if (is.null(x$norm.spikes)) { stop("no 'norm.spikes' are present in 'colData', run normalize() first") }
        cur.assay <- do.call(cbind, x$norm.spikes)
    }
    return(cur.assay)
})
