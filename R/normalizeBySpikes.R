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

