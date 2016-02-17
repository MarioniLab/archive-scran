setGeneric("normalizeBySpikes", function(x, ...) { standardGeneric("normalizeBySpikes") })

setMethod("normalizeBySpikes", "ANY", function(x, is.spike) 
# Uses the total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
{
    if (!missing(is.spike)) { x <- x[is.spike,] }
    colSums(x)
})

setMethod("normalizeBySpikes", "SummarizedExperiment0", function(x, i="counts") {
    normalizeBySpikes(assay(x, i=i), is.spike=mcols(x)$spike)
})
