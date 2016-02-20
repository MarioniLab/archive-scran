countsToSE <- function(counts, spikes)
# Sets up the SummarizedExperiment(0) object from the raw count data.
# A bit easier to manage than ExpressionSet objects, I think.
#
# written by Aaron Lun
# created 17 February 2016
{
    libs <- colSums(counts)
    llibs <- log(libs)
    colData <- DataFrame(lib.size=libs, size.factor=exp(llibs - mean(llibs)))
    if (!missing(spikes)) {
        colData <- DataFrame(colData, .breakToList(spikes, "spikes"))
    } 

    out <- SummarizedExperiment(List(counts=as.matrix(counts)), colData=colData)
    return(out)
}

.breakToList <- function(y, colname) {
    ncells <- ncol(y)
    y <- as.matrix(y)
    y <- split(y, col(y))
    dim(y) <- c(ncells, 1)
    colnames(y) <- colname
    return(y)
}

setMethod("normalize", "ANY", function(object, size.factor=NULL, log=TRUE, prior.count=1) 
# Computes the normalized log-expression values.
# 
# written by Aaron Lun
# created 17 February 2016
# last modified 19 February 2016
{
    object <- as.matrix(object)
    if (is.null(size.factor)) { size.factor <- colSums(object) } 
    lsf <- log(size.factor) # Mean-centered size factors, for valid comparisons between size factor sets.
    size.factor <- exp(lsf - mean(lsf))
    cpm.default(object, lib.size=size.factor, prior.count=prior.count, log=log)
})

setMethod("normalize", "SummarizedExperiment0", function(object, ..., separate.spikes=TRUE) {
    out <- normalize(assay(object, "counts"), size.factor=object$size.factor, ...)
    assay(object, "exprs") <- out

    if (!is.null(object$spikes)) {
        if (separate.spikes) { 
            sf <- normalizeBySpikes(object)
        } else {
            sf <- object$size.factor 
        }
        out <- normalize(spikes(object, type="counts"), size.factor=sf, ...)
        if (!is.null(object$norm.spikes)) { object$norm.spikes <- NULL }
        colData(object) <- DataFrame(colData(object), .breakToList(out, "norm.spikes"))
    } 
    return(object)
})
