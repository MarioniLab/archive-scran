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

setMethod("normalize", "SCESet", function(object, ..., separate.spikes=TRUE) {
    out <- normalize(assayDataElement(object, "counts"), size.factor=object$size.factor, ...) # Normalizing everything, not just spikes.

    if (separate.spikes) { 
        sf <- normalizeBySpikes(object)
        out2 <- normalize(spikes(object, type="counts"), size.factor=sf, ...)
        out[is.spike(object),] <- out2
    } 
    
    assayDataElement(object, "exprs") <- out
    return(object)
})
