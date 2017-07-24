setGeneric("computeSpikeFactors", function(x, ...) { standardGeneric("computeSpikeFactors") })

setMethod("computeSpikeFactors", "SingleCellExperiment", function(x, type=NULL, sf.out=FALSE, general.use=TRUE) 
# Uses the mean-centred total of spike-in transcripts as the size factor.
#
# written by Aaron Lun
# created 17 February 2016
# last modified 24 July 2017
{
    out <- colSums(spikes(x, type=type))
    if (any(out < 1e-8)) { 
        warning("zero spike-in counts during spike-in normalization")
    } 
    sf <- out/mean(out)

    # Returning size factors directly.
    if (sf.out) {
        return(sf)
    }

    # Saving size factors for general use, or for specific use by one (or all) of the spike-in sets.
    if (general.use) {
        sizeFactors(x) <- sf
    } 
    if (is.null(type)) {
        type <- spikeNames(x)
    }
    for (f in type) {
        sizeFactors(x, type=f) <- sf
    }        
    x
})

