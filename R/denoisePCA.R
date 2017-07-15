.denoisePCA <- function(x, technical, design=NULL, subset.row=NULL)
# Performs PCA and chooses the number of PCs to keep based on the technical noise.
# This is done on the residuals if a design matrix is supplied.
#
# written by Aaron Lun
# created 13 March 2017    
# last modified 15 July 2017
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    x <- x[subset.row,] # Might as well, need to do PCA on the subsetted matrix anyway.
    subset.row <- seq_len(nrow(x))
    all.means <- rowMeans(x)

    if (!is.null(design)) { 
        checked <- .make_var_defaults(x, fit=NULL, design=design)
        design <- checked$design
        QR <- .ranksafe_qr(design)
        
        # Computing residuals; don't set a lower bound.
        rx <- .calc_residuals_wt_zeroes(x, QR=QR, subset.row=subset.row, lower.bound=NA) 

        # Rescaling residuals so that the variance is unbiased.
        # This is necessary because variance of residuals is underestimated.
        xout <- .Call(cxx_estimate_variance, QR$qr, QR$qraux, x, subset.row - 1L)
        xvar <- xout[[2]]
        rvar <- apply(rx, 1, var)
       
        # Replacing 'x' with the scaled residuals (these shoud have a mean of zero,
        # see http://math.stackexchange.com/questions/494181/ for a good explanation).
        x <- rx * sqrt(xvar/rvar)
    }

    # Computing the technical variance sum.
    if (is.function(technical)) { 
        technical <- sum(technical(all.means))
    } else if (is.numeric(technical)) { 
        if (is.null(rownames(x))) { 
            stop("rows of 'x' should be named with gene names")
        }
        technical <- sum(technical[rownames(x)])
        if (is.na(technical)) {
            stop("missing gene names in 'technical'")
        }
    } else {
        stop("'technical' should be a function or a scalar")
    }

    # Performing PCA and discarding later PCs that add up to the technical sum.
    pcout <- prcomp(t(x)) # no scaling, otherwise technical sum isn't comparable.
    npcs <- ncol(pcout$x)
    prog.var <- pcout$sdev^2
    flipped.prog.var <- rev(prog.var)
    estimated.contrib <- cumsum(flipped.prog.var) + flipped.prog.var * (npcs:1 - 1L)
    estimated.contrib <- rev(estimated.contrib)

    below.noise <- technical > estimated.contrib
    if (any(below.noise)) { 
        to.keep <- min(which(below.noise))
    } else {
        to.keep <- npcs
    }

    return(pcout$x[,seq_len(to.keep),drop=FALSE])
} 

setGeneric("denoisePCA", function(x, ...) standardGeneric("denoisePCA"))

setMethod("denoisePCA", "matrix", .denoisePCA)

setMethod("denoisePCA", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) {
        subset.row <- .spike_subset(x, get.spikes)
    }
    out <- .denoisePCA(assayDataElement(x, assay), ..., subset.row=subset.row)
    reducedDimension(x) <- out
    return(x)
})

