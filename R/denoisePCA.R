.denoisePCA <- function(x, technical, design=NULL, subset.row=NULL,
                        value=c("pca", "n", "lowrank"), min.rank=5, max.rank=100, 
                        preserve.dim=FALSE)
# Performs PCA and chooses the number of PCs to keep based on the technical noise.
# This is done on the residuals if a design matrix is supplied.
#
# written by Aaron Lun
# created 13 March 2017    
# last modified 15 July 2017
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    checked <- .make_var_defaults(x, fit=NULL, design=design)
    QR <- .ranksafe_qr(checked$design)
    stats <- .Call(cxx_fit_linear_model, QR$qr, QR$qraux, x, subset.row - 1L, FALSE)
    all.means <- stats[[1]]
    all.var <- stats[[2]]

    # Filtering out genes with negative biological components.
    tech.var <- technical(all.means)
    keep <- all.var > tech.var
    use.rows <- subset.row[keep]
    all.means <- all.means[keep]
    all.var <- all.var[keep]
    tech.var <- tech.var[keep]
    technical <- sum(tech.var)

    if (!is.null(design)) {  
        # Computing residuals; don't set a lower bound.
        # Note that this function implicitly subsets by subset.row.
        rx <- .calc_residuals_wt_zeroes(x, QR=QR, subset.row=use.rows, lower.bound=NA) 

        # Rescaling residuals so that the variance is unbiased.
        # This is necessary because variance of residuals is underestimated.
        rvar <- apply(rx, 1, var)
       
        # Replacing 'x' with the scaled residuals (these should already have a mean of zero,
        # see http://math.stackexchange.com/questions/494181/ for a good explanation).
        y <- rx * sqrt(all.var/rvar)
    } else {
        y <- x[use.rows,,drop=FALSE] - all.means
    }

    # Performing SVD to get the variance of each PC, and choosing the number of PCs to keep.
    y <- t(y)
    svd.out <- svd(y, nu=0, nv=0)
    var.exp <- svd.out$d^2/(ncol(x) - 1)
    to.keep <- .get_npcs_to_keep(var.exp, technical)
    to.keep <- min(max(to.keep, min.rank), max.rank)
    
    # Figuring out what value to return; the number of PCs, the PCs themselves, or a denoised low-rank matrix.
    value <- match.arg(value)
    if (value=="n") {
        return(to.keep)
    } else if (value=="pca") {
        pc.out <- prcomp(y, rank.=to.keep, scale.=FALSE, center=FALSE)
        return(pc.out$x)
    } else if (value=="lowrank") {
        more.svd <- La.svd(y, nu=to.keep, nv=to.keep)
        denoised <- more.svd$u %*% (more.svd$d[seq_len(to.keep)] * more.svd$vt) 
        denoised <- t(denoised) + all.means

        # Returning as a the full matrix (or that subsetted with subset.row)
        # where the discarded genes are set to zero.
        output <- x
        output[] <- 0
        if (preserve.dim) { 
            output[use.rows,] <- denoised
        } else {
            output <- output[subset.row,,drop=FALSE]
            output[match(use.rows, subset.row),] <- denoised
        }
        return(output)
    }
} 

.get_npcs_to_keep <- function(var.exp, tech.var) 
# Discarding PCs until we get rid of as much technical noise as possible
# while preserving the biological signal. This is done by assuming that 
# the biological signal is fully contained in earlier PCs, such that we 
# discard the later PCs until we account for 'tech.var'.
{
    npcs <- length(var.exp)
    flipped.var.exp <- rev(var.exp)
    estimated.contrib <- cumsum(flipped.var.exp) + flipped.var.exp * (npcs:1 - 1L)

    below.noise <- tech.var > estimated.contrib
    if (any(below.noise)) { 
        to.keep <- npcs - max(which(below.noise))
    } else {
        to.keep <- npcs
    }
    return(to.keep)
}

setGeneric("denoisePCA", function(x, ...) standardGeneric("denoisePCA"))

setMethod("denoisePCA", "ANY", .denoisePCA)

setMethod("denoisePCA", "SingleCellExperiment", 
          function(x, ..., subset.row=NULL, value=c("pca", "n", "lowrank"), 
                   assay.type="exprs", get.spikes=FALSE) {

    subset.row <- .SCE_subset_genes(subset.row=subset.row, x=x, get.spikes=get.spikes)
    out <- .denoisePCA(assay(x, i=assay.type), ..., value=value, subset.row=subset.row, preserve.dim=TRUE)

    value <- match.arg(value) 
    if (value=="pca"){ 
        reducedDim(x, "PCA") <- out
    } else if (value=="n") {
        metadata(x)$denoised.npcs <- out
    } else if (value=="lowrank") {
        assay(x, i="lowrank") <- out
    }
    return(x)
})

