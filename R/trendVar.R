setGeneric("trendVar", function(x, ...) standardGeneric("trendVar"))

setMethod("trendVar", "ANY", function(x, trend=c("poly", "loess"), df=5, span=0.3, prior.count=1, design=NULL, weight.by.bin=5) 
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 21 February 2016
{
    x <- as.matrix(x)
    lmeans <- rowMeans(x)
    if (is.null(design)) { design <- .interceptModel(ncol(x)) }
    lvar <- .estimateVariance(design, x)

    is.okay <- lvar > 1e-8
    kept.means <- lmeans[is.okay]
    llvar <- log2(lvar)[is.okay]
    trend <- match.arg(trend)

    # Fitting the trend with weights.
    bins <- seq(from=min(kept.means), to=max(kept.means), length.out=weight.by.bin+1L)[-1]
    bin.id <- findInterval(kept.means, bins, rightmost.closed=TRUE) + 1
    bin.pop <- tabulate(bin.id)
    weights <- 1/bin.pop[bin.id]
    if (trend=="loess") { 
        fit <- loess(llvar ~ kept.means, span=span, degree=1, weights=weights)
    } else if (trend=="poly") {
        fit <- lm(llvar ~ poly(kept.means, df=df), weights=weights)
    } 

    left.edge <- which.min(kept.means)
    right.edge <- which.max(kept.means)
    FUN <- function(x) {
        out <- predict(fit, data.frame(kept.means=x))
        out[x < kept.means[left.edge]] <- fitted(fit)[left.edge]
        out[x > kept.means[right.edge]] <- fitted(fit)[right.edge]
        return(2^out)
    }
    return(list(mean=lmeans, var=lvar, trend=FUN, prior.count=prior.count, design=design))
})

.interceptModel <- function(ncells) {
    as.matrix(rep(1, ncells)) 
}

setMethod("trendVar", "SCESet", function(x, ..., use.spikes=TRUE) {
    if (use.spikes) {
        cur.assay <- spikes(x, "norm_exprs")
    } else {
        cur.assay <- .getUsedMatrix(x, "norm_exprs")
    }
    out <- trendVar(cur.assay, ...)
    return(out)
})

.estimateVariance <- function(X, y) {
    fit <- lm.fit(x=X, y=t(y))
    return(colMeans(fit$effects[-seq_len(fit$rank),]^2))
}

