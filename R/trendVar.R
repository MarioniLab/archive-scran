setGeneric("trendVar", function(x, ...) { standardGeneric("trendVar") })

setMethod("trendVar", "ANY", function(x, trend=c("poly", "loess"), df=5, span=0.3, prior.count=1, design=NULL) 
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 17 February 2016
{
    x <- as.matrix(x)
    if (is.null(design)) { design <- as.matrix(rep(1, ncol(x))) }
    lmeans <- rowMeans(x)
    lvar <- .estimateVariance(design, x)

    is.okay <- lvar > 1e-8
    kept.means <- lmeans[is.okay]
    llvar <- log2(lvar)[is.okay]
    trend <- match.arg(trend)
    if (trend=="loess") { 
        fit <- loess(llvar ~ kept.means, span=span, degree=1)
    } else {
        fit <- lm(llvar ~ poly(kept.means, df=df))
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

setMethod("trendVar", "SummarizedExperiment0", function(x, ..., use.spikes=TRUE, i="exprs") {
    if (use.spikes) {
        if (is.null(x$norm.spikes)) { stop("no 'norm.spikes' are present in 'colData'") }
        cur.assay <- do.call(cbind, x$norm.spikes)
    } else {
        cur.assay <- assay(x, i=i)
    }
    out <- trendVar(cur.assay, ...)
    out$assay <- i
    return(out)
})

.estimateVariance <- function(X, y) {
    fit <- lm.fit(x=X, y=t(y))
    return(colMeans(fit$effects[-seq_len(fit$rank),]^2))
}

