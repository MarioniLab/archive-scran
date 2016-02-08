fitTechTrend <- function(counts, size.factor=NULL, trend=c("poly", "loess"), df=5, span=0.3, prior.count=1, design=NULL) 
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 9 February 2016
{
    counts <- as.matrix(counts)
    if (is.null(size.factor)) { size.factor <- colSums(counts) } 
    if (is.null(design)) { design <- as.matrix(rep(1, ncol(counts))) }
    adjc <- cpm.default(counts, lib.size=size.factor, prior.count=prior.count, log=TRUE)
    lmeans <- rowMeans(adjc)
    lvar <- .estimateVariance(design, adjc)

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
    return(list(mean=lmeans, var=lvar, trend=FUN, 
                size.factor=size.factor, prior.count=prior.count, design=design))
}

getBioVar <- function(counts, tech.fit, size.factor=NULL, design=NULL)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
# last modified 9 February 2016
{
    if (is.null(size.factor)) { size.factor <- tech.fit$size.factor }
    if (is.null(design)) { design <- tech.fit$design }
    adjc <- cpm.default(counts, lib.size=size.factor, prior.count=tech.fit$prior.count, log=TRUE) 
    lmeans <- rowMeans(adjc)
    lvar <- .estimateVariance(design, adjc)
    tech.var <- tech.fit$trend(lmeans)
    bio.var <- lvar - tech.var
    return(list(mean=lmeans, total=lvar, bio=bio.var, tech=tech.var))
}

.estimateVariance <- function(X, y) {
    fit <- lm.fit(x=X, y=t(y))
    return(colMeans(fit$effects[-seq_len(fit$rank),]^2))
}
