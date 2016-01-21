fitTechTrend <- function(spikes, size.factor=NULL, df=5, prior.count=1) 
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016    
{
    spikes <- as.matrix(spikes)
    if (is.null(size.factor)) { size.factor <- colSums(spikes) } 
    adjc <- cpm.default(spikes, lib.size=size.factor, prior.count=prior.count, log=TRUE)
    lmeans <- rowMeans(adjc)
    lvar <- apply(adjc, 1, var)

    is.okay <- lvar > 0
    kept.means <- lmeans[is.okay]
    fit <- lm(log2(lvar)[is.okay] ~ poly(kept.means, df=df))
    left.edge <- which.min(kept.means)
    right.edge <- which.max(kept.means)

    FUN <- function(x) {
        out <- predict(fit, data.frame(kept.means=x))
        out[x < kept.means[left.edge]] <- fitted(fit)[left.edge]
        out[x > kept.means[right.edge]] <- fitted(fit)[right.edge]
        return(2^out)
    }
    return(list(mean=lmeans, var=lvar, trend=FUN, 
                size.factor=size.factor, prior.count=prior.count))
}

getBioVar <- function(counts, tech.fit, size.factor=NULL)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
{
    if (is.null(size.factor)) { size.factor <- tech.fit$size.factor }
    adjc <- cpm.default(counts, lib.size=size.factor, prior.count=tech.fit$prior.count, log=TRUE) 
    lmeans <- rowMeans(adjc)
    lvar <- apply(adjc, 1, var)
    tech.var <- tech.fit$trend(lmeans)
    bio.var <- lvar - tech.var
    return(list(mean=lmeans, total=lvar, bio=bio.var, tech=tech.var))
}
