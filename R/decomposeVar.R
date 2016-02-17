setGeneric("decomposeVar", function(x, fit, ...) { standardGeneric("decomposeVar") })

setMethod("decomposeVar", c("ANY", "list"), function(x, fit, design=NULL)
# Computes the biological variability of the log-CPMs by subtracting the
# inferred technical variance from the total variance.
#
# written by Aaron Lun
# created 21 January 2016 
# last modified 17 February 2016
{
    x <- as.matrix(x)
    if (is.null(design)) { design <- fit$design }
    else if (is.na(design)) { design <- as.matrix(rep(1, ncol(x))) }
    lmeans <- rowMeans(x)
    lvar <- .estimateVariance(design, x)
    tech.var <- fit$trend(lmeans)
    bio.var <- lvar - tech.var
    return(data.frame(mean=lmeans, total=lvar, bio=bio.var, tech=tech.var))
})

setMethod("decomposeVar", c("SummarizedExperiment0", "list"), function(x, fit, ..., i="exprs") {
    if (is.null(fit$assay)) {
        fit$assay <- i
    }
    out <- decomposeVar(assay(x, i=fit$assay), fit, ...)
    return(out)   
})

testVar <- function(total, null, df, design=NULL) 
# Tests that total > null given variances estimated on 'df' degrees of freedom.
# You can also give it the design matrix directly if you can't be bothered estimating 'df'.
# Obviously there's an assumption of normality here, regarding the observations from which estimation was performed.
#
# written by Aaron Lun
# created 9 February 2016    
{
    if (missing(df)) { df <- nrow(design) - qr(design)$rank }
    pchisq(total/null*df, df=df, lower.tail=FALSE)
}

