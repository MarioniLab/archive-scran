.trend_var <- function(x, parametric=TRUE, method=c("loess", "spline", "semiloess"), 
                       span=0.3, family="symmetric", degree=1, df=4, start=NULL, 
                       var.model=c("f", "normal"), design=NULL, subset.row=NULL)
# Fits a polynomial trend to the technical variability of the log-CPMs,
# against their abundance (i.e., average log-CPM).
# 
# written by Aaron Lun
# created 21 January 2016
# last modified 11 July 2017
{
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    checked <- .make_var_defaults(x, fit=NULL, design=design)
    design <- checked$design
    QR <- .ranksafe_qr(design)

    lout <- .Call(cxx_estimate_variance, QR$qr, QR$qraux, x, subset.row - 1L)
    means <- lout[[1]]
    vars <- lout[[2]]
    names(means) <- names(vars) <- rownames(x)[subset.row]

    is.okay <- vars > 1e-8
    kept.vars <- vars[is.okay]
    kept.means <- means[is.okay]

    # Switching from old-school 'semiloess'.
    method <- match.arg(method) 
    if (method=="semiloess") {
        warning("'semiloess' is deprecated, use parametric=TRUE and method='loess' instead")
        method <- "loess"
        parametric <- TRUE
    }
   
    # Fitting a parametric curve to try to flatten the shape.
    if (parametric) { 
        if (length(kept.vars) <= 3L) {
            stop("need at least 4 values for non-linear curve fitting")
        } 

        if (is.null(start)) start <- .get_nls_starts(kept.vars, kept.means)
        init.fit <- nls(kept.vars ~ (a*kept.means)/(kept.means^n + b), start=start,
                        control=nls.control(warnOnly=TRUE, maxiter=500), algorithm="port",
                        lower=list(a=0, b=0, n=1))

        to.fit <- log(kept.vars/fitted(init.fit))
        SUBSUBFUN <- function(x) { predict(init.fit, data.frame(kept.means=x)) }
    } else {
        to.fit <- log(kept.vars)
        SUBSUBFUN <- function(x) { 1 }
    } 
    
    # Fitting loess or splines to what's left after removing the shape.
    if (method=="loess") { 
        after.fit <- loess(to.fit ~ kept.means, degree=degree, family=family, span=span)
        trend.df <- after.fit$enp
    } else {
        after.fit <- MASS::rlm(to.fit ~ ns(kept.means, df=df))
        trend.df <- df+1
    }
    SUBFUN <- function(x) { exp(predict(after.fit, data.frame(kept.means=x))) * SUBSUBFUN(x) }

    # Modelling the variability around the trend, including the scale shift from estimating the mean of logs.
    leftovers <- kept.vars/SUBFUN(kept.means)
    var.model <- match.arg(var.model)
    if (var.model=="f") { 
        f.fit <- fitFDistRobustly(leftovers, df1=nrow(design) - ncol(design))
        f.df2 <- f.fit$df2
     
        # We don't just want the scaling factor, we want the scaled mean of the F-distribution (see explanation below).
        rescale <- f.fit$scale
        if (is.infinite(f.df2)) { 
            # do nothing. 
        } else if (f.df2 > 2) {
            rescale <- rescale * f.df2/(f.df2 - 2) 
        } else {
            warning("undefined expectation for small df2")
        }
        extra <- list(df=f.df2)
    } else {
        V <- var(leftovers)
        soln <- uniroot(function(x) exp(2*x)*x - V, interval=c(0, V))
        rescale <- exp(soln$root)
        extra <- list(s2=soln$root, df=trend.df)
    }
    
    # Creating a predictive function, with special behaviour at the ends.
    left.edge <- min(kept.means)
    left.val <- SUBFUN(left.edge) * rescale
    right.edge <- max(kept.means)
    right.val <- SUBFUN(right.edge) * rescale

    FUN <- function(x) {
        out <- SUBFUN(x) * rescale
        lower <- x < left.edge
        out[lower] <- x[lower] * (left.val/left.edge) # Gradient towards zero.
        out[x > right.edge] <- right.val
        return(out)
    }
    return(list(mean=means, var=vars, trend=FUN, 
                .internal=c(list(start=start, design=design, var.model=var.model), extra)))
}

.get_nls_starts <- function(vars, means, grad.prop=0.5, grid.length=100, grid.max=10) {
    lvars <- log2(vars)
    fit <- loess(lvars ~ means)

    # Getting a rough peak location from the fitted values.
    toppt <- which.max(fitted(fit))
    top.x <- means[toppt]
    top.y <- 2^fitted(fit)[toppt]

    # Getting the initial gradient, from all points halfway to the peak.
    min.x <- min(means)
    lower.set <- means < min.x + (top.x - min.x)*grad.prop
    lower.vars <- vars[lower.set]
    lower.means <- means[lower.set]
    line.fit <- lm(lower.vars ~ lower.means)
    grad <- coef(line.fit)[2]

    # If 'n' is known, the others can be solved by knowing the gradient near zero (grad=a/b), 
    # knowing the location of the peak (top.x^n=b/(n-1)) and the height of the peak (top.y=a*top.x*(n-1)/nb).
    # The last could be used to solve for 'n', but this seems unstable, so we do a grid search instead.
    getB <- function(n) (n-1)*top.x^n
    getA <- function(n) getB(n) * grad

    # Getting a rough 'n' estimate with a grid search.
    grid.pts <- seq(from=0, to=grid.max, length.out=grid.length)
    grid.ss <- sapply(grid.pts, FUN=function(i) {
        n <- 2^i
        resid <- vars - (getA(n)*means)/(means^n + getB(n))
        sum(resid^2)
    })

    n <- 2^grid.pts[which.min(grid.ss)]
    b <- unname(getB(n))
    a <- unname(getA(n))
    list(n=n, b=b, a=a)
}

setGeneric("trendVar", function(x, ...) standardGeneric("trendVar"))

setMethod("trendVar", "ANY", .trend_var)

setMethod("trendVar", "SCESet", function(x, subset.row=NULL, ..., assay="exprs", use.spikes=TRUE) {
    .check_centered_SF(x, assay=assay)
    if (is.null(subset.row)) {
        if (is.na(use.spikes)) {   
            subset.row <- NULL
        } else if (use.spikes) {
            subset.row <- isSpike(x, warning=FALSE)
            if (is.null(subset.row)) { 
                subset.row <- logical(nrow(x)) # no spikes at all.
            }
        } else {
            subset.row <- .spike_subset(x, get.spikes=FALSE)
        }
    }
    out <- .trend_var(assayDataElement(x, assay), ..., subset.row=subset.row)
    return(out)
})

# EXPLANATION OF TREND:
# To wit, the technical variance is modelled as:
#    var ~ s0*F(df1, df2)
# When we log it, we get:
#    log(var) ~ log(s0) + log(F(df1, df2))
# When we fit a trend to the log-variances, we get E(log(var)) for any mu.
# This can also be expressed as:
#    E(log(var)) = log(s0) + E(log(F(df1, df2)))
# So, dividing var by exp[E(log(var))] would give us:
#    var/exp[E(log(var))] ~ F(df1, df2)/exp[E(log(F(df1, df2)))]
# The estimated scale from fitFDistRobustly represents 1/exp[E(log(F(df1, df2)))].
# Thus, if we want E(var), we should be getting:
#    E(var) = s0 * df2/(df2 - 2)
#           = exp[E(log(var))] * 1/exp[E(log(F(df1, df2)))] * df2/(df2 - 2)
# ... the first term of which is 2^predict, and the product of the latter two terms is 'f.scale'.

