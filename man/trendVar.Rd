\name{trendVar}
\alias{trendVar}
\alias{trendVar,ANY-method}
\alias{trendVar,SCESet-method}

\title{Get the biological variability}
\description{Compute the biological and technical components of the gene-specific variance in single-cell RNA-seq data.}

\usage{
\S4method{trendVar}{ANY}(x, trend=c("poly", "loess"), df=5, span=0.3, prior.count=1, design=NULL) 
\S4method{trendVar}{SCESet}(x, ..., use.spikes=TRUE)
}

\arguments{
\item{x}{
    A numeric matrix of normalized expression values, where each column corresponds to a cell and each row corresponds to a spike-in transcript.
    Alternatively, a SCESet object that contains such values.
}
\item{trend}{A string indicating whether the trend should be polynomial or loess-based.}
\item{df}{An integer scalar specifying the degrees of freedom for polynomial fitting.}
\item{span}{An numeric scalar specifying the span for loess fitting.}
\item{prior.count}{A numeric scalar specifying the prior count to add before log-transformation.}
\item{design}{A numeric matrix describing the systematic factors contributing to expression in each cell.}
\item{...}{Additional arguments to pass to \code{trendVar,ANY,list-method}.}
\item{use.spikes}{A logical scalar specifying whether the trend should be fitted to variances for spike-in or endogenous genes.}
}

\details{
The strategy is to fit an abundance-dependent trend to the variance of the log-CPMs for the spike-in genes, using \code{trendVar}.
Log-CPMs are computed from the spike-in counts using \code{\link[edgeR]{cpm}} with the library sizes set to \code{size.factor} and the specified \code{prior.count}.
Normalization should be performed based on the spike-in totals using \code{\link{normalizeBySpikes}} followed by \code{\link{normalize}}.
Log-transformed values are used as these tend to be more robust to genes with strong expression in only one or two outlier cells.

The mean and variance of the log-CPMs is calculated for each spike-in gene.
A polynomial or loess trend is then fitted to the variance against the mean for all genes.
This represents technical variability due to sequencing, drop-outs during capture, etc.
Variance decomposition to biological and technical components for endogenous genes can then be performed 

The design matrix can be set if there are factors that should be blocked, e.g., batch effects, known (and uninteresting) clusters.
Otherwise, it will default to an all-ones matrix, effectively treating all cells as part of the same group.
}

\section{Additional notes}{
By default, a polynomial trend with \code{df} degrees of freedom is fitted to the spike-in variances as it is more precise than the loess curve.
Note that this method is rather dependent on the quality of the spike-ins -- the fit will obviously be poor if the coverage of all spike-ins is low.
In some data sets, a loess curve may yield better results, though this may require some fiddling with the \code{span}.

When spike-ins are not available, \code{trendVar} can also be applied directly to the counts for endogenous genes by setting \code{use.spikes=FALSE} 
    (or by manually supplying a matrix of normalized expression for endogenous genes, for \code{trendVar,ANY-method}).
In this case, the trend is fitted to variances for endogenous genes and represents the uninteresting mean-variance relationship.
Subtracting the trend from the variance of each gene will eliminate the mean-variance dependency, yielding a residual that can be ranked to identify highly variable genes. 
A loess curve is recommended here as it is more robust to a subset of genes with very large or very small variances.
}

\value{
A named list is returned, containing:
\describe{
\item{\code{mean}:}{A numeric vector of mean log-CPMs for all spike-in genes.}
\item{\code{var}:}{A numeric vector of the variances of log-CPMs for all spike-in genes.}
\item{\code{trend}:}{A function that returns the fitted value of the trend at any mean log-CPM.}
\item{\code{prior.count}:}{A numeric scalar of the prior count used.}
\item{\code{design}:}{A numeric matrix, containing the design matrix that was used.}
}
}

\seealso{
\code{\link{poly}},
\code{\link{loess}},
\code{\link{decomposeVar}},
\code{\link{normalizeBySpikes}},
\code{\link{normalize}}
}

\author{
Aaron Lun
}

\examples{
set.seed(100)

nspikes <- ncells <- 100
spike.means <- 2^runif(nspikes, 3, 8)
spike.disp <- 100/spike.means + 0.5
spike.data <- matrix(rnbinom(nspikes*ncells, mu=spike.means, size=1/spike.disp), ncol=ncells)

ngenes <- 10000
cell.means <- 2^runif(ngenes, 2, 10)
cell.disp <- 100/cell.means + 0.5
cell.data <- matrix(rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)

# Normalizing.
y <- countsToSE(cell.data, spike.data)
y$size.factor <- normalizeBySpikes(y)
y <- normalize(y)

# Fitting a trend to the spike-ins.
fit <- trendVar(y)
plot(fit$mean, fit$var)
x <- sort(fit$mean)
lines(x, fit$trend(x), col="red", lwd=2)

# Fitting a trend to the endogenous genes. 
fit <- trendVar(y, use.spikes=FALSE)
plot(fit$mean, fit$var)
x <- sort(fit$mean)
lines(x, fit$trend(x), col="red", lwd=2)
}

\keyword{variance}

