\name{Normalized Expression}
\alias{normalize,ANY-method}
\alias{normalize,SummarizedExperiment0-method}

\title{Compute normalized expression values}
\description{Compute (log-)expression values from read counts, using pre-calculated size factors for normalization.}

\usage{
\S4method{normalize}{ANY}(object, size.factor=NULL, log=TRUE, prior.count=1)
\S4method{normalize}{SummarizedExperiment0}(object, ...)
}

\arguments{
\item{object}{A matrix of read counts, or a SummarizedExperiment0 object with an assay named \code{"counts"}.}
\item{size.factor}{A numeric vector of size factors.}
\item{log}{A logical scalar specifying whether the expression should be log-transformed.}
\item{prior.count}{A numeric scalar indicating the prior count to add prior to log-transformation, to avoid undefined values from zero counts.}
\item{...}{Additional arguments to pass to \code{normalize,ANY-method}.}
}

\details{
This function computes normalized log-expression values by adding \code{prior.count} to each count, dividing by the \code{size.factor} for that cell, and log-transforming.
If \code{object} is a SummarizedExperiment0 object, size factors are taken from the appropriate field in the \code{colData}.
These size factors can be computed with a number of functions like \code{\link{normalizeBySums}} or \code{\link{normalizeBySpikes}}.
If spike-in counts are present, these will also be converted into normalized values using the same parameters.
}

\value{
For \code{normalize,SummarizedExperiment0-method}, a SummarizedExperiment0 object is returned with an additional assay named \code{"exprs"}.
This contains normalized log-expression values for the endogenous genes.
If spike-ins are present, normalized values for the spike-in transcripts are stored in the \code{norm.spikes} field of the \code{colData}.
}

\seealso{
\code{\link{cpm}},
\code{\link{SummarizedExperiment}},
\code{\link{normalizeBySums}},
\code{\link{normalizeBySpikes}}
}

\author{
Aaron Lun
}

\examples{
set.seed(100)
popsize <- 10
ngenes <- 1000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=10*all.facs, size=1), ncol=popsize, byrow=TRUE)
spikes <- matrix(rnbinom(100*popsize, mu=10*all.facs, size=0.5), ncol=popsize, byrow=TRUE)

y <- countsToSE(counts, spikes)
assay(normalize(y), "exprs")[1:10,]

y$size.factor <- normalizeBySpikes(y)
assay(normalize(y), "exprs")[1:10,]
}
