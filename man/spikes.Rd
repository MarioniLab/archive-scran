\name{Get spikes}
\alias{spikes}
\alias{spikes,SummarizedExperiment0-method}

\title{Construct spike-in matrix}
\description{Reconstruct the matrix of counts or normalized expression values for the spike-in transcripts.}

\usage{
\S4method{spikes}{SummarizedExperiment0}(x, type=c("counts", "exprs"))
}

\arguments{
\item{x}{A SummarizedExperiment0 object with spike-in data in the \code{colData}.}
\item{type}{A string specifying whether counts or normalized expression values are to be extracted.}
}

\details{
This function extracts the spike-in data from \code{x}, into a numeric matrix that can be used for downstream analyses.
Users should construct \code{x} via \code{\link{countsToSE}} with a spike-in matrix supplied.
If \code{type="exprs"}, users should have run \code{x} through \code{\link{normalize}}.
}

\value{
A numeric matrix of counts or normalized expression values, with one column per cell and one row per spike-in transcript.
}

\seealso{
\code{\link{SummarizedExperiment}},
\code{\link{normalize}},
\code{\link{countsToSE}}
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
y$size.factor <- normalizeBySpikes(y)
y <- normalize(y)
spikes(y)[1:10,]
spikes(y, type="exprs")[1:10,]
}
