\name{Deconvolution Methods}
\alias{normalizeBySpikes}
\alias{normalizeBySpikes,SummarizedExperiment0-method}

\title{Normalization by spikes}
\description{Compute size factors based on the coverage of spike-in transcripts.}

\usage{
\S4method{normalizeBySpikes}{SummarizedExperiment0}(x)
}

\arguments{
\item{x}{A SummarizedExperiment0 object containing spike-in counts in the \code{colData}.}
}

\details{
The size factor for each cell is defined as the sum of all spike-in counts in each cell.
This is equivalent to normalizing to equalize spike-in coverage between cells.
Spike-in counts are assumed to be stored in \code{colData(x)$spikes} -- see \code{\link{countsToSE}} for more details.
Note that the geometric mean of all size factors is set to unity, for standardization purposes if one were to compare different sets of size factors.
}

\value{
A numeric vector of size factors for all cells in \code{x} is returned.
}

\author{
Aaron Lun
}

\examples{
set.seed(100)
popsize <- 200
ngenes <- 1000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs*10, size=1), ncol=popsize, byrow=TRUE)
spikes <- matrix(rnbinom(100*popsize, mu=all.facs*10, size=0.5), ncol=popsize, byrow=TRUE)

y <- countsToSE(counts, spikes)
out.facs <- normalizeBySpikes(y)
}

\keyword{normalization}
