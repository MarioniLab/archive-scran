\name{Get spikes}
\alias{spikes}
\alias{spikes,SCESet-method}

\title{Construct spike-in matrix}
\description{Reconstruct the matrix of counts or normalized expression values for the spike-in transcripts.}

\usage{
\S4method{spikes}{SCESet}(x, type=c("counts", "norm_exprs"))
}

\arguments{
\item{x}{A SCESet object with spike-in data in the \code{colData}.}
\item{type}{A string specifying whether counts or normalized expression values are to be extracted.}
}

\details{
This function extracts the spike-in data from \code{x}, into a numeric matrix that can be used for downstream analyses.
Users should put a logical vector in \code{fData(x)$is_feature_spike} indicating if each row of \code{x} the spike-in controls
If \code{type="norm_exprs"}, users should have run \code{x} through \code{\link{normalize}}.
}

\value{
A numeric matrix of counts or normalized expression values, with one column per cell and one row per spike-in transcript.
}

\seealso{
\code{\link{normalize}},
\code{\link{SCESet}}
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

combined <- rbind(counts, spikes)
colnames(combined) <- seq_len(popsize)
rownames(combined) <- seq_len(nrow(combined))
y <- newSCESet(countData=combined)
fData(y)$is_feature_spike <- rep(c(FALSE, TRUE), c(ngenes, 100))

y$size.factor <- normalizeBySpikes(y)
y <- normalize(y)
spikes(y)[1:10,]
spikes(y, type="norm_exprs")[1:10,]
}
