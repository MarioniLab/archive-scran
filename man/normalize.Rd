\name{Normalized Expression}
\alias{normalize,ANY-method}
\alias{normalize,SCESet-method}

\title{Compute normalized expression values}
\description{Compute (log-)expression values from read counts, using pre-calculated size factors for normalization.}

\usage{
\S4method{normalize}{ANY}(object, size.factor=NULL, log=TRUE, prior.count=1)
\S4method{normalize}{SCESet}(object, ..., separate.spikes=TRUE)
}

\arguments{
\item{object}{A matrix of read counts, or a SCESet object with an assay named \code{"counts"}.}
\item{size.factor}{A numeric vector of size factors.}
\item{log}{A logical scalar specifying whether the expression should be log-transformed.}
\item{prior.count}{A numeric scalar indicating the prior count to add prior to log-transformation, to avoid undefined values from zero counts.}
\item{...}{Additional arguments to pass to \code{normalize,ANY-method}.}
\item{separate.spikes}{A logical scalar indicating whether spike-in counts should be normalized separately.}
}

\details{
This function computes normalized log-expression values by adding \code{prior.count} to each count, dividing by the \code{size.factor} for that cell, and log-transforming.
If \code{object} is a SCESet object, size factors are taken from the appropriate field in the \code{colData}.
These size factors can be computed with a number of functions like \code{\link{normalizeBySums}} or \code{\link{normalizeBySpikes}}.

If spike-in counts are present in the SCESet object, these will also be converted into normalized values.
If \code{separate.spikes=FALSE}, this is done with the same set of size factors that was used for the endogenous genes.
Otherwise, a separate set of spike-in size factors will be used instead -- these are defined by calling \code{\link{normalizeBySpikes}}.

All size factors are mean-centered so that their geometric mean is equal to unity prior to computing normalized expression values.
This ensures that expression values are roughly comparable when different sets of size factors are used.
}

\value{
For \code{normalize,SCESet-method}, a SCESet object is returned with an additional assay named \code{"norm_exprs"}.
This contains normalized log-expression values for the endogenous genes.
If spike-ins are present, normalized values for the spike-in transcripts are stored in the \code{norm.spikes} field of the \code{colData}.
}

\seealso{
\code{\link{cpm}},
\code{\link{SCESet}},
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

combined <- rbind(counts, spikes)
colnames(combined) <- seq_len(popsize)
rownames(combined) <- seq_len(nrow(combined))
y <- newSCESet(countData=combined)
y <- calculateQCMetrics(y, feature_controls=list(Spike=rep(c(FALSE, TRUE), c(ngenes, 100))))

assayDataElement(normalize(y), "norm_exprs")[1:10,] # Library size normalization, basically.
y$size.factor <- normalizeBySpikes(y)
assayDataElement(normalize(y), "norm_exprs")[1:10,]
}
