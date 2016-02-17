\name{countsToSE}
\alias{countsToSE}

\title{Counts to SummarizedExperiment0}
\description{Construct a SummarizedExperiment0 object from read count data for a scRNA-seq experiment.}

\usage{
countsToSE(counts, spikes)
}

\arguments{
\item{counts}{A matrix of read counts for endogenous genes.}
\item{spikes}{An optional matrix of read counts for spike-in transcripts, with the same number of columns as \code{counts}.}
}

\details{
This function constructs a SummarizedExperiment0 object with a read count matrix for endogenous genes as the assay.
Library sizes and size factors are included in the \code{colData} (by default, the mean-centered library sizes are used as the size factors).
If spike-in counts are supplied, this are also placed in the \code{colData}.
They are not stored with the counts for the endogenous genes, to avoid having to continually remove spike-in transcripts from the results for downstream analyses.
}

\value{
A SummarizedExperiment0 object.
}

\author{
Aaron Lun
}

\seealso{
\code{\link[SummarizedExperiment]{SummarizedExperiment}}
}

\examples{
set.seed(100)
popsize <- 20
ngenes <- 1000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=10*all.facs, size=1), ncol=popsize, byrow=TRUE)
spikes <- matrix(rnbinom(100*popsize, mu=10*all.facs, size=0.5), ncol=popsize, byrow=TRUE)
y <- countsToSE(counts, spikes)
}
