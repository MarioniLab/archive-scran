\name{Deconvolution Methods}
\alias{quickCluster}
\alias{quickCluster,ANY-method}
\alias{quickCluster,SummarizedExperiment0-method}

\title{Quick clustering of cells}
\description{Cluster similar cells based on rank correlations in their gene expression profiles.}

\usage{
\S4method{quickCluster}{ANY}(x, min.size=200, ...)
\S4method{quickCluster}{SummarizedExperiment0}(x, ..., i="counts")
}

\arguments{
\item{x}{
    A numeric count matrix where rows are genes and columns are cells.
    Alternatively, a SummarizedExperiment0 object containing such a matrix.
}
\item{min.size}{An integer scalar specifying the minimum size of each cluster.}
\item{...}{
    For \code{quickCluster,ANY-method}, additional arguments to be passed to \code{\link{cutreeDynamic}}.
    For \code{quickCluster,SummarizedExperiment0-method}, additional arguments to pass to \code{quickCluster,ANY-method}.
}
\item{i}{A numeric or character scalar indicating which matrix to use in the SummarizedExperiment0 object.}
}

\details{
This function provides a correlation-based approach to quickly define clusters of a minimum size \code{min.size}.
A distance matrix is constructed using Spearman's correlation on the counts between cells.
Hierarchical clustering is performed and a dynamic tree cut is used to define clusters of cells.
A correlation-based approach is preferred here as it is invariant to scaling normalization.
This avoids circularity between normalization and clustering.

Note that some cells may not be assigned to any cluster.
In most cases, this is because those cells belong in a separate cluster with fewer than \code{min.size} cells.
The function will not be able to call this as a cluster as the minimum threshold on the number of cells has not been passed.
Users are advised to check that the unassigned cells do indeed form their own cluster.
If so, it is generally safe to ignore this warning and to treat all unassigned cells as a single cluster.
Otherwise, it may be necessary to use a custom clustering algorithm.
}

\value{
A vector of cluster identities for each cell in \code{counts}. 
Values of \code{"0"} are used to indicate cells that are not assigned to any cluster.
}

\author{
Aaron Lun and Karsten Bach
}

\seealso{
\code{\link{cutreeDynamic}},
\code{\link{normalizeBySums}}
}

\examples{
set.seed(100)
popsize <- 200
ngenes <- 10000
all.facs <- 2^rnorm(popsize, sd=0.5)
counts <- matrix(rnbinom(ngenes*popsize, mu=all.facs, size=1), ncol=popsize, byrow=TRUE)
clusters <- quickCluster(counts, min.size=20)
}

\keyword{normalization}