\name{genePairRho}
\alias{genePairRho}
\alias{buildNullRho}

\title{Test for significant correlations}
\description{Identify pairs of genes that are significantly correlated based on a modified Spearman's rho.}

\usage{
buildNullRho(ncells, iters=1e6) 
genePairRho(exprs, null.dist=NULL, BPPARAM=bpparam(), use.names=TRUE)
}

\arguments{
\item{ncells}{An integer scalar indicating the number of cells in the data set.}
\item{iters}{An integer scalar specifying the number of observations in the null distribution.}
\item{exprs}{A numeric matrix of (normalized) expression values, where rows are genes and columns are cells.}
\item{null.dist}{A numeric vector of rho values under the null hypothesis.}
\item{BPPARAM}{A BiocParallelParam object to use in \code{bplapply} for parallel processing.}
\item{use.names}{A logical scalar specifying whether the row names of \code{exprs} should be used in the output.}
}

\details{
The aim of the \code{genePairRho} function is to identify significant correlations between all pairs of genes in \code{exprs}.
This allows prioritization of genes that are driving systematic substructure in the data set.
By definition, such genes should be correlated as they are behaving in the same manner across cells.
In contrast, genes driven by random noise should not exhibit any correlations with other genes.

An approximation of Spearman's rho is used to quantify correlations robustly based on ranks.
The significance of non-zero correlations is assessed using a permutation test.
The \code{buildNullRho} function constructs an empirical null distribution for the correlations with \code{iters} observations.
The p-value for each gene pair is defined as the tail probability of this distribution at the observed correlation (with some adjustment to avoid zero p-values).
Correction for multiple testing is done using the BH method.

Note that the approximation is that any tied values are randomly assigned untied ranks.
This means that a common empirical distribution can be used for all gene pairs, rather than having to do new permutations for every pair to account for the different pattern of ties.
Generally, this modification has little effect on the results - some correlations may end up being spuriously large, but this is handled by the error control machinery.

Users should select their genes in \code{exprs} with some care.
Using a top set of 100-200 highly variable genes (HVGs) is recommended.
This will focus on genes contributing to cell-to-cell heterogeneity (and thus more likely to be involved in driving substructure).
There is no need to account for HVG pre-selection in multiple testing, because rank correlations are unaffected by the variance.
For more genes, set \code{BPPARAM} to use more workers.

For \code{genePairRho}, a pre-computed empirical distribution can be supplied as \code{null.dist} if available.
Otherwise, it will be automatically constructed via \code{buildNullRho} with \code{ncells} set to the number of columns in \code{exprs}.
}

\value{
For \code{buildNullRho}, a numeric vector of length \code{iters} is returned containing the sorted correlations under the null hypothesis of no correlations.

For \code{genePairRho}, a dataframe is returned with one row per gene pair and the following fields:
\describe{
\item{\code{gene1, gene2}:}{Character or integer fields (depending on \code{use.names}) specifying the genes in the pair.}
\item{\code{rho}:}{A numeric field containing the approximate Spearman's rho.}
\item{\code{p.value, FDR}:}{Numeric fields containing the permutation p-value and its BH-corrected equivalent.}
} 
Rows are sorted by increasing \code{p.value} and, if tied, decreasing absolute size of \code{rho}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{bpparam}},
\code{\link{cor}}
}

\references{
Phipson B and Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Stat. Appl. Genet. Mol. Biol.} 9:Article 39.
}

\examples{
set.seed(0)
ncells <- 100
null.dist <- buildNullRho(ncells, iters=100000)

exprs <- matrix(rpois(ncells*100, lambda=10), ncol=ncells)
out <- genePairRho(exprs, null.dist=null.dist)
hist(out$p.value) 
}

\keyword{
correlation
}
