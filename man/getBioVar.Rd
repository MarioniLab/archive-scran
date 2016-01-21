\name{getBioVar}
\alias{getBioVar}
\alias{fitTechTrend}

\title{Get the biological variability}
\description{Compute the biological and technical components of the gene-specific variance in single-cell RNA-seq data.}

\usage{
fitTechTrend(spikes, size.factor=NULL, df=5, prior.count=1) 
getBioVar(counts, tech.fit, size.factor=NULL)
}

\arguments{
\item{spikes}{A numeric matrix of spike-in counts, where each column corresponds to a cell and each row corresponds to a spike-in gene.}
\item{size.factor}{A numeric vector of size factors for all cells, defaults to the sum of spike-in counts for each cell.}
\item{df}{An integer scalar specifying the degrees of freedom for polynomial fitting.}
\item{prior.count}{A numeric scalar specifying the prior count to add before log-transformation.}
\item{counts}{A numeric matrix of cellular counts, where each column corresponds to a cell and each row corresponds to a cellular gene.}
\item{tech.fit}{A list containing the output of \code{fitTechTrend}.}
}

\details{
The strategy is to fit an abundance-dependent trend to the variance of the log-CPMs for the spike-in genes, using \code{fitTechTrend}.
Log-CPMs are computed from the spike-in counts using \code{\link[edgeR]{cpm}} with the library sizes set to \code{size.factor} and the specified \code{prior.count}.
The mean and variance of the log-CPMs is calculated for each spike-in gene.
A polynomial with \code{df} degrees of freedom is then fitted to the variance against the mean for all genes.

The variance of the log-CPMs for each cellular gene is similarly computed in \code{getBioVar}.
The technical component of the variance for each gene is determined by interpolating the fitted trend at the mean log-CPM for that gene.
This represents variability due to sequencing, drop-outs during capture, etc.
The biological component is determined by subtracting the technical component from the total variance.

The size factors used to normalize the cellular counts should be the same as those used to normalize the spike-in counts.
This ensures that the average abundances computed for the former are comparable to the abundance-dependent trend fitted to the variances of the latter.
In \code{getBioVar}, \code{size.factor} is automatically taken from \code{tech.fit}.
Users should only set it manually if, e.g., only subset of cells are used in \code{counts} and \code{tech.fit$size.factor} needs to be subsetted accordingly.

Unlike other methods for decomposition, this approach estimates the variance of the log-CPMs rather than of the counts themselves.
The log-transformation blunts the impact of large positive outliers and ensures that the HVG list is not dominated by genes with expression in just one cell.
Moreover, it does not compromise interpretation -- highly variable genes will still be so, regardless of whether counts or log-CPMs are considered.

% # Compare the pair - 'a' has higher CV2 and lower var(log), while 'b' has lower CV2 and higher var(log).
% a <- integer(100)
% a[1] <- 10000
% mean(a)
% var(a)/mean(a)^2
% var(log(a+1))
% b <- rep(c(0, 200), each=50)
% mean(b)
% var(b)/mean(b)^2
% var(log(b+1))

Note that this method is rather dependent on the quality of the spike-ins -- the fit will obviously be poor if the coverage of all spike-ins is low.
}

\value{
\code{fitTechTrend} will return a named list:
\describe{
\item{\code{mean}:}{A numeric vector of mean log-CPMs for all spike-in genes.}
\item{\code{var}:}{A numeric vector of the variances of log-CPMs for all spike-in genes.}
\item{\code{trend}:}{A function that returns the fitted value of the trend at any mean log-CPM.}
\item{\code{size.factor}:}{A numeric vector of the size factors used for normalization.}
\item{\code{prior.count}:}{A numeric scalar of the prior count used.}
}

\code{getBioVar} will also return a named list:
\describe{
\item{\code{mean}:}{A numeric vector of mean log-CPMs for all cellular genes.}
\item{\code{total}:}{A numeric vector of the variances of log-CPMs for all cellular genes.}
\item{\code{bio}:}{A numeric vector containing the biological component of the variance for all genes.}
\item{\code{tech}:}{A numeric vector containing the technical component of the variance for all genes.}
}
}

\seealso{
\code{\link[edgeR]{cpm}},
\code{\link{poly}}
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
fit <- fitTechTrend(spike.data)

# Making a pretty plot.
plot(fit$mean, fit$var)
x <- sort(fit$mean)
lines(x, fit$trend(x), col="red", lwd=2)

# Getting cellular data.
ngenes <- 10000
cell.means <- 2^runif(ngenes, 2, 10)
cell.disp <- 100/cell.means +0.5 #+ 5/rchisq(ngenes, df=5)
cell.data <- matrix(rnbinom(ngenes*ncells, mu=cell.means, size=1/cell.disp), ncol=ncells)
results <- getBioVar(cell.data, fit)

plot(results$mean, results$total)
o <- order(results$mean)
lines(results$mean[o], results$tech[o], col="red", lwd=2)
plot(results$mean, results$bio)
}



\keyword{variance}
