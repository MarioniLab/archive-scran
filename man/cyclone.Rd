\name{cyclone}
\alias{cyclone}
\alias{sandbag}

\title{Cell cycle phase classification}
\description{Classify single cells into their cell cycle phases based on gene expression data.}

\usage{
sandbag(is.G1, is.S, is.G2M, training.data, 
    gene.names=rownames(training.data), fraction=0.5)

cyclone(test.data, pairs, gene.names=rownames(test.data), 
    iter=1000, min.iter=100, min.pairs=50, verbose=FALSE)
}

\arguments{
\item{is.GI, is.S, is.G2M}{A vector indicating which cells are in each phase of the cell cycle.}
\item{training.data, test.data}{A numeric matrix of gene expression values where rows are genes and columns are cells.}
\item{pairs}{A list of data.frames produced by \code{sandbag}.}
\item{gene.names}{A character vector of gene names.}
\item{fraction}{A numeric scalar specifying the minimum fraction to define a marker gene pair.}
\item{iter}{An integer scalar specifying the number of iterations for random sampling to obtain a cycle score.}
\item{min.iter}{An integer scalar specifying the minimum number of iterations for score estimation.}
\item{min.pairs}{An integer scalar specifying the minimum number of pairs for cycle estimation.}
\item{verbose}{A logical scalar specifying whether diagnostics should be printed to screen.}
}

\details{
This function implements the pair-based prediction method described by Scialdone \emph{et al.} (2015).
In the \code{sandbag} function, pairs of genes (A, B) are identified from a training data set where in each pair,
    the fraction of cells in phase G1 with expression of A > B (based on expression values in \code{training.data}) 
    and the fraction with B > A in each other phase exceeds \code{fraction}.
These pairs are defined as the marker pairs for G1.
This is repeated for each phase to obtain a separate marker set.

Phase assignment is then performed using the \code{cyclone} function.
For each phase and each cell, the function calculates the proportion of all marker pairs where the expression of the first gene is greater than the second
    (pairs with the same expression are ignored).
A distribution of proportions is constructed by shuffling the expression values within the cell and recalculating the proportion at each iteration.
The phase score for that cell is then defined as the lower tail probability of this distribution.

By default, shuffling is performed \code{iter} times 
    -- however, some iterations may not be used if there are fewer than \code{min.pairs} pairs with different expression,
       such that the proportion cannot be calculated precisely.
A score is only returned if the distribution is large enough for stable calculation of the tail probability, i.e., consists of results from at least \code{min.iter} iterations.

Cells with G1 and G2M scores above 0.5 should be assigned to the G1 and G2M phases, respectively.
This is based on the interpretation of the score as 1 minus the p-value, under the null hypothesis of independent gene expression within each cell.
Cells can be assigned to S phase based on the S phase score, but a more reliable approach is to define S phase cells based on those cells with G1 and G2M scores below 0.5.

% The default implementation shuffles expression values from all genes in the training set, rather than only genes in the marker set for each phase.
% In theory, this shouldn't matter, as the score calculation depends on the ranking of each gene.
% That should be the same regardless of the distribution of expression values -- each set of rankings is equally likely, no matter what.
% In practice, the number of ties will differ between different set of genes.
% The most appropriate comparison would involve the same number of ties as that used to calculate the observed score.
% Thus, I would set the marker set genes as the empirical distribution, rather than the entire training set of genes.
}

\value{
For \code{sandbag}, a named list of data.frames. 
Each data frame corresponds to a cell cycle phase and contains the names of the genes in each marker pair.

For \code{cyclone}, a list of two data.frames -- 
    \code{scores}, containing the phase scores for each phase and cell; 
    and \code{normalized.scores}, containing the row-normalized scores (i.e., where the row sum for each cell is equal to 1).
}

\author{
Antonio Scialdone,
with modifications by Aaron Lun
}

\examples{
ncells <- 50
ngenes <- 20
training <- matrix(rnorm(ncells*ngenes), ncol=ncells)
rownames(training) <- paste0("X", seq_len(ngenes))

is.G1 <- 1:20
is.S <- 21:30
is.G2M <- 31:50

# Training
out <- sandbag(is.G1, is.S, is.G2M, training) 

# Classifying (obviously, you would use something else in real data)
test <- training 
assignments <- cyclone(test, out)

# Visualizing
col <- character(ncells)
col[is.G1] <- "red"
col[is.G2M] <- "blue"
col[is.S] <- "darkgreen"
plot(assignments$score$G1, assignments$score$G2M, col=col, pch=16)
}

\references{
Scialdone A, Natarajana KN, Saraiva LR \emph{et al.} (2015). 
Computational assignment of cell-cycle stage from single-cell transcriptome data.
\emph{Methods} 85:54--61
}

\keyword{clustering}
