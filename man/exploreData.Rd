\name{Explore Data}
\alias{exploreData}

\title{Shiny app for explorative data analysis}
\description{Generate an interactive Shiny app to explore gene expression patterns in single-cell data}

\usage{
exploreData(m, pD, fD, reDim, run=TRUE) {
}

\arguments{
\item{m}{numeric gene*cell count matrix with normalized expression values}
\item{pD}{data frame of cell features with rows corresponding to the columns of \code{m}}
\item{fD}{data frame of gene features with rows corresponding to the rows of \code{m}}
\item{reDim}{numeric matrix with two colums corressponding to the dimensions used for plotting and rows corresponding to the cells}
\item{run}{logical scalar indicating whether to run the app immediately or to return the app object}
}

\details{
This function will return a Shiny app object that can be run with \code{\link{runApp}}.
The app allows the user to interactively explore gene expression patterns in the data. 
Explorative analysis is focused on comparing gene exression between different groups of cells as defined by the covariates of \code{pD}.

Three plots are shown: two scatterplots based on the \code{reDim} coordinates colored by expression of a gene or a covariate and
a boxplot of log-transformed gene-expression grouped by the covariate.
The \code{fD} data frame is rendered as an interactive table using the JavaScript library DataTable.
This allows the user to subsection and search the feature data and select a gene by clicking on a row of interest.


Options within the app:
\describe{
\item{\dQuote{Color By}:}{Covariate to be used for coloring the first scatter plot.}
\item{\dQuote{Group Boxplot By}:}{Covariate to group boxplot by.}
\item{\dQuote{Click on rows}:}{Selection of a gene of interest}
}
}

\value{
If \code{run=FALSE}, a Shiny app object is returned, which can be run with \code{\link{runApp}}.
If \code{run=TRUE}, a Shiny app object is created and run.
}

\seealso{
\code{\link{runApp}}
}

\author{
Karsten Bach
}

\examples{
# Set up example data
require(viridis)
require(DT)
require(Rtsne)
example(newSCESet)
m <- counts(example_sceset)
pD <- pData(example_sceset)
fD <- data.frame("meanExp"=rowMeans(m))
rownames(fD) <- rownames(m)
tsn <- Rtsne(t(log(m+1)),perplexity=10)
reDim <- tsn$Y[,1:2]


# Creating the app object.
app <- exploreData(m, pD, fD, reDim, run=FALSE)

\dontrun{# Running the app object.
saved <- shiny::runApp(app) 

# Running the app directly from the function.
saved <- selectorPlot(x, y)
}

}
