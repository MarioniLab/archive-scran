.buildSNNGraph <- function(x, k=10, subset.row=NULL) 
# Builds a shared nearest-neighbor graph, where edges are present between each 
# cell and its 'k' nearest neighbours. Edges are weighted based on the ranks of 
# the shared nearest neighbours of the two cells, as described in the SNN-Cliq paper.
#
# written by Aaron Lun
# created 3 April 2017    
{ 
    if (!is.null(subset.row)) {
        x <- x[.subset_to_index(subset.row, x, byrow=TRUE),]
    }
    nn.out <- get.knn(t(x), k=k, algorithm="cover_tree") 
    nn.index <- nn.out$nn.index
    all.nodes <- row(nn.out$nn.index)
    edges <- as.vector(rbind(as.integer(all.nodes), as.integer(nn.index)))

    weights <- matrix(0, nrow=nrow(all.nodes), ncol=ncol(all.nodes))
    for (i in seq_len(nrow(all.nodes))) {
        i.nn <- c(i, nn.index[i,])
        for (x in seq_len(k)) {
            j <- i.nn[x+1]
            j.nn <- c(j, nn.index[j,])
            shared <- intersect(i.nn, j.nn) 
            weights[i, x] <- k - 0.5 * min(match(shared, i.nn) + match(shared, j.nn))
        }
    }

    g <- make_graph(edges, directed=FALSE)
    E(g)$weight <- weights
    g <- simplify(g, edge.attr.comb="first") # symmetric, so doesn't really matter.
    return(g)
}

setGeneric("buildSNNGraph", function(x, ...) standardGeneric("buildSNNGraph"))

setMethod("buildSNNGraph", "matrix", .buildSNNGraph)

setMethod("buildSNNGraph", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { 
        subset.row <- .spikeSubset(x, get.spikes)
    }
    .buildSNNGraph(assayDataElement(x, assay), ..., subset.row=subset.row)
})
