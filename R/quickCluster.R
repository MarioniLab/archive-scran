quickCluster <- function(counts, min.size=200, ...)
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach
# created 1 December 2015
# last modified 15 January 2016
{   
    if (ncol(counts) <= min.size){
        stop('fewer cells than the mininimum cluster size')
    }
    if (ncol(counts) < 2*min.size) {
        min.size <- as.integer((ncol(counts)/5L))
        warning(paste("min.size scaled down to", min.size))
    }

    distM <- as.dist( 1 - cor(counts, method='spearman'))
    htree <- hclust(distM, method='ward.D2')
    clusters <- unname(cutreeDynamic(htree, minClusterSize=min.size, distM=as.matrix(distM), verbose=0, ...))

    unassigned <- clusters==0L
    if (any(unassigned)) { 
        clusters[unassigned] <- NA_integer_
        warning(paste(sum(unassigned), "cells were not assigned to any cluster"))
    }
    clusters <- factor(clusters)
    return(clusters)
}
