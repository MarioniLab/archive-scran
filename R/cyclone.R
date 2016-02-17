setGeneric("cyclone", function(x, ...) { standardGeneric("cyclone") })

setMethod("cyclone", "ANY", function(x, pairs, gene.names=rownames(x), iter=1000, min.iter=100, min.pairs=50, verbose=FALSE)
# Takes trained pairs and test data, and predicts the cell cycle phase from that. 
#
# written by Antonio Scialdone
# with modifications by Aaron Lun
# created 22 January 2016    
# last modified 17 February 2016
{ 
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    Ngenes <- nrow(x)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' must be equal to 'x' nrows")
    }
    iter <- as.integer(iter)
    min.iter <- as.integer(min.iter)
    min.pairs <- as.integer(min.pairs)

    # Only keeping training pairs where both genes are in the test data;
    # and subsetting the test data to only the genes in the training pairs.
    chosen.x <- list()
    for (p in names(pairs)) {
        curp <- pairs[[p]]
        m1 <- match(curp$first, gene.names)
        m2 <- match(curp$second, gene.names)
        keep <- !is.na(m1) & !is.na(m2)
        m1 <- m1[keep]
        m2 <- m2[keep]
       
        all.present <- union(m1, m2)
        chosen.x[[p]] <- x[all.present,,drop=FALSE]
        pairs[[p]] <- data.frame(first=match(m1, all.present),
                                 second=match(m2, all.present))
    }

    if (verbose) { 
        cat(sprintf("Number of G1 pairs: %d\n", nrow(pairs$G1)))
        cat(sprintf("Number of S pairs: %d\n", nrow(pairs$S)))
        cat(sprintf("Number of G2M pairs: %d\n", nrow(pairs$G2)))
    }
  
    # Run the allocation algorithm
    ncells <- ncol(x)
    score.G1 <- .Call("shuffle_scores", ncells, nrow(chosen.x$G1), chosen.x$G1, pairs$G1[,1], pairs$G1[,2], iter, min.iter, min.pairs) 
    score.S <- .Call("shuffle_scores", ncells, nrow(chosen.x$S), chosen.x$S, pairs$S[,1], pairs$S[,2], iter, min.iter, min.pairs) 
    score.G2M <- .Call("shuffle_scores", ncells, nrow(chosen.x$G2M), chosen.x$G2M, pairs$G2M[,1], pairs$G2M[,2], iter, min.iter, min.pairs) 
    
    scores <- data.frame(G1=score.G1, S=score.S, G2M=score.G2M)
    scores.normalised <- scores/rowSums(scores)
    return(list(scores=scores, normalized.scores=scores.normalised))  
})

setMethod("cyclone", "SummarizedExperiment0", function(x, ..., i="counts") {
    cyclone(assay(x, i=i), ...)          
})

