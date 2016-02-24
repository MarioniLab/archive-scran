setGeneric("cyclone", function(x, ...) standardGeneric("cyclone"))

setMethod("cyclone", "ANY", function(x, pairs, gene.names=rownames(x), iter=1000, min.iter=100, min.pairs=50, BPPARAM=bpparam(), verbose=FALSE)
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
    workass <- .workerAssign(ncells, BPPARAM)
    out <- bplapply(seq_along(workass$start), FUN=function(core) {
        to.use <- workass$start[core]:workass$end[core]
        cur.ncells <- length(to.use)
        G1 <- .Call("shuffle_scores", cur.ncells, nrow(chosen.x$G1), chosen.x$G1[,to.use], pairs$G1[,1], pairs$G1[,2], iter, min.iter, min.pairs) 
        S <- .Call("shuffle_scores", cur.ncells, nrow(chosen.x$S), chosen.x$S[,to.use], pairs$S[,1], pairs$S[,2], iter, min.iter, min.pairs) 
        G2M <- .Call("shuffle_scores", cur.ncells, nrow(chosen.x$G2M), chosen.x$G2M[,to.use], pairs$G2M[,1], pairs$G2M[,2], iter, min.iter, min.pairs) 
        return(list(G1, S, G2M))
    }, BPPARAM=BPPARAM)

    # Assembling the output.
    score.G1 <- score.S <- score.G2M <- list()
    for (x in seq_along(out)) {
        current <- out[[x]]
        lapply(current, FUN=function(y) { if (is.character(y)) { stop(y) } })
        score.G1[[x]] <- current[[1]]
        score.S[[x]] <- current[[2]]
        score.G2M[[x]] <- current[[3]]
    }
    score.G1 <- unlist(score.G1)
    score.S <- unlist(score.S)
    score.G2M <- unlist(score.G2M)
    
    scores <- data.frame(G1=score.G1, S=score.S, G2M=score.G2M)
    scores.normalised <- scores/rowSums(scores)
    return(list(scores=scores, normalized.scores=scores.normalised))  
})

setMethod("cyclone", "SCESet", function(x, ...) {
    cyclone(.getUsedMatrix(x, "counts"), ...)          
})

