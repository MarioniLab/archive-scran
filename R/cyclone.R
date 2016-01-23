find.markers <- function(id1, id2, id3, training.data, fraction=0.5, gene.names=rownames(training.data))
# This identifies pairs of genes whose relative expression is > 0 in 
# at least a 'fraction' of cells in one phase is < 0 in at least 
# 'fraction' of the cells in each of the other phases.
{
    Ngenes <- nrow(training.data)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' vector must be equal to 'training.data' nrows")
    }

    data1 <- t(training.data[,id1])
    data2 <- t(training.data[,id2])
    data3 <- t(training.data[,id3])  

    Nthr1 <- ceiling(nrow(data1) * fraction)
    Nthr2 <- ceiling(nrow(data2) * fraction)
    Nthr3 <- ceiling(nrow(data3) * fraction)

    collected <- list()
    counter <- 1L
    if (Ngenes) { 
        for (i in seq_len(Ngenes-1L)) { 
            others <- (i+1):Ngenes
            diff1 <- data1[,i] - data1[,others,drop=FALSE]
            diff2 <- data2[,i] - data2[,others,drop=FALSE]
            diff3 <- data3[,i] - data3[,others,drop=FALSE]

            npos1 <- colSums(diff1 > 0)
            npos2 <- colSums(diff2 > 0)
            npos3 <- colSums(diff3 > 0)
            nneg1 <- colSums(diff1 < 0)
            nneg2 <- colSums(diff2 < 0)
            nneg3 <- colSums(diff3 < 0)

            chosen <- others[npos1 >= Nthr1 & nneg2 >= Nthr2 & nneg3 >= Nthr3]
            if (length(chosen)) { 
                collected[[counter]] <- cbind(i, chosen)
                counter <- counter + 1L
            }
            chosen.flip <- others[nneg1 >= Nthr1 & npos2 >= Nthr2 & npos3 >= Nthr3]
            if (length(chosen.flip)) { 
                collected[[counter]] <- cbind(chosen.flip, i)
                counter <- counter + 1L
            }
        }
    }

    collected <- do.call(rbind, collected)
    return(data.frame(first=gene.names[collected[,1]], 
                      second=gene.names[collected[,2]]))
}

classify.single <- function(cell, markers, Nmin.couples) 
# Count number of successes for a single set of markers, given a cell's expression values.
# Successes are defined as those pairs where the first gene is more highly expressed.
{
    test <- cell[markers[,1]] - cell[markers[,2]]
    t1 <- sum(test>0)
    tot <- sum(test!=0)
    if (tot < Nmin.couples) { return(NA) }
    return(t1/tot)
}

random.success <- function(cell, markers, N, Nmin, Nmin.couples)
# Given a set of markers, a number N of trials and a cell, find the number of hits in each of the N randomised set of markers.
# This returns the probability of randoml obtaining a fraction of hits lower than the observed fraction.
# The null hypothesis is that the expression of each gene is sampled from the empirical distribution in 'cell',
# in a manner that is independent of the pairings between genes. We then calculate the classification based on those pairings.
{
    success <- sapply(seq_len(N), function(x) {    
                      cell.random <- cell[sample(length(cell))]
                      classify.single(cell.random, markers, Nmin.couples)
    })
  
    success <- success[!is.na(success)]
    test <- classify.single(cell, markers, Nmin.couples)

    if(length(success) < Nmin || is.na(test)) { 
        warning("not enough gene pairs with different expression values")
        return(NA) 
    }
    return(sum(success<test)/length(success))
}

sandbag <- function(is.G1, is.S, is.G2M, training.data, gene.names=rownames(training.data), fraction=0.5) 
# Identifies the relevant pairs before running 'cyclone'.
# Basically runs through all combinations of 'find.markers' for each phase. 
#
# written by Aaron Lun
# based on code by Antonio Scialdone
# created 22 January 2016 
{
    G1.marker.pairs <- find.markers(id1=is.G1, id2=is.S, id3=is.G2M, training.data=training.data, fraction=fraction, gene.names=gene.names)
    S.marker.pairs <- find.markers(id1=is.S, id2=is.G1, id3=is.G2M, training.data=training.data, fraction=fraction, gene.names=gene.names)
    G2M.marker.pairs <- find.markers(id1=is.G2M, id2=is.G1, id3=is.S, training.data=training.data, fraction=fraction, gene.names=gene.names)
    return(list(G1=G1.marker.pairs, S=S.marker.pairs, G2M=G2M.marker.pairs))
}

cyclone <- function(test.data, pairs, gene.names=rownames(test.data), iter=1000, min.iter=100, min.pairs=50, verbose=FALSE)
# Takes trained pairs and test data, and predicts the cell cycle phase from that. 
#
# written by Antonio Scialdone
# with modifications by Aaron Lun
# created 22 January 2016    
{ 
    Ngenes <- nrow(test.data)
    if (length(gene.names)!=Ngenes) {
        stop("length of 'gene.names' vector must be equal to 'test.data' nrows")
    }
    iter <- as.integer(iter)
    min.iter <- as.integer(min.iter)
    min.pairs <- as.integer(min.pairs)
    test.data <- as.matrix(test.data)
    storage.mode(test.data) <- "double"

    # Only keeping training pairs where both genes are in the test data;
    # and subsetting the test data to only the genes in the training pairs.
    best.test.data <- list()
    for (x in names(pairs)) {
        curx <- pairs[[x]]
        m1 <- match(curx$first, gene.names)
        m2 <- match(curx$second, gene.names)
        keep <- !is.na(m1) & !is.na(m2)
        m1 <- m1[keep]
        m2 <- m2[keep]
       
        all.present <- union(m1, m2)
        best.test.data[[x]] <- test.data[all.present,,drop=FALSE]
        pairs[[x]] <- data.frame(first=match(m1, all.present),
                                 second=match(m2, all.present))
    }

    if (verbose) { 
        cat(sprintf("Number of G1 pairs: %d\n", nrow(pairs$G1)))
        cat(sprintf("Number of S pairs: %d\n", nrow(pairs$S)))
        cat(sprintf("Number of G2M pairs: %d\n", nrow(pairs$G2)))
    }
  
    # Run the allocation algorithm
    ncells <- ncol(test.data)
    score.G1 <- .Call("shuffle_scores", ncells, nrow(best.test.data$G1), best.test.data$G1, pairs$G1[,1], pairs$G1[,2], iter, min.iter, min.pairs) 
    score.S <- .Call("shuffle_scores", ncells, nrow(best.test.data$S), best.test.data$S, pairs$S[,1], pairs$S[,2], iter, min.iter, min.pairs) 
    score.G2M <- .Call("shuffle_scores", ncells, nrow(best.test.data$G2M), best.test.data$G2M, pairs$G2M[,1], pairs$G2M[,2], iter, min.iter, min.pairs) 
#    score.G1 <- apply(best.test.data$G1, 2, FUN=random.success, markers=pairs$G1, N=iter, Nmin=min.iter, Nmin.couples=min.pairs)
#    score.S <- apply(best.test.data$S, 2, FUN=random.success, markers=pairs$S, N=iter, Nmin=min.iter, Nmin.couples=min.pairs)
#    score.G2M <- apply(best.test.data$G2M, 2, FUN=random.success, markers=pairs$G2M, N=iter, Nmin=min.iter, Nmin.couples=min.pairs)
    
    scores <- data.frame(G1=score.G1, S=score.S, G2M=score.G2M)
    scores.normalised <- scores/rowSums(scores)
    return(list(scores=scores, normalized.scores=scores.normalised))  
}

