correctMNN <- function(..., k=20, sigma=1, cos.norm=TRUE, svd.dim=20, order=NULL) 
# Performs correction based on the batches specified in the ellipsis.
#    
# written by Laleh Haghverdi
# with modifications by Aaron Lun
# created 7 April 2017
# last modified 9 April 2017
{ 
    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { stop("at least two batches must be specified") }

    # Setting up the order.
    if (is.null(order)) {
        order <- seq_len(nbatches)
    } else {
        order <- as.integer(order)
        if (!identical(seq_len(nbatches), sort(order))) { 
            stop(sprintf("'order' should contain values in 1:%i", nbatches))
        }
    }

    # Setting up the variables.
    ref <- order[1]
    ref.batch <- batches[[ref]]
    if (cos.norm) { ref.batch <- cosine.norm(ref.batch) }
    num.mnn <- matrix(NA_integer_, nbatches, 2)
    output <- vector("list", nbatches)
    output[[ref]] <- ref.batch

    for (b in 2:nbatches) { 
        other.batch <- batches[[order[b]]]
        if (cos.norm) { other.batch <- cosine.norm(other.batch) } 

        # Finding pairs of mutual nearest neighbours.
        sets <- find.mutual.nn(t(ref.batch), t(other.batch), k1=k, k2=k, sigma=sigma, svd.dim=svd.dim)
        s1 <- sets$set1
        s2 <- sets$set2

        # Computing the biological subspace in both batches.
        ndim <- min(c(svd.dim, dim(ref.batch), dim(other.batch)))
        span1 <- get.bio.span(ref.batch[,s1,drop=FALSE], min(ndim, length(s1)))
        span2 <- get.bio.span(other.batch[,s2,drop=FALSE], min(ndim, length(s2)))
        bio.span <- cbind(span1, span2)
  
        # Identifying the biological component of the batch correction vector 
        # (i.e., the part that is parallel to the biological subspace) and removing it.
        bv <- sets$batchvect.all
        bio.comp <- bv %*% bio.span %*% t(bio.span)
        correction <- t(bv) - t(bio.comp) 

        # Applying the correction and storing the numbers of nearest neighbors.
        other.batch <- other.batch + correction
        num.mnn[b,] <- c(length(s1), length(s2))
        output[[b]] <- other.batch

        # Expanding the reference batch to include the new, corrected data.
        ref.batch <- cbind(ref.batch, other.batch)
    }

    # Formatting output to be consistent with input.
    names(output) <- names(batches)
    list(corrected=output, num.mnn=num.mnn)
}

# Internal functions:
find.mutual.nn <- function(data1, data2, k1, k2, sigma=1, svd.dim=20) 
# Finds mutal neighbors between data1 and data2.
# Computes the batch correction vector for each cell in data2.
{
    n1 <- nrow(data1)
    n2 <- nrow(data2)
    n.total <- n1 + n2
    W <- matrix(0, n.total, n.total)

    W21 <- FNN::get.knnx(data2, query=data1, k=k1)
    js1 <- matrix(seq_len(n1), n1, k1)
    is1 <- n1 + W21$nn.index
    indices1 <- cbind(as.vector(js1), as.vector(is1))
    W[indices1] <- 1

    W12 <- FNN::get.knnx(data1, query=data2, k=k2)
    js2 <- matrix(n1 + seq_len(n2), n2, k2)
    is2 <- W12$nn.index
    indices2 <- cbind(as.vector(js2), as.vector(is2))
    W[indices2] <- 1

    W <- W * t(W)  #elementwise multiplication to keep mutual nns only
    A <- which(W>0, arr.ind=TRUE) # row/col indices of mutual NNs
    set <- A

    # Computing the batch correction vector between MNN pairs.
    A1 <- A[,1]
    A1 <- A1[A1 <= n1]
    A2 <- A[,2] - n1
    A2 <- A2[A2 > 0]
    vect <- data1[A1,] - data2[A2,]    
    #row.names(vect) <- as.character(A2)

    # Checking if batches share any subspace.
    exprs1 <- t(data1[unique(A1),])
    exprs2 <- t(data2[unique(A2),])
    ndim <- min(c(svd.dim, dim(exprs1), dim(exprs2)))
    span1 <- get.bio.span(exprs1, ndim)
    span2 <- get.bio.span(exprs2, ndim)
    nshared <- find.shared.subspace(span1, span2)$nshared
    if (nshared==0L) { warning("batches not sufficiently related") }

    # Gaussian smoothing of individual correction vectors for MNN pairs.
    dd2 <- as.matrix(dist(data2))
    if (sigma==0) {
        G <- matrix(1, nrow(data2), ncol(data2))
    } else {
        G <- exp(-dd2^2/sigma)  
    }
    nA2<-table(A2)
    
    D <- colSums(G)
    batchvect <- matrix(0, nrow(data2), ncol(data2))   
    
    F1 <- matrix(D[A2], nrow=nrow(data2), ncol=length(A2), byrow=TRUE)
    F2<-matrix(rep(nA2[as.character(A2)],each=dim(data2)[1]),nrow=dim(data2)[1],ncol=length(A2))

    batchvect<-batchvect+(G[,A2]/ (F1*F2)) %*% vect#[A2,] #[as.character(A2),]  #density normalized for cancelling strong effect from dense parts
    partitionf<-rowSums(G[,A2]/(F1*F2))
    
    
    batchvect <- batchvect/partitionf

    # Report cells that are MNNs, and the correction vector per cell in data2.
    set1 <- set[set<(n1+1)]
    set2 <- set[set>n1]-n1
    list(set1=unique(set1), set2=unique(set2), batchvect.all=batchvect, nshared=nshared) 
}

get.bio.span <- function(exprs, ndim) 
# Computes the biological span within a data set for a given number of SVD dimensions.
# Avoids extra dimensions dominated by technical noise, which will result in both 
# trivially large and small angles when using find.shared.subspace().
{ 
    S <- svd(exprs)
    used.dim <- seq_len(ndim)
    S$u[,used.dim,drop=FALSE] %*% S$v[used.dim,used.dim,drop=FALSE]
}

find.shared.subspace <- function(A, B, sin.threshold=0.85, cos.threshold=1/sqrt(2)) 
# Computes the maximum angle between subspaces, to determine if spaces are orthogonal.
# Also identifies if there are any shared subspaces. 
{
    A <- pracma::orth(A)
    B <- pracma::orth(B)
     
    # Singular values close to 1 indicate shared subspace A \invertsymbol{U} B
    # Otherwise A and B are completely orthogonal, i.e., all angles=90.
    S <- svd(t(A) %*% B)
    shared <- sum(S$d > sin.threshold)

    # Computing the angle from singular values; using cosine for large angles,
    # sine for small angles (due to differences in relative accuracy).
    costheta <- min(S$d) 
    if (costheta < cos.threshold){ 
        theta <- acos(min(1, costheta))
    } else {
        if (ncol(A) < ncol(B)){ 
            sintheta <- svd(t(A) - (t(A) %*% B) %*% t(B))$d[1]
        } else {
            sintheta <- svd(t(B) - (t(B) %*% A) %*% t(A))$d[1]
        }
        theta <- asin(min(1, sintheta)) 
    }
    
    list(angle=180*theta/pi, nshared=shared)
}

cosine.norm <- function(X)
# Computes the cosine norm, with some protection from zero-length norms.
{
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
}

