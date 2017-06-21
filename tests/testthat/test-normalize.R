# This tests out the normalization methods in scran - specifically, compute*Factors and normalize().
# require(scran); require(testthat); source("test-normalize.R")

set.seed(20000)
ncells <- 200
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)

out <- computeSumFactors(dummy)
expect_equal(out, count.sizes/mean(count.sizes))

# Adding some DE genes.

count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
is.de <- sample(ngenes, 100)
dummy[is.de,] <- rnbinom(ncells*length(is.de), mu=100, size=1)
out <- computeSumFactors(dummy)
expect_equal(out, count.sizes/mean(count.sizes))

count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
is.de <- sample(ngenes, 400)
dummy[is.de,] <- rnbinom(ncells*length(is.de), mu=100, size=1)
out <- computeSumFactors(dummy)
expect_equal(out, count.sizes/mean(count.sizes))

# Checking the ring construction.

lib.sizes <- runif(100)
out <- scran:::.generateSphere(lib.sizes)
r <- rank(lib.sizes)
expect_identical(r[out][1:50], 1:50*2-1L) # All odd ranks
expect_identical(r[out][51:100], 50:1*2) # All even ranks
expect_identical(r[out][1:100], r[out][101:200]) # Repeated for easy windowing

lib.sizes <- runif(101)
out <- scran:::.generateSphere(lib.sizes)
r <- rank(lib.sizes)
expect_identical(r[out][1:51], 1:51*2-1L) # All odd ranks
expect_identical(r[out][52:101], 50:1*2) # All even ranks
expect_identical(r[out][1:101], r[out][102:202]) # Repeated for easy windowing

# Creating an R-only implementation for comparison.

sumInR <- function(x, sizes, center=TRUE) {
    lib.sizes <- colSums(x)
    x <- t(t(x)/lib.sizes)
    ref <- rowMeans(x)

    keep <- ref > 0
    ref <- ref[keep]
    x <- x[keep,,drop=FALSE]

    ncells <- length(lib.sizes)
    o <- scran:::.generateSphere(lib.sizes)
    all.mat <- all.vals <- vector("list", sum(sizes)*ncells) 
    i <- 1L

    for (s in sizes) {
        for (w in seq_len(ncells)) {
            chosen <- o[w+seq_len(s)-1L]

            current <- integer(ncells)
            current[chosen] <- 1L
            all.mat[[i]] <- current

            ratios <- rowSums(x[,chosen,drop=FALSE])/ref
            all.vals[[i]] <- median(ratios)
            i <- i+1L
        }
    }

    # Adding the low weight additional equations.
    extra.mat <- diag(ncells)*sqrt(scran:::LOWWEIGHT)
    extra.val <- apply(x/ref, 2, median)*sqrt(scran:::LOWWEIGHT)
    final.mat <- rbind(do.call(rbind, all.mat), extra.mat)
    final.val <- c(unlist(all.vals), extra.val)

    nf <- solve(qr(final.mat), final.val)
    sf <- nf * lib.sizes
    if (center) {
        sf <- sf/mean(sf)
    }
    return(sf)
}

ngenes2 <- 200
x <- matrix(rpois(ngenes2*ncells, lambda=10), nrow=ngenes2, ncol=ncells)
sizes <- seq(20, 100, 5)
ref <- sumInR(x, sizes)
obs <- computeSumFactors(x, sizes=sizes)
expect_equal(ref, obs)

x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
x[sample(nrow(x), 100),] <- 0L # Throwing in some zeroes.
ref <- sumInR(x, sizes)
obs <- computeSumFactors(x, sizes=sizes)
expect_equal(ref, obs)

x <- matrix(rpois(ncells*ngenes2, lambda=10), nrow=ngenes2, ncol=ncells)
subset.row <- sample(nrow(x), 100)
ref <- sumInR(x[subset.row,,drop=FALSE], sizes)
obs <- computeSumFactors(x, subset.row=subset.row, sizes=sizes)
expect_equal(ref, obs)

####################################################################################################

# Trying it out with other options.

dummy <- matrix(rpois(ncells*ngenes, lambda=10), nrow=ngenes, ncol=ncells)
out <- computeSumFactors(dummy)
if (.Platform$OS.type!="windows") { # Because limSolve doesn't build on Windows, apparently.
outx <- computeSumFactors(dummy, positive=TRUE)
expect_true(all(abs(outx -  out) < 1e-3)) # need to be a bit generous here, the solution code is different.
}
expect_warning(outx <- computeSumFactors(dummy, errors=TRUE), "errors=TRUE is no longer supported")
expect_equal(as.numeric(outx), out)

# Checking the the clustering works as expected.

clusters <- rep(1:2, 100)
sizes <- seq(20, 100, 5)
obs <- computeSumFactors(dummy, sizes=sizes, cluster=clusters)
ref1 <- sumInR(dummy[,clusters==1], sizes, center=FALSE) # Avoid centering, as this destroys relative information.
ref2 <- sumInR(dummy[,clusters==2], sizes, center=FALSE)

adj <- t(t(dummy)/colSums(dummy))
pseudo1 <- rowMeans(adj[,clusters==1])
pseudo2 <- rowMeans(adj[,clusters==2])
rescale2 <- median(pseudo2/pseudo1)

ref <- numeric(ncells)
ref[clusters==1] <- ref1
ref[clusters==2] <- ref2*rescale2
ref <- ref/mean(ref)
expect_equal(ref, obs)

# Trying it out on a SCESet object.

count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
out <- computeSumFactors(X)
expect_equal(unname(sizeFactors(out)), computeSumFactors(dummy))

# Throwing in some silly inputs.

expect_error(computeSumFactors(dummy[,0,drop=FALSE]), "not enough cells in each cluster")
expect_error(computeSumFactors(dummy[0,,drop=FALSE]), "cells should have non-zero library sizes")
expect_error(computeSumFactors(dummy, sizes=c(10, 10, 20)), "'sizes' is not unique")
expect_error(computeSumFactors(dummy, clusters=integer(0)), "'x' ncols is not equal to 'clusters' length")

####################################################################################################

# Checking out what happens with clustering.

set.seed(20001)
ncells <- 700
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- seq_len(ngenes)/100
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0

emp.clusters <- quickCluster(dummy)
expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)
shuffled <- c(1:50, 301:350, 601:650)
expect_identical(quickCluster(dummy, subset.row=shuffled), emp.clusters)

# Checking out the ranks.

emp.ranks <- quickCluster(dummy, get.ranks=TRUE)
ref <- apply(dummy, 2, FUN=function(y) {
    r <- rank(y)
    r <- r - mean(r)
    r/sqrt(sum(r^2))/2
})
expect_equal(emp.ranks, ref)

emp.ranks <- quickCluster(dummy, get.ranks=TRUE, subset.row=shuffled)
ref <- apply(dummy, 2, FUN=function(y) {
    r <- rank(y[shuffled])
    r <- r - mean(r)
    r/sqrt(sum(r^2))/2
})
expect_equal(emp.ranks, ref)

# Checking out that clustering is consistent with that based on correlations.

set.seed(200011)
mat <- matrix(rpois(10000, lambda=5), nrow=20)
obs <- quickCluster(mat)

refM <- sqrt(0.5*(1 - cor(mat, method="spearman")))
distM <- as.dist(refM) 
htree <- hclust(distM, method='ward.D2')
clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=200, distM=refM, verbose=0))
expect_identical(clusters, as.integer(obs))

obs <- quickCluster(mat, min.size=50)
clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=refM, verbose=0))
expect_identical(clusters, as.integer(obs))

mat <- matrix(rpois(10000, lambda=5), nrow=20)
subset.row <- 15:1 # With subsetting
refM <- sqrt(0.5*(1 - cor(mat[subset.row,], method="spearman")))
distM <- as.dist(refM) 
htree <- hclust(distM, method='ward.D2')
obs <- quickCluster(mat, min.size=50, subset.row=subset.row)
clusters <- unname(dynamicTreeCut::cutreeDynamic(htree, minClusterSize=50, distM=refM, verbose=0))
expect_identical(clusters, as.integer(obs))

# Other checks

expect_identical(length(quickCluster(mat, method="igraph", d=NA)), ncol(mat)) # Checking that the dimensions are correct for igraph.
suppressWarnings(expect_false(identical(quickCluster(mat), quickCluster(mat[subset.row,])))) # Checking that subsetting gets different results.
suppressWarnings(expect_identical(quickCluster(mat, subset.row=subset.row), quickCluster(mat[subset.row,]))) # Checking that subset.row works.

# Checking out what happens with silly inputs.

expect_error(quickCluster(dummy[0,]), "rank variances of zero detected for a cell")
expect_error(quickCluster(dummy[,0]), "fewer cells than the minimum cluster size")

leftovers <- 100
expect_warning(forced <- quickCluster(dummy[,c(which(known.clusters==1), which(known.clusters==2), which(known.clusters==3)[1:leftovers])]), 
               sprintf("%i cells were not assigned to any cluster", leftovers))
expect_identical(as.character(tail(forced, leftovers)), rep("0", leftovers))

# Seeing how it interacts with the normalization method.

out <- computeSumFactors(dummy, cluster=known.clusters)
expect_equal(out, count.sizes/mean(count.sizes)) # Even though there is a majority of DE, each pair of clusters is still okay.

out1 <- computeSumFactors(dummy, cluster=known.clusters, ref=1)
expect_equal(out, out1)
out2 <- computeSumFactors(dummy, cluster=known.clusters, ref=2)
expect_equal(out, out2)
out3 <- computeSumFactors(dummy, cluster=known.clusters, ref=3)
expect_equal(out, out3)

expect_error(computeSumFactors(dummy, cluster=known.clusters, ref=0), "'ref.clust' value not in 'clusters'")

# Trying it out on a SCESet object.

set.seed(20002)
count.sizes <- rnbinom(ncells, mu=100, size=5)
multiplier <- sample(seq_len(ngenes)/100)
dummy <- outer(multiplier, count.sizes)

known.clusters <- sample(3, ncells, replace=TRUE)
dummy[1:300,known.clusters==1L] <- 0
dummy[301:600,known.clusters==2L] <- 0  
dummy[601:900,known.clusters==3L] <- 0

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
emp.clusters <- quickCluster(X)
expect_true(length(unique(paste0(known.clusters, emp.clusters)))==3L)

####################################################################################################

# Checking out the behaviour of the computeSpikeFactors function.

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
is.spike <- rbinom(ngenes, 1, 0.7)==0L
dummy[is.spike,] <- matrix(rnbinom(sum(is.spike)*ncells, mu=20, size=5), ncol=ncells, nrow=sum(is.spike), byrow=TRUE)

rownames(dummy) <- paste0("X", seq_len(ngenes))
X <- newSCESet(countData=data.frame(dummy))
X <- calculateQCMetrics(X, list(MySpike=is.spike))
setSpike(X) <- "MySpike"
out <- computeSpikeFactors(X)
ref <- colSums(dummy[is.spike,])
expect_equal(unname(sizeFactors(out)), ref/mean(ref))
expect_equal(sizeFactors(out), sizeFactors(out, type="MySpike"))

# Checking out what happens when you have multiple spike-ins supplied.
X2 <- newSCESet(countData=data.frame(dummy))
subset <- split(which(is.spike), rep(1:2, length.out=sum(is.spike)))
X2 <- calculateQCMetrics(X2, list(MySpike=subset[[1]], SecondSpike=subset[[2]]))
setSpike(X2) <- c("MySpike", "SecondSpike")

out.sub <- computeSpikeFactors(X2, type="MySpike") # Sanity check, to make sure that it's calculating it differently for each spike-in.
subref <- colSums(dummy[subset[[1]],])
expect_equal(unname(sizeFactors(out.sub)), subref/mean(subref))
expect_equal(sizeFactors(out.sub), sizeFactors(out.sub, type="MySpike"))
expect_warning(sizeFactors(out.sub, type="SecondSpike"), "'sizeFactors' have not been set for 'SecondSpike'")

out2 <- computeSpikeFactors(X2)
expect_equal(sizeFactors(out), sizeFactors(out2))
expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))

out2 <- computeSpikeFactors(X2, type=c("MySpike", "SecondSpike"))
expect_equal(sizeFactors(out), sizeFactors(out2))
expect_equal(sizeFactors(out), sizeFactors(out2, type="MySpike"))
expect_equal(sizeFactors(out), sizeFactors(out2, type="SecondSpike"))

# Checking out the general use function.
sizeFactors(X) <- 1
out <- computeSpikeFactors(X, general.use=FALSE)
expect_equal(unname(sizeFactors(out)), rep(1, ncells))
expect_equal(unname(sizeFactors(out, type="MySpike")), ref/mean(ref))

# Breaks if you try to feed it silly inputs.
expect_warning(out <- computeSpikeFactors(X[0,]), "zero spike-in counts during spike-in normalization")
expect_identical(unname(sizeFactors(out)), rep(NaN, ncol(out)))
out <- computeSpikeFactors(X[,0])
expect_identical(unname(sizeFactors(out)), numeric(0))

####################################################################################################

# Checking out the behaviour of the normalize() function.

set.seed(20003)
ncells <- 200
ngenes <- 1000
dummy <- matrix(rnbinom(ncells*ngenes, mu=100, size=5), ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))
colnames(dummy) <- paste0("Y", seq_len(ncells))
X <- newSCESet(countData=dummy)

ref <- colSums(dummy)
sizeFactors(X) <- ref
out <- normalize(X)
sf <- ref/mean(ref)
expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1))
out <- normalize(X, logExprsOffset=3)
expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+3))

ref <- runif(ncells, 10, 20)
sizeFactors(X) <- ref
out <- normalize(X)
sf <- ref/mean(ref)
expect_equivalent(exprs(out), log2(t(t(dummy)/sf)+1)) 

expect_equivalent(sf, sizeFactors(out))
Xb <- X
sizeFactors(Xb) <- ref
outb <- normalize(Xb, centre_size_factors=FALSE)
expect_equivalent(ref, sizeFactors(outb))
expect_equivalent(exprs(out), exprs(outb))

# Now adding some controls.

chosen <- rbinom(ngenes, 1, 0.7)==0L
X <- calculateQCMetrics(X, feature_controls=list(whee=chosen))
X3 <- normalize(X)
expect_equal(exprs(out), exprs(X3))

Xb <- X
setSpike(Xb) <- "whee"
expect_warning(X3b <- normalize(Xb), "spike-in transcripts in 'whee'")
expect_equal(exprs(X3b), exprs(X3))

sizeFactors(X, type="whee") <- colSums(counts(X)[chosen,])
expect_warning(X4 <- normalize(X), NA) # i.e., no warning.
expect_equivalent(exprs(out)[!chosen,], exprs(X4)[!chosen,])
ref <- sizeFactors(X, type="whee")
sf <- ref/mean(ref)
expect_equivalent(exprs(X4)[chosen,], log2(t(t(dummy[chosen,])/sf)+1))

expect_equivalent(sizeFactors(X4, type="whee"), sf)
X4b <- normalize(X, centre_size_factors=FALSE)
expect_equivalent(sizeFactors(X4b, type="whee"), sizeFactors(X, type="whee"))
expect_equivalent(exprs(X4), exprs(X4b))

# Checking out silly inputs.

expect_equal(unname(dim(normalize(X[,0,drop=FALSE]))), c(ngenes, 0L))
expect_equal(unname(dim(normalize(X[0,,drop=FALSE]))), c(0L, ncells)) 


