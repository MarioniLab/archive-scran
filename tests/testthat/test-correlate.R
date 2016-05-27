# This checks the correlateNull function.

# require(scran); require(testthat); source("test-correlate.R")

refnull <- function(niters, ncells, resort=TRUE) {
    rankings <- as.double(seq_len(ncells))
    shuffled <- .Call(scran:::cxx_auto_shuffle, rankings, as.integer(niters))
    out <- cor(shuffled, rankings, method="spearman")
    if (resort) { out <- sort(out) }
    out
}

set.seed(100)
ref <- refnull(1e3, 121)
set.seed(100)
out <- correlateNull(121, iters=1e3)
expect_equal(ref, as.double(out))

set.seed(100)
ref <- refnull(1e3, 12)
set.seed(100)
out <- correlateNull(12, iters=1e3)
expect_equal(ref, as.double(out))

# Checking with design matrix.

design <- model.matrix(~factor(rep(c(1,2), each=10)))
Q <- qr.Q(qr(design), complete=TRUE)
df <- nrow(design)-ncol(design)

set.seed(100)
collected <- list()
for (x in seq_len(1e3)) {
    first.half <- Q %*% c(0,0, rnorm(df))
    second.half <- Q %*% c(0, 0, rnorm(df)) 
    collected[[x]] <- cor(first.half, second.half, method="spearman")
}
out1 <- sort(unlist(collected))

set.seed(100)
out2 <- correlateNull(design=design, iters=1e3, simulate=TRUE)
expect_equal(out1, as.double(out2))
expect_equal(attr(out2, "design"), design)
expect_equal(attr(out2, "simulate"), TRUE)

# A more complicated design.

design <- model.matrix(~seq_len(10))
Q <- qr.Q(qr(design), complete=TRUE)
df <- nrow(design)-ncol(design)

set.seed(100)
collected <- list()
for (x in seq_len(1e3)) {
    first.half <- Q %*% c(0,0, rnorm(df))
    second.half <- Q %*% c(0, 0, rnorm(df)) 
    collected[[x]] <- cor(first.half, second.half, method="spearman")
}
out1 <- sort(unlist(collected))

set.seed(100)
out2 <- correlateNull(design=design, iters=1e3, simulate=TRUE)
expect_equal(out1, as.double(out2))
expect_equal(attr(out2, "design"), design)
expect_equal(attr(out2, "simulate"), TRUE)

# A one-way layout without simulating the residuals.

grouping <- rep(1:5, each=3)
design <- model.matrix(~factor(grouping))

set.seed(100)
out1 <- 0L
for (gl in table(grouping)) { 
    out1 <- out1 + refnull(1e3, gl, resort=FALSE) * gl
}
out1 <- out1/length(grouping)
out1 <- sort(out1)

set.seed(100)
out2 <- correlateNull(design=design, iters=1e3)
expect_equal(out1, as.double(out2))
expect_equal(attr(out2, "design"), design)
expect_equal(attr(out2, "simulate"), FALSE)

# Checking nonsense inputs.

expect_error(correlateNull(ncells=100, iters=0), "number of iterations should be positive")
expect_error(correlateNull(ncells=0), "number of cells should be greater than 2")
expect_error(correlateNull(ncells=100, design=design), "cannot specify both 'ncells' and 'design'")

####################################################################################################
# Checking what happens for the error-tolerant ranking.

whee <- runif(100, -1e-16, 1e-16)
set.seed(100)
r <- scran:::.tolerant_rank(whee)
set.seed(100)
r2 <- rank(integer(100), ties.method="random")
expect_identical(r, r2)

set.seed(200)
extra <- sample(10, 100, replace=TRUE)
set.seed(100)
r <- scran:::.tolerant_rank(whee + extra)
set.seed(100)
r2 <- rank(extra, ties.method="random")
expect_identical(r, r2)

####################################################################################################

checkCorrelations <- function(out, exprs, null.dist) {
    ranked.exprs <- apply(exprs, 1, FUN=scran:::.tolerant_rank)
    colnames(ranked.exprs) <- rownames(exprs)
    assembled.pval <- assembled.rho <- numeric(nrow(out))
    for (p in seq_along(assembled.rho)) { 
        assembled.rho[p] <- cor(ranked.exprs[,out$gene1[p]], ranked.exprs[,out$gene2[p]], method="spearman")        
        assembled.pval[p] <- min(sum(null.dist <= assembled.rho[p] + 1e-8), sum(null.dist >= assembled.rho[p] - 1e-8))
    }
    assembled.pval <- 2*(assembled.pval + 1)/(length(null.dist)+1)
    return(data.frame(rho=assembled.rho, pvalue=assembled.pval, FDR=p.adjust(assembled.pval, method="BH")))
}

####################################################################################################

# This checks the correlatePairs function.

set.seed(10000)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

# With a pre-specified distribution.

nulls <- sort(runif(1e6, -1, 1))

set.seed(100)
out <- correlatePairs(X, null.dist=nulls)
set.seed(100)
ref <- checkCorrelations(out, X, null.dist=nulls)

expect_equal(out$rho, ref$rho)
expect_equal(out$p.value, ref$pvalue)
expect_equal(out$FDR, ref$FDR)

# Without a pre-specified distribution.

set.seed(10001)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))

set.seed(100)
out <- correlatePairs(X)
set.seed(100)
nulls <- correlateNull(Ncells)
ref <- checkCorrelations(out, X, null.dist=nulls)

expect_equal(out$rho, ref$rho)
expect_equal(out$p.value, ref$pvalue)
expect_equal(out$FDR, ref$FDR)

# With a design matrix.

set.seed(10002)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))
grouping <- factor(rep(c(1,2), each=50))
design <- model.matrix(~grouping)

set.seed(200)
nulls <- correlateNull(design=design, iter=1e3, simulate=TRUE)
set.seed(100) # Need because of random ranking.
out <- correlatePairs(X, design=design, null=nulls, simulate=TRUE)
fit <- lm.fit(x=design, y=t(X))
exprs <- t(fit$residual)
set.seed(100)
ref <- checkCorrelations(out, exprs, null.dist=nulls)

expect_equal(out$rho, ref$rho)
expect_equal(out$p.value, ref$pvalue)
expect_equal(out$FDR, ref$FDR)

# Repeating without simulation (need to use a normal matrix to avoid ties;
# bplapply mucks up the seed).

set.seed(200)
X[] <- rnorm(length(X))
nulls <- correlateNull(design=design, iter=1e3)

out <- correlatePairs(X, design=design, null=nulls)
collected.rho <- 0L
for (group in split(seq_along(grouping), grouping)) { 
    ref <- checkCorrelations(out, X[,group], null.dist=nulls)
    collected.rho <- collected.rho + length(group)/length(grouping) * ref$rho
}

expect_equal(out$rho, collected.rho)

collected.p <- numeric(length(collected.rho)) 
for (x in seq_along(collected.rho)) { 
    collected.p[x] <- min(sum(nulls <= collected.rho[x] + 1e-8), sum(nulls >= collected.rho[x] - 1e-8))
}
collected.p <- 2*(collected.p + 1)/(length(nulls)+1)
expect_equal(out$p.value, collected.p)

# Checking that it works with a SCESet object.

set.seed(10003)
Ngenes <- 20
Ncells <- 100
X <- log(matrix(rpois(Ngenes*Ncells, lambda=10), nrow=Ngenes) + 1)
rownames(X) <- paste0("X", seq_len(Ngenes))
nulls <- sort(runif(1e6, -1, 1))

set.seed(100)
ref <- correlatePairs(X, null.dist=nulls)
set.seed(100)
X2 <- newSCESet(exprsData=data.frame(X))
out <- correlatePairs(X2, null.dist=nulls)
expect_equal(out, ref)

# With spikes.

isSpike(X2) <- rbinom(Ngenes, 1, 0.6)==0L
set.seed(100)
ref <- correlatePairs(exprs(X2)[!isSpike(X2),,drop=FALSE], null.dist=nulls)
set.seed(100)
out <- correlatePairs(X2, null.dist=nulls)
expect_equal(out, ref)

# Checking nonsense inputs.

expect_error(correlatePairs(X[0,], nulls), "need at least two genes to compute correlations")
expect_error(correlatePairs(X[,0], nulls), "number of cells should be greater than 2")
out <- correlatePairs(X, numeric(0))
expect_equal(out$p.value, rep(1, nrow(out)))

####################################################################################################
# A high-level test, to make sure that our stuff gives a uniform distribution of p-values.
#
# design <- model.matrix(~factor(rep(1:5, 2)))
# y <- matrix(rnorm(1000, mean=rep(1:5, 5), sd=2), ncol=10, byrow=TRUE)
# null <- correlateNull(ncol(y))
# out <- correlatePairs(y, design=design, null=null)
# plot(log10(sort(out$p.value)/1:nrow(out)*nrow(out))) # wrong
# null <- correlateNull(design=design, simulate=TRUE)
# out <- correlatePairs(y, design=design, null=null)
# plot(log10(sort(out$p.value)/1:nrow(out)*nrow(out))) # right

