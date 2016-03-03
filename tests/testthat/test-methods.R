# This tests the various SCESet methods in scran.

require(scran); require(testthat);

set.seed(20000)
ncells <- 200
ngenes <- 1000
count.sizes <- rnbinom(ncells, mu=100, size=5)
dummy <- matrix(count.sizes, ncol=ncells, nrow=ngenes, byrow=TRUE)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- newSCESet(countData=data.frame(dummy))
is.spike <- rbinom(ngenes, 1, 0.5)==0L
isSpike(X) <- is.spike
expect_identical(isSpike(X), is.spike)

sf <- runif(ncells, 0.5, 1.5)
sizeFactors(X) <- sf
expect_identical(sf, unname(sizeFactors(X)))
expect_identical(colnames(X), names(sizeFactors(X)))

# Checking silly inputs

expect_error(isSpike(X) <- "whee", "must be a logical vector")
expect_error(sizeFactors(X) <- "whee", "size factors should be numeric")
expect_identical(isSpike(X[0,]), logical(0))
expect_identical(unname(sizeFactors(X[,0])), numeric(0))

