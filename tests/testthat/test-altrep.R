# This tests that various functions are applicable with alternative matrix representations.
# library(scran); library(testthat); source("test-altrep.R")

set.seed(99999)
library(Matrix)
X <- as(matrix(rpois(100000, lambda=1), ncol=100), "dgCMatrix")
X_ <- as.matrix(X)

library(HDF5Array)
Y <- as(matrix(rpois(100000, lambda=5), ncol=100), "HDF5Array")
Y_ <- as.matrix(Y)

test_that("cyclone runs properly", {
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    rownames(X) <- rownames(X_) <- sample(mm.pairs$G1[,1], nrow(X))
    rownames(Y) <- rownames(Y_) <- sample(mm.pairs$G1[,1], nrow(Y))

    set.seed(100)
    assignments1 <- cyclone(X[,1:10], mm.pairs)
    set.seed(100)
    assignments2 <- cyclone(X_[,1:10], mm.pairs)
    expect_identical(assignments1, assignments2)

    set.seed(100)
    assignments1 <- cyclone(Y[,1:10], mm.pairs)
    set.seed(100)
    assignments2 <- cyclone(Y_[,1:10], mm.pairs)
    expect_identical(assignments1, assignments2)
})

test_that("computeSumFactors runs properly", {
    sf1 <- computeSumFactors(X_)
    sf2 <- computeSumFactors(X)
    expect_identical(sf1, sf2)

    sf1 <- computeSumFactors(Y_)
    sf2 <- computeSumFactors(Y)
    expect_identical(sf1, sf2)
})

test_that("Variance estimation runs properly", {
    fit1 <- trendVar(X, parametric=FALSE) # because parametric=TRUE doesn't work properly with non-log values.
    dec1 <- decomposeVar(X, fit1)
    fit2 <- trendVar(X_, parametric=FALSE)
    dec2 <- decomposeVar(X_, fit1)
    expect_equal(fit1, fit2)
    expect_equal(dec1, dec2)

    fit1 <- trendVar(Y, parametric=FALSE)
    dec1 <- decomposeVar(Y, fit1)
    fit2 <- trendVar(Y_, parametric=FALSE)
    dec2 <- decomposeVar(Y_, fit1)
    expect_equal(fit1, fit2)
    expect_equal(dec1, dec2)
})

test_that("correlatePairs runs properly", {
    set.seed(100) 
    ref <- correlatePairs(X[1:10,])
    set.seed(100) 
    alt <- correlatePairs(X_[1:10,])
    expect_equal(ref, alt)

    set.seed(100) 
    ref <- correlatePairs(Y[20:50,])
    set.seed(100) 
    alt <- correlatePairs(Y_[20:50,])
    expect_equal(ref, alt)
})
