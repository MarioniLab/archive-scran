correlateNull <- function(ncells, iters=1e6, design=NULL, simulate=FALSE) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 27 May 2016
{
    if (!is.null(design)) { 
        if (!missing(ncells)) { 
            stop("cannot specify both 'ncells' and 'design'")
        }

        groupings <- .isOneWay(design)
        if (is.null(groupings) || simulate) { 
            # Using simulated residual effects if the design matrix is not a one-way layout (or if forced by simulate=TRUE).
            out <- .Call(cxx_get_null_rho_design, design, ncol(design), nrow(design), as.integer(iters))
            if (is.character(out)) { 
                stop(out)
            }
        } else {
            # Otherwise, estimating the correlation as a weighted mean of the correlations in each group.
            # This avoids the need for the normality assumption in the residual effect simulation.
            out <- 0
            for (gr in groupings) {
                out.g <- .Call(cxx_get_null_rho, length(gr), as.integer(iters))
                if (is.character(out.g)) { 
                    stop(out.g)
                }
                out <- out + out.g * length(gr)
            }
            out <- out/nrow(design)
        }
        attrib <- list(design=design, simulate=simulate)

    } else {
        out <- .Call(cxx_get_null_rho, as.integer(ncells), as.integer(iters))
        if (is.character(out)) { 
            stop(out)
        }
        attrib <- NULL
    }

    # Storing attributes, to make sure it matches up.
    out <- sort(out)
    attributes(out) <- attrib
    return(out)  
}

.isOneWay <- function(design) {
    if (nrow(design) <= ncol(design)) {
        stop("design matrix has no residual degrees of freedom")
    }
    group <- designAsFactor(design)
    if (nlevels(group) == ncol(design)) {
        # Stripping out groups with only one level.
        groupings <- split(seq_len(nrow(design)), group)
        groupings[lengths(groupings)==1L] <- NULL
        return(groupings)
    } 
    return(NULL)
}

setGeneric("correlatePairs", function(x, ...) standardGeneric("correlatePairs"))

setMethod("correlatePairs", "matrix", function(x, null.dist=NULL, design=NULL, BPPARAM=bpparam(), use.names=TRUE, tol=1e-8, simulate=FALSE)
# This calculates a (modified) Spearman's rho for each pair of genes.
#
# written by Aaron Lun
# created 10 February 2016
# last modified 27 May 2016
{
    if (!is.null(design)) { 
        groupings <- .isOneWay(design)
        if (is.null(groupings) || simulate) { 
            fit <- lm.fit(y=t(x), x=design)
            exprs.list <- list(t(fit$residuals))
        } else {
            exprs.list <- list()
            for (g in seq_along(groupings)) {
                exprs.list[[g]] <- x[,groupings[[g]],drop=FALSE]
            }
        }
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(design=design, simulate=simulate)
        }
    } else {
        exprs.list <- list(x)
        if (is.null(null.dist)) { 
            null.dist <- correlateNull(ncol(x))
        } 
    }

    # Checking that the null distribution is sensible.
    if (!identical(design, attr(null.dist, "design"))) { 
        stop("'design' is not the same as that used to generate 'null.dist'")
    }
    if (!is.null(design)) { 
        if (!identical(simulate, attr(null.dist, "simulate"))) {
            stop("'simulate' is not the same as that used to generate 'null.dist'")
        }
    }
    null.dist <- as.double(null.dist)
    if (is.unsorted(null.dist)) { 
        null.dist <- sort(null.dist)
    }

    # Generating all pairs of genes
    ngenes <- nrow(x)
    if (ngenes < 2L) { stop("need at least two genes to compute correlations") }
    all.pairs <- combn(ngenes, 2L)
    gene1 <- all.pairs[1,]
    gene2 <- all.pairs[2,]

    # Iterating through all subgroups (for one-way layouts; otherwise, this is a loop of length 1).
    all.rho <- 0L
    for (exprs in exprs.list) { 
        ncells <- ncol(exprs)

        # Ranking genes, in an error-tolerant way. This avoids getting untied rankings for zeroes
        # (which should have the same value +/- precision, as the prior count scaling cancels out).
        ranked.exprs <- apply(exprs, 1, FUN=.tolerant_rank, tol=tol)

        # Running through each set of jobs 
        workass <- .workerAssign(length(gene1), BPPARAM)
        out <- bplapply(seq_along(workass$start), FUN=.get_correlation,
            work.start=workass$start, work.end=workass$end,
            gene1=gene1, gene2=gene2, ncells=ncells, ranked.exprs=ranked.exprs, 
            BPPARAM=BPPARAM)

        # Peeling apart the output.
        current.rho <- list()
        for (i in seq_along(out)) {
            current <- out[[i]]
            if (is.character(current)) { stop(current) }
            current.rho[[i]] <- current
        }
        current.rho <- unlist(current.rho)

        # Adding a weighted value to the final.
        all.rho <- all.rho + current.rho * (ncells/ncol(x))
    }

    # Estimating the p-values (need to shift values to break ties conservatively by increasing the p-value).
    left <- findInterval(all.rho + 1e-8, null.dist)
    right <- length(null.dist) - findInterval(all.rho - 1e-8, null.dist)
    all.pval <- (pmin(left, right)+1)*2/(length(null.dist)+1)
    all.pval <- pmin(all.pval, 1)

    # Returning some useful output
    newnames <- NULL
    if (is.logical(use.names)) {
        if (use.names) {
            newnames <- rownames(x)
        }
    } else if (is.character(use.names)) {
        if (length(use.names)!=nrow(x)) {
            stop("length of 'use.names' does not match 'x' nrow")
        }
        newnames <- use.names
    }
    if (!is.null(newnames)) {
        gene1 <- newnames[gene1]
        gene2 <- newnames[gene2]
    }

    out <- data.frame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, 
                      FDR=p.adjust(all.pval, method="BH"), stringsAsFactors=FALSE)
    out <- out[order(out$p.value, -abs(out$rho)),]
    rownames(out) <- NULL
    return(out)
})

.workerAssign <- function(njobs, BPPARAM) {
    ncores <- bpworkers(BPPARAM)
    starting <- as.integer(seq(from=1, to=njobs+1, length.out=ncores+1))
    starting <- unique(starting[seq_len(ncores)])
    ending <- c((starting - 1L)[-1], njobs)
    return(list(start=starting, end=ending))
}

.get_correlation <- function(core, work.start, work.end, gene1, gene2, ncells, ranked.exprs) {
    to.use <- work.start[core]:work.end[core]
    .Call(cxx_compute_rho, gene1[to.use], gene2[to.use], ncells, ranked.exprs)
}

.tolerant_rank <- function(y, tol=1e-6) {
    if (!length(y)) { return(integer(0)) }
    o <- order(y)                          
    rle.out <- rle(y[o])
    okay <- c(TRUE, diff(rle.out$values) > tol)
    to.use <- cumsum(okay)
    rle.out$values <- rle.out$values[okay][to.use]
    y[o] <- inverse.rle(rle.out)
    rank(y, ties.method="random")
}

setMethod("correlatePairs", "SCESet", function(x, ..., assay="exprs", get.spikes=FALSE) {
    correlatePairs(.getUsedMatrix(x, assay, get.spikes), ...)             
})

