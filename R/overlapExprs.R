.overlapExprs <- function(x, groups, design=NULL, residuals=FALSE, tol=1e-8, subset.row=NULL)
# Computes the gene-specific overlap in expression profiles between two groups of cells.
# This aims to determine whether two distributions of expression values are well-separated.    
# 
# written by Aaron Lun
# created 17 April 2017    
{
    compute.residuals <- FALSE
    if (!is.null(design)) { 
        QR <- qr(design, LAPACK=TRUE)
        groupings <- .isOneWay(design)
        if (is.null(groupings) || residuals) { 
            compute.residuals <- TRUE
            groupings <- list(seq_len(ncol(x)))
        } 
    } else {
        groupings <- list(seq_len(ncol(x)))
    }

    # Checking dimensions.
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    ncells <- ncol(x)
    if (length(groups)!=ncells) { 
        stop("length of 'groups' not equal to number of cells")
    }
    if (!is.null(design) && nrow(design)!=ncells) { 
        stop("'nrow(design)' not equal to number of cells")
    }

    # Computing residuals; also replacing the subset vector, as it'll already be subsetted.
    if (compute.residuals) { 
        use.x <- .Call(cxx_get_residuals, x, QR$qr, QR$qraux, subset.row - 1L)
        if (is.character(use.x)) { stop(use.x) }
        use.subset.row <- seq_len(nrow(use.x)) - 1L
    } else {
        use.x <- x[subset.row,,drop=FALSE]
        use.subset.row <- subset.row - 1L
    }

    # Setting up the output matrices.
    unique.groups <- sort(unique(groups))
    ngroups <- length(unique(groups))
    output <- used.cells <- vector("list", ngroups)
    ngenes <- length(subset.row)
    for (g in seq_len(ngroups)) { 
        temp <- matrix(0, ngenes, ngroups-1L) 
        colnames(temp) <- unique.groups[-g]
        rownames(temp) <- rownames(x)[subset.row]
        output[[g]] <- temp
        temp.n <- integer(ngroups-1)
        names(temp.n) <- colnames(temp)
        used.cells[[g]] <- temp.n
    }
    names(output) <- names(used.cells) <- unique.groups

    # Running through each blocking level and computing the proportions.
    for (subset.col in groupings) { 
        cur.groups <- groups[subset.col]
        by.group <- split(subset.col, cur.groups)
        if (length(by.group)==1L) { next }

        for (b1 in seq_along(by.group)) {
            g1 <- by.group[[b1]] 
            for (b2 in seq_len(b1-1)) {  
                g2 <- by.group[[b2]]
                cur.output <- numeric(ngenes)

                for (i in seq_len(ngenes)) {
                    e1 <- use.x[i,g1]
                    e2 <- use.x[i,g2]
                    de <- outer(e1, e2, `-`)
                    neq <- sum(abs(de) <= tol)
                    ngr <- sum(de > tol)
                    cur.output[i] <- neq*0.5 + ngr
                }
                cur.output <- cur.output / (length(g1) * length(g2))

                # Storing in output.
                cur.used <- length(g1) + length(g2)
                name1 <- names(by.group)[b1] 
                name2 <- names(by.group)[b2] 
                output[[name1]][,name2] <- output[[name1]][,name2] + cur.output * cur.used
                output[[name2]][,name1] <- output[[name2]][,name1] + (1 - cur.output) * cur.used
                used.cells[[name1]][[name2]] <- used.cells[[name1]][[name2]] + cur.used
                used.cells[[name2]][[name1]] <- used.cells[[name2]][[name1]] + cur.used
            }
        }
    }

    # Normalizing the output matrices.
    for (g in names(output)) { 
        output[[g]] <- t(t(output[[g]])/used.cells[[g]])
    }
    return(output)
}

setGeneric("overlapExprs", function(x, ...) standardGeneric("overlapExprs"))

setMethod("overlapExprs", "matrix", .overlapExprs)

setMethod("overlapExprs", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { subset.row <- .spikeSubset(x, get.spikes) }
    .overlapExprs(assayDataElement(x, assay), ..., subset.row=subset.row)
})                                 

