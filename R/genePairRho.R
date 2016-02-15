buildNullRho <- function(ncells, iters=1e6) 
# This builds a null distribution for the modified Spearman's rho.
#
# written by Aaron Lun
# created 10 February 2016    
{
    out <- .Call("get_null_rho", as.integer(ncells), as.integer(iters), PACKAGE="scran") 
    if (is.character(out)) { 
        stop(out)
    }
    out <- sort(out)
    return(out)  
}

genePairRho <- function(exprs, null.dist=NULL, BPPARAM=bpparam(), use.names=TRUE)
# This calculates a (modified) Spearman's rho for each pair of genes.
#
# written by Aaron Lun
# created 10 February 2016    
{
    exprs <- as.matrix(exprs)
    ncells <- ncol(exprs)
    if (is.null(null.dist)) { 
        null.dist <- buildNullRho(ncells)
    } else {
        null.dist <- as.double(null.dist)
    }
    if (is.unsorted(null.dist)) { 
        null.dist <- sort(null.dist)
    }
    ranked.exprs <- apply(exprs, 1, FUN=rank, ties.method="random")

    # Generating all pairs of genes
    ngenes <- nrow(exprs)
    all.pairs <- combn(ngenes, 2L)
    gene1 <- all.pairs[1,]
    gene2 <- all.pairs[2,]

    # Assigning to cores 
    ncores <- bpworkers(BPPARAM)
    starting <- seq(from=1, to=length(gene1)+1, length.out=ncores+1)
    starting <- unique(starting[seq_len(ncores)])
    ending <- c((starting - 1)[-1], length(gene1))

    # Running through each set of jobs 
    out <- bplapply(seq_len(ncores), FUN=function(core) {
        to.use <- starting[core]:ending[core]
        .Call("compute_rho", gene1[to.use], gene2[to.use], ncells, ranked.exprs, null.dist, PACKAGE="scran")
    }, BPPARAM=BPPARAM)

    # Peeling apart the output
    all.rho <- all.pval <- list()
    for (x in seq_along(out)) {
        current <- out[[x]]
        if (is.character(current)) { stop(current) }
        all.rho[[x]] <- current[[1]]
        all.pval[[x]] <- current[[2]]
    }
    all.pval <- unlist(all.pval)
    all.rho <- unlist(all.rho)

    # Returning some useful output
    if (use.names) {
        newnames <- rownames(exprs)
        if (!is.null(newnames)) {
            gene1 <- newnames[gene1]
            gene2 <- newnames[gene2]
        }
    }

    out <- data.frame(gene1=gene1, gene2=gene2, rho=all.rho, p.value=all.pval, 
                      FDR=p.adjust(all.pval, method="BH"), stringsAsFactors=FALSE)
    out <- out[order(out$p.value, -abs(out$rho)),]
    return(out)
}
