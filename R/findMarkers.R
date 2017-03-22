.findMarkers <- function(x, clusters, block=NULL, subset.row=NULL)
# Uses limma to find the markers that are differentially expressed between clusters,
# given a log-expression matrix and some blocking factors.
#
# written by Aaron Lun
# created 22 March 2017    
{
    # Creating a design matrix.
    clusters <- as.factor(clusters)
    design <- model.matrix(~0 + clusters)
    colnames(design) <- clust.vals <- levels(clusters)

    if (!is.null(block)) {
        # Removing the intercept, if it exists. 
        out <- qr.solve(block, cbind(rep(1, nrow(block))))
        to.drop <- abs(out) > 1e-8
        if (any(to.drop)) {
            block <- block[,-which(to.drop)[1],drop=FALSE]
        }
        design <- cbind(design, block) # Other dependencies will trigger warnings.
    }
    
    subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
    lfit <- lmFit(x[subset.row,,drop=FALSE], design)
    output <- vector("list", length(clust.vals))
    names(output) <- clust.vals  

    for (host in clust.vals) { 
        not.host <- clust.vals!=host
        targets <- clust.vals[not.host]
        all.p <- all.lfc <- vector("list", length(targets))
        names(all.p) <- names(all.lfc) <- targets
       
        con <- -diag(length(clust.vals))
        rownames(con) <- colnames(con) <- clust.vals
        con[host,] <- 1
        con <- con[,not.host,drop=FALSE]

        fit2 <- contrasts.fit(lfit, con)
        fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
        
        for (target in targets) { 
            res <- topTable(fit2, n=Inf, sort.by="none", coef=target)
            all.p[[target]] <- res$P.Value
            all.lfc[[target]] <- res$logFC
        }

        collected.ranks <- lapply(all.p, rank, ties="first")
        min.rank <- do.call(pmin, collected.ranks)
        adj.min.p <- p.adjust(do.call(pmin, all.p), n=sum(lengths(all.p)), method="BH")
        marker.set <- data.frame(Top=min.rank, Gene=rownames(x)[subset.row], FDR=adj.min.p,
                                 logFC=do.call(cbind, all.lfc), 
                                 stringsAsFactors=FALSE, check.names=FALSE)
        marker.set <- marker.set[order(marker.set$Top),]
        output[[host]] <- marker.set
    }

    return(output)
}

setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

setMethod("findMarkers", "matrix", function(x, ...) .findMarkers) 

setMethod("findMarkers", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
    if (is.null(subset.row)) { subset.row <- .spikeSubset(x, get.spikes) }
    .findMarkers(assayDataElement(x, assay), ..., subset.row=subset.row)
})                                 


