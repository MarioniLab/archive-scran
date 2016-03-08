# A function to convert between formats

setGeneric("convertTo", function(x, ...) standardGeneric("convertTo"))

setMethod("convertTo", "SCESet", function(x, type=c("edgeR", "DESeq2", "monocle"),
    fData.col=NULL, pData.col=NULL, ..., assay, logged.exprs=TRUE, get.spikes=FALSE) {

    type <- match.arg(type)
    if (type=="edgeR" || type=="DESeq2") { 
        sf <- suppressWarnings(sizeFactors(x))
        fd <- fData(x)[,fData.col,drop=FALSE]
        pd <- pData(x)[,pData.col,drop=FALSE] 
    } else if (type=="monocle") {
        fd <- featureData(x)[,fData.col,drop=FALSE]
        pd <- phenoData(x)[,pData.col,drop=FALSE] 
    }

    if (type=="edgeR") {
        if (missing(assay)) { assay <- "counts" }
        y <- DGEList(assayDataElement(x, assay), ...)
        if (ncol(fd)) { y$genes <- fd }
        if (!is.null(sf)) { 
            nf <- log(sf/y$samples$lib.size)
            nf <- exp(nf - mean(nf))
            y$samples$norm.factors <- nf
        }
        if (ncol(pd)) { y$samples <- cbind(y$samples, pd) }
        if (!get.spikes) { y <- y[!isSpike(x),] }
        return(y)

    } else if (type=="DESeq2") {
        if (missing(assay)) { assay <- "counts" }
        dds <- DESeq2::DESeqDataSetFromMatrix(assayDataElement(x, assay), pd, ~1, ...)
        S4Vectors::mcols(dds) <- fd
        sizeFactors(dds) <- sf
        if (!get.spikes) { dds <- dds[!isSpike(x),] }
        return(dds)

    } else if (type=="monocle") {
        if (missing(assay)) { assay <- "exprs" }
        cur.exprs <- assayDataElement(x, assay)
        if (logged.exprs) { cur.exprs <- 2^cur.exprs }
        out <- monocle::newCellDataSet(cur.exprs, phenoData=pd, featureData=fd, ...)
        if (!get.spikes) { out <- out[!isSpike(x),] }
        return(out)

    }
})
