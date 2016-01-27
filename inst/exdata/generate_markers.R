library(scran)

# Need to replace with actual code from stable public sources.
download.file("https://github.com/PMBio/cyclone/raw/master/R/pairs_method/core/pairs_functions.RData", 
              "pairs_functions.RData", quiet=TRUE)
load("pairs_functions.RData")

all.pairs <- sandbag(is.G1=id.G1, is.S=id.S, is.G2M=id.G2M, training.data[genes.training,], fraction=0.5)
saveRDS(file="mouse_cycle_markers.rds", all.pairs)

# HUMAN:

count.file <- "GSE64016.csv.gz"
download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE64016&format=file&file=GSE64016%5FH1andFUCCI%5Fnormalized%5FEC%2Ecsv%2Egz", count.file)
hs.counts <- read.csv(count.file, header=TRUE, row.names=1)
hs.G1 <- grepl("G1", colnames(hs.counts))
hs.S <- grepl("S", colnames(hs.counts))
hs.G2 <- grepl("G2", colnames(hs.counts))

library(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, keytype="SYMBOL", key=rownames(hs.counts), column="ENSEMBL")
m <- match(anno$SYMBOL, rownames(hs.counts))
hs.counts2 <- as.matrix(hs.counts[m,])
rownames(hs.counts2) <- anno$ENSEMBL
hs.counts2 <- hs.counts2[!is.na(anno$ENSEMBL),]
hs.cycle <- select(org.Hs.eg.db, keytype="GOALL", key="GO:0007049", column="ENSEMBL")
hs.training <- rownames(hs.counts2) %in% hs.cycle$ENSEMBL

all.pairs <- sandbag(is.G1=hs.G1, is.S=hs.S, is.G2M=hs.G2, hs.counts2[hs.training,], fraction=0.5)
saveRDS(file="human_cycle_markers.rds", all.pairs)
