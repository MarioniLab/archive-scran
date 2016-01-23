library(scran)

# Need to replace with actual code from stable public sources.
download.file("https://github.com/PMBio/cyclone/blob/master/R/pairs_method/core/pairs_functions.RData", "pairs_functions.RData")
load("pairs_functions.RData")

all.pairs <- sandbag(is.G1=id.G1, is.S=id.S, is.G2M=id.G2M, training.data[genes.training,], fraction=0.5)
saveRDS(file="mouse_cycle_markers.rds", all.pairs)
