library(vissE)
library(msigdb)
library(GSEABase)
library(igraph)
library(RCy3)
library(ExperimentHub)

# msigdb = msigdb.hs.SYM()
load('../msigdb/msigdb.hs.SYM.rda')
# load('../msigdb/msigdb.mm.SYM.rda')
gsc = msigdb.mm.SYM

#add KEGG
gsc = appendKEGG(gsc)

#msignetwork
msigOverlap = computeMsigOverlap(gsc, thresh = 0.15)
msigOverlapNetwork = computeMsigNetwork(msigOverlap, gsc)
save(msigOverlapNetwork, file = 'msigOverlapNetwork.rda')

#layout computation using cytoscape
cytoscapePing()
createNetworkFromIgraph(msigOverlapNetwork)
RCy3::layoutNetwork(layout.name = 'force-directed-cl')

#get layout
lyt = as.matrix(RCy3::getNodePosition())
lyt = t(apply(lyt, 1, as.numeric))
lyt = apply(lyt, 2, function(x) (x - mean(x)) / sd(x))

V(msigOverlapNetwork)$x = lyt[V(msigOverlapNetwork)$name, 1]
V(msigOverlapNetwork)$y = -lyt[V(msigOverlapNetwork)$name, 2]
save(msigOverlapNetwork, file = 'msigOverlapNetwork.rda')
