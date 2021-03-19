library(msigdb)
library(GSEABase)

#----load collections----
msigdb.v7.2.hs.SYM = msigdb.v7.2.hs.SYM()
msigdb.v7.2.mm.SYM = msigdb.v7.2.mm.SYM()

#----append KEGG----
msigdb.v7.2.hs.EZID = appendKEGG(msigdb.v7.2.hs.EZID)
msigdb.v7.2.mm.EZID = appendKEGG(msigdb.v7.2.mm.EZID)

computeMemMatrix <- function(msigGsc) {
  msigGsc = GeneSetCollection(msigGsc)
  x = geneIds(msigGsc)
  vals = unique(unlist(x))
  
  # #create membership matrix for x - on desktop
  # matx = Matrix::Matrix(
  #   0,
  #   nrow = length(x),
  #   ncol = length(vals),
  #   dimnames = list(names(x), vals),
  #   sparse = TRUE
  # )
  # for (i in 1:length(x)) {
  #   matx[i, ] = as.numeric(vals %in% x[[i]])
  # }
  
  # #create membership matrix for x - on servers
  matx = plyr::laply(x, function(gs) {
    as.numeric(vals %in% gs)
  })
  matx = Matrix::Matrix(
    matx,
    dimnames = list(names(x), vals),
    sparse = TRUE
  )
  
  return(matx)
}

mem_mat_hs = computeMemMatrix(msigdb.v7.2.hs.EZID)
mem_mat_mm = computeMemMatrix(msigdb.v7.2.mm.EZID)

usethis::use_data(mem_mat_hs)
usethis::use_data(mem_mat_mm)
