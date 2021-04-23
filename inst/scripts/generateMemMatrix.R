library(msigdb)
library(GSEABase)

#----load collections----
msigdb.hs = msigdb.v7.2.hs.EZID()
msigdb.mm = msigdb.v7.2.mm.EZID()

#----append KEGG----
msigdb.hs = appendKEGG(msigdb.hs)
msigdb.mm = appendKEGG(msigdb.mm)

computeMemMatrix <- function(msigGsc) {
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

selcolc = c('h', 'c5')
mem_mat_hs = computeMemMatrix(subsetCollection(msigdb.hs, selcolc))
mem_mat_mm = computeMemMatrix(subsetCollection(msigdb.mm, selcolc))

usethis::use_data(mem_mat_hs)
usethis::use_data(mem_mat_mm)
