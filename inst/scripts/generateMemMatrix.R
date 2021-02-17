library(msigdb)
library(GSEABase)

#----load collections----
# msigdb.hs.SYM = msigdb.hs.SYM()
# msigdb.mm.SYM = msigdb.mm.SYM()

e = new.env()
load('../msigdb/msigdb.hs.EZID.rda', envir = e)
load('../msigdb/msigdb.mm.EZID.rda', envir = e)
msigdb.hs.EZID = e$msigdb.hs.EZID
msigdb.mm.EZID = e$msigdb.mm.EZID

#----append KEGG----
msigdb.hs.EZID = appendKEGG(msigdb.hs.EZID)
msigdb.mm.EZID = appendKEGG(msigdb.mm.EZID)

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

mem_mat_hs = computeMemMatrix(msigdb.hs.EZID)
mem_mat_mm = computeMemMatrix(msigdb.mm.EZID)

usethis::use_data(mem_mat_hs, mem_mat_mm, internal = TRUE)
