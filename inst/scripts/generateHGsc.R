library(msigdb)

#----load collections----
msigdb.hs = getMsigdb()
hgsc = subsetCollection(msigdb.hs, 'h')

usethis::use_data(hgsc)