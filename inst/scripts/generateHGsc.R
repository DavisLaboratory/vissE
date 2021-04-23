library(msigdb)

#----load collections----
msigdb.hs = msigdb.v7.2.hs.SYM()
hgsc = subsetCollection(msigdb.hs, 'h')

usethis::use_data(hgsc)