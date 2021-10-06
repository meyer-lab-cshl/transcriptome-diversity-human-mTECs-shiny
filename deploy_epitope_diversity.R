library(rsconnect)

options(repos = BiocManager::repositories())
getOption("repos")
options(rsconnect.max.bundle.size = 314572800000000)

rsconnect::setAccountInfo(name='chapin',
                          token='1D2291431F22163DC0BAF790A564E998',
                          secret='pgC+KHHrX7KZEnsEuNH1iMtW6s/8eDp1RLaq5D+k')

rsconnect::deployApp("epitope_diversity/")
