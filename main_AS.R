library('xts')
library('hydromad')
library('ggplot2')
library('stringr')
library(gridExtra)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('utils_AS.R')


est_J_H(t_catchment = 'Gingera', sub_catchment = '10-samples', t_year = '70s', no_seed = 2025, perts='1e-04',
        Nsamp =10, method='latin.hypercube')


# seeds=c(2025,2026)
# years=c('80s','90s','00s')
# for (s in seeds){
#   est_J_H(t_catchment = 'Gingera', sub_catchment = '1000-samples', t_year = '70s', no_seed = s, perts='1e-04',
#           Nsamp =100)
# }

plot_eigen(t_catchment = 'Jacobian-Hessian/Gingera', sub_catchment = '500-samples', additional_dir = 'hess/')

  construct_RSM(t_catchment = 'Jacobian-Hessian/Gingera', sub_catchment = '500-samples', seperation.n = 1, plot = T,
                additional_dir = 'hess/')

