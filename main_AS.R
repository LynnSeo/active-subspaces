source('utils_AS.R')

est_J_H(t_catchment = 'Gingera', sub_catchment = '500 samples', t_year = '70s', no_seed = 2025, perts='1e-04',
        Nsamp =10)

plot_eigen(t_catchment = 'Jacobian-Hessian/Gingera', sub_catchment = '500 samples')

construct_RSM(sub_catchment = '500 samples', t_year = '70s', perts = '1e-04', seperation.n = 3)