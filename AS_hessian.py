# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 15:18:47 2017

@author: seol
"""
import numpy as np
import os
import sys
import pandas as pd


def active_subspaces_from_hessian(hessian, build_samples, build_values,
                                  test_samples, test_values, active_subspaces_dim,
                                  plot=False):
    import active_subspaces as asub

    # sys.path.append('C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/')

    # initiate subspaces
    ss = asub.subspaces.Subspaces()
    # ss.compute based on sstype
    ss.compute(df=hessian, sstype=-1, nboot=2000)

    if plot:
        asub.utils.plotters.eigenvalues(ss.eigenvalues, e_br=ss.e_br)
        asub.utils.plotters.subspace_errors(ss.sub_br, out_label='errors')

    # Define the active subspace to be one-dimensional.
    ss.partition(active_subspaces_dim)

    # Instantiate the domain of functions of the active variables
    avdom = asub.domains.BoundedActiveVariableDomain(ss)
    avmap = asub.domains.BoundedActiveVariableMap(avdom)
    print 'done building active subspace map'

    # Instantiate a map between the active variables and the original
    # variables
    avmap = asub.domains.BoundedActiveVariableMap(avdom)

    # Map the build samples into the active subspace
    active_build_samples = avmap.forward(build_samples)[0].T

    # Map the test samples into the active subspace
    active_test_samples = avmap.forward(test_samples)[0].T

    # Plot the outputs versus the value of the active variables
    if plot:
        asub.utils.plotters.sufficient_summary(active_test_samples.T,
                                               test_values)

    # Instantiate a PolynomialApproximation object of degree 2
    pr = asub.utils.response_surfaces.PolynomialApproximation(N=2)
    # Fit the response surface on the active variables
    pr.train(active_build_samples.T, build_values.reshape(build_values.shape[0], 1))
    # Evaluate the response surface at test points
    predicted_values = pr.predict(active_test_samples.T)[0]
    print 'Abs. response error', np.linalg.norm(predicted_values.squeeze() - test_values.squeeze()) / np.sqrt(predicted_values.shape[0])
    print 'Rel. response error', np.linalg.norm(predicted_values.squeeze() - test_values.squeeze()) / np.std(test_values)

    return ss


def eig_hessian(catchment='Hessian-based/Gingera/',
                t_year='70s', no_seed='2025', size_pert='1e-06'):

    data_dir = catchment + '/' + t_year + '/seed-2025/' + size_pert + '/'
    hessian_filename = size_pert + '_Hessian.csv'
    hessian_filename = os.path.join(data_dir, hessian_filename)

    # read hessian.csv
    hessian = np.loadtxt(hessian_filename, skiprows=1, delimiter=',')
    # model dimension : it's opposite to Jacobian-based
    Npar = hessian.shape[0]
    Nsamp = hessian.shape[1] / hessian.shape[0]

    # convert to 3d array
    hess_list = np.empty(shape=(Nsamp, Npar, Npar))
    hess_list[:] = np.hsplit(hessian, Nsamp)

    c_h = hess_list[0]

    for i in xrange(1, Nsamp):
        c_h += hess_list[i]

    c_h /= Nsamp
    # average of Hessians
    hess_list = c_h

    # read coordinates
    build_samples_filename = 'xy.csv'
    build_samples_filename = os.path.join(data_dir, build_samples_filename)
    data = np.loadtxt(build_samples_filename, skiprows=1, delimiter=',')

    # check size
    assert data.shape[1] == Npar + 1
    build_samples = data[:, :-1]
    build_values = data[:, -1]
    # print 'structure of build_sample', len(build_samples.shape)
    # print build_samples.shape

    test_samples_filename = 'xy.csv'
    test_data_dir = catchment + '/' + t_year + '/seed-2026/' + size_pert + '/'
    test_samples_filename = os.path.join(test_data_dir, test_samples_filename)
    data = np.loadtxt(test_samples_filename, skiprows=1, delimiter=',')

    # check size
    assert data.shape[1] == Npar + 1
    test_samples = data[:, :-1]
    test_values = data[:, -1]
    
    ss = active_subspaces_from_hessian(hess_list, build_samples, build_values,
                                       test_samples, test_values, 1, plot=False)

    np.savetxt(os.path.join(data_dir, 'eigenvalues.csv'), ss.eigenvalues, delimiter=',')
    np.savetxt(os.path.join(data_dir, 'eigenvectors.csv'), ss.eigenvectors, delimiter=',')
    print ss.eigenvalues
    print ss.eigenvectors

if __name__ == '__main__':

    # set seed
    np.random.seed(3)
    # quadratic_study()
    eig_hessian(catchment='Hessian-based/Gingera/constrained-range-1',
                t_year='70s')
