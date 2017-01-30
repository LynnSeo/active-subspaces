# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 15:18:47 2017

@author: seol
"""
import numpy as np
import os
import sys
import pandas as pd

def print_AS_result (ss):
    """
    ss = output from subspaces.py
    """
    df = pd.DataFrame(ss.eigenvectors)
    
    ss.eigenvalues
    ss.eigenvectors
    ss.W1
    ss.W2
    ss.e_br
    ss.sub_br
    ss.partition
#    f = open('stat_table/eigenvectors.txt','w')
#    f.write(str(df))
    return(df)    
    
def active_subspaces_from_hessian(hessian,build_samples,build_values,
                                  test_values,test_samples,active_subspaces_dim,
                                  plot=False):
    sys.path.append('C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/')
    import active_subspaces as asub
    
    #initiate subspaces
    ss=asub.subspaces.Subspaces()
    #ss.compute based on sstype
    ss.compute(df=hessian, sstype=-1, nboot=2000)
    return ss
    if plot:
        asub.utils.plotters.eigenvalues(ss.eigenvalues, e_br=ss.e_br)
        asub.utils.plotters.subspace_errors(ss.sub_br,out_label='errors')
        
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
    active_build_samples = avmap.forward(build_samples.T)[0].T

    # Map the test samples into the active subspace
    active_test_samples = avmap.forward(test_samples.T)[0].T

    # Plot the outputs versus the value of the active variables
    if plot:
        asub.utils.plotters.sufficient_summary(active_test_samples.T,
                                               test_values)

    # Instantiate a PolynomialApproximation object of degree 2
    pr = asub.utils.response_surfaces.PolynomialApproximation(N=2)
    # Fit the response surface on the active variables
    pr.train(active_build_samples.T, build_values.reshape(build_values.shape[0],1))
    # Evaluate the response surface at test points
    predicted_values = pr.predict(active_test_samples.T)[0]
    print 'Abs. response error',np.linalg.norm(predicted_values.squeeze()-test_values.squeeze())/np.sqrt(predicted_values.shape[0])
    print 'Rel. response error',np.linalg.norm(predicted_values.squeeze()-test_values.squeeze())/np.std(test_values)
    
    
def lynn_study(catchment = 'Hessian-based/Hessian matrix/Gingera/',
               t_year = '70s', size_pert = '1e-06'):
    
    data_dir = catchment+'/'+t_year+'/seed-2025/'+size_pert+'/'
    hessian_filename = size_pert+'_Hessian.csv'
    hessian_filename = os.path.join(data_dir,hessian_filename)

    with open(hessian_filename) as f:
        ncols = len(f.readline().split(','))
        
    hessian = np.loadtxt(hessian_filename,skiprows=1, usecols=range(0,ncols),delimiter =',')
    Nsamp = 1000    
    m = hessian.shape[0]
    hess_list = np.empty(shape=(Nsamp, m, m))
    hess_list[:] = np.hsplit(hessian, Nsamp)
                         
    num_vars = m
    
    build_samples_filename = 'xy.csv'
    build_samples_filename=os.path.join(data_dir,build_samples_filename)
    
#    data=np.loadtxt(build_samples_filename,skiprows=1,
#                          usecols=range(1,num_vars+2),delimiter=',')
    data=np.loadtxt(build_samples_filename,skiprows=1,delimiter=',')
    print(data.shape[1])    
    
    #check size
    assert data.shape[1] == num_vars+1
    build_samples = data[:,:-1].T
    build_values = data[:,-1]

    test_samples_filename = 'xy.csv'
    test_data_dir = catchment+'/'+t_year+'/seed-2026/'+size_pert+'/'
    test_samples_filename=os.path.join(test_data_dir,test_samples_filename)

    data=np.loadtxt(test_samples_filename,skiprows=1,delimiter=',')
    #check size    
    assert data.shape[1] == num_vars+1
    test_samples = data[:,:-1].T
    test_values = data[:,-1]

    ss = active_subspaces_from_hessian(hess_list,build_samples,build_values,
                          test_samples,test_values,2,plot=False)
    np.savetxt(os.path.join(data_dir,'eigenvalues.csv'), ss.eigenvalues, delimiter =',')
    np.savetxt(os.path.join(data_dir,'eigenvectors.csv'), ss.eigenvectors, delimiter =',')

if __name__ == '__main__':
    
    # set seed
    np.random.seed(3)
    #quadratic_study()
    lynn_study(catchment = 'Hessian-based//Gingera/constrained-range-1',
               t_year='70s')
    


    
                                  
