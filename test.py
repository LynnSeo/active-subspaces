# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 16:01:05 2017
conduct eigen decomposition based on Hessian

@author: seol
"""
import sys
import os 
import numpy as np
import active_subspaces as asub

#import pandas as pd
#import csv

def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue
        del globals()[var]
        

def sorted_eigh(C):
    """
    TODO: docs
    """
    e, W = np.linalg.eigh(C)
    e = abs(e)
    ind = np.argsort(e)
    e = e[ind[::-1]]
    W = W[:,ind[::-1]]
    s = np.sign(W[0,:])
    s[s==0] = 1
    W = W*s
    return e.reshape((e.size,1)), W
    
def active_subspace_from_hessian(hessian):
    """
    hessian : list of hessians
    """
    #compute the matrix
    c_h = hessian[0]

    nsamp = len(hessian)
    for i in xrange(1, nsamp):
        c_h += hessian[i]

    c_h /= nsamp

    return c_h, sorted_eigh(c_h)

def eig_hessian(catchment = 'Hessian-based/Gingera/',
               t_year = '70s', no_seed='2025', size_pert = '1e-06'):
#    sys.path.append('C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/')
    
    data_dir = catchment+'/'+t_year+'/seed-'+no_seed+'/'+size_pert+'/'
    hessian_filename = size_pert+'_Hessian.csv'
    hessian_filename = os.path.join(data_dir,hessian_filename)
        
    #model dimension
    hessian = np.loadtxt(hessian_filename,skiprows=1, delimiter =',')
    m = hessian.shape[0]
    Nsamp = hessian.shape[1]/hessian.shape[0]
    
    hess_list = np.empty(shape=(Nsamp, m, m))
    hess_list[:] = np.hsplit(hessian, 1000)
     
    #hess_list[:] = np.array_split(hessian, Nsamp, axis=1)
    #print "Match?"
    #print np.allclose(hess_list[0], hessian[:,0:6])
    
    #eigen decomposition for list of hessians
    c_h, eig_dcomp = active_subspace_from_hessian(hess_list)
    
    print c_h[3,1]

    eig_dcomp = list(eig_dcomp)
    ##print '$$$$$$$$$$$$'
#    print(eig_dcomp[0])
#    print(eig_dcomp[1])
    
    #write eigenvalues and eigenvectors as .csv       
    np.savetxt(os.path.join(data_dir,'eigenvalues.csv'), eig_dcomp[0], delimiter =',')
    np.savetxt(os.path.join(data_dir,'eigenvectors.csv'), eig_dcomp[1], delimiter =',')

    
#nrow, ncol
#print ('hessian size=',hessian.shape[0], hessian.shape[1])
if __name__ == "__main__":
    clear_all()
    eig_hessian(catchment = 'Hessian-based/Gingera/constrained-range-1',t_year= '70s',
                no_seed='2025')
    
    # asub.utils.plotters.eigenvalues(ss.eigenvalues, e_br=ss.e_br)
    # asub.utils.plotters.subspace_errors(ss.sub_br,out_label='errors')
    



    