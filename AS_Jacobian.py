import numpy as np
import os
import sys
#import pylab as plt

def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue
        del globals()[var]

def active_subspace_study(gradients,build_samples,build_values,
                          test_samples,test_values,active_subspace_dim,
                          plot=False):
    sys.path.append('C:/UserData/seol/Sensitivity Analyses/IHACRES/AS/')
    import active_subspaces as asub

    # Instantiate the Subspaces object
    ss = asub.subspaces.Subspaces()

    # Use the gradient samples to compute the eigenvectors that define
    # the active subspace
    # gradients is num_dimsxnum_pts
    # ss.compute(gradients.T)
    ss.compute(df=gradients.T)

    # Plot the eigenvectors. The magnitude of each entry in an eigenvector
    # indicates large sensitivity,   e.g. that parameter is important. Two or
    # more large values in a given eigen vector indicate strong correlation.
    #for i in xrange(ss.eigenvectors.shape[1]):
    #    plt.plot(range(ss.eigenvectors.shape[1]),ss.eigenvectors[:,i],'o')
    #    plt.title("%d"%(i+1))
    #    plt.show()

    # Plot the estimated eigenvalues and errors in the estimated active
    # subspace for varying dimension./n",
    if plot:
        asub.utils.plotters.eigenvalues(ss.eigenvalues, e_br=ss.e_br)
        asub.utils.plotters.subspace_errors(ss.sub_br,out_label='errors')

    # return ss
    
    # Define the active subspace to be one-dimensional.
    ss.partition(active_subspace_dim)
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
    pr = asub.utils.response_surfaces.PolynomialApproximation(N=4)
    # Fit the response surface on the active variables
    pr.train(active_build_samples.T, build_values.reshape(build_values.shape[0],1))
    # Evaluate the response surface at test points
    predicted_values = pr.predict(active_test_samples.T)[0]
    print 'Abs. response error',np.linalg.norm(predicted_values.squeeze()-test_values.squeeze())/np.sqrt(predicted_values.shape[0])
    print 'Rel. response error',np.linalg.norm(predicted_values.squeeze()-test_values.squeeze())/np.std(test_values)
    return (ss)

def quadratic_study():
    # total number of random variables
    num_vars = 6
    # dimension of active subspace
    rank = 3
    # Generate a matrix with normally distributed entries
    Amatrix = np.random.normal(0,1,(num_vars,num_vars))
    # Make A symmetric positive definite
    Amatrix = np.dot( Amatrix.T, Amatrix )
    # Construct low rank approximation of A
    eigvals, eigvecs = np.linalg.eigh( Amatrix.copy() )
    # Set smallest eigenvalues to zero. Note eigenvals are in
    # ascending order
    eigvals[:(num_vars-rank)] = 0.
    # reconstruct Amatrix which will now have an active subspace of
    # dimension=rank

    Amatrix = np.dot( eigvecs, np.dot( np.diag(eigvals),eigvecs.T))

    # A quadratic function of 3 variables.
    def quad_fun(x):
        if x.ndim==1:
            x = x.reshape((x.shape[0],1))
            return 0.5*np.dot(x.T,np.dot(Amatrix,x))
        else:
            vals = np.empty((x.shape[1]),float)
            for i in xrange(x.shape[1]):
                xi = x[:,i]
                vals[i] = 0.5*np.dot(xi.T,np.dot(Amatrix,xi))
            return vals


    # The gradient of the quadratic function.
    def quad_dfun(x):
        if x.ndim==1:
            x = x.reshape((x.shape[0],1))
        return np.dot(Amatrix,x)

    num_build_samples=100
    build_samples=np.random.uniform(-1.,1.,(num_vars,num_build_samples))
    build_values=quad_fun(build_samples)
    gradients = quad_dfun(build_samples)
    print build_values.shape

    num_test_samples=900
    test_samples = np.random.uniform(-1.,1.,(num_vars,num_test_samples))
    test_values = quad_fun(test_samples)
    active_subspace_study(gradients,build_samples,build_values,
                          test_samples,test_values,2,plot=True)

def lynn_study(catchment = 'Jacobian-based/Gingera/Constrained range',
               t_year = '70s',no_seed = '2025', size_pert = '1e-06'):
    # load gradient data
               
    data_dir = catchment+'/'+t_year+'/seed-'+no_seed+'/'+size_pert+'/'
    print(os.getcwd()+data_dir)    
    
    gradient_filename = size_pert+'.csv'
    gradient_filename = os.path.join(data_dir,gradient_filename)

    with open(gradient_filename) as f:
        ncols = len(f.readline().split(','))

#    gradients = np.loadtxt(gradient_filename,skiprows=1,
#                              usecols=range(1,ncols),delimiter=',')
    gradients = np.loadtxt(gradient_filename,skiprows=1,delimiter=',')
    gradients = gradients[:1000,:].T
    num_vars = gradients.shape[0]

    build_samples_filename = 'xy.csv'
    build_samples_filename=os.path.join(data_dir,build_samples_filename)
    
#    data=np.loadtxt(build_samples_filename,skiprows=1,
#                          usecols=range(1,num_vars+2),delimiter=',')
    data=np.loadtxt(build_samples_filename,skiprows=1,delimiter=',')
    assert data.shape[1] == num_vars+1
    build_samples = data[:,:-1].T
    build_values = data[:,-1]

    test_samples_filename = 'xy.csv'
    test_data_dir = catchment+'/'+t_year+'/seed-2026/'+size_pert+'/'
    test_samples_filename=os.path.join(test_data_dir,test_samples_filename)

#    data=np.loadtxt(test_samples_filename,skiprows=1,
#                       usecols=range(1,num_vars+2),delimiter=',')
    data=np.loadtxt(test_samples_filename,skiprows=1,delimiter=',')
    assert data.shape[1] == num_vars+1
    test_samples = data[:,:-1].T
    test_values = data[:,-1]

    ss = active_subspace_study(gradients,build_samples,build_values,
                          test_samples,test_values,2,plot=False)
    
    np.savetxt(os.path.join(data_dir,'eigenvalues.csv'), ss.eigenvalues, delimiter =',')
    np.savetxt(os.path.join(data_dir,'eigenvectors.csv'), ss.eigenvectors, delimiter =',')
    # np.savetxt(os.path.join(data_dir,'e_br.csv'), ss.e_br, delimiter =',')
    # np.savetxt(os.path.join(data_dir,'sub_br.csv'), ss.sub_br, delimiter =',')
        
if __name__ == '__main__':
    clear_all()
    # set seed
    np.random.seed(3)
    #quadratic_study()
    lynn_study(catchment = 'Jacobian-Hessian/Gingera/1000-samples',
               t_year = '70s',size_pert = '1e-04')
    
    
