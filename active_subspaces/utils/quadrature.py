"""
Gaussian quadrature utilities for use with the Python Active-subspaces Utility
Library.
"""

import numpy as np
import misc as mi
import logging

def r_hermite(N):
    """Recurrence coefficients for the Hermite orthogonal polynomials.

    :param int N: The number of recurrence coefficients.

    :return: ab, An `N`-by-2 array of the recurrence coefficients.
    :rtype: ndarray

    **See Also**

    utils.quadrature.jacobi_matrix
    utils.quadrature.gauss_hermite

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    if not isinstance(N, int):
        raise TypeError('N must be an int')

    if N <= 0:
        raise ValueError('Parameters out of range.')

    if N == 1:
        return np.array([[0.0, 1.0]])
    else:
        n = np.array(range(1, N+1))
        b = np.vstack((1.0, n.reshape((N, 1))))
        a = np.zeros(b.shape)
        ab = np.hstack((a, b))
        return ab
        
def r_jacobi(N,l,r,a,b):
    """Recurrence coefficients for the Legendre orthogonal polynomials.

    :param int N: The number of recurrence coefficients.
    :param float l: The left endpoint of the interval.
    :param float r: The right endpoint of the interval.
    :param float a: Jacobi weight parameter
    :param float b: Jacobi weight parameter

    :return: ab, An `N`-by-2 array of the recurrence coefficients.
    :rtype: ndarray

    **See Also**

    utils.quadrature.jacobi_matrix
    utils.quadrature.gauss_legendre

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    if not isinstance(N, int):
        raise TypeError('N must be an int')

    if N <= 0:
        raise ValueError('Parameters out of range.')
    
    a0 = (b-a)/(a+b+2.0)
    ab = np.zeros((N+1,2))
    b2a2 = b**2 - a**2
    s, o = (r-l)/2.0, l + (r-l)/2.0
    
    # first row
    ab[0,0] = s*a0 + o
    ab[0,1] = 1
    
    for k in range(1,N+1):
        ab[k,0] = s*b2a2/((2*(k)+a+b)*(2*(k+1) + a+b)) + o
        if k==1:
            ab[k,1] = ((r-l)**2*(k)*(k+a)*(k+b)) / ((2.0*(k)+a+b)**2*(2.0*(k)+a+b+1))
        else:
            ab[k,1] = ((r-l)**2*(k)*(k+a)*(k+b)*(k+a+b)) / ((2.0*(k)+a+b)**2*(2.0*(k)+a+b+1)*(2.0*(k)+a+b-1))
    return ab

def jacobi_matrix(ab):
    """
    Tri-diagonal Jacobi matrix of recurrence coefficients.

    :param ndarray ab: N-by-2 array of recurrence coefficients

    :return: J, (N-1)-by-(N-1) symmetric, tridiagonal Jacobi matrix associated
        with the orthogonal polynomials.
    :rtype: ndarray

    **See Also**

    utils.quadrature.r_hermite
    utils.quadrature.gauss_hermite

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    if len(ab.shape) != 2:
        raise ValueError('ab must be 2 dimensional')

    if ab.shape[1] != 2:
        raise ValueError('ab must have two columns')

    n = ab.shape[0] - 1
    if n == 0:
        return ab[0,0]
    else:
        J = np.zeros((n, n))
        J[0,0] = ab[0,0]
        J[0,1] = np.sqrt(ab[1,1])
        for i in range(1,n-1):
            J[i,i] = ab[i,0]
            J[i,i-1] = np.sqrt(ab[i,1])
            J[i,i+1] = np.sqrt(ab[i+1,1])
        J[n-1,n-1] = ab[n-1,0]
        J[n-1,n-2] = np.sqrt(ab[n-1,1])
        return J

def gl1d(N):
    """
    One-dimensional Gauss-Legendre quadrature rule.

    :param int N: Number of nodes in the quadrature rule

    :return: x, N-by-1 array of quadrature nodes
    :rtype: ndarray

    :return: w, N-by-1 array of quadrature weights
    :rtype: ndarray

    **See Also**

    utils.quadrature.gauss_legendre

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    return g1d(N, 'Legendre')
    
def gh1d(N):
    """
    One-dimensional Gauss-Hermite quadrature rule.

    :param int N: Number of nodes in the quadrature rule

    :return: x, N-by-1 array of quadrature nodes
    :rtype: ndarray

    :return: w, N-by-1 array of quadrature weights
    :rtype: ndarray

    **See Also**

    utils.quadrature.gauss_hermite

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    return g1d(N, 'Hermite')

def g1d(N, quadtype):
    """
    One-dimensional Gaussian quadrature rule.

    :param int N: Number of nodes in the quadrature rule
    :param str quadtype: Type of quadrature rule {'Legendre', 'Hermite'}

    :return: x, N-by-1 array of quadrature nodes
    :rtype: ndarray

    :return: w, N-by-1 array of quadrature weights
    :rtype: ndarray

    **See Also**

    utils.quadrature.gauss_hermite

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    if N > 1:
        if quadtype == 'Hermite':
            ab = r_hermite(N)
        elif quadtype == 'Legendre':
            ab = r_jacobi(N, -1, 1, 0, 0)
        else:
            raise ValueError('quadtype must be Legendre or Hermite')
        
        J = jacobi_matrix(ab)
        e, V = np.linalg.eig(J)
        ind = np.argsort(e)
        x = e[ind].reshape((N, 1))
        x[np.fabs(x) < 1e-12] = 0.0
        w = (V[0,ind]*V[0,ind]).reshape((N, 1))
    else:
        x, w = np.array([[0.0]]),np.array([[1.0]])
    return x, w

def gauss_hermite(N):
    """
    Tensor product Gauss-Hermite quadrature rule.

    :param int[] N: Number of nodes in each dimension of the quadrature rule

    :return: x, N-by-1 array of quadrature nodes
    :rtype: ndarray

    :return: w, N-by-1 array of quadrature weights
    :rtype: ndarray

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    if isinstance(N, int):
        N = [N]

    if type(N) is not list:
        raise TypeError('N must be a list.')

    Npts = int(np.prod(np.array(N)))
    logging.getLogger(__name__).debug('Making a tensor product Gauss-Hermite rule with {:d} points in {:d} dimensions.'.format(Npts, len(N)))

    if len(N) == 1:
        x, w = gh1d(N[0])
    else:
        x = np.array([[1.0]])
        w = np.array([[1.0]])
        for n in N:
            xi, wi = gh1d(n)

            xL = np.kron(x.copy(), np.ones(xi.shape))
            xU = np.kron(np.ones((x.shape[0],1)), xi)
            x = np.hstack((xL, xU))
            w = np.kron(w.copy(), wi)
        x, w = np.atleast_2d(x[:,1:]), mi.atleast_2d_col(w)

    return x, w

def gauss_legendre(N):
    """
    Tensor product Gauss-Legendre quadrature rule.

    :param int[] N: Number of nodes in each dimension of the quadrature rule

    :return: x, N-by-1 array of quadrature nodes
    :rtype: ndarray

    :return: w, N-by-1 array of quadrature weights
    :rtype: ndarray

    **Notes**

    This computation is inspired by Walter Gautschi's code at
    https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html.
    """

    if isinstance(N, int):
        N = [N]

    if type(N) is not list:
        raise TypeError('N must be a list.')

    Npts = int(np.prod(np.array(N)))
    logging.getLogger(__name__).debug('Making a tensor product Gauss-Legendre rule with {:d} points in {:d} dimensions.'.format(Npts, len(N)))

    if len(N) == 1:
        x, w = gl1d(N[0])
    else:
        x = np.array([[1.0]])
        w = np.array([[1.0]])
        for n in N:
            xi, wi = gl1d(n)

            xL = np.kron(x.copy(), np.ones(xi.shape))
            xU = np.kron(np.ones((x.shape[0],1)), xi)
            x = np.hstack((xL, xU))
            w = np.kron(w.copy(), wi)
        x, w = np.atleast_2d(x[:,1:]), mi.atleast_2d_col(w)

    return x, w