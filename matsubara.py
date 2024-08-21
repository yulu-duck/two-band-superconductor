# -*- mode:python; coding: utf-8; eval: (blacken-mode) -*-
"""
Efficient Matsubara summation
"""
# Author: Pauli Virtanen
# License: GNU Affero General Public License, see LICENSE.txt for details

#
# Common utility code: matsubara
#
from __future__ import division, print_function, absolute_import

import numpy as np
import operator
import functools
from numpy import pi
from scipy import linalg


@functools.lru_cache(4096)
def _gauss_series_weights(n):
    r"""
    Compute Gauss quadrature weights for an infinite series.

    .. math::

       \sum_{j=1}^\infty f(j) \simeq \sum_{k=1}^n w_k f(x_k)

    Uses the method in [1]_.

    Parameters
    ----------
    n : int
        Number of quadrature points to return.

    Returns
    -------

    References
    ----------
    .. [1] H. Monien, "Gaussian quadrature for sums: a rapidly
       convergent summation scheme", Math. Comp. 79, 857 (2010).
       doi:10.1090/S0025-5718-09-02289-3

    """
    n = operator.index(n)
    if n <= 0:
        raise ValueError("n must be positive")
    elif n == 1:
        return np.array([1.0]), np.array([1.0])

    j = np.arange(n)
    j1 = j[1:]

    # Method in [1]: gauss quadrature weights for
    #
    #    sum_{k=1}^oo 1/k^2 f(1/k^2)
    #
    # The abscissas are zeros of `p_n(lmb)` where:
    #
    # lmb = 1/x**2
    # p_n = lambda n, lmb: bessel_poly_yn(2*n+1, -1j*np.sqrt(lmb)/np.pi).real
    # bessel_poly_yn = lambda n, x: np.sqrt(2/(np.pi*x)) * np.exp(1/x) * special.kv(n+0.5, 1/x)

    # Recurrence formula from [1]
    a = 2 * np.pi ** 2 / (4 * j + 1) / (4 * j + 5)
    a[0] = np.pi ** 2 / 15
    sqrt_b = np.pi ** 2 / np.sqrt((4 * j1 - 1) * (4 * j1 + 3)) / (4 * j1 + 1)

    m0 = pi ** 2 / 6  # = sum_{k=1}^oo 1/k^2

    # Compute eigenvalues and first elements of eigenvectors
    #
    # What follows is equivalent to::
    #
    #     z, w = linalg.eigh_tridiagonal(a, sqrt_b)
    #     w = w[0,:]
    #
    # but uses less memory (we only need the first element of each
    # eigenvector).

    m, z, iblock, isplit, info = linalg.lapack.dstebz(a, sqrt_b, 0, 0, 0, 0, 0, 0, "B")
    if info != 0:
        raise ValueError("dstebz failed")

    block_size = max(10, 6250000 // n)

    w = np.empty_like(z)
    iblock2 = np.empty(n, dtype=int)
    for k in range(0, n, block_size):
        zb = z[k : k + block_size]
        iblock2[: len(zb)] = iblock[k : k + block_size]
        wb, info = linalg.lapack.dstein(a, sqrt_b, zb, iblock2, isplit)
        if info != 0:
            raise ValueError("dstein failed at k={}".format(k))
        w[k : k + block_size] = wb[0, :]

    # Then find eigenvectors, one by one

    # Convert back to unweighted sum (z = 1/x^2)
    w = m0 * w ** 2 / z
    x = z ** (-0.5)

    # Swap order
    x = x[::-1]
    w = w[::-1]

    x.flags.writeable = False
    w.flags.writeable = False

    return x, w


def get_matsubara_sum(T, E_typical=1.0, max_ne=2000, steps=1):
    r"""
    Get Matsubara sum/quadrature.

    Gives the approximation:

    .. math::

       T \sum_{n=-\infty}^\infty f(\omega_n) \simeq \sum_{j=0}^{M-1} a_j f(w_j)

    where :math:`f(\omega)` is an analytic function decaying sufficiently
    fast at :math:`\omega\to\pm\infty`.

    The quadrature applies the method of [1]_, which is based on applying
    the idea of Gaussian quadrature to infinite summation of polynomials
    in :math:`z = n^{-2}`.

    Parameters
    ----------
    T : float
        Temperature
    E_typical : float
        Expected typical energy scale of the summand
    max_ne : int
        Max number of frequencies in summation

    Returns
    -------
    w : ndarray
        Energies to evaluate at
    a : ndarray
        Quadrature weights

    Examples
    --------
    >>> import numpy as np
    >>> from nfiscomp.matsubara import get_matsubara_sum
    >>> w, a = get_matsubara_sum(T=1, E_typical=1)
    >>> n = w / (2*np.pi) + 0.5; n
    array([-55.01190294, -17.89542335, -10.6225115 ,  -7.63686822,
            -6.07573899,  -5.00201745,  -4.00000629,  -3.        ,
            -2.        ,  -1.        ,   0.        ,   1.        ,
             2.        ,   3.        ,   4.        ,   5.00000629,
             6.00201745,   7.07573899,   8.63686822,  11.6225115 ,
            18.89542335,  56.01190294])
    >>> (a * 1/w**2).sum(), 1/4
    (0.25000271006449587, 0.25)

    The above Gaussian quadrature based sum converges much faster
    than the naive one:

    >>> w_naive = 2*np.pi*(np.arange(-60, 61) + 0.5)
    >>> (1/w_naive**2).sum()
    0.2491625967192363

    The quadrature is exact for :math:`\sum_{n=1}^\infty p(n^{-2})`
    where `p` is a polynomial up to high order:

    >>> a = a[n>0]; n = n[n>0]
    >>> (a * (1/n**2 + 2/n**4 + 4/n**16)).sum()
    7.809641663308138
    >>> import sympy; m = sympy.symbols('m')
    >>> r = sympy.Sum(1/m**2 + 2/m**4 + 4/m**16, (m, 1, sympy.oo)).doit()
    >>> r
    pi**2/6 + pi**4/45 + 7234*pi**16/162820783125
    >>> r.n()
    7.80964166330814

    References
    ----------
    .. [1] H. Monien, "Gaussian quadrature for sums: a rapidly
       convergent summation scheme", Math. Comp. 79, 857 (2010).
       doi:10.1090/S0025-5718-09-02289-3

    """

    E_max = 100 * abs(E_typical)
    try:
        ne = 5 + int(np.sqrt(4 * E_max / (pi ** 2 * abs(T))))
        k, rem = divmod(ne, steps)
        if rem != 0:
            k += 1
        ne = k * steps
    except (ZeroDivisionError, OverflowError):
        ne = np.inf

    if not (ne < max_ne):
        raise ValueError("Too many matsubara frequencies required")

    # Get weights for sum_{n=1}^oo f(n) =~ sum(weights * f(x))
    x, weights = _gauss_series_weights(ne)

    E = 2 * pi * T * (x - 0.5)
    weights = weights * T

    assert E[-1] >= E_max

    # Get two-tailed
    E = np.r_[-E[::-1], E]
    weights = np.r_[weights[::-1], weights]

    return E, weights
