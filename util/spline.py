import numpy as np

def derivatives(xs, ys):
    """ Calculate the derivatives used in the Spline-Interpolation"""
    num = xs.size
    u = np.empty_like(xs)
    v = np.empty_like(ys)
    z = np.empty_like(ys)
    h = np.diff(xs)
    b = (np.diff(ys, axis=0).T / h).T
    #Gaussian Elimination
    u[1] = 2. * (h[0] + h[1])
    v[1, ...] = 6. * (b[1, ...] - b[0, ...])
    for i in xrange(2, num - 1):
        u[i] = 2. * (h[i-1] + h[i]) - ((h[i-1] * h[i-1]) / u[i-1])
        v[i, ...] = 6. * (b[i, ...] - b[i-1, ...]) - ((h[i-1] * v[i-1, ...].T) / u[i-1]).T
    #Back-substitution
    z[num-1, ...] = 0.
    for i in xrange(num-2,0,-1):
        z[i, ...] = ((v[i, ...].T - h[i] * z[i+1, ...].T) / u[i]).T
    z[0, ...] = 0.
    return z, h

def spline_coeffs(xs,ys):
    """ Calculate the Spline-Interpolation through the points given by ys at the times given by xs
     Returns the coefficients for the function in the given intervals
    """
    num = xs.size
    #Spline-coefficients
    coeffs_shape = (num-1,) + ys.shape[1:] + (4,)
    coeffs = np.empty(coeffs_shape, dtype=ys.dtype)
    #Coefficients used to compute the Spline-coefficients
    z, h = derivatives(xs,ys)
    coeffs[..., 0] = ys[:-1, ...].copy()
    coeffs[..., 1] = (- (h / 6.) * z[1:, ...].T - (h / 3.) * z[:-1, ...].T + np.diff(ys, axis=0).T / h).T
    coeffs[..., 2] = z[:-1, ...] / 2.
    coeffs[..., 3] = (np.diff(z, axis=0).T / (6. * h)).T
    return coeffs

def S(x,xs,coeffs):
    """ Calculate the interpolated Spline at the time given by x
     Uses the Coefficients the Spline-Interpolation gives and the original time-intervals
    """
    num = xs.size
    i = 0
    while ((x > xs[i+1]) and (i < num - 2)):
        i += 1
    a = coeffs[i, ..., 0]
    b = coeffs[i, ..., 1]
    c = coeffs[i, ..., 2]
    d = coeffs[i, ..., 3]
    diff = x - xs[i]
    ret = a + diff * (b + diff * (c + diff * d))
    return ret
