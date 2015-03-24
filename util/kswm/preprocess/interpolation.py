import numpy as np
import spline as sp

def average(number, xs, coeffs, units):
    """ Calculate the average over one time-interval of the cubic spline
     Using:
       number: the interval in which the average is to be computed (half between number-1 and number to half between number and number+1)
       units: the number of points over which the average is to be computed
    """
    if number == 0:
        left = xs[number] - (xs[number+1] - xs[number])/2.
    else:
        left = xs[number] - (xs[number] - xs[number-1])/2.
    if number == xs.size-1:
        right = xs[number] + (xs[number] - xs[number-1])/2.
    else:
        right = xs[number] + (xs[number+1] - xs[number])/2.
    diff = right - left
    avg = np.zeros(coeffs.shape[1:-1])
    for i in xrange(units):
        avg += sp.S(left + i * diff/float(units),xs,coeffs)
    return avg/float(units)

def interpolate(xs,ys,points,epsilon):
    """ Interpolate the points given by ys, at the dates given by xs, using the algorithm proposed in
     "A. Harzallah, 'The Interpolation of Data Series Using a Constrained Iterating Technique'"
     Using:
       xs: Time-Values of the points
       ys: Values of the points at the given times
       points: Number of time-points in the interval around xs
       epsilon: Level of precision
     Returns:
       coeffs_sum: Coefficients for the time-intervals
    """
    avgs = np.empty_like(ys)
    coeffs = sp.spline_coeffs(xs,ys)
    coeffs_sum = np.copy(coeffs)
    diffs = np.copy(ys)
    while True:
        for i in xrange(xs.size):
            avgs[i, ...] = average(i,xs,coeffs,points[i])
        diffs = diffs - avgs
        coeffs = sp.spline_coeffs(xs,diffs)
        coeffs_sum += coeffs
        if np.all(np.sum(diffs**2, axis=0) < epsilon**2):
            break
    return coeffs_sum
