from netCDF4 import Dataset
import interpolation as ip
import numpy as np


def test(xs_name, ys_name, terrain_names):
    tauy = Dataset('tauy_sc.nc', 'r', format='NETCDF3')
    units = "empty"
    xs = tauy.variables[xs_name][:]
    points = np.empty(xs.size, dtype=int)
    units = tauy.variables[xs_name].units
    if units != "empty":
        if units[0:5] == "years":
            for i in xrange(points.size - 1):
                points[i] = (xs[i+1] - xs[i]) * 6.
            points[points.size-1] = points[points.size-2]
        elif units[0:6] == "months":  # Monate mit mehr/weniger als 30 Tage?
            for i in xrange(points.size - 1):
                points[i] = (xs[i+1] - xs[i]) * 15.
            points[points.size-1] = points[points.size-2]
        elif units[0:4] == "days":
            for i in xrange(points.size - 1):
                points[i] = (xs[i+1] - xs[i]) * 12.
            points[points.size-1] = points[points.size-2]
        elif units[0:5] == "hours":
            for i in xrange(points.size - 1):
                points[i] = (xs[i+1] - xs[i]) * 30.
            points[points.size-1] = points[points.size-2]
        elif units[0:7] == "minutes":
            for i in xrange(points.size - 1):
                points[i] = (xs[i+1] - xs[i]) * 30.
            points[points.size-1] = points[points.size-2]
        elif units[0:7] == "seconds":
            for i in xrange(points.size):
                points[i] = 1
    terrains = len(terrain_names)
    sizes = np.empty(terrains + 2, dtype=int)
    for i in xrange(terrains):
        sizes[i] = tauy.variables[terrain_names[i]].size
    sizes[terrains] = 4
    sizes[terrains+1] = xs.size
    ys = tauy.variables[ys_name][:]
    coefficients = ip.interpolate(xs, ys, points, 0.1)
    tauy.close()
    return coefficients

if __name__ == "__main__":
    coeffs = test('time', 'TAUY', ['LONGITUDE', 'LATITUDE'])
