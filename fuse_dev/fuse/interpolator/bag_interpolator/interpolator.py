# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:31:40 2019

@author: Casiano.Koprowski
"""
import numpy as _np
import scipy as _scipy
import astropy.convolution as _apc
from datetime import datetime as _dt

def tupleGrid(grid, nodata):
    """Takes an input matrix and an assumed nodata value. The function iterates
    through the matrix and compiles a list of 'edge' points [[x, y, z], ...]
    where:

    1. the current value is not a nodata value and previous value was a nodata value.
        - sets io True, indicating that the next value to compare against should be a nodata value.
    2. the current value is a nodata value and the previous value was not a nodata value.
        - sets io False, indicating that the next value to compare against should not be a nodata value.
    3. returns a list of the points found.

    Parameters
    ----------
    grid : numpy.array
        An input array
    nodata : int
        The array's nodata value

    Returns
    -------
    np.array
        Array of indecies where nodata values meet data values
        in order x, y, z

    """
    points = []
    a = 0
    for x in range(grid.shape[1]):
        io = False
        for y in range(grid.shape[0]):
            if grid[y,x] == nodata:
                if grid[y-1,x] != nodata:
                    val = grid[y-1,x]
                    point = [x, y-1, val]
                    if a == 1:
                        a += 1
                    points.append(point)
                    io = False
                else:
                    pass
            else:
                if io == False:
                    val = grid[y,x]
                    point = [x, y, val]

                    if a == 0:
                        a += 1
                    points.append(point)
                    io = True
    return _np.array(points)

def concatGrid(arr_1, arr_2, nodata):
    """Takes an input of an array of grid objects and the assumed nodata value
    Passes the assumed nodata value and the arrays held within each of the
    listed grid objects to :func:`tupleGrid` for a return of an array of
    edge points for each grid [[x, y, z], ...]

    Takes the results of both tupleGrid calls and combines them into a single
    array of edge points.  Then uses .view().sort() to sort the combined
    products [[x, y, z], ...] based on ascending x, then column y.

    The results of the sorted point array are then split into a list of points
    [[x, y], ...] and values [[z], ...] and returned

    Parameters
    ----------
    grids : list
        The BAG and GeoTiff objects
    maxVal : int
        The BAG data's nodata value

    Returns
    -------
    xy : numpy.array
        Array of indecies where nodata values meet data values
        in order x, y
    z : numpy.array
        Array of values where nodata values meet data values
        in order z

    """
    points_1 = tupleGrid(arr_1, nodata)
    points_2 = tupleGrid(arr_2, nodata)
    if len(points_1) == 0 or len(points_2) == 0:
        xy, z = [], []
    else:
        comb = _np.concatenate([points_1, points_2])
        comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
        grid = _np.hsplit(comb, [2, 4])
        z = grid[1].squeeze()
        xy = grid[0]
    return xy, z

class linear:
    """Interpolates input data and convolves the ouput of the interpolation, if
    applicable.

    Takes input bathy and coverage arrays (tile or complete data) as well as
    the uncertainty array.  This data is used to inform the shape/size of the
    resulting output of the interpolation function. 
    
    The input bathy and covrg arrays are passed to :func:`concatGrid` in order
    to gather the edges of the each of the arrays and combine them into a
    single list of combined xy and z egde values. If xy and z are empty, the 
    interpolation process is skipped and the original tile or complete data is 
    passed back from the function. Otherwise, the lists xy, z, and the shape of
    the BAG data input are passed to :func:`scipy.interpolate.griddata` an
    output to the variable "grid_pre". The output of this function is also 
    saved to a seperate variable "grid" while the original output is left 
    untouched. The seperate output is then passed to 
    :func:`astropy.convolution.convolve` to smooth the output.

    The uncertainty layer is calculated using catzoc and grid::

        >>> m, b = catzoc
        >>> uncrt = (grid*m)+b
    
    Parameters
    ----------
    bathy : numpy.array
        The input bathemetry data
    uncrt : numpy.array
        The input uncertainty data
    covrg : numpy.array
        The input coverage data
    catzoc : tuple
        The input values for uncertainty calculation
    nodata : float, optional
        The default value is 1000000.0, the nodata value associated with the
        BAG format
    
    Returns
    -------
    bathy : numpy.array
        The interpolated and convolved bathemetry
    uncrt : numpy.array
        The calculated uncertainty from the interpolated and convolved
        bathemetry
    unint : numpy.array
        The interpolated bathemetry before convolution is performed
    
    """
    def __init__(self, bathy, uncrt, covrg, catzoc, nodata=1000000.0):
        x, y = _np.arange(bathy.shape[1]), _np.arange(bathy.shape[0])
        xi, yi = _np.meshgrid(x, y)
        xy, z = concatGrid(bathy, covrg, nodata)
        if len(xy) != 0:
            bathy, uncrt, unint = self.interpolate(xy, z, xi, yi, catzoc, nodata)
            return bathy, uncrt, unint
        else:
            return bathy, uncrt, bathy

    def _interpolate(self, xy, z, xi, yi, uval, nodata):
        m, b = uval
        grid_pre = _scipy.interpolate.griddata(xy, z, (xi, yi),
                                            method='linear', fill_value=nodata)
        grid = grid_pre
        grid = _np.asarray(grid, dtype='float64')
        grid[grid>0] = _np.nan
        kernel = _apc.Gaussian2DKernel(3)
        grid = _apc.convolve(grid,kernel)
        grid[_np.isnan(grid)]=nodata
        grid[grid>=0] = nodata
        uncr = (grid*m)+b
        return grid, uncr, grid_pre


