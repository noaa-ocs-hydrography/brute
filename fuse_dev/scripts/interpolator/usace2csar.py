# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:30:01 2019

@author: Casiano.Koprowski
"""


import os
import sys
import re as _re
import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.mlab as _mlab
import matplotlib
import scipy as _scipy

print(matplotlib.__version__)

#import fuse.proc_io.proc_io.proc_io as _io
#import fuse.interpolator.point_interpolator.point_interpolator as _pi

#print(_io.__version__)

progLoc = os.getcwd()
_ussft2m = 0.30480060960121924 # US survey feet to meters

#path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\GR_LD_GR1_20180817_CS_15_16_SORT.DAT'
#path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\GR_LD_GR1_20180817_CS_15_16_SORT_A.XYZ'
#path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\MR_54_NO1_20190108_CS_10X10_A.xyz'
path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\MR_54_NO1_20190108_CS_10X10.dat'

def _start_xyz(infilename):
    """
    from cemvn.py, slightly modified

    looks for the first line of the xyz data after the header
    returns the row number of first line of data
    """
    x = 0
    pattern_coordinates = '[0-9]{6}'#at least six digits# should be seven then . plus two digits
    with open(infilename, 'r') as infile:
        for (index1, line) in enumerate (infile):
            print(index1, line)
            if _re.match(pattern_coordinates, line) is not None:
                break
            x += 1
#    print (x)
    return x

def read_bathymetry_dat(infilename):
    """
    from cemvn.py, slightly modified

    Read the bathymetry from the .dat file. The dat file is less precise,
    but had no header and is in a standardized format
    """
    # get the dat file for CEMVN
    stub, ext = os.path.splitext(infilename)
    bathyfilename = stub + '.dat'
    xyz = _np.loadtxt(bathyfilename, usecols=(0,1,2))
    return xyz

def read_bathymetry_xyz(infilename):
    """
    from cemvn.py, slightly modified

    Read the bathymetry from the xyz files, this tells it to not include
    the header when reading the file

    Note: The high resolution multibeam files are available as .xyz on E-Hydro
    """
    first_instance = _start_xyz(infilename)
    if first_instance != 0:
        xyz = _np.loadtxt(infilename, skiprows = first_instance, usecols=(0,1,2))
    else:
        xyz = _np.loadtxt(infilename, usecols=(0,1,2))
    return xyz

def tupleGrid(grid, maxVal):
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
    maxVal : int
        The array's nodata value

    Returns
    -------
    np.array
        Array of indecies where nodata values meet data values
        in order x, y, z

    """
    print ('tupleGrid')
    points = []
    a = 0
    for x in range(grid.shape[1]):
        io = False
        for y in range(grid.shape[0]):
            if grid[y,x] == maxVal:
                if grid[y-1,x] != maxVal:
                    val = grid[y-1,x]
                    point = [x, y-1, val]
                    if a == 1:
                        print (point, val)
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
                        print (point, val)
                        a += 1
                    points.append(point)
                    io = True
    return _np.array(points)

def make_grid(data):
    maxVal = 1000000.0
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    e, n = _np.amax(x), _np.amax(y)
    w, s = _np.amin(x), _np.amin(y)
    shape = (int(_np.round(n-s)), int(_np.round(e-w)))
    extents = (n, w), (s, e)
    mdx = _np.min(abs(_np.diff(x)))
    mdy = _np.min(abs(_np.diff(y)))
#    dx = _np.median(mdx[_np.where(mdx>0.0)[0]]) * _ussft2m
#    dy = _np.median(mdy[_np.where(mdy>0.0)[0]]) * _ussft2m
    dx = 1
    dy = 1
    res = _np.round(dx), _np.round(dy)
#    res = 1, 1
    print (res, (mdx, mdy))
    xi= _np.arange(w,e+1,1)
    yi= _np.arange(s,n+1,1)
    zi= _np.ones((len(yi),len(xi)))*maxVal
    print (_np.shape(zi), zi)
    
    # calculate indices in full grid (zi) to stick the input z values
    ix = _np.round((x-w)/dx).astype(int)
    iy = _np.round((y-s)/dy).astype(int)
    zi[iy,ix] = z
    
    zi = _np.flipud(zi)
    ziShape = zi.shape
    
    _plt.figure()
    _plt.imshow(zi)
    _plt.show()
    
    points = tupleGrid(zi, maxVal)
    
    comb = points
    comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
    grid = _np.hsplit(comb, [2, 4])
    vals = grid[1].squeeze()
    grid = grid[0]
    
    print (grid, vals)
    
    xyz_data = xyz_grid(grid, vals, extents, ziShape, res)
    return xyz_data

def natInterp(xy, z, shape):
    print ('natInterp')
    print (shape)
    maxVal = 1000000.0
    a = xy[:,0]
    b = xy[:,1]
#    shape = xy.shape, z.shape
#    print (shape)
    x, y = _np.arange(shape[1]), _np.arange(shape[0])
    print (a, b, z)
    xi, yi = _np.meshgrid(x, y)
#    interp = _mlab.griddata(a, b, z, xi, yi)
    interp = _scipy.interpolate.griddata(xy, z, (xi, yi),
                                            method='linear', fill_value=_np.nan)
    
    print ('interp done')
    
    _plt.figure()
    _plt.imshow(interp)
    _plt.show()
    
    return interp

class xyz_grid():
    
    def __init__(self, grid, vals, extents, shape, res):
        self.grid = grid
        self.vals = vals
        self.extents = extents
        self.shape = shape
        self.res = res

ext = os.path.splitext(path)[1].lower()
if ext == '.xyz':
    data = read_bathymetry_xyz(path)
elif ext == '.dat':
    data = read_bathymetry_dat(path)
d = make_grid(data)
interp = natInterp(d.grid,d.vals,d.shape)
#print (data)
