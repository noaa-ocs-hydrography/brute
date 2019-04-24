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

#import fuse.proc_io.proc_io.proc_io as _io
#import fuse.interpolator.point_interpolator.point_interpolator as _pi

#print(_io.__version__)

progLoc = os.getcwd()
_ussft2m = 0.30480060960121924 # US survey feet to meters

path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\GR_LD_GR1_20180817_CS_15_16_SORT.DAT'
#path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\GR_LD_GR1_20180817_CS_15_16_SORT_A.XYZ'
#path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\MR_54_NO1_20190108_CS_10X10_A.xyz'
#path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\MR_54_NO1_20190108_CS_10X10.dat'

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

def make_grid(data):
    comb = data
    comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
    grid = _np.hsplit(comb, [2, 4])
    vals = grid[1].squeeze()
    grid = grid[0]
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    e, n = _np.amax(x), _np.amax(y)
    w, s = _np.amin(x), _np.amin(y)
    shape = (int(_np.round(n-s)), int(_np.round(e-w)))
    extents = (n, w), (s, e)
    mdx = abs(_np.diff(x))
    mdy = abs(_np.diff(y))
#    dx = _np.median(mdx[_np.where(mdx>0.0)[0]]) * _ussft2m
#    dy = _np.median(mdy[_np.where(mdy>0.0)[0]]) * _ussft2m
    dx = 1
    dy = 1
    res = _np.round(dx), _np.round(dy)
#    res = 1, 1
    print (res)
    xi= _np.arange(w,e+1,1)
    yi= _np.arange(s,n+1,1)
    zi= _np.ones((len(yi),len(xi)))*_np.nan
    print (_np.shape(zi), zi)
    
    # calculate indices in full grid (zi) to stick the input z values
    ix = _np.round((x-w)/dx).astype(int)
    iy = _np.round((y-s)/dy).astype(int)
    zi[iy,ix] = z
    
    zi = _np.flipud(zi)
    
    _plt.figure()
    _plt.imshow(zi)
    _plt.show()
    
    xyz_data = xyz_grid(grid, vals, extents, shape, res)
    return xyz_data

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
gridObj = make_grid(data)
print (data)

