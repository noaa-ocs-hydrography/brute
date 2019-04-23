# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:30:01 2019

@author: Casiano.Koprowski
"""


import os
import sys
import re as _re
import numpy as _np

import fuse.proc_io.proc_io.proc_io as _io
import fuse.interpolator.point_interpolator.point_interpolator as _pi

print(_io.__version__)

progLoc = os.getcwd()

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
    print (x)
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

ext = os.path.splitext(path)[1].lower()
if ext == '.xyz':
    data = read_bathymetry_xyz(path)
elif ext == '.dat':
    data = read_bathymetry_dat(path)
print (data)

#    indata = np.loadtxt(path, usecols=(0,1,2))
#    print (indata)