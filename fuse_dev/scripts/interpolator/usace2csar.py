# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:30:01 2019

@author: Casiano.Koprowski
"""

import datetime as _dt
import os
import re as _re
from typing import Union

import astropy.convolution as _apc
import matplotlib.pyplot as _plt
import numpy as _np
import scipy as _scipy
from osgeo import gdal as _gdal
from osgeo import ogr as _ogr
from osgeo import osr as _osr

# import fuse.proc_io.proc_io.proc_io as _io
# import fuse.interpolator.point_interpolator.point_interpolator as _pi

# print(_io.__version__)

progLoc = os.getcwd()
_ussft2m = 0.30480060960121924  # US survey feet to meters

# path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\GR_LD_GR1_20180817_CS_15_16_SORT.DAT'
# path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\GR_LD_GR1_20180817_CS_15_16_SORT_A.XYZ'
path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\MR_54_NO1_20190108_CS_10X10_A.xyz'


# path = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\scripts\\interpolator\\MR_54_NO1_20190108_CS_10X10.dat'

def write_raster(raster_array: _np.array, gt, data_obj, outputpath: str, dtype=_gdal.GDT_UInt32,
                 options: int = 0, color_table: int = 0, nbands: int = 1, nodata: Union[int, bool] = False):
    """
    Directly From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
    (https://gis.stackexchange.com/a/278965)]

    Parameters
    ----------
    raster_array :
        param gt:
    data_obj :
        param outputpath:
    dtype :
        Default value = _gdal.GDT_UInt32)
    options :
        Default value = 0)
    color_table :
        Default value = 0)
    nbands :
        Default value = 1)
    nodata :
        Default value = False)
    raster_array: _np.array :
        
    gt :
        
    outputpath: str :
        
    options: int :
         (Default value = 0)
    color_table: int :
         (Default value = 0)
    nbands: int :
         (Default value = 1)
    nodata: Union[int :
        
    bool] :
         (Default value = False)

    Returns
    -------

    """

    print('write_raster')

    height, width = raster_array.shape

    # Prepare destination file
    driver = _gdal.GetDriverByName("GTiff")
    if options != 0:
        dest = driver.Create(outputpath, width, height, nbands, dtype, options)
    else:
        dest = driver.Create(outputpath, width, height, nbands, dtype)

    # Write output raster
    if color_table != 0:
        dest.GetRasterBand(1).SetColorTable(color_table)

    dest.GetRasterBand(1).WriteArray(raster_array)

    if nodata is not False:
        dest.GetRasterBand(1).SetNoDataValue(nodata)

    # Set transform and projection
    dest.SetGeoTransform(gt)
    wkt = data_obj.GetProjection()
    srs = _osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dest.SetProjection(srs.ExportToWkt())

    # Close output raster dataset
    dest = None


def _dat2gdal(self, infilename: str, dest_epsg: int):
    """
    Read a dat file and turn it into a GDAL point cloud.

    Parameters
    ----------
    infilename :
        param dest_epsg:
    infilename: str :
        
    dest_epsg: int :
        

    Returns
    -------

    """

    points = _np.loadtxt(infilename, delimiter=' ')
    # datum and unit conversion function declaration
    dest = _osr.SpatialReference()
    dest.ImportFromEPSG(dest_epsg)
    # turn numpy points into ogr points in a gdal dataset
    dataset = _gdal.GetDriverByName('Memory').Create('', 0, 0, 0, _gdal.GDT_Unknown)
    layer = dataset.CreateLayer('pts', dest, geom_type=_ogr.wkbPoint)
    for p in points:
        newp = _ogr.Geometry(_ogr.wkbPoint)
        newp.AddPoint(p[0], p[1], p[2])
        feature = _ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(newp)
        layer.CreateFeature(feature)
    return dataset


def _zone2epsg(self, zone: int):
    """
    Assume the EPSG code for a UTM zone is 3707 + the zone number.  This should
    work for zones 1 through 19 for NSRS2007.
    http://spatialreference.org/ref/?search=nad83+utm+zone

    Parameters
    ----------
    zone :
        
    zone: int :
        

    Returns
    -------

    """

    return int(zone) + 3707


def _start_xyz(infilename: str):
    """
    from cemvn.py, slightly modified
    
    looks for the first line of the xyz data after the header
    returns the row number of first line of data

    Parameters
    ----------
    infilename :
        
    infilename: str :
        

    Returns
    -------

    """

    x = 0
    pattern_coordinates = '[0-9]{6}'  # at least six digits# should be seven then . plus two digits
    with open(infilename, 'r') as infile:
        for (index1, line) in enumerate(infile):
            print(index1, line)
            if _re.match(pattern_coordinates, line) is not None:
                break
            x += 1
    #    print (x)
    return x


def read_bathymetry_dat(infilename: str):
    """
    from cemvn.py, slightly modified
    
    Read the bathymetry from the .dat file. The dat file is less precise,
    but had no header and is in a standardized format

    Parameters
    ----------
    infilename :
        
    infilename: str :
        

    Returns
    -------

    """

    # get the dat file for CEMVN
    stub, ext = os.path.splitext(infilename)
    bathyfilename = f'{stub}.dat'
    xyz = _np.loadtxt(bathyfilename, usecols=(0, 1, 2))
    return xyz


def read_bathymetry_xyz(infilename: str):
    """
    from cemvn.py, slightly modified
    
    Read the bathymetry from the xyz files, this tells it to not include
    the header when reading the file
    
    Note: The high resolution multibeam files are available as .xyz on E-Hydro

    Parameters
    ----------
    infilename :
        
    infilename: str :
        

    Returns
    -------

    """

    first_instance = _start_xyz(infilename)
    if first_instance != 0:
        xyz = _np.loadtxt(infilename, skiprows=first_instance, usecols=(0, 1, 2))
    else:
        xyz = _np.loadtxt(infilename, usecols=(0, 1, 2))
    return xyz


def tupleGrid(grid: _np.array, maxVal: int):
    """
    Takes an input matrix and an assumed nodata value. The function iterates
    through the matrix and compiles a list of 'edge' points [[x, y, z], ...]
    where:
    
    1. the current value is not a nodata value and previous value was a nodata value.
        - sets io True, indicating that the next value to compare against should be a nodata value.
    2. the current value is a nodata value and the previous value was not a nodata value.
        - sets io False, indicating that the next value to compare against should not be a nodata value.
    3. returns a list of the points found.

    Parameters
    ----------
    grid :
        An input array
    maxVal :
        The array's nodata value
    grid: _np.array :
        
    maxVal: int :
        

    Returns
    -------

    """

    print('tupleGrid')
    points = []
    a = 0
    for x in range(grid.shape[1]):
        io = False
        for y in range(grid.shape[0]):
            if grid[y, x] == maxVal:
                if grid[y - 1, x] != maxVal:
                    val = grid[y - 1, x]
                    point = [x, y - 1, val]
                    if a == 1:
                        print(point, val)
                        a += 1
                    points.append(point)
                    io = False
                else:
                    pass
            elif not io:
                val = grid[y, x]
                point = [x, y, val]

                if a == 0:
                    print(point, val)
                    a += 1
                points.append(point)
                io = True
    return _np.array(points)


def make_grid(data):
    """
    TODO write description

    Parameters
    ----------
    data :
        

    Returns
    -------

    """

    maxVal = 1000000.0
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2] * _ussft2m
    e, n = _np.amax(x), _np.amax(y)
    w, s = _np.amin(x), _np.amin(y)
    shape = (int(_np.round(n - s)), int(_np.round(e - w)))
    extents = (n, w), (s, e)
    mdx = abs(_np.diff(x))
    mdy = abs(_np.diff(y))
    diffx = _np.median(mdx[_np.where(mdx > 0.0)[0]])
    diffy = _np.median(mdy[_np.where(mdy > 0.0)[0]])
    diff = _np.round(_np.amax([diffx, diffy]))
    print((diffx, diffy), diff)
    dx = 5
    dy = 5
    res = _np.round(dx), _np.round(dy)
    #    res = 1, 1
    #    print (res)
    #    xi= _np.arange(w,e+1,1)
    #    yi= _np.arange(s,n+1,1)
    #    zi= _np.ones((len(yi),len(xi)))*maxVal
    #    print (_np.shape(zi), zi)

    # calculate indices in full grid (zi) to stick the input z values
    #    ix = _np.round((x-w)/dx).astype(int)
    #    iy = _np.round((y-s)/dy).astype(int)
    #    zi[iy,ix] = z
    #
    #    zi = _np.flipud(zi)
    zi = ''
    #    ziShape = zi.shape
    ziShape = int(_np.round((n - s) / 5)), int(_np.round((e - w) / 5))
    #    _plt.figure()
    #    _plt.imshow(zi)
    #    _plt.show()
    #
    #    points = tupleGrid(zi, maxVal)

    comb = data
    comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
    grid = _np.hsplit(comb, [2, 4])
    vals = grid[1].squeeze()
    grid = grid[0]

    xyz_data = xyz_grid(data, vals, extents, ziShape, res, diff)
    return zi, xyz_data


def natInterp(grid: _np.array, xy: _np.array, z: _np.array, shape: tuple, diff: int):
    """
    Applies natural neighbor interpolation to data and then trims it based
    on a mask genereated via the convolution of a binary grid derived from
    original data

    Parameters
    ----------
    grid :
        The complete grid of original data
    xy :
        x and y points for data within the grid
    z :
        z values of the provided data points
    shape :
        Shape of the grid
    diff :
        Median of the difference of distance between points in the grid
    grid: _np.array :
        
    xy: _np.array :
        
    z: _np.array :
        
    shape: tuple :
        
    diff: int :
        

    Returns
    -------

    """

    print('natInterp')
    print(shape, diff)
    maxVal = 1000000.0
    #    a = xy[:,0]
    #    b = xy[:,1]
    #    shape = xy.shape, z.shape
    #    print (shape)
    x, y = _np.arange(shape[1]), _np.arange(shape[0])
    xi, yi = _np.meshgrid(x, y)
    print('try tri', _dt.datetime.now())
    #    interp = _mlab.griddata(a, b, z, xi, yi)
    interp = _scipy.interpolate.griddata(xy, z, (xi, yi),
                                         method='linear', fill_value=_np.nan)
    _plt.figure()
    _plt.imshow(interp)
    _plt.show()
    print('done')
    print('interp', _dt.datetime.now())
    #    kernel = _apc.Gaussian2DKernel(3)
    kernel = _apc.Tophat2DKernel(diff)
    mask_bin = (interp < maxVal).astype(_np.int)
    mask_con = _apc.convolve(mask_bin, kernel)
    mask = (mask_con > 0).astype(_np.int)
    mask_int = _np.where(mask, interp, _np.nan)
    _plt.figure()
    _plt.imshow(mask_int)
    _plt.show()
    print('interp done', _dt.datetime.now())

    return mask_int


class xyz_grid():
    """TODO write description"""

    def __init__(self, grid, vals, extents, shape, res, diff):
        self.grid = grid
        self.vals = vals
        self.extents = extents
        self.shape = shape
        self.res = res
        self.diff = diff


if __name__ == '__main__':
    ext = os.path.splitext(path)[1].lower()
    if ext == '.xyz':
        data = read_bathymetry_xyz(path)
    elif ext == '.dat':
        data = read_bathymetry_dat(path)
    print(ext, data)
    # grid, d = make_grid(data)
    # interp = natInterp(grid, d.grid, d.vals, d.shape, d.diff)
