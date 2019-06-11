# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:05 2019

@author: Casiano.Koprowski
"""
import os as _os
import numpy as _np

from hyo2 import bag as _bag
from osgeo import gdal as _gdal

_gdal.UseExceptions()

class bag_file:
    """This class serves as the main container for BAG data.

    In order to read and assign the data relevant to the interpolation process,
    this class heavily depends on HydrOffice's hyo2.bag module to open and
    populate the data needed

    Attributes
    ----------
    nodata : float
        1000000.0 is the nodata value associated with the BAG format
    elevation : numpy.array
        The elevation layer of the BAG
    uncertainty : numpy.array
        The uncertainty layer of the BAG
    shape : tuple
        The (y, x) dimensions of the layer data
    bounds : touple
        The geographic bounds of the data in order (nw, se)=([sx,ny],[nx,sy])
    wkt : str
        The WKT representation of the data CRS

    """
    def __init__(self):
        self.name = None
        self.nodata = 1000000.0
        self.elevation = None
        self.uncertainty = None
        self.shape = None
        self.bounds = None
        self.resolution = None
        self.wkt = None
        self.size = None
        self.outfilename = None

    def open_file(self, filepath):
        """
        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """
        _fName = _os.path.split(filepath)[-1]
        self.name = _os.path.splitext(_fName)[0]
        self.nodata = 1000000.0
        _bag_obj = _bag.BAGFile(filepath)
        _bag_obj.populate_metadata()
        self.elevation = self._nan2ndv(_bag_obj.elevation(), self.nodata)
        self.uncertainty = self._nan2ndv(_bag_obj.uncertainty(), self.nodata)
        self.shape = _bag_obj.elevation_shape()
        self.bounds = self._meta2bounds(_bag_obj.meta)
        self.resolution = (_bag_obj.meta.res_x, _bag_obj.meta.res_y)
        self.wkt = _bag_obj.meta.wkt_srs
        self.size = self._size_finder(filepath)

    def _nan2ndv(self, arr, nodata):
        arr[_np.isnan(arr)] = nodata
        arr = _np.flipud(arr)
        return arr

    def _meta2bounds(self,meta):
        sx,sy = meta.sw
        nx,ny = meta.ne
        bounds = ([sx,ny],[nx,sy])
        return bounds

    def _size_finder(self, filepath):
        size = int(_np.round(_os.path.getsize(filepath)/1000))
        return size

    def generate_name(self, outlocation, io):
        if io:
            ext = '_INTERP_ONLY.bag'
        else:
            ext = '_INTERP_FULL.bag'
        name = self.name + ext
        self.outfilename = _os.path.join(outlocation, name)

class gdal_create:
    def __init__(self):
        self.dataset = None

    def bag2gdal(self, bag):
        arrays = [bag.elevation, bag.uncertainty]
        bands = len(arrays)
        x_orig, y_orig = bag.bounds[0]
        y_cols, x_cols = bag.shape
        res_x, res_y = bag.resolution[0], bag.resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols,
                                                       bands, _gdal.GDT_Float32)
        target_gt = (x_orig, res_x, 0, y_orig, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        target_ds.SetProjection(bag.wkt)
        x = 1
        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetNoDataValue(bag.nodata)
            band.WriteArray(_np.flipud(item))
            band = None
            x+=1
        self.dataset = target_ds
        target_ds = None

    def components2gdal(self, arrays, shape, bounds, resolution, wkt,
                 nodata):
        bands = len(arrays)
        x_orig, y_orig = bounds[0]
        y_cols, x_cols = shape
        res_x, res_y = resolution[0], resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols,
                                                       bands, _gdal.GDT_Float32)
        target_gt = (x_orig, res_x, 0, y_orig, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        target_ds.SetProjection(wkt)
        x = 1
        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetNoDataValue(nodata)
            band.WriteArray(_np.flipud(item))
            band = None
            x+=1
        self.dataset = target_ds
        target_ds = None
