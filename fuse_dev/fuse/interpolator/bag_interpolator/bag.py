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


    def open_file(self, filepath, method):
        if method == 'gdal':
            self._file_gdal(filepath)
        elif method == 'hyo':
            self._file_hyo(filepath)
        else:
            raise ValueError('Open method not implemented.')

    def _file_hyo(self, filepath):
        """Used to read a BAG file using HydrOffice's hyo2.bag module.

        This function reads and populates this object's attributes

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """
        self._known_data(filepath)
        bag_obj = _bag.BAGFile(filepath)
        bag_obj.populate_metadata()
        self.elevation = self._nan2ndv(bag_obj.elevation(), self.nodata)
        self.uncertainty = self._nan2ndv(bag_obj.uncertainty(), self.nodata)
        self.shape = bag_obj.elevation_shape()
        self.bounds = self._meta2bounds(bag_obj.meta)
        self.resolution = (bag_obj.meta.res_x, bag_obj.meta.res_y)
        self.wkt = bag_obj.meta.wkt_srs

        bag_obj = None

    def _file_gdal(self, filepath):
        self._known_data(filepath)
        bag_obj = _gdal.Open(filepath)
        self.elevation = self._npflip(self._gdalreadarray(bag_obj, 1))
        self.uncertainty = self._npflip(self._gdalreadarray(bag_obj, 2))
        self.shape = self.elevation.shape
        self.bounds, self.resolution = self._gt2bounds(bag_obj.GetGeoTransform(),
                                                       self.shape)
        self.wkt = bag_obj.GetProjectionRef()

        bag_obj = None

    def _known_data(self, filepath):
        _fName = _os.path.split(filepath)[-1]
        self.name = _os.path.splitext(_fName)[0]
        self.nodata = 1000000.0
        self.size = self._size_finder(filepath)

    def _nan2ndv(self, arr, nodata):
        arr[_np.isnan(arr)] = nodata
        return _np.flipud(arr)

    def _npflip(self, arr):
        return _np.flipud(arr)

    def _meta2bounds(self,meta):
        sx,sy = meta.sw
        nx,ny = meta.ne
        return ([sx,ny],[nx,sy])

    def _gt2bounds(self,meta, shape):
        # (605260.0, 4.0, 0.0, 4505852.0, 0.0, -4.0)
        # ([605260.0, 4494888.0], [620528.0, 4483924.0])
        y, x = shape
        res = (_np.round(meta[1]), _np.round(meta[5]))
        sx, sy = meta[0], meta[3]
        nx = sx + (x * res[0])
        ny = sy + (y * res[1])
        return ([sx,sy],[nx,ny]), res

    def _gdalreadarray(self, bag_obj, band):
        return bag_obj.GetRasterBand(band).ReadAsArray()

    def _size_finder(self, filepath):
        return int(_np.round(_os.path.getsize(filepath)/1000))

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
