# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:05 2019

@author: Casiano.Koprowski
"""
import os as _os

import numpy as _np
from hyo2 import bag as _bag
from osgeo import gdal as _gdal
from osgeo import osr as _osr

_gdal.UseExceptions()

class bag_file:
    """This class serves as the main container for BAG data."""
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
        """Used to read a BAG file using the method determined by input.
        
        The current methods available are 'gdal' and 'hyo'

        :param filepath: 
        :param method: 

        
        """
        if method == 'gdal':
            self._file_gdal(filepath)
        elif method == 'hyo':
            self._file_hyo(filepath)
        else:
            raise ValueError('Open method not implemented.')

    def _file_hyo(self, filepath):
        """Used to read a BAG file using HydrOffice's hyo2.bag module.
        
        This function reads and populates this object's attributes

        :param filepath: 

        
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
        """Used to read a BAG file using OSGEO's GDAL module.
        
        This function reads and populates this object's attributes

        :param filepath: 

        
        """
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
        """

        :param filepath: 

        """
        _fName = _os.path.split(filepath)[-1]
        self.name = _os.path.splitext(_fName)[0]
        self.nodata = 1000000.0
        self.size = self._size_finder(filepath)

    def _nan2ndv(self, arr, nodata):
        """

        :param arr: 
        :param nodata: 

        """
        arr[_np.isnan(arr)] = nodata
        return _np.flipud(arr)

    def _npflip(self, arr):
        """

        :param arr: 

        """
        return _np.flipud(arr)

    def _meta2bounds(self,meta):
        """

        :param meta: 

        """
        sx,sy = meta.sw
        nx,ny = meta.ne
        return ([sx,ny],[nx,sy])

    def _gt2bounds(self ,meta, shape):
        """

        :param meta: 
        :param shape: 

        """
        y, x = shape
        res = (_np.round(meta[1]), _np.round(meta[5]))
        sx, sy = meta[0], meta[3]
        nx = sx + (x * res[0])
        ny = sy + (y * res[1])
        return ([sx,sy],[nx,ny]), res

    def _gdalreadarray(self, bag_obj, band):
        """

        :param bag_obj: 
        :param band: 

        """
        return bag_obj.GetRasterBand(band).ReadAsArray()

    def _size_finder(self, filepath):
        """

        :param filepath: 

        """
        return int(_np.round(_os.path.getsize(filepath)/1000))

    def generate_name(self, outlocation, io):
        """

        :param outlocation: 
        :param io: 

        """
        if io:
            ext = '_INTERP_ONLY.bag'
        else:
            ext = '_INTERP_FULL.bag'
        name = self.name + ext
        self.outfilename = _os.path.join(outlocation, name)

class gdal_create:
    """ """
    _descriptions  = ['Elevation', 'Uncertainty', 'Interpolated']

    def __init__(self, out_verdat=None):
        self.dataset = None
        self.out_verdat = out_verdat

    def bag2gdal(self, bag):
        """

        :param bag: 

        """
        arrays = [bag.elevation, bag.uncertainty]
        bands = len(arrays)
        nw, se = bag.bounds
        nwx, nwy = nw
        scx, scy = se
        y_cols, x_cols = bag.shape
        res_x, res_y = bag.resolution[0], bag.resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols,
                                                       bands, _gdal.GDT_Float32)
        target_gt = (nwx, res_x, 0, scy, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        srs = _osr.SpatialReference(wkt=bag.wkt)
        srs.SetVertCS(self.out_verdat, self.out_verdat, 2000)
        wkt = srs.ExportToWkt()
        target_ds.SetProjection(wkt)
        x = 1
        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetDescription(self._descriptions[x])
            band.SetNoDataValue(bag.nodata)
            band.WriteArray(_np.flipud(item))
            band = None
            x+=1
        self.dataset = target_ds
        target_ds = None

    def components2gdal(self, arrays, shape, bounds, resolution, prj,
                 nodata):
        """

        :param arrays: 
        :param shape: 
        :param bounds: 
        :param resolution: 
        :param prj: 
        :param nodata: 

        """
        bands = len(arrays)
        nw, se = bounds
        nwx, nwy = nw
        scx, scy = se
        y_cols, x_cols = shape
        res_x, res_y = resolution[0], resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols,
                                                       bands, _gdal.GDT_Float32)
        target_gt = (nwx, res_x, 0, scy, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        srs = _osr.SpatialReference(wkt=prj)
        srs.SetVertCS(self.out_verdat, self.out_verdat, 2000)
        wkt = srs.ExportToWkt()
        target_ds.SetProjection(wkt)
        x = 1
        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetDescription(self._descriptions[x-1])
            band.SetNoDataValue(nodata)
            band.WriteArray(_np.flipud(item))
            band = None
            x+=1
        self.dataset = target_ds
        target_ds = None
