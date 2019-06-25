# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:05 2019

@author: Casiano.Koprowski
"""

import os as _os
from typing import Tuple, List

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

    def open_file(self, filepath: str, method: str):
        """
        Used to read a BAG file using the method determined by input.

        The current methods available are 'gdal' and 'hyo'

        Parameters
        ----------
        filepath :
            param method:
        filepath: str :

        method: str :


        Returns
        -------

        """

        if method == 'gdal':
            self._file_gdal(filepath)
        elif method == 'hyo':
            self._file_hyo(filepath)
        else:
            raise ValueError('Open method not implemented.')

    def _file_hyo(self, filepath: str):
        """
        Used to read a BAG file using HydrOffice's hyo2.bag module.

        This function reads and populates this object's attributes

        Parameters
        ----------
        filepath :

        filepath: str :


        Returns
        -------

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

        print (self.bounds)
        bag_obj = None

    def _file_gdal(self, filepath: str):
        """
        Used to read a BAG file using OSGEO's GDAL module.

        This function reads and populates this object's attributes

        Parameters
        ----------
        filepath :

        filepath: str :


        Returns
        -------

        """

        self._known_data(filepath)
        bag_obj = _gdal.Open(filepath)
        self.elevation = self._gdalreadarray(bag_obj, 1)
        self.uncertainty = self._gdalreadarray(bag_obj, 2)
        self.shape = self.elevation.shape
        print (bag_obj.GetGeoTransform())
        self.bounds, self.resolution = self._gt2bounds(bag_obj.GetGeoTransform(),
                                                       self.shape)
        self.wkt = bag_obj.GetProjectionRef()

        print (self.bounds)
        bag_obj = None

    def _known_data(self, filepath: str):
        """
        TODO write description

        Parameters
        ----------
        filepath :

        filepath: str :


        Returns
        -------

        """

        _fName = _os.path.split(filepath)[-1]
        self.name = _os.path.splitext(_fName)[0]
        self.nodata = 1000000.0
        self.size = self._size_finder(filepath)

    def _nan2ndv(self, arr: _np.array, nodata: float) -> _np.array:
        """
        TODO write description

        Parameters
        ----------
        arr :
            param nodata:
        arr: _np.array :

        nodata: float :


        Returns
        -------
        type
            array

        """

        arr[_np.isnan(arr)] = nodata
        return _np.flipud(arr)

    def _npflip(self, arr: _np.array) -> _np.array:
        """
        TODO write description

        Parameters
        ----------
        arr :
            returns: array
        arr: _np.array :


        Returns
        -------
        type
            array

        """

        return _np.flipud(arr)

    def _meta2bounds(self, meta) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """
        TODO write description

        Parameters
        ----------
        meta :
            returns: tuple of bounds

        Returns
        -------
        type
            tuple of bounds

        """

        sx, sy = meta.sw
        nx, ny = meta.ne
        return ((sx, ny), (nx, sy))

    def _gt2bounds(self, meta, shape: Tuple[int, int]) -> Tuple[Tuple[Tuple[float, float], Tuple[float, float]], float]:
        """
        TODO write description

        Parameters
        ----------
        meta :
            param shape:
        shape:


        Returns
        -------

        """
        y, x = shape
        res = (meta[1], meta[5])
        sx, sy = _np.round(meta[0]), _np.round(meta[3])
        nx = sx + (x * res[0])
        ny = sy + (y * res[1])
        print ([sx,sy],[nx,ny])
        res = (_np.round(meta[1], 2), _np.round(meta[5], 2))
        return ((sx, ny), (nx, sy)), res

    def _gdalreadarray(self, bag_obj, band: int) -> _np.array:
        """
        TODO write description

        Parameters
        ----------
        bag_obj :
            param band:
        band: int :


        Returns
        -------
        type
            array

        """

        return bag_obj.GetRasterBand(band).ReadAsArray()

    def _size_finder(self, filepath: str) -> int:
        """
        TODO write description

        Parameters
        ----------
        filepath :
            returns: size
        filepath: str :


        Returns
        -------
        type
            size

        """

        return int(_np.round(_os.path.getsize(filepath) / 1000))

    def generate_name(self, outlocation: str, io: bool):
        """
        TODO write description

        Parameters
        ----------
        outlocation :
            param io:
        outlocation: str :

        io: bool :


        Returns
        -------

        """

        name = self.name + '_INTERP_ONLY.bag' if io else '_INTERP_FULL.bag'
        self.outfilename = _os.path.join(outlocation, name)


class gdal_create:
    """TODO write description"""
    _descriptions = ['Elevation', 'Uncertainty', 'Interpolated']

    def __init__(self, out_verdat: str = None):
        """

        Parameters
        ----------
        out_verdat
        """

        self.dataset = None
        self.out_verdat = out_verdat

    def bag2gdal(self, bag):
        """
        TODO write description

        Parameters
        ----------
        bag :


        Returns
        -------

        """

        arrays = [bag.elevation, bag.uncertainty]
        bands = len(arrays)
        sw, ne = bag.bounds
        scx, scy = sw
        nex, ney = ne
        y_cols, x_cols = bag.shape
        res_x, res_y = bag.resolution[0], bag.resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols,
                                                        bands, _gdal.GDT_Float32)
        target_gt = (scx, res_x, 0, scy, 0, res_y)
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
            band.WriteArray(item)
            band = None
            x += 1
        self.dataset = target_ds
        target_ds = None

    def components2gdal(self, arrays: List[_np.array], shape: Tuple[int, int], bounds: Tuple[float, float],
                        resolution: Tuple[float, float], prj: str, nodata: float):
        """
        TODO write description

        Parameters
        ----------
        arrays :
            param shape:
        bounds :
            param resolution:
        prj :
            param nodata:
        arrays: List[_np.array] :

        shape: Tuple[int :

        int] :

        bounds: Tuple[float :

        float] :

        resolution: Tuple[float :

        prj: str :

        nodata: float :


        Returns
        -------

        """

        bands = len(arrays)
        sw, ne = bounds
        scx, scy = sw
        nex, ney = ne
        y_cols, x_cols = shape
        res_x, res_y = resolution[0], resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols,
                                                        bands, _gdal.GDT_Float32)
        target_gt = (scx, res_x, 0, scy, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        srs = _osr.SpatialReference(wkt=prj)
        srs.SetVertCS(self.out_verdat, self.out_verdat, 2000)
        wkt = srs.ExportToWkt()
        target_ds.SetProjection(wkt)
        x = 1
        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetDescription(self._descriptions[x - 1])
            band.SetNoDataValue(nodata)
            band.WriteArray(item)
            band = None
            x += 1
        self.dataset = target_ds
        target_ds = None
