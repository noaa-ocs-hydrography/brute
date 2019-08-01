# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:05 2019

@author: Casiano.Koprowski
"""

import os as _os
from typing import Tuple, List

import numpy as _np
import tables as _tb
from hyo2 import bag as _bag
from osgeo import gdal as _gdal
from osgeo import osr as _osr

from . import bag_hack as _bh

_gdal.UseExceptions()


class BagFile:
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
        self.version = None

    def open_file(self, filepath: str, method: str):
        """
        Used to read a BAG file using the method determined by input.

        The current methods available are 'gdal' and 'hyo'

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file
        method : str
            The method used to open the file

        """

        if method == 'gdal':
            self._file_gdal(filepath)
        elif method == 'hyo':
            self._file_hyo(filepath)
        elif method == 'hack':
            self._file_hack(filepath)
        else:
            raise ValueError('Open method not implemented.')

    def _file_hyo(self, filepath: str):
        """
        Used to read a BAG file using HydrOffice's hyo2.bag module.

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

        print(self.bounds)
        del bag_obj

    def _file_gdal(self, filepath: str):
        """
        Used to read a BAG file using OSGEO's GDAL module.

        This function reads and populates this object's attributes

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """

        self._known_data(filepath)
        bag_obj = _gdal.Open(filepath)
        self.elevation = self._gdalreadarray(bag_obj, 1)
        self.uncertainty = self._gdalreadarray(bag_obj, 2)
        self.shape = self.elevation.shape
        print(bag_obj.GetGeoTransform())
        self.bounds, self.resolution = self._gt2bounds(bag_obj.GetGeoTransform(), self.shape)
        self.wkt = bag_obj.GetProjectionRef()
        self.version = bag_obj.GetMetadata()

        print(self.bounds)
        del bag_obj

    def _file_hack(self, filepath: str):
        """
        Used to read a BAG file using pytables and HDF5.

        This function reads and populates this object's attributes

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """

        self._known_data(filepath)

        with _tb.open_file(filepath, mode='r') as bagfile:
            self.elevation = _np.flipud(bagfile.root.BAG_root.elevation.read())
            self.uncertainty = _np.flipud(bagfile.root.BAG_root.uncertainty.read())
            self.shape = self.elevation.shape
            meta_read = [str(x, 'utf-8', 'ignore') for x in bagfile.root.BAG_root.metadata.read()]
            # print (meta_read)
            meta_xml = ''.join(meta_read)
            # print (meta_xml)
            encodeVal = 0

            for x in meta_xml:
                if meta_xml[encodeVal] == '>':
                    meta_xml = meta_xml[encodeVal:]
                    break
                else:
                    encodeVal += 1
            startVal = 0

            for x in meta_xml:
                if meta_xml[startVal] == '<':
                    meta_xml = meta_xml[startVal:]
                    break
                else:
                    startVal += 1

            xml_tree = _bh._et.XML(meta_xml)
            self.wkt = _bh.read_wkt_prj(xml_tree)
            self.resolution = _bh.read_res_x_and_y(xml_tree)
            sw, ne = _bh.read_corners_sw_and_ne(xml_tree)
            sx, sy = sw
            nx = (sx + (self.resolution[0] * self.shape[1]))
            ny = (sy + (self.resolution[0] * self.shape[0]))
            print(ne, (nx, ny))
            self.bounds = ([sx, ny], [nx, sy])

    def _known_data(self, filepath: str):
        """
        Assigns class attributes that are extension agnostic

        :attr:`name`, :attr:`nodata`, and :attr:`size` can be determined
        without using a specific library or method to open file contents, so
        they are assigned by using this function

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """

        _fName = _os.path.split(filepath)[-1]
        self.name = _os.path.splitext(_fName)[0]
        self.nodata = 1000000.0
        self.size = self._size_finder(filepath)

    def _nan2ndv(self, arr: _np.array, nodata: float) -> _np.array:
        """
        Reassigns numpy.nan values with a replacement no data value.

        Parameters
        ----------
        arr : numpy.array
            TODO write description
        nodata : float
            The no data value to be assigned to the numpy.array

        Returns
        -------
        type
            numpy.array

        """

        arr[_np.isnan(arr)] = nodata
        return _np.flipud(arr)

    def _npflip(self, arr: _np.array) -> _np.array:
        """
        Performs :func:`numpy.flipud` on the input numpy.array

        Parameters
        ----------
        arr : numpy.array
            TODO write description

        Returns
        -------
        type
            numpy.array

        """

        return _np.flipud(arr)

    def meta2bounds(self, meta) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        """
        Breaks up and assigns the NW and SE corners from the NE and SW corners
        of a :obj:`hyo2.bag.meta` object

        Parameters
        ----------
        meta : hyo2.bag.meta
            TODO write description

        Returns
        -------
        type
            tuple of bounds

        """

        sx, sy = meta.sw
        nx, ny = meta.ne
        return (sx, ny), (nx, sy)

    def gt2bounds(self, meta, shape: Tuple[int, int]) -> Tuple[Tuple[Tuple[float, float], Tuple[float, float]], float]:
        """
        Formats and returns the bounds and resolution

        This function takes a GeoTransform object and array shape and
        calculates the NW and SE corners.

        Parameters
        ----------
        meta : gdal.GetGeoTransform
            TODO write description
        shape : tuple of int
            (y, x) shape of the bag object

        Returns
        -------
        type
            tuple of bounds, resolution

        """

        y, x = shape
        res = (meta[1], meta[5])
        # ulx, uly = _np.round(meta[0]), _np.round(meta[3])
        ulx, uly = meta[0], meta[3]
        lrx = ulx + (x * res[0])
        lry = uly + (y * res[1])
        print([ulx, uly], [lrx, lry])
        # res = (_np.round(meta[1], 2), _np.round(meta[5], 2))
        return ([ulx, uly], [lrx, lry]), res

    def _gdalreadarray(self, bag_obj, band: int) -> _np.array:
        """
        Returns a numpy.array of a GDAL object

        This function uses a :obj:`gdal.Dataset` object and band number to pull
        and read the appropriate array from the object

        Parameters
        ----------
        bag_obj : gdal.Dataset
            TODO write description
        band : int
            raster band number

        Returns
        -------
        type
            numpy.array

        """

        return bag_obj.GetRasterBand(band).ReadAsArray()

    def _size_finder(self, filepath: str) -> int:
        """
        Returns the rounded size of a file in MB as an integer

        Parameters
        ----------
        filepath : str, os.Pathlike
            TODO write description

        Returns
        -------
        int

        """

        return int(_np.round(_os.path.getsize(filepath) / 1000))

    def generate_name(self, outlocation: str, io: bool):
        """
        Assigns the :attr:`outfilename` of the bag object.

        This function uses the input bool ``io`` to determine the naming
        convention of the output file name. This requires an assigned
        :attr:`name` for this object that is not ``None``

        Parameters
        ----------
        outlocation : str, os.Pathlike
            Folder path of the output file
        io : bool
            Boolean determination of interpolated only or full bag data

        """

        if self.name is not None:
            name = f'{self.name}_INTERP_{"ONLY" if io else "FULL"}.bag'
            self.outfilename = _os.path.join(outlocation, name)


class BagToGDALConverter:
    """This class serves as the main container for converting processed bag
    data into a :obj:`gdal.Dataset` object.

    """
    _descriptions = ['Elevation', 'Uncertainty', 'Interpolated']

    def __init__(self, out_verdat: str = None, flip: bool = False):
        """

        Parameters
        ----------
        out_verdat : str
            Output Vertical Coordinate System (ie. 'MLLW')
        """

        self.dataset = None
        self.out_verdat = out_verdat
        self.flip = flip

    def bag2gdal(self, bag):
        """
        Converts a :obj:`bag` object into a :obj:`gdal.Dataset`

        Parameters
        ----------
        bag : :obj:`bag`
            TODO write description
        """

        arrays = [bag.elevation, bag.uncertainty]
        bands = len(arrays)
        nw, se = bag.bounds
        nwx, nwy = nw
        scx, scy = se
        y_cols, x_cols = bag.shape
        res_x, res_y = bag.resolution[0], bag.resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols, bands, _gdal.GDT_Float32)
        target_gt = (nwx, res_x, 0, nwy, 0, res_y)
        target_gt = self.translate_bag2gdal_extents(target_gt)
        target_ds.SetGeoTransform(target_gt)
        srs = _osr.SpatialReference(wkt=bag.wkt)
        srs.SetVertCS(self.out_verdat, self.out_verdat, 2000)
        wkt = srs.ExportToWkt()
        target_ds.SetProjection(wkt)
        x = 1

        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetDescription(self._descriptions[x - 1])
            band.SetNoDataValue(bag.nodata)

            if self.flip:
                item = _np.flipud(item)

            band.WriteArray(item)
            del band
            x += 1

        self.dataset = target_ds
        del target_ds

    def components2gdal(self, arrays: List[_np.array], shape: Tuple[int, int],
                        bounds: Tuple[Tuple[float, float], Tuple[float, float]], resolution: Tuple[float, float],
                        prj: str, nodata: float = 1000000.0):
        """
        Converts raw dataset components into a :obj:`gdal.Dataset` object

        Parameters
        ----------
        arrays : list of numpy.array
            Arrays [elevation, uncertainty]
        shape : tuple of int
            (y, x) shape of the arrays
        bounds : tuple of tuple of float
            (NW, SE) corners of the data
        resolution : tuple of float
            (x_res, y_res) of the data
        prj : str
            WKT string of the projection of the data
        nodata : float, optional
            no data value of the arrays, default is 1000000.0


        """

        bands = len(arrays)
        nw, se = bounds
        nwx, nwy = nw
        scx, scy = se
        y_cols, x_cols = shape
        res_x, res_y = resolution[0], resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols, bands, _gdal.GDT_Float32)
        target_gt = (nwx, res_x, 0, nwy, 0, res_y)
        target_gt = self.translate_bag2gdal_extents(target_gt)
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

            if self.flip:
                item = _np.flipud(item)

            band.WriteArray(item)
            del band
            x += 1

        self.dataset = target_ds
        del target_ds

    def translate_bag2gdal_extents(self, geotransform: Tuple[float, float, float, float, float, float]):
        orig_x, res_x, skew_x, orig_y, skew_y, res_y = geotransform
        new_x = orig_x - (res_x / 2)
        new_y = orig_y + (res_y / 2)
        return new_x, res_x, skew_x, new_y, skew_y, res_y
