# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:05 2019

@author: Casiano.Koprowski
"""
import os as _os
import numpy as _np
from hyo2 import bag as _bag


class bag_file:
    """This class serves as the main container for BAG data.

    In order to read and assign the data relevant to the interpolation process,
    this class heavily depends on HydrOffice's hyo2.bag module to open and
    populate the data needed

    Parameters
    ----------
    filepath : str
        The complete file path of the input BAG file

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
    def __init__(self, filepath):
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

    def _nan2ndv(self, arr, nodata):
        arr[_np.isnan(arr)] = nodata
        return arr

    def _meta2bounds(self,meta):
        sx,sy = meta.sw
        nx,ny = meta.ne
        bounds = ([sx,ny],[nx,sy])
        return bounds