# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

import fuse.interpolator.point_interpolator as pinterp

class interpolator:
    """
    An abstraction for data interpolation.
    """
    def __init__(self, interpolation_engine, interp_type, resolution):
        """
        Set the interpolation method.
        """
        self._interp_engine = interpolation_engine
        self._interp_type = interp_type
        self._resolution = resolution
        self._setup()
        
    def _setup(self):
        """
        Set up and configure the interpolation tools.
        """
        if self._interp_engine == 'point':
            self._engine = pinterp.point_interpolator()
        else:
            raise ValueError('No interpolation engine type specified')
        
    def interpolate(self, dataset, shapefile=None):
        """
        Take a gdal dataset and run the interpolation, returning a gdal raster.
        """
        return self._engine.interpolate(dataset, self._interp_type, self._resolution)
