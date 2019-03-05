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
    def __init__(self, interp_type, resolution):
        """
        Set the interpolation method.
        """
        self._interp_type = interp_type
        self._resolution = resolution
        self._setup()
        
    def _setup(self):
        """
        Set up and configure the interpolation tools.
        """
        if 'vdatum_path' in self._config:
            self._engine = pinterp.point_interpolator()
        else:
            raise ValueError('No java path provided')
        
    def interpolate(self, dataset):
        """
        Take a gdal dataset and run the interpolation, returning a gdal raster.
        """
        return self._engine.interpolate(dataset, self._interp_type, self._resolution)
