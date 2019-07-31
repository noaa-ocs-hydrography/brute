# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

import fuse.interpolator.bag_interpolator as binterp
import fuse.interpolator.point_interpolator as pinterp
from osgeo import gdal


class Interpolator:
    """An abstraction for data interpolation."""

    def __init__(self, interpolation_engine: str, interp_type: str, resolution: float):
        """
        Set the interpolation method.
        """

        self._interp_engine = interpolation_engine
        self._interp_type = interp_type
        self._resolution = resolution
        self._setup()

    def _setup(self):
        """Set up and configure the interpolation tools."""

        if self._interp_engine == 'point':
            self._engine = pinterp.PointInterpolator()
        elif self._interp_engine == 'bag':
            self._engine = binterp.bag_interpolator()
        else:
            raise ValueError('No interpolation engine type specified')

    def interpolate(self, dataset: gdal.Dataset, shapefile: str = None) -> gdal.Dataset:
        """
        Take a gdal dataset and run the interpolation, returning a gdal raster.

        Parameters
        ----------
        dataset :
            param shapefile:  (Default value = None)
        dataset: gdal.Dataset :

        shapefile: str :
             (Default value = None)

        Returns
        -------

        """

        if self._interp_engine == 'point':
            if shapefile is None:
                return self._engine.interpolate(dataset, self._interp_type, self._resolution)
            else:
                return self._engine.interpolate(dataset, self._interp_type, self._resolution, shapefile)
