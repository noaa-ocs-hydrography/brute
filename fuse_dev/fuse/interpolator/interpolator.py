# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

import os as _os
import sys as _sys
import logging as _logging

import fuse.interpolator.bag_interpolator as binterp
import fuse.interpolator.point_interpolator as pinterp
from osgeo import gdal


class Interpolator:
    """An abstraction for data interpolation."""

    def __init__(self, interpolation_engine: str, interp_type: str, resolution: float):
        """
        Set the interpolation method.
        """
        self._logger = _logging.getLogger('data_log')
        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            self._logger.addHandler(ch)
        self._interp_engine = interpolation_engine
        self._interp_type = interp_type
        self._resolution = resolution
        self._setup()

    def _setup(self):
        """Set up and configure the interpolation tools."""

        if self._interp_engine == 'point':
            self._engine = pinterp
        elif self._interp_engine == 'raster':
            self._engine = binterp.process.RasterInterpolator()
        else:
            raise ValueError('No interpolation engine type specified')

    def interpolate(self, dataset: gdal.Dataset, metadata: dict) -> (gdal.Dataset, dict):
        """
        Take a gdal dataset and run the interpolation, returning a gdal raster.

        Parameters
        ----------
        dataset
            GDAL dataset
        metadata
            dictionary of metadata

        Returns
        -------
            interpolated GDAL dataset and dictionary of metadata
        """

        if 'support_files' in metadata:
            support_files = metadata['support_files']
        else:
            support_files = []

        if 'file_size' in metadata:
            file_size = metadata['file_size']
        else:
            file_size = None

        root, filename = _os.path.split(metadata['outpath'])
        base, ext = _os.path.splitext(filename)

        # Point Interpolation
        if self._interp_engine == 'point':
            if len(support_files) == 0:
                interpolated_dataset = self._engine.interpolate(dataset, self._interp_type, self._resolution)
                metadata['from_filename'] = self.gettag(base)
                metadata['interpolated'] = True
            else:
                interpolated_dataset = self._engine.interpolate(dataset, self._resolution, support_files)

        # Raster Interpolation
        elif self._interp_engine == 'raster':
            if len(support_files) == 0:
                raise ValueError("No coverage files provided; no interpolation can occur")
            else:
                interpolated_dataset = self._engine.interpolate(dataset, self._interp_type, support_files, file_size)
                metadata['from_filename'] = self.gettag(base)
                metadata['interpolated'] = True

        resolution = interpolated_dataset.GetGeoTransform()[1]
        metadata['to_filename'] = f"{_os.path.join(root, base)}_{int(resolution if resolution >= 1 else resolution * 100)}" + \
                                  f"{'m' if resolution >= 1 else 'cm'}_interp.{metadata['new_ext']}"

        return interpolated_dataset, metadata

    def gettag(self, from_name: str) -> str:
        """
        Return the tag for the interpolated dataset given the original filename.
        """
        return f"{from_name}.interpolated"
