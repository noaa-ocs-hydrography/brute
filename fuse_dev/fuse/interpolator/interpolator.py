# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

import os as _os

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
            GDAL point cloud dataset
        metadata
            dictionary of metadata

        Returns
        -------

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
        metadata['from_filename'] = f"{base}.interpolated"

        # Point Interpolation
        if self._interp_engine == 'point':
            if not support_files:
                interpolated_dataset = self._engine.interpolate(dataset, self._interp_type, self._resolution)
            else:
                interpolated_dataset = self._engine.interpolate(dataset, self._interp_type, self._resolution,
                                                                support_files[0])

            dataset_resolution = interpolated_dataset.GetGeoTransform()[1]
            if dataset_resolution < 1:
                resolution = f'{int(dataset_resolution * 100)}cm'
            elif dataset_resolution >= 1:
                resolution = f'{int(dataset_resolution)}m'

            metadata['to_filename'] = f"{_os.path.join(root, base)}_{resolution}_interp.{metadata['new_ext']}"

        # Raster Interpolation
        elif self._interp_engine == 'raster':
            if not support_files:
                raise ValueError("No coverage files provided; no interpolation can occur")
            else:
                interpolated_dataset = self._engine.interpolate(dataset, self._interp_type, support_files, file_size)

            metadata['to_filename'] = f"{_os.path.join(root, base)}_interp.{metadata['new_ext']}"

        metadata['interpolated'] = True

        return interpolated_dataset, metadata
