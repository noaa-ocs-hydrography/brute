# -*- coding: utf-8 -*-
"""
transform.py

@author: grice

Abstract datum transformation.
"""

import gdal
from fuse.datum_transform import use_vdatum as uv


class DatumTransformer:
    """
    An object for abstracting the datum transformation API.  This should allow
    for different transformation machines and versions.
    """

    def __init__(self, config: dict, reader):
        """
        Set up and configure the transformation tools based on the information provided in the configruation file.
        """

        self._reader = reader

        if 'vdatum_path' in config:
            self._engine = uv.VDatum(config, self._reader)
        else:
            raise ValueError('no vdatum path provided')

    def translate(self, filename: str, metadata: dict) -> (gdal.Dataset, bool):
        """
        Run the specified transformation engine to translate the provided
        dataset.

        Parameters
        ----------
        filename
            filename of data file
        metadata
            dictionary of metadata

        Returns
        -------
            GDAL point cloud and whether data was reprojected (`True`) or projected (`False`)
        """

        if metadata['from_horiz_type'].lower() != metadata['to_horiz_type'].lower() or \
                metadata['from_vert_key'].lower() != metadata['to_vert_key'].lower():
            return self._engine.translate(filename, metadata), True
        else:
            return self._engine.create(filename, metadata), False
