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
    _from_datum_info = [
        'from_horiz_frame',
        'from_horiz_type',
        'from_horiz_units',
        'from_horiz_key',
        'from_vert_key',
        'from_vert_units',
        'from_vert_direction',
    ]
    _to_datum_info = [
        'to_horiz_frame',
        'to_horiz_type',
        'to_horiz_units',
        'to_horiz_key',
        'to_vert_key',
        'to_vert_units',
        'to_vert_direction',
    ]

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

        if False in [metadata[self._from_datum_info[key]] != metadata[self._to_datum_info[key]] for key in range(len(self._from_datum_info))]:
            if self._reader.version == 'BAG':
                metadata['interpolate'] == 'False'
                return self._engine.create(filename, metadata), metadata, False
            else:
                return self._engine.translate(filename, metadata), metadata, True
        else:
            return self._engine.create(filename, metadata), metadata, False
