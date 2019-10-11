# -*- coding: utf-8 -*-
"""
transform.py

@author: grice

Abstract datum transformation.
"""

import gdal
from fuse.datum_transform import use_vdatum as uv
from fuse.raw_read.noaa.bag import BAGRawReader


class DatumTransformer:
    """ An object for abstracting the datum transformation API.  This should allow for different transformation machines and versions. """

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

    def __init__(self, vdatum_path: str, java_path: str, reader):
        """
        Set up and configure the transformation tools based on the information provided in the configruation file.
        """

        self._reader = reader
        self._engine = uv.VDatum(vdatum_path, java_path, self._reader)

    def translate(self, filename: str, metadata: dict) -> (gdal.Dataset, dict, bool):
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
        gdal.Dataset, bool
            GDAL point cloud, metadata, and boolean value of whether data was reprojected
        """

        if any(metadata[self._from_datum_info[index]].lower() != metadata[self._to_datum_info[index]].lower() for index in
               range(len(self._from_datum_info))):
            if type(self._reader) is BAGRawReader:
                metadata['interpolate'] = 'False'
                return 'None', metadata, False
            else:
                return self._engine.translate(filename, metadata), metadata, True
        elif metadata['interpolate'].lower() != 'false':
            return self._engine.create(filename, metadata), metadata, False
        else:
            if type(self._reader) is BAGRawReader:
                return 'None', metadata, False
            else:
                return self._engine.create(filename, metadata), metadata, False

    def translate_support_files(self, metadata: dict):
        """
        Check the horizontal georeferencing for the support files.  If they are
        not in the same datum as the output datum they are translated and
        written back to disk.
        
        Parameters
        ----------
        metadata
            A dictionary of the metadata associated with a survey.  The support
            files referenced in the dictionary will be updated if their
            horizontal datum does not match the output datum.
        
        Returns
        -------
        None
        """
        # get the support file types
        # build transform
        # transform files
        # write back out
        pass