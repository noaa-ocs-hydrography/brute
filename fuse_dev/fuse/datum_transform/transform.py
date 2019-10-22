# -*- coding: utf-8 -*-
"""
transform.py

@author: grice

Abstract datum transformation.
"""

import os
import gdal
from tempfile import TemporaryDirectory as tempdir

from fuse.datum_transform import use_vdatum as uv
from fuse.datum_transform import use_gdal as ug
from fuse.raw_read.noaa.bag import BAGRawReader
gdal.UseExceptions()

class DatumTransformer:
    """ An object for abstracting the datum transformation API.  This should allow for different transformation machines and versions. """

    _from_horiz_datum_info = [
        'from_horiz_frame',
        'from_horiz_type',
        'from_horiz_units',
        'from_horiz_key',
    ]
    _to_horiz_datum_info = [
        'to_horiz_frame',
        'to_horiz_type',
        'to_horiz_units',
        'to_horiz_key',
    ]
    _from_vert_datum_info = [
        'from_vert_key',
        'from_vert_units',
        'from_vert_direction',
    ]
    _to_vert_datum_info = [
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
        gdal.Dataset, dict, bool
            GDAL point cloud, metadata, and boolean value of whether data was reprojected
        """
        not_same_horiz = any(metadata[self._from_horiz_datum_info[index]].lower() != metadata[self._to_horiz_datum_info[index]].lower() for index in
               range(len(self._from_horiz_datum_info)))
        not_same_vert = any(metadata[self._from_vert_datum_info[index]].lower() != metadata[self._to_vert_datum_info[index]].lower() for index in
               range(len(self._from_vert_datum_info)))
        is_bag = type(self._reader) is BAGRawReader
        # VDatum and rasters are giving us trouble, so this is a temp workaround
        if is_bag:
            # we can't deal with a change in BAG vertical datum, so punt
            if not_same_vert:
                # avoid interpolating files in the wrong vertical datum.
                metadata['interpolate'] = 'False'
                transformed = False
                data_obj = None
            # if just the horizontal needs updating, warp it
            elif not_same_horiz:
                data_obj = ug._reproject_via_geotransform(filename, metadata, self._reader)
                transformed = True
            # otherwirse, no datume change needed, hooray!
            else:
                metadata['to_filename'] = filename
                transformed = False
                data_obj = self.create(filename, metadata)
            return data_obj, metadata, transformed
        else:
            if not_same_horiz or not_same_vert:
                return self._engine.translate(filename, metadata), metadata, True
            else:
              return self.create(filename, metadata), metadata, False


    def create(self, filename: str, instructions: dict) -> gdal.Dataset:
        """
        Get a GDAL dataset from the given file.

        Parameters
        ----------
        filename
            filename of data file
        instructions
            dictionary of metadata

        Returns
        -------
            GDAL dataset
        """

        if instructions['read_type'] == 'ehydro':
            points = self._reader.read_bathymetry(filename)

            if 'to_horiz_key' in instructions:
                utm_zone = int(instructions['from_horiz_key'])
            if 'from_vert_key' in instructions:
                vertical_datum = instructions['from_vert_key']

            return ug.__xyz2gdal(points, utm_zone, vertical_datum)
        elif instructions['read_type'] == 'bag':
            return self._reader.read_bathymetry(filename, instructions['to_vert_key'])
        else:
            raise ValueError('Reader type not implemented')


    def translate_support_files(metadata: dict, dest_dir: str):
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

        dest_dir
            The path to the directory where updated files should be stored.

        Returns
        -------
        dict
            The updated metadata dictionary with reference to the transformed
            files.
        """

        return ug.translate_support_files(metadata, dest_dir)
