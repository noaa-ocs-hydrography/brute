# -*- coding: utf-8 -*-
"""
transform.py

@author: grice

Abstract datum transformation.
"""

import os
import gdal
from fuse.datum_transform import use_vdatum as uv
from fuse.raw_read.noaa.bag import BAGRawReader
gdal.UseExceptions()

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

    def translate_support_files(self, metadata: dict, dest_dir: str):
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
        if 'support_files' in metadata:
            sf = metadata['support_files']
            t = []
            srs = self._build_srs(metadata)
            for f in sf:
                root, ext = os.path.splitext(f)
                if ext == '.tiff' or ext == '.tif':
                    src = gdal.Open(f)
                    if src is not None:
                        prj = src.GetProjection()
                        if prj != srs.ExportToWkt():
                            path, base = os.path.split(root)
                            newf = os.path.join(dest_dir, base + ext)
                            options = gdal.WarpOptions(dstSRS=srs, format='GTiff')
                            gdal.Warp(newf, src, options = options)
                            t.append(newf)
                        else:
                            t.append(f)
                    else:
                        t.append(f)
                        print(f'{f} failed to open with gdal')
                else:
                    t.append(f)
                    print('Only the translation of GeoTiffs is currently supported.')
            metadata['support_files'] = t
        return metadata
    
    def _build_srs(self, metadata):
        """
        Build a gdal osr spatial reference object from the provided metadata 
        dictionary to_horiz_* fields.  If the 'to_horiz_frame','to_horiz_type'
        and 'to_horiz_key' fields are not populated a value error is raised.
        
        Parameters
        ----------
        metadata
            A dictionary of the metadata associated with a survey.
        
        Returns
        -------
        osr object for the target horizontal reference system.
        """
        req = ['to_horiz_frame','to_horiz_type','to_horiz_key']
        if all(key in metadata for key in req):
            if metadata['to_horiz_type'] == 'utm':
                proj4_str = f'+proj=utm +zone={metadata["to_horiz_key"]} +datum={metadata["to_horiz_frame"]}'
                srs = gdal.osr.SpatialReference()
                srs.ImportFromProj4(proj4_str)
                return srs
            else:
                raise ValueError('We still need to sort our when we are not working in utm...')
        else:
            raise ValueError('Not all metadata is available to build the proj4 string')
            