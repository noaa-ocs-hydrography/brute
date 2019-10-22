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
                outdir = tempdir()
                outtif = os.path.join(outdir, 'tmp.tif')
                outbag = os.path.join(outdir, 'tmp.bag')
                dest_srs = self._build_srs(metadata)
                source = self._engine.create(filename, metadata)
                tif_obj = gdal.Warp(outtif, source, dstSRS=dest_srs)
                data_obj = gdal.Translate(outbag, tif_obj)
                transformed = True
            # otherwirse, no datume change needed, hooray! 
            else:
                metadata['to_filename'] = filename
                transformed = False
                data_obj = self._engine.create(filename, metadata)
            return data_obj, metadata, transformed
        else:
            if not_same_horiz or not_same_vert:
                return self._engine.translate(filename, metadata), metadata, True
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
            dest_srs = self._build_gdal_horiz_srs(metadata)
            for f in sf:
                newf, ext = self._dest_filename(f,dest_dir)
                if ext == '.tiff' or ext == '.tif':
                    src = gdal.Open(f)
                    if src is not None:
                        source_prj = gdal.osr.SpatialReference(wkt=src.GetProjectionRef())
                        if dest_srs.IsSame(source_prj):
                            options = gdal.WarpOptions(dstSRS=dest_srs, format='GTiff')
                            gdal.Warp(newf, src, options=options)
                            t.append(newf)
                        else:
                            t.append(f)
                    else:
                        t.append(f)
                        print(f'{f} failed to open with gdal')
                elif ext == '.gpkg':
                    resulting_file = self._reproject_geopackage(f, newf, dest_srs)
                    t.append(resulting_file)
                else:
                    t.append(f)
                    print(f'Unsupported support file format: {ext}')
            metadata['support_files'] = t
        return metadata

    def _build_gdal_horiz_srs(self, metadata):
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

    def _reproject_geopackage(self, fromfilename: str, tofilename: str, dest_srs: str):
        """
        Convert a Geopackage to a provided reference frame.  A single layer is
        assumed.

        Parameters
        ----------
        filename: str :
            The complete file path of the input coverage file
        to_crs: str :
            WKT object with destination spatial reference system

        Returns
        -------
        str
            The file name for the resulting file.
        """

        fName = os.path.split(fromfilename)[-1]
        splits = os.path.splitext(fName)
        name = splits[0]

        # Open the data source and read in the extent
        source_ds = gdal.ogr.Open(fromfilename)
        source_layer = source_ds.GetLayer()
        source_srs = source_layer.GetSpatialRef()
        # check to see if the file is already projected as it should be
        if dest_srs.IsSame(source_srs):
            outname = fromfilename
        else:
            coordTrans = gdal.osr.CoordinateTransformation(source_srs, dest_srs)
            # build new file
            driver = gdal.ogr.GetDriverByName('gpkg')
            dest_ds = driver.CreateDataSource(tofilename)
            outname = tofilename
            layer = dest_ds.CreateLayer(name, dest_srs, gdal.ogr.wkbMultiPolygon)

            for feature in source_layer:
                if feature is not None:
                    # get the data and transform
                    geom = feature.GetGeometryRef()
                    ds_geom = gdal.ogr.CreateGeometryFromWkt(geom.ExportToWkt())
                    ds_geom.Transform(coordTrans)
                    # Add one attribute
                    layer.CreateField(gdal.ogr.FieldDefn('Survey', gdal.ogr.OFTString))
                    defn = layer.GetLayerDefn()
                    # Create a new feature (attribute and geometry)
                    feat = gdal.ogr.Feature(defn)
                    feat.SetField('Survey', name)
                    feat.SetGeometry(ds_geom)
                    layer.CreateFeature(feat)
            dest_ds = None
        source_ds = None
        return outname
    
    def _dest_filename(self, filename: str, dest_dir: str):
        """
        Build the filename for the transformed file.
        
        Parameters
        ----------
        
        filename : str
            The filename of the source file
            
        dest_dir : str
            The path to the target directory for a tranlated file.
            
        Returns
        -------
        str
            The resulting path to the transformed file.
        """
        root, ext = os.path.splitext(filename)
        path, base = os.path.split(root)
        newfname = os.path.join(dest_dir, base + ext)
        return newfname, ext