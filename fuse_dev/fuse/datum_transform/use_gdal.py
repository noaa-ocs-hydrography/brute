# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:01:28 2019

@author: Casiano.Koprowski
"""

import gdal
import logging
import os

from fuse.raw_read.noaa.bag import BAGRawReader


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
    if 'support_files' in metadata:
        sf = metadata['support_files']
        t = []
        dest_srs = _build_gdal_horiz_srs(metadata)
        for f in sf:
            newf, ext = _dest_filename(f,dest_dir)
            if ext in ('.tif', '.tiff'):
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
                    logging.warning(f'{f} failed to open with gdal')
            elif ext == '.gpkg':
                resulting_file = _reproject_geopackage(f, newf, dest_srs)
                t.append(resulting_file)
            else:
                t.append(f)
                if ext != '.tfw':
                    logging.warning(f'Unsupported support file format: {ext}')
        metadata['support_files'] = t
    return metadata

def _build_gdal_horiz_srs(metadata):
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


def _reproject_via_geotransform(filename: str, instructions: dict, reader: BAGRawReader) -> gdal.Dataset:
    """
    Get a GDAL reprojected dataset from the given file.

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

    source_dataset = reader.read_bathymetry(filename, None)
    source_geotransform = source_dataset.GetGeoTransform()
    source_spatialref = gdal.osr.SpatialReference(wkt=source_dataset.GetProjectionRef())
    dest_spatialref =  _build_gdal_horiz_srs(instructions)

    transform = gdal.osr.CoordinateTransformation(source_spatialref, dest_spatialref)
    dest_point = gdal.ogr.CreateGeometryFromWkt(f"POINT ({source_geotransform[0]} {source_geotransform[3]})")
    dest_point.Transform(transform)

    target_geotransform = (dest_point.GetX(), source_geotransform[1], 0, dest_point.GetY(), 0, source_geotransform[5])

    source_dataset.SetGeoTransform(target_geotransform)
    source_dataset.SetProjection(dest_spatialref.ExportToWkt())

    return source_dataset


def _reproject_geopackage(fromfilename: str, tofilename: str, dest_srs: str):
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


def _dest_filename(filename: str, dest_dir: str):
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


def __xyz2gdal(points: [(float, float, float)], utm_zone: int, vertical_datum: str) -> gdal.Dataset:
        """
        Get a GDAL point cloud dataset from XYZ points.

        Parameters
        ----------
        points
            XYZ points
        utm_zone
            desired UTM zone
        vertical_datum
            name of desired vertical datum

        Returns
        -------
            GDAL point cloud dataset
        """

        # setup the gdal bucket
        spatial_reference = osr.SpatialReference()
        spatial_reference.SetWellKnownGeogCS('NAD83')
        # positive UTM zone is in the northern hemisphere
        spatial_reference.SetUTM(abs(utm_zone), 1 if utm_zone > 0 else 0)
        spatial_reference.SetVertCS(vertical_datum, vertical_datum, 2000)
        dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = dataset.CreateLayer('pts', spatial_reference, geom_type=gdal.ogr.wkbPoint)
        for point in points:
            geometry = gdal.ogr.Geometry(gdal.ogr.wkbPoint)
            geometry.AddPoint(*point)
            feature = gdal.ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(geometry)
            layer.CreateFeature(feature)
        return dataset