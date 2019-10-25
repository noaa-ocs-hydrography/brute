# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:01:28 2019

@author: Casiano.Koprowski
"""

import logging
import os

from fuse.raw_read.noaa.bag import BAGRawReader
from osgeo import osr, gdal, ogr


def reproject_support_files(metadata: dict, output_directory: str) -> dict:
    """
    Horizontally transform the given support files, writing reprojected files to the given output directory.

    Parameters
    ----------
    metadata
        dictionary of metadata
    output_directory
        path to directory for transformed files

    Returns
    -------
    dict
        dictionary of metadata with updated list of support files
    """

    if 'support_files' in metadata:
        support_filenames = metadata['support_files']
        reprojected_filenames = []
        output_spatial_reference = _spatial_reference_from_metadata(metadata)
        for input_filename in support_filenames:
            basename, extension = os.path.splitext(input_filename)
            basename = os.path.basename(basename)
            output_filename = os.path.join(output_directory, basename + extension)
            if '.tif' in extension:
                input_dataset = gdal.Open(input_filename)
                if input_dataset is not None:
                    input_spatial_reference = osr.SpatialReference(wkt=input_dataset.GetProjectionRef())
                    if output_spatial_reference.IsSame(input_spatial_reference):
                        gdal.Warp(output_filename, input_dataset,
                                  options=gdal.WarpOptions(dstSRS=output_spatial_reference, format='GTiff'))

                        if not os.path.exists(output_filename):
                            logging.warning(f'file not created: {output_filename}')

                        reprojected_filenames.append(output_filename)
                    else:
                        reprojected_filenames.append(input_filename)
                else:
                    print(f'unable to open file with GDAL: {input_filename}')
                    reprojected_filenames.append(input_filename)
            elif extension == '.gpkg':
                output_filename = _reproject_geopackage(input_filename, output_filename, output_spatial_reference)
                reprojected_filenames.append(output_filename)
            else:
                if extension != '.tfw':
                    logging.warning(f'unsupported file format "{extension}"')

                reprojected_filenames.append(input_filename)
        metadata['support_files'] = reprojected_filenames

    return metadata


def _spatial_reference_from_metadata(metadata: dict) -> osr.SpatialReference:
    """
    Build an OSR spatial reference from the `to_horiz_*` fields in the provided metadata.
    If 'to_horiz_frame','to_horiz_type' and 'to_horiz_key' are not populated, raise an error.

    Parameters
    ----------
    metadata
        dictionary of metadata

    Returns
    -------
    osr.SpatialReference
        target horizontal reference system
    """

    if all(key in metadata for key in ('to_horiz_frame', 'to_horiz_type', 'to_horiz_key')):
        if metadata['to_horiz_type'] == 'utm':
            spatial_reference = osr.SpatialReference()
            spatial_reference.ImportFromProj4(f'+proj=utm +zone={metadata["to_horiz_key"]} +datum={metadata["to_horiz_frame"]}')
            return spatial_reference
        else:
            # TODO We still need to sort our when we are not working in utm...
            raise NotImplementedError('working outside of UTM is not implemented')
    else:
        raise ValueError('Not all metadata is available to build the PROJ4 string')


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
    gdal.Dataset
        reprojected GDAL dataset
    """

    source_dataset = reader.read_bathymetry(filename, None)
    source_geotransform = source_dataset.GetGeoTransform()
    source_spatialref = gdal.osr.SpatialReference(wkt=source_dataset.GetProjectionRef())
    dest_spatialref = _spatial_reference_from_metadata(instructions)

    coordinate_transform = gdal.osr.CoordinateTransformation(source_spatialref, dest_spatialref)
    dest_point = gdal.ogr.CreateGeometryFromWkt(f"POINT ({source_geotransform[0]} {source_geotransform[3]})")
    dest_point.Transform(coordinate_transform)

    target_geotransform = (dest_point.GetX(), source_geotransform[1], 0, dest_point.GetY(), 0, source_geotransform[5])

    source_dataset.SetGeoTransform(target_geotransform)
    source_dataset.SetProjection(dest_spatialref.ExportToWkt())

    return source_dataset


def _reproject_geopackage(input_filename: str, output_filename: str, spatial_reference: osr.SpatialReference) -> str:
    """
    Convert a GeoPackage to the provided reference frame, assuming a single layer.

    Parameters
    ----------
    input_filename
        file path of the input coverage file
    output_filename
        file path to which to save reprojected file
    spatial_reference
        OSR spatial reference of output

    Returns
    -------
    str
        file name of the resulting file
    """

    survey_name = os.path.splitext(os.path.basename(input_filename))[0]

    # Open the data source and read in the extent
    input_dataset = ogr.Open(input_filename)
    input_layer = input_dataset.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()

    # check to see if the file is already projected in the given SRS
    if spatial_reference.IsSame(input_spatial_reference):
        output_filename = input_filename
    else:
        transform = osr.CoordinateTransformation(input_spatial_reference, spatial_reference)

        output_dataset = ogr.GetDriverByName('gpkg').CreateDataSource(output_filename)
        output_layer = output_dataset.CreateLayer(survey_name, spatial_reference, ogr.wkbMultiPolygon)
        output_layer.CreateField(ogr.FieldDefn('Survey', ogr.OFTString))
        output_layer_definition = output_layer.GetLayerDefn()

        for feature in input_layer:
            if feature is not None:
                geometry = ogr.CreateGeometryFromWkt(feature.GetGeometryRef().ExportToWkt())
                geometry.Transform(transform)

                feature = ogr.Feature(output_layer_definition)
                feature.SetField('Survey', survey_name)
                feature.SetGeometry(geometry)
                output_layer.CreateFeature(feature)
        del output_dataset
    del input_dataset

    return output_filename


def _filename_in_other_directory(filename: str, directory: str) -> (str, str):
    """
    Get the given filename in a given directory.

    Parameters
    ----------
    filename
        input filename
    directory
        other directory

    Returns
    -------
    str, str
        path to output file in the given directory, and extension
    """

    return os.path.join(directory, os.path.basename(filename)), os.path.splitext(filename)[-1]


def _xyz_to_gdal(points: [(float, float, float)], utm_zone: int, vertical_datum: str) -> gdal.Dataset:
    """
    Create a GDAL point cloud from the given XYZ points.

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
    gdal.Dataset
        GDAL point cloud dataset
    """

    # positive UTM zone is in the northern hemisphere
    spatial_reference = osr.SpatialReference()
    spatial_reference.SetWellKnownGeogCS('NAD83')
    spatial_reference.SetUTM(abs(utm_zone), 1 if utm_zone > 0 else 0)
    spatial_reference.SetVertCS(vertical_datum, vertical_datum, 2000)

    dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
    layer = dataset.CreateLayer('pts', spatial_reference, geom_type=gdal.ogr.wkbPoint)
    feature_definition = layer.GetLayerDefn()

    for point in points:
        geometry = gdal.ogr.Geometry(gdal.ogr.wkbPoint)
        geometry.AddPoint(*point)
        feature = gdal.ogr.Feature(feature_definition)
        feature.SetGeometry(geometry)
        layer.CreateFeature(feature)

    return dataset
