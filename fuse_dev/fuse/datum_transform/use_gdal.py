# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:01:28 2019

@author: Casiano.Koprowski
"""

import logging
import os
from typing import Union

import fiona
import numpy
import rasterio
from fiona.transform import transform_geom
from fuse.raw_read.noaa.bag import BAGRawReader
from fuse.utilities import bounds_from_opposite_corners
from osgeo import osr, gdal
from rasterio.crs import CRS
from rasterio.warp import reproject, calculate_default_transform, Resampling


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

        output_crs = spatial_reference_from_metadata(metadata, fiona_crs=True)
        for input_filename in support_filenames:
            basename, extension = os.path.splitext(input_filename)
            basename = os.path.basename(basename)
            output_filename = os.path.join(output_directory, basename + extension)
            if '.tif' in extension:
                with rasterio.open(input_filename) as input_raster:
                    input_crs = input_raster.crs
                    if output_crs != input_crs:
                        input_shape = input_raster.height, input_raster.width
                        input_transform = input_raster.transform
                        input_resolution = numpy.array((input_transform.a, input_transform.e))
                        input_origin = numpy.array((input_transform.c, input_transform.f))

                        left, bottom, right, top = bounds_from_opposite_corners(input_origin,
                                                                                input_origin + numpy.flip(input_shape) * input_resolution)
                        output_transform, output_width, output_height = calculate_default_transform(input_crs, output_crs,
                                                                                                    width=input_shape[1],
                                                                                                    height=input_shape[0], left=left,
                                                                                                    bottom=bottom, right=right, top=top)

                        with rasterio.open(output_filename, 'w', 'GTiff', width=output_width, height=output_height,
                                           count=input_raster.count, crs=output_crs, transform=output_transform,
                                           dtype=rasterio.float32, nodata=input_raster.nodata) as output_raster:
                            for band_index in range(1, input_raster.count + 1):
                                reproject(rasterio.band(input_raster, band_index), rasterio.band(output_raster, band_index),
                                          resampling=Resampling.min)

                        if not os.path.exists(output_filename):
                            logging.warning(f'file not created: {output_filename}')

                        reprojected_filenames.append(output_filename)
                    else:
                        reprojected_filenames.append(input_filename)
            elif extension == '.gpkg':
                output_filename = _reproject_geopackage(input_filename, output_filename, output_crs)
                reprojected_filenames.append(output_filename)
            else:
                if extension != '.tfw':
                    logging.warning(f'unsupported file format "{extension}"')

                reprojected_filenames.append(input_filename)
        metadata['support_files'] = reprojected_filenames

    return metadata


def _reproject_via_geotransform(filename: str, metadata: dict, reader: BAGRawReader) -> gdal.Dataset:
    """
    Get a GDAL reprojected dataset from the given file.

    Parameters
    ----------
    filename
        filename of data file
    metadata
        dictionary of metadata

    Returns
    -------
    gdal.Dataset
        reprojected GDAL dataset
    """

    dataset = reader.read_bathymetry(filename, None)

    input_geotransform = dataset.GetGeoTransform()
    input_origin = numpy.array((input_geotransform[0], input_geotransform[3]))
    input_resolution = numpy.array((input_geotransform[1], input_geotransform[5]))
    input_shape = dataset.RasterYSize, dataset.RasterXSize

    input_crs = CRS.from_string(dataset.GetProjectionRef())
    output_crs = spatial_reference_from_metadata(metadata, fiona_crs=True)

    left, bottom, right, top = bounds_from_opposite_corners(input_origin, input_origin + numpy.flip(input_shape) * input_resolution)
    output_transform, output_width, output_height = calculate_default_transform(input_crs, output_crs, width=input_shape[1],
                                                                                height=input_shape[0], left=left, bottom=bottom,
                                                                                right=right, top=top)

    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromProj4(output_crs.to_string())

    dataset.SetProjection(spatial_reference.ExportToWkt())
    dataset.SetGeoTransform(output_transform.to_gdal())

    return dataset


def _reproject_geopackage(input_filename: str, output_filename: str, output_crs: CRS, input_layer=None) -> str:
    """
    Convert a GeoPackage to the provided reference frame, assuming a single layer.

    Parameters
    ----------
    input_filename
        file path of the input coverage file
    output_filename
        file path to which to save reprojected file
    output_crs
        OSR spatial reference of output
    input_layer
        name of layer in input file

    Returns
    -------
    str
        file name of the resulting file
    """

    # TODO implement handling for rasters
    output_layer = os.path.splitext(os.path.basename(input_filename))[0]

    with fiona.open(input_filename, layer=input_layer) as input_layer:
        input_crs = input_layer.crs

        # check to see if the file is already projected in the given CRS
        if output_crs != input_crs:
            records = [{'id': record['id'], 'geometry': transform_geom(input_crs, output_crs, record['geometry']),
                        'properties': record['properties']} for record in input_layer]

            with fiona.open(output_filename, 'w', 'GPKG', layer=output_layer, schema=input_layer.schema, crs=output_crs) as output_layer:
                output_layer.writerecords(records)
        else:
            output_filename = input_filename

    return output_filename


def spatial_reference_from_metadata(metadata: dict, fiona_crs: bool = False) -> Union[CRS, osr.SpatialReference]:
    """
    Build an OSR spatial reference from the `to_horiz_*` fields in the provided metadata.
    If 'to_horiz_frame','to_horiz_type' and 'to_horiz_key' are not populated, raise an error.

    Parameters
    ----------
    metadata
        dictionary of metadata

    Returns
    -------
    CRS
        horizontal coordinate reference system
    """

    if all(key in metadata for key in ('to_horiz_frame', 'to_horiz_type', 'to_horiz_key')):
        if metadata['to_horiz_type'] == 'utm':
            proj4_string = f'+proj=utm +zone={metadata["to_horiz_key"]} +datum={metadata["to_horiz_frame"]}'
            if fiona_crs:
                return CRS.from_string(proj4_string)
            else:
                spatial_reference = osr.SpatialReference()
                spatial_reference.ImportFromProj4(proj4_string)
                return spatial_reference
        else:
            # TODO We still need to sort our when we are not working in utm...
            raise NotImplementedError('working outside of UTM is not implemented')
    else:
        raise ValueError('Not all metadata is available to build the PROJ4 string')


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
