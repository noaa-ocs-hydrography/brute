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
from affine import Affine
from fiona._crs import CRSError
from fiona.transform import transform_geom
from fuse.raw_read.noaa.bag import BAGRawReader
from osgeo import osr, gdal
from rasterio.crs import CRS
from rasterio.warp import transform as transform_points


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
                        input_transform = input_raster.transform
                        input_shape = input_raster.height, input_raster.width

                        # # TODO uncomment if we ever need to reproject between systems with different units
                        # left, bottom, right, top = bounds_from_opposite_corners(input_origin,
                        #                                                         input_origin + numpy.flip(input_shape) * input_resolution)
                        # output_transform, output_width, output_height = calculate_default_transform(input_crs, output_crs,
                        #                                                                             width=input_shape[1],
                        #                                                                             height=input_shape[0], left=left,
                        #                                                                             bottom=bottom, right=right, top=top)

                        output_transform = reproject_transform(input_transform, input_crs, output_crs)

                        input_data = input_raster.read()
                        with rasterio.open(output_filename, 'w', 'GTiff', width=input_shape[1], height=input_shape[0],
                                           count=input_raster.count, crs=output_crs, transform=output_transform, dtype=input_data.dtype,
                                           nodata=input_raster.nodata) as output_raster:
                            output_raster.write(input_data)
                            # for band_index in range(1, input_raster.count + 1):
                            #     reproject(rasterio.band(input_raster, band_index), rasterio.band(output_raster, band_index),
                            #               resampling=Resampling.min)

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
    Get a reprojected GDAL dataset from the given file.

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

    input_crs = CRS.from_string(dataset.GetProjectionRef())
    spatial_reference = spatial_reference_from_metadata(metadata)
    output_geotransform = reproject_transform(dataset.GetGeoTransform(), input_crs, spatial_reference.ExportToWkt())

    dataset.SetProjection(spatial_reference.ExportToWkt())
    dataset.SetGeoTransform(output_geotransform)

    return dataset


def reproject_transform(transform: Union[tuple, Affine], input_crs: CRS, output_crs: CRS) -> Union[tuple, Affine]:
    """
    Update the origin of the given transform to the given CRS while keeping the resolution.

    Parameters
    ----------
    transform
        GDAL or affine transform
    input_crs
        CRS of input
    output_crs
        desired CRS

    Returns
    -------
    Union[tuple, Affine]
        GDAL or affine transform in new projection
    """

    try:
        input_crs = CRS.from_user_input(input_crs)
    except CRSError:
        raise NotImplementedError(f'could not parse input CRS of type "{type(input_crs)}"')

    try:
        output_crs = CRS.from_user_input(output_crs)
    except CRSError:
        raise NotImplementedError(f'could not parse output CRS of type "{type(output_crs)}"')

    if type(transform) is tuple:
        output_origin = numpy.ravel(transform_points(input_crs, output_crs, [transform[0]], [transform[3]]))
        return output_origin[0], transform[1], transform[2], output_origin[1], transform[4], transform[5]
    elif type(transform) is Affine:
        output_origin = numpy.ravel(transform_points(input_crs, output_crs, [transform.c], [transform.f]))
        return rasterio.transform.from_origin(*output_origin, transform.a, transform.e * -1)
    else:
        raise NotImplementedError(f'could not parse transform of type "{type(transform)}"')


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
    Union[CRS, osr.SpatialReference]
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
