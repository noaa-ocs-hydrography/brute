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
from fuse.raw_read.noaa.bag import BAGRawReader, Open, BagToGDALConverter
from fuse.utilities import bounds_from_opposite_corners
from osgeo import osr, gdal
from rasterio.crs import CRS
from rasterio.warp import transform as transform_points, calculate_default_transform


def _maxValue(arr: numpy.array):
    """
    Returns the most used value in the array as an integer

    Takes an input array and finds the most used value in the array, this
    value is used by the program to assume the array's nodata value

    Parameters
    ----------
    arr: _np.array :
        An input array

    Returns
    -------

    """

    nums, counts = numpy.unique(arr, return_counts=True)
    index = numpy.where(counts == numpy.amax(counts))
    return int(nums[index])


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
        for input_filename in support_filenames:
            basename, extension = os.path.splitext(input_filename)
            basename = os.path.basename(basename)
            output_filename = os.path.join(output_directory, basename + extension)
            if '.tif' in extension:
                input_dataset = gdal.Open(input_filename)
                input_crs = osr.SpatialReference(wkt=input_dataset.GetProjectionRef())
                nodata = input_dataset.GetRasterBand(1).GetNoDataValue()
                if nodata is None:
                    nodata =_maxValue(input_dataset.GetRasterBand(1).ReadAsArray())
                del input_dataset
                output_crs = spatial_reference_from_metadata(metadata)
                if not input_crs.IsSame(output_crs):
                    options = gdal.WarpOptions(format='GTiff', srcSRS=input_crs, dstSRS=output_crs, srcNodata=nodata, dstNodata=nodata)
                    gdal.Warp(output_filename, input_filename, options=options)
                    if not os.path.exists(output_filename):
                        logging.warning(f'file not created: {output_filename}')
                    else:
                        reprojected_filenames.append(output_filename)
            elif extension == '.gpkg':
                output_crs = spatial_reference_from_metadata(metadata, fiona_crs=True)
                output_filename = _reproject_geopackage(input_filename, output_filename, output_crs)
                reprojected_filenames.append(output_filename)
            else:
                if extension != '.tfw':
                    logging.warning(f'unsupported file format "{extension}"')
                reprojected_filenames.append(input_filename)
        metadata['support_files'] = reprojected_filenames

    return metadata


def _reproject_via_reprojectimage(filename: str, metadata: dict, reader: BAGRawReader) -> gdal.Dataset:
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

    bag = Open(filename)
    input_crs = osr.SpatialReference(wkt=bag.wkt)
    output_crs = spatial_reference_from_metadata(metadata)
    coord_transform = osr.CoordinateTransformation(input_crs, output_crs)
    input_shape = bag.shape
    input_ul, input_lr = bag.bounds
    input_resolution = bag.resolution

    ulx, uly, ulz = coord_transform.TransformPoint(input_ul[0], input_ul[1])
    llx, lly, llz = coord_transform.TransformPoint(input_ul[0], input_lr[1])
    lrx, lry, lrz = coord_transform.TransformPoint(input_lr[0], input_lr[1])
    urx, ury, urz = coord_transform.TransformPoint(input_lr[0], input_ul[1])

    x_vals = [ulx, llx, lrx, urx]
    y_vals = [uly, lly, lry, ury]

    min_x, max_x = min(x_vals), max(x_vals)
    min_y, max_y = min(y_vals), max(y_vals)

    reprojected_shape = [int((max_y - min_y)/input_resolution[0]), int((max_x - min_x)/input_resolution[0])]

    reprojected_geotransform = (min_x, input_resolution[0], 0, max_y, 0, input_resolution[1])

    memory_driver = gdal.GetDriverByName('MEM')
    reprojected_dataset = memory_driver.Create('', reprojected_shape[1], reprojected_shape[0], 2, gdal.GDT_Float32)
    for band_index in (1, 2):
        band_out = reprojected_dataset.GetRasterBand(band_index)
        band_out.SetNoDataValue(1000000.0)
        band_out.WriteArray(numpy.full(reprojected_shape, 1000000.0))
        del band_out
    reprojected_dataset.SetProjection(output_crs.ExportToWkt())
    reprojected_dataset.SetGeoTransform(reprojected_geotransform)

    build_dataset = BagToGDALConverter()
    build_dataset.bag2gdal(bag)
    del bag

    dataset = build_dataset.dataset
    gdal.ReprojectImage(dataset, reprojected_dataset, input_crs.ExportToWkt(), output_crs.ExportToWkt(), gdal.GRA_Min)
    del build_dataset, dataset

    return reprojected_dataset


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
