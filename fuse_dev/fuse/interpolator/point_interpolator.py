# -*- coding: utf-8 -*-
"""
point_interpolator.py

grice 20180625

V0.0.4 201901115

This is a collection of functions for creating interpolated bathymetry for the
national bathymetry.

This is a rework of the original point interpolator.
"""

from datetime import datetime
from re import match
from types import GeneratorType
from typing import Union, Tuple

import fiona
import numpy as np
from affine import Affine
from osgeo import gdal, osr
from rasterio.features import geometry_mask
from scipy.spatial.ckdtree import cKDTree
from scipy.spatial.qhull import Delaunay
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, MultiPoint, Point, MultiLineString, JOIN_STYLE
from shapely.geometry import shape as shapely_shape
from shapely.ops import unary_union, polygonize

CATZOC = {
    'A1': [.01, .5],
    'A2': [.02, 1],
    'B': [.02, 1],
    'C': [.05, 2]
}


# ============================================================================

def interpolate(dataset: gdal.Dataset, resolution: float, ancillary_coverage_files: [str] = None,
                catzoc: str = 'B', nodata=1000000) -> gdal.Dataset:
    """
    Generate a raster from the input data using the given interpolation method.

    Parameters
    ----------
    dataset
        GDAL dataset (point cloud)
    ancillary_coverage_files
        filenames of USACE survey frameworks
    resolution
        output raster resolution 
    catzoc
        CATZOC score ('A1', 'A2'. 'B', or 'C') for the output uncertainty layer

    Returns
    -------
    gdal.Dataset
        interpolated GDAL raster dataset
    """

    if type(resolution) is not float:
        resolution = float(resolution)

    catzoc = catzoc.upper()
    assert catzoc in CATZOC, f'CATZOC score "{catzoc}" not supported'

    points = _gdal_to_xyz(dataset)[:, :2]

    # get the data bounds as (ulx,uly,lrx,lry)
    tmp = np.concatenate((np.min(points, axis=0), np.max(points, axis=0)))
    input_bounds = tmp[[0, 3, 2, 1]]

    # set the size of the interpolation window size to be the maximum of all nearest neighbor distances
    max_neighbor_distance = maximum_nearest_neighbor_distance(points)
    # but if the number is stupid small (single beam), default to resolution for interpolation
    atomic_length = max((max_neighbor_distance, resolution))

    start_time = datetime.now()
    # get the concave hull of the survey points
    approximate_alpha_hull = alpha_hull(points, max_length=atomic_length, tolerance=resolution / 2,
                                        max_nn=max_neighbor_distance, iterations=1)

    # if there is one reulting polygon we think we know where to interpolate
    if type(approximate_alpha_hull) is Polygon:
        interpolation_region = approximate_alpha_hull
    else:
        # but if there are many, see if we can sort out if this is a set line spacing survey
        uniform_density = ((np.sqrt(approximate_alpha_hull.area) / max_neighbor_distance) - 1) ** 2 / approximate_alpha_hull.area
        actual_density = len(points) / approximate_alpha_hull.area
        # if the point density is less than if uniformly distributed, this is set line spacing
        if actual_density < max((uniform_density, max_neighbor_distance)):
            interpolation_region = buffer_hull(points, buffer=atomic_length, tolerance=resolution / 2,
                                               max_nn=max_neighbor_distance)
        else:
            interpolation_region = alpha_hull(points, max_length=atomic_length, tolerance=resolution / 2,
                                              max_nn=max_neighbor_distance)

    print(f'concave hull calculation took {datetime.now() - start_time}')

    # Make sure the channel has no bathymetry holes
    if ancillary_coverage_files is not None and len(ancillary_coverage_files) > 0:
        ancillary_coverage = []
        for ancillary_coverage_filename in ancillary_coverage_files:
            with fiona.open(ancillary_coverage_filename) as ancillary_coverage_file:
                ancillary_coverage.extend(shapely_shape(feature['geometry']) for feature in ancillary_coverage_file)

    interpolation_region = unary_union([interpolation_region] + ancillary_coverage)

    output_shape, output_bounds = _get_raster_properties(resolution, input_bounds)

    interpolated_raster = gdal.Grid('', dataset, format='MEM', width=output_shape[1], height=output_shape[0], outputBounds=output_bounds,
                                    algorithm=f'linear:radius=0:nodata={int(nodata)}')
    interpolated_band = interpolated_raster.GetRasterBand(1)
    interpolated_values = interpolated_band.ReadAsArray()
    geotransform = interpolated_raster.GetGeoTransform()

    uncertainty = _uncertainty(np.where(interpolated_values != nodata, interpolated_values, np.nan), catzoc)

    uncertainty[np.isnan(uncertainty)] = nodata

    output_raster = gdal.GetDriverByName('MEM').Create('', int(output_shape[1]), int(output_shape[0]), 2, gdal.GDT_Float32)
    output_raster.SetGeoTransform(geotransform)
    output_raster.SetProjection(crs_wkt_from_gdal(dataset))

    band_1 = output_raster.GetRasterBand(1)
    band_1.SetNoDataValue(nodata)
    band_1.WriteArray(interpolated_values)
    del band_1

    band_2 = output_raster.GetRasterBand(2)
    band_2.SetNoDataValue(nodata)
    band_2.WriteArray(uncertainty)
    del band_2

    interpolation_mask = raster_mask_like(interpolation_region, output_raster)
    interpolated_dataset = apply_raster_mask(output_raster, interpolation_mask)

    return interpolated_dataset


def _gdal_to_xyz(dataset: gdal.Dataset) -> np.array:
    """
    Take a gdal vector xyz point cloud and return a numpy array.

    Parameters
    ----------
    dataset: gdal.Dataset :
        TODO write description

    Returns
    -------
    type
        vector dataset

    """

    # get the data out of the gdal data structure
    lyr = dataset.GetLayerByIndex(0)
    count = lyr.GetFeatureCount()
    data = np.zeros((count, 3))

    for n in np.arange(count):
        f = lyr.GetFeature(n)
        data[n, :] = f.geometry().GetPoint()

    return data


def nearest_neighbors(points: np.array) -> (np.array, np.array):
    """
    Get the nearest neighbors to each of the given points.

    Parameters
    ----------
    points
        N x M array of points

    Returns
    -------
    numpy.array, numpy.array
        distances to and indices of each nearest neighbor
    """

    if type(points) is GeneratorType:
        points = MultiPoint(list(points))

    return cKDTree(points).query(points, k=[2])


def maximum_nearest_neighbor_distance(points: np.array) -> float:
    """
    Get the maximum of all closest neighbor distances within the given set of points.

    Parameters
    ----------
    points
        N x M array of points

    Returns
    -------
    float
        maximum nearest neighbor distance
    """

    return np.max(nearest_neighbors(points)[0])


def _get_raster_properties(resolution: float, data_bounds: Tuple[float, float, float, float]) -> Tuple[
    int, int, Tuple[float, float, float, float]]:
    """
    Get the bounds and number of rows and columns for a raster given the
    desired resolution and the data max and min in X and Y values.  This algorithm use
    the data min as the anchor for the rows and columns, but floored to the
    nearest multiple of the resolution.

    Parameters
    ----------
    resolution:
        Resolution of the desired raster grid.

    data_bounds:
        Data bounds as (ulx,uly,lrx,lry)


    Returns
    -------
    type
        output shape as (number of rows, number of columns), and bounds as
        (ulx,uly,lrx,lry)

    """

    xmin, ymax, xmax, ymin = data_bounds
    numcolumns = int(np.ceil((1. * (xmax - xmin) / resolution)))
    xrem = xmin % resolution
    xnmin = xmin - xrem
    xnmax = xnmin + 1. * numcolumns * resolution
    numrows = int(np.ceil((1. * (ymax - ymin) / resolution)))
    yrem = ymin % resolution
    ynmin = ymin - yrem
    ynmax = ynmin + 1. * numrows * resolution
    return (numrows, numcolumns), (xnmin, ynmax, xnmax, ynmin)


def alpha_hull(points: np.array, max_length: float = None, tolerance: float = None, max_nn: float = None, iterations: int = 50,
               triangles: Delaunay = None) -> MultiPolygon:
    """
    Calculate the concave hull of the given points using triangulation.
    inspired by https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html

    Parameters
    ----------
    points
        N x M array of points
    max_length
        maximum length to include in the convex boundary
    tolerance
        distance within generated hull to dissolve points
    max_nn
        maximum of nearest neighbor distances in the given points
    iterations
        number of iterations to complete (default will continue until all points are within the hull)
    triangles
        Delaunay triangulation of points

    Returns
    ----------
    MultiPolygon
        alpha shape (concave hull) of points
    """

    if len(points) < 4:
        return MultiPoint(points).convex_hull

    if max_nn is None:
        max_nn = maximum_nearest_neighbor_distance(points)

    if max_length is None:
        max_length = max_nn

    if tolerance is None:
        tolerance = max_nn / 2

    if triangles is None:
        triangles = Delaunay(points)

    indices = np.stack(((triangles.simplices[:, 0], triangles.simplices[:, 1]),
                        (triangles.simplices[:, 1], triangles.simplices[:, 2]),
                        (triangles.simplices[:, 2], triangles.simplices[:, 0])), axis=1).T
    edges = points[indices]
    vectors = np.squeeze(np.diff(edges, axis=2))
    lengths = np.hypot(vectors[:, :, 0], vectors[:, :, 1])
    semiperimeters = np.sum(lengths, axis=1) / 2
    areas = np.sqrt(semiperimeters * np.product(np.expand_dims(semiperimeters, axis=1) - lengths, axis=1))
    circumcircle_radii = np.product(lengths, axis=1) / (4 * areas)
    indices = indices[circumcircle_radii <= max_length]

    if len(indices) > 0:
        indices = np.sort(np.reshape(indices, (indices.shape[0] * indices.shape[1], indices.shape[2])), axis=1)
        boundary_edge_indices, counts = np.unique(indices, axis=0, return_counts=True)
        boundary_edge_indices = boundary_edge_indices[counts == 1]

        polygons = remove_interior_duplicates(polygonize(MultiLineString(points[boundary_edge_indices].tolist())))
        output_hull = MultiPolygon(polygons).buffer(0)
        if output_hull.area == 0:
            output_hull = unary_union(polygons)

        points_geometry = MultiPoint(points)

        points_outside_hull = points_geometry.difference(output_hull)
        if type(points_outside_hull) is Point:
            points_outside_hull = GeometryCollection([points_outside_hull])

        if len(points_outside_hull) > 0:
            segments_before_dissolve = len(output_hull) if type(output_hull) is MultiPolygon else 1
            output_hull = GeometryCollection((output_hull, points_outside_hull)).buffer(max_nn)

            if iterations > 1:
                output_hull = output_hull.buffer(-max_nn, join_style=JOIN_STYLE.mitre)
                segments_after_dissolve = len(output_hull) if type(output_hull) is MultiPolygon else 1

                points_outside_hull = points_geometry.difference(output_hull.buffer(tolerance))
                if type(points_outside_hull) is Point:
                    points_outside_hull = GeometryCollection([points_outside_hull])

                if segments_after_dissolve > segments_before_dissolve or len(points_outside_hull) > 0:
                    print(f'alpha hull created from length threshold {max_length:.5} contains ' +
                          f'{(len(points) - len(points_outside_hull)) / len(points) * 100:.5}% of points ({len(points_outside_hull)} exist outside of the current hull)')
                    output_hull = alpha_hull(points, max_length + max_nn, tolerance, max_nn, iterations - 1, triangles)
    else:
        print(f'no edges were found to be shorter than the length threshold ({max_length:.5})')
        output_hull = alpha_hull(points, max_length + max_nn, tolerance, max_nn, iterations, triangles)

    return output_hull


def remove_interior_duplicates(polygons: [Polygon]) -> [Polygon]:
    """
    Given a list of polygons, return a list without polygons whose exteriors represent the interior of another polygon.

    Parameters
    ----------
    polygons
        list of polygons

    Returns
    -------
    [Polygon]
        polygons without duplicates of interiors
    """

    if type(polygons) is not list:
        polygons = list(polygons)

    for polygon_with_interiors in (polygon for polygon in polygons if len(polygon.interiors) > 0):
        for interior in polygon_with_interiors.interiors:
            interior_polygon = Polygon(interior)
            duplicate_polygons = []

            for polygon_index in range(len(polygons)):
                polygon = polygons[polygon_index]
                if polygon.equals(interior_polygon):
                    duplicate_polygons.append(polygon)

            for polygon in duplicate_polygons:
                polygons.remove(polygon)

    return polygons


def buffer_hull(points: Union[MultiPoint, np.array], buffer: float = None, tolerance: float = None, iterations: int = 50,
                max_nn: float = None):
    """
    Calculate the concave hull of the given points by buffering from points.

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    buffer
        length by which to expand
    tolerance
        distance within generated hull to dissolve points
    iterations
        number of iterations to complete (default will continue until all points are within the hull)
    max_nn
        maximum of nearest neighbor distances in the given points

    Returns
    ----------
    MultiPolygon
        alpha shape (concave hull) of points
    """

    if type(points) is not MultiPoint:
        points = MultiPoint(points)

    if max_nn is None:
        max_nn = maximum_nearest_neighbor_distance(points)

    if buffer is None:
        buffer = max_nn

    if tolerance is None:
        tolerance = max_nn / 2

    max_cluster_distance = max_nn * 10

    output_hull = points.buffer(buffer)

    if type(output_hull) is MultiPolygon:
        cluster_distances = [distance for distance in polygon_nearest_neighbor_distances(output_hull) if distance < max_cluster_distance]
        if len(cluster_distances) > 0:
            cluster_max_nn = np.max(cluster_distances)
            output_hull = output_hull.buffer(cluster_max_nn / 2).buffer(-cluster_max_nn / 2)

    output_hull = output_hull.buffer(-buffer, join_style=JOIN_STYLE.mitre)

    points_outside_hull = points.difference(output_hull.buffer(tolerance))
    if type(points_outside_hull) is Point:
        points_outside_hull = GeometryCollection([points_outside_hull])

    if len(points_outside_hull) > 0:
        output_hull = GeometryCollection((output_hull, points_outside_hull)).buffer(max_nn)

        if iterations > 1:
            output_hull.buffer(-max_nn, join_style=JOIN_STYLE.mitre)

            points_outside_hull = points.difference(output_hull.buffer(tolerance))
            if type(points_outside_hull) is Point:
                points_outside_hull = GeometryCollection([points_outside_hull])

            if len(points_outside_hull) > 0:
                print(f'buffer hull created by expanding and shrinking by {buffer:.5} contains ' +
                      f'{(len(points) - len(points_outside_hull)) / len(points) * 100:.5}% of points ({len(points_outside_hull)} exist outside of the current hull)')
                output_hull = buffer_hull(points, buffer + max_nn, tolerance, iterations - 1, max_nn)

    if type(output_hull) is MultiPolygon:
        cluster_distances = [distance for distance in polygon_nearest_neighbor_distances(output_hull) if distance < max_cluster_distance]
        if len(cluster_distances) > 0:
            cluster_max_nn = max(np.ravel(cluster_distances)) / 2
            output_hull = output_hull.buffer(cluster_max_nn).buffer(-cluster_max_nn)

    return output_hull


def polygon_nearest_neighbor_distances(polygons: [Polygon]) -> [float]:
    """
    Get nearest neighbor distances within the given polygons.

    Parameters
    ----------
    polygons
        list of Shapely polygons

    Returns
    -------
    [float]
        distances to the nearest neighbor
    """

    # return [min(polygon.distance(other_polygon) for other_polygon in polygons if not other_polygon.equals(polygon)) for polygon in polygons]
    return nearest_neighbors(polygon.centroid for polygon in polygons)[0]


def _uncertainty(values: np.array, catzoc: str) -> np.array:
    """
    Calculate the uncertainty from the given interpolated values.

    Parameters
    ----------
    values
        array of data
    catzoc
        CATZOC score

    Returns
    -------
    numpy.array
        array of uncertainty values
    """

    m, b = CATZOC[catzoc]
    return (values * m) + b


def raster_mask(polygon: Union[Polygon, MultiPolygon], shape: (float, float), origin: (float, float), resolution: (float, float),
                nodata: float, crs_wkt: str) -> gdal.Dataset:
    """
    Convert the given Shapely polygon into a GDAL raster dataset mask (with values outside of the polygon set to `nodata`).

    Parameters
    ----------
    polygon
        Shapely polygon (or multipolygon) with which to create mask
    shape
        shape of mask
    resolution
        cell size of mask
    origin
        raster origin
    crs_wkt
        spatial reference of mask
    nodata
        value for no data in mask

    Returns
    -------
    gdal.Dataset
        GDAL raster dataset
    """

    mask, transform = rasterize_polygon(polygon, shape, origin, resolution)

    output_dataset = gdal.GetDriverByName('MEM').Create('', int(shape[1]), int(shape[0]), 1, gdal.GDT_Float32)
    output_dataset.SetProjection(crs_wkt)
    output_dataset.SetGeoTransform(transform.to_gdal())
    output_band = output_dataset.GetRasterBand(1)
    output_band.SetNoDataValue(nodata)
    output_band.WriteArray(np.where(mask, 1, nodata))

    return output_dataset


def raster_mask_like(polygon: Union[Polygon, MultiPolygon], like_raster: gdal.Dataset, band: int = 0) -> gdal.Dataset:
    """
    Get a mask of the survey area with the same properties as the given raster.

    Parameters
    ----------
    polygon
        Shapely polygon with which to create mask
    like_raster
        raster to copy
    band
        zero-based index of raster band

    Returns
    -------
    gdal.Dataset
        mask
    """

    geotransform = like_raster.GetGeoTransform()
    raster_band = like_raster.GetRasterBand(band + 1)
    return raster_mask(polygon, shape=(like_raster.RasterYSize, like_raster.RasterXSize), origin=(geotransform[0], geotransform[3]),
                       resolution=(geotransform[1], geotransform[5]), nodata=raster_band.GetNoDataValue(),
                       crs_wkt=crs_wkt_from_gdal(like_raster))


def apply_raster_mask(raster: gdal.Dataset, mask: gdal.Dataset, mask_value: float = None, band: int = 0) -> gdal.Dataset:
    """
    Mask a raster using the given mask values in the given mask band.
    Both rasters are assumed to be collocated, as all operations are conducted on the pixel level.

    Parameters
    ----------
    raster
        GDAL raster dataset
    mask
        GDAL raster dataset of mask to apply
    mask_value
        value to use as mask
    band
        zero-based index of raster band of mask

    Returns
    -------
    gdal.Dataset
        GDAL raster dataset with mask applied
    """

    band = mask.GetRasterBand(band + 1)
    mask_array = band.ReadAsArray()

    if mask_value is None:
        mask_value = band.GetNoDataValue()

    for band_index in range(1, raster.RasterCount + 1):
        raster_band = raster.GetRasterBand(band_index)
        raster_array = raster_band.ReadAsArray()
        raster_array[mask_array == mask_value] = raster_band.GetNoDataValue()
        raster_band.WriteArray(raster_array)

    return raster


def crs_wkt_from_gdal(dataset: gdal.Dataset, layer: int = 0) -> str:
    """
    Extract the well-known text of the CRS from the given GDAL dataset.

    Parameters
    ----------
    dataset
        GDAL dataset (vector or raster)
    layer
        index of vector layer

    Returns
    -------
    str
        well-known text of CRS
    """

    crs_wkt = dataset.GetProjectionRef()

    if crs_wkt == '':
        vector_layer = dataset.GetLayerByIndex(layer)
        if vector_layer is not None:
            crs_wkt = vector_layer.GetSpatialRef().ExportToWkt()
        else:
            raise ValueError("No Spatial Ref Found?")
    elif match('^EPSG:[0-9]+$', crs_wkt):
        crs_wkt = crs_wkt_from_epsg(int(crs_wkt[5:]))

    return crs_wkt


def rasterize_polygon(polygon: Union[Polygon, MultiPolygon], shape: (float, float), origin: (float, float),
                      resolution: (float, float)) -> (np.array, Affine):
    """
    Convert the given Shapely polygon into a boolean mask with its requisite transform.

    Parameters
    ----------
    polygon
        Shapely polygon (or multipolygon) with which to create mask
    shape
        shape of output mask
    resolution
        cell size of output mask
    origin
        origin of raster

    Returns
    -------
    numpy.array, Affine
        boolean array mask (True is within the polygon, False outside) and its affine transform
    """

    if type(polygon) is Polygon:
        polygon = MultiPolygon([polygon])

    transform = affine_from_georeference(origin, resolution)
    return ~geometry_mask(polygon, shape, transform), transform


def crs_wkt_from_epsg(crs: Union[int, str]) -> str:
    """
    Get the well-known text of a CRS from EPSG code, SPCS string, Esri string, or well-known text of CRS.

    Parameters
    ----------
    crs
        EPSG code, SPCS string, Esri string, or well-known text of CRS

    Returns
    -------
    str
        well-known text of CRS
    """

    if type(crs) is str:
        spatial_reference = osr.SpatialReference(wkt=crs)
        spatial_reference.MorphFromESRI()
    else:
        spatial_reference = osr.SpatialReference()
        spatial_reference.ImportFromEPSG(crs)

    return spatial_reference.ExportToWkt()


def affine_from_georeference(origin: (float, float), resolution: (float, float), rotation: (float, float) = None) -> Affine:
    """
    Calculate the affine transformation from the given geographic parameters.

    Parameters
    ----------
    origin
        XY coordinates of origin
    resolution
        XY cell size

    Returns
    -------
        Affine transform
    """

    if type(resolution) is float:
        resolution = resolution, resolution

    if rotation is None:
        rotation = 0, 0

    return Affine.translation(*origin) * Affine.scale(*resolution) * Affine.rotation(rotation[0])
