from re import match
from typing import Union

import fiona
import fiona.crs
import numpy
import rasterio
import rasterio.crs
from affine import Affine
from matplotlib import pyplot
from matplotlib.cm import get_cmap
from osgeo import gdal, osr
from rasterio import MemoryFile
from rasterio.features import shapes as rasterio_shapes
from rasterio.mask import mask
from scipy.spatial import Delaunay, cKDTree
from shapely.geometry import Polygon, MultiLineString, MultiPolygon, shape as shapely_shape, mapping
from shapely.ops import unary_union, polygonize


def gdal_to_xyz(dataset: gdal.Dataset, index: int = 0, nodata: float = None) -> numpy.array:
    """
    Extract XYZ points from a GDAL dataset.

    Parameters
    ----------
    dataset
        GDAL dataset
    index
        index of layer (for vector, 0-indexed) or band (for raster, 1-indexed) to read in given dataset
    nodata
        value to exclude from point creation

    Returns
    -------
    numpy.array
        N x 3 array of XYZ points
    """

    if dataset.GetProjectionRef() == '':
        if nodata is None:
            nodata = numpy.nan

        point_layer = dataset.GetLayerByIndex(index)
        num_points = point_layer.GetFeatureCount()

        output_points = numpy.empty((num_points, 3))
        for point_index in range(num_points):
            feature = point_layer.GetFeature(point_index)
            output_point = feature.geometry().GetPoint()

            if output_point[2] != nodata:
                output_points[point_index, :] = output_point

        return output_points
    else:
        if index < 1:
            index += 1

        raster_band = dataset.GetRasterBand(index)
        geotransform = dataset.GetGeoTransform()

        if nodata is None:
            nodata = raster_band.GetNoDataValue()

        return geoarray_to_points(raster_band.ReadAsArray(), origin=(geotransform[0], geotransform[3]),
                                  resolution=(geotransform[1], geotransform[5]), nodata=nodata)


def geoarray_to_points(array: numpy.array, origin: (float, float), resolution: (float, float), nodata: float = None) -> numpy.array:
    """
    Extract XYZ points from an array of data using the given raster-like georeference (origin  and resolution).

    Parameters
    ----------
    array
        2D array of gridded data
    origin
        X, Y coordinates of northwest corner
    resolution
        cell size
    nodata
        value to exclude from point creation from the input grid

    Returns
    -------
    numpy.array
        N x 3 array of XYZ points
    """

    x_values, y_values = numpy.meshgrid(numpy.linspace(origin[0], origin[0] + resolution[0] * array.shape[1], array.shape[1]),
                                        numpy.linspace(origin[1], origin[1] + resolution[1] * array.shape[0], array.shape[0]))

    if nodata is not None:
        x_values = x_values[array != nodata]
        y_values = y_values[array != nodata]
        array = array[array != nodata]

    return numpy.stack((x_values, y_values, array), axis=1)


def raster_edge_points(data: numpy.array, origin: (float, float), resolution: (float, float), nodata: float) -> numpy.array:
    """
    Get the edge points of the array of data using the given raster-like georeference (origin  and resolution).

    Parameters
    ----------
    data
        array of raster data
    origin
        geographic origin of data
    resolution
        geographic resolution of data
    nodata
        value for no data in raster

    Returns
    -------
    numpy.array
        N x 3 array of points
    """

    elevation_coverage = array_coverage(data, nodata)

    horizontal_difference = numpy.concatenate((numpy.full((elevation_coverage.shape[0], 1), 0),
                                               numpy.diff(numpy.where(elevation_coverage, 1, 0), axis=1)), axis=1)
    vertical_difference = numpy.concatenate((numpy.full((1, elevation_coverage.shape[1]), 0),
                                             numpy.diff(numpy.where(elevation_coverage, 1, 0), axis=0)), axis=0)

    horizontal_edges = (horizontal_difference == 1) | numpy.roll(horizontal_difference == -1, -1, axis=1)
    vertical_edges = (vertical_difference == 1) | numpy.roll(vertical_difference == -1, -1, axis=0)

    return geoarray_to_points(numpy.where(horizontal_edges | vertical_edges, data, nodata), origin, resolution, nodata)


def vertices(region: MultiPolygon):
    if type(region) is Polygon:
        return numpy.stack(region.exterior.xy, axis=1)
    else:
        points = []
        for polygon in region:
            points.append(numpy.stack(polygon.exterior.xy, axis=1))

        return numpy.concatenate(points, axis=0)


def alpha_hull(points: numpy.array, max_length: float = None) -> MultiPolygon:
    """
    Calculate the alpha shape (concave hull) of the given points.
    inspired by https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    max_length
        maximum length to include in the convex boundary

    Returns
    ----------
    MultiPolygon
        multipolygon of alpha hull
    """

    if points.shape[0] < 4:
        raise ValueError('need at least 4 points to perform triangulation')

    if max_length is None:
        max_length = numpy.max(cKDTree(points).query(points, k=2)[0][:, 1])

    triangles = Delaunay(points)

    indices = numpy.stack(((triangles.simplices[:, 0], triangles.simplices[:, 1]),
                           (triangles.simplices[:, 1], triangles.simplices[:, 2]),
                           (triangles.simplices[:, 2], triangles.simplices[:, 0])), axis=1).T
    edges = points[indices]
    vectors = numpy.squeeze(numpy.diff(edges, axis=2))
    lengths = numpy.hypot(vectors[:, :, 0], vectors[:, :, 1])
    indices = indices[lengths < max_length]

    if indices.shape[0] > 0:
        boundary_edge_indices = numpy.unique(numpy.sort(indices, axis=1), axis=0)
        return unary_union(list(polygonize(MultiLineString(points[boundary_edge_indices].tolist()))))
    else:
        print(f'no edges were found to be shorter than the specified length {max_length}; reverting to maximum nearest-neighbor distance')
        return alpha_hull(points)


def raster_mask(region: Polygon, shape: (float, float), resolution: (float, float), origin: (float, float), nodata: float,
                crs_wkt: str) -> gdal.Dataset:
    """
    Convert the given Shapely polygon into a GDAL raster dataset mask (with values outside of the polygon set to `nodata`).

    Parameters
    ----------
    region
        Shapely polygon with which to create mask
    shape
        shape of mask
    resolution
        cell size of mask
    origin
        northwest corner
    crs_wkt
        spatial reference of mask
    nodata
        value for no data in mask

    Returns
    -------
    gdal.Dataset
        GDAL raster dataset
    """

    if type(origin) is not numpy.array:
        origin = numpy.array(origin)

    with MemoryFile() as rasterio_memory_file:
        with rasterio_memory_file.open(driver='GTiff', width=shape[1], height=shape[0], count=1, crs=rasterio.crs.CRS.from_wkt(crs_wkt),
                                       transform=georeference_to_affine(origin, resolution), dtype=rasterio.float64,
                                       nodata=numpy.array([nodata]).astype(rasterio.float64).item()) as memory_raster:
            memory_raster.write(numpy.full(shape, 1, dtype=rasterio.float64), 1)

        with rasterio_memory_file.open() as memory_raster:
            masked_data, masked_transform = mask(memory_raster, [region])

    masked_data = masked_data[0, :, :]

    output_dataset = gdal.GetDriverByName('MEM').Create('', int(shape[1]), int(shape[0]), 1, gdal.GDT_Float32)
    output_dataset.SetProjection(crs_wkt)
    output_dataset.SetGeoTransform(masked_transform.to_gdal())
    output_band = output_dataset.GetRasterBand(1)
    output_band.SetNoDataValue(nodata)
    output_band.WriteArray(masked_data)

    return output_dataset


def raster_mask_like(region: Polygon, like_raster: gdal.Dataset, band: int = 1) -> gdal.Dataset:
    """
    Get a mask of the survey area with the same properties as the given raster.

    Parameters
    ----------
    region
        Shapely polygon with which to create mask
    like_raster
        raster to copy
    band
        raster band (1-indexed)

    Returns
    -------
    gdal.Dataset
        mask
    """

    geotransform = like_raster.GetGeoTransform()
    raster_band = like_raster.GetRasterBand(band)
    return raster_mask(region, shape=(like_raster.RasterYSize, like_raster.RasterXSize), resolution=(geotransform[1], geotransform[5]),
                       origin=(geotransform[0], geotransform[3]), nodata=raster_band.GetNoDataValue(), crs_wkt=gdal_crs_wkt(like_raster))


def apply_raster_mask(raster: gdal.Dataset, mask: gdal.Dataset, mask_value: float = None, band: int = 1) -> gdal.Dataset:
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
        raster band of mask (1-indexed)

    Returns
    -------
    gdal.Dataset
        GDAL raster dataset with mask applied
    """

    band = mask.GetRasterBand(band)
    mask_array = band.ReadAsArray()

    if mask_value is None:
        mask_value = band.GetNoDataValue()

    for band_index in range(1, raster.RasterCount + 1):
        raster_band = raster.GetRasterBand(band_index)
        raster_array = raster_band.ReadAsArray()
        raster_array[mask_array == mask_value] = raster_band.GetNoDataValue()
        raster_band.WriteArray(raster_array)

    return raster


def overwrite_raster(from_raster: gdal.Dataset, onto_raster: gdal.Dataset):
    assert (from_raster.RasterYSize, from_raster.RasterXSize) == \
           (onto_raster.RasterYSize, onto_raster.RasterXSize), 'rasters must be the same shape'

    for band_index in range(1, from_raster.RasterCount + 1):
        from_band = from_raster.GetRasterBand(band_index)
        onto_band = onto_raster.GetRasterBand(band_index)

        from_nodata = from_band.GetNoDataValue()
        from_values = from_band.ReadAsArray()

        onto_band.WriteArray(numpy.where(from_values != from_nodata, from_values, onto_band.ReadAsArray()))
        del from_band, onto_band

    return onto_raster


def shape_from_cell_size(resolution: (float, float), bounds: (float, float, float, float)) -> ((int, int), (float, float, float, float)):
    """
    Given the cell size and bounds of a raster, calculate the shape and bounds.
    This algorithm uses the minimum X and Y point as the anchor for rows and columns.
    Output bounds are floored to the nearest cell.

    Parameters
    ----------
    resolution
        cell size
    bounds
        min X, min Y, max X, max Y

    Returns
    -------
    (int, int), (float, float, float, float)
        shape and bounds (rounded to cell)
    """

    if type(resolution) is not numpy.array:
        resolution = numpy.array(resolution)

    if type(bounds) is not numpy.array:
        bounds = numpy.array(bounds)

    sw_corner = bounds[:2]
    ne_corner = bounds[2:]

    cell_remainder = (ne_corner - sw_corner) % resolution
    ne_corner -= cell_remainder

    shape = numpy.flip((ne_corner - sw_corner) / resolution).astype(int)

    return shape, numpy.concatenate((sw_corner, ne_corner))


def epsg_to_wkt(epsg: int) -> str:
    """
    Use OSR to get the well-known text of a CRS from its EPSG code.

    Parameters
    ----------
    epsg
        EPSG code of CRS

    Returns
    -------
    str
        well-known text of CRS
    """

    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromEPSG(epsg)
    return spatial_reference.ExportToWkt()


def gdal_crs_wkt(dataset: gdal.Dataset, layer: int = 0) -> str:
    """
    extract the well-known text of the CRS from the given GDAL dataset

    Parameters
    ----------
    dataset
        GDAL dataset (either raster or vector)
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
            crs_wkt = epsg_to_wkt(4326)
    elif match('^EPSG:[0-9]+$', crs_wkt):
        crs_wkt = epsg_to_wkt(int(crs_wkt[5:]))

    return crs_wkt


def array_coverage(array: numpy.array, nodata: float = None) -> numpy.array:
    """
    Return a boolean array of where data exists in the given array.

    Parameters
    ----------
    array
        array of gridded data with dimensions (Z)YX
    nodata
        value where there is no data in the given array

    Returns
    -------
    numpy.array
        array of booleans indicating where data exists
    """

    if len(array.shape) > 2:
        array = numpy.squeeze(array)

    if nodata is None:
        nodata = numpy.nan

    # TODO find reduced generalization of multiple bands
    if array.shape[0] == 3:
        coverage = (array[0, :, :] != nodata) | (array[1, :, :] != nodata) | (array[2, :, :] != nodata)
    else:
        coverage = array != nodata

    return coverage


def vectorize_geoarray(array: numpy.array, transform: Affine, nodata: float = None) -> MultiPolygon:
    """
    Vectorize the extent of the given georeferenced array where data exists.

    Parameters
    ----------

    array
        array of gridded data
    transform
        Affine transform of geoarray
    nodata
        value where there is no data

    Returns
    -------
    MultiPolygon
        Shapely polygon or multipolygon of coverage extent
    """

    raster_mask = array_coverage(array, nodata)
    return MultiPolygon(shapely_shape(shape[0]) for shape in rasterio_shapes(numpy.where(raster_mask, 1, 0), mask=raster_mask,
                                                                             transform=transform))


def vectorize_raster(raster: Union[gdal.Dataset, str], band: int = 1) -> MultiPolygon:
    """
    Vectorize the extent of the given raster where data exists.

    Parameters
    ----------
    raster
        GDAL raster dataset or filename of raster
    band
        raster band (1-indexed)

    Returns
    -------
    MultiPolygon
        Shapely multipolygon of coverage extent
    """

    if type(raster) is gdal.Dataset:
        raster_band = raster.GetRasterBand(band)
        raster_array = raster_band.ReadAsArray()
        del raster_band
        transform = Affine.from_gdal(*raster.GetGeoTransform())
    elif type(raster) is str:
        with rasterio.open(raster) as raster:
            raster_array = raster.read()
            transform = raster.transform
    else:
        raise ValueError(f'unsupported input type {type(raster)}')

    return vectorize_geoarray(raster_array, transform)


def write_geometry(filename: str, geometry: Union[Polygon, MultiPolygon], crs_wkt: str = None, name: str = None, layer: str = None):
    """
    Write the given Shapely geometry to the given vector file.

    Parameters
    ----------
    filename
        file path to vector file
    geometry
        Shapely geometry
    crs_wkt
        well-known text of CRS
    name
        value to write to the `name` attribute
    layer
        name of layer to write to
    """

    write_geojson(filename, mapping(geometry), crs_wkt, name, layer)


def write_geojson(filename: str, geojson: dict, crs_wkt: str = None, name: str = None, layer: str = None):
    """
    Write the given GeoJSON dictionary to the given vector file.

    Parameters
    ----------
    filename
        file path to vector file
    geojson
        dictionary with GeoJSON mappings
    crs_wkt
        well-known text of CRS
    name
        value to write to the `name` attribute
    layer
        name of layer to write to
    """

    if crs_wkt is None:
        crs_wkt = fiona.crs.to_string(fiona.crs.from_epsg(4326))

    if name is None:
        name = geojson['type']

    with fiona.open(filename, 'w', 'GPKG', schema={'geometry': geojson['type'], 'properties': {'name': 'str'}}, crs_wkt=crs_wkt,
                    layer=layer) as output_vector_file:
        output_vector_file.write({'geometry': geojson, 'properties': {'name': name}})


def georeference_to_affine(origin: (float, float), resolution: (float, float), rotation: (float, float) = None) -> Affine:
    """
    Calculate the affine transformation from the given geographic parameters.

    Parameters
    ----------
    origin
        XY coordinates of the northwest corner
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


def bounds_from_opposite_corners(corner_1: (float, float), corner_2: (float, float)) -> (float, float, float, float):
    """
    Get bounds from two opposite XY points.

    Parameters
    ----------
    corner_1
        XY point
    corner_2
        XY point

    Returns
    -------
    float, float, float, float
        min X, min Y, max X, max Y
    """

    return numpy.ravel(numpy.sort(numpy.stack((corner_1, corner_2), axis=0), axis=0))


def raster_bounds(raster: gdal.Dataset) -> (float, float, float, float):
    """
    Get the bounds (grouped by dimension) of the given unrotated raster.

    Parameters
    ----------
    raster
        GDAL raster dataset

    Returns
    -------
    float, float, float, float
        min X, min Y, max X, max Y
    """

    geotransform = raster.GetGeoTransform()
    origin = numpy.array((geotransform[0], geotransform[3]))
    resolution = numpy.array((geotransform[1], geotransform[5]))
    rotation = numpy.array((geotransform[2], geotransform[4]))
    shape = raster.RasterYSize, raster.RasterXSize

    if numpy.any(rotation != 0):
        raise NotImplementedError('rotated rasters not supported')

    return bounds_from_opposite_corners(origin, origin + numpy.flip(shape) * resolution)


def extent_from_bounds(bounds: (float, float, float, float)) -> (float, float, float, float):
    """
    Get the extent (grouped by min / max) of the given unrotated raster.

    Parameters
    ----------
    bounds
        min X, min Y, max X, max Y

    Returns
    -------
    float, float, float, float
        min X, max X, min Y, max Y
    """

    if type(bounds) is not numpy.array:
        bounds = numpy.array(bounds)

    return bounds[[0, 2, 1, 3]]


def plot_region(region: Union[Polygon, MultiPolygon], axis: pyplot.Axes = None, show: bool = False, **kwargs):
    """
    PLot the given region (a Shapely polygon or multipolygon).

    Parameters
    ----------
    region
        shapely polygon or multipolygon
    axis
        `pyplot` axis to plot to
    show
        whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    if type(region) is Polygon:
        axis.plot(*region.exterior.xy, **kwargs)
    else:
        if 'c' not in kwargs:
            try:
                color = next(axis._get_lines.color_cycle)
            except AttributeError:
                color = 'r'
            kwargs['c'] = color

        for geometry in region:
            axis.plot(*geometry.exterior.xy, **kwargs)

    if show:
        pyplot.show()


def plot_regions(regions: [Polygon], colors: [str] = None, axis: pyplot.Axes = None, show: bool = False, **kwargs):
    """
    PLot the given regions using the given colors.

    Parameters
    ----------
    regions
        list of shapely polygons or multipolygons
    colors
        colors to plot each region
    axis
        `pyplot` axis to plot to
    show
        whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    if colors is None:
        colors = [get_cmap('gist_rainbow')(color_index / len(regions)) for color_index in range(len(regions))]

    for color_index, geometry in enumerate(regions):
        color = colors[color_index]
        if type(geometry) is Polygon:
            axis.plot(*geometry.exterior.xy, c=color, **kwargs)
        else:
            for polygon in geometry:
                axis.plot(*polygon.exterior.xy, c=color, **kwargs)

    if show:
        pyplot.show()


def plot_bounding_box(sw: (float, float), ne: (float, float), axis: pyplot.Axes = None, show: bool = False, **kwargs):
    """
    Plot the bounding box of the given extent.

    Parameters
    ----------
    sw
        XY coordinates of southwest corner
    ne
        XY coordinates of northeast corner
    axis
        `pyplot` axis to plot to
    show
        whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    corner_points = numpy.array([sw, (ne[0], sw[1]), ne, (sw[0], ne[1]), sw])

    axis.plot(corner_points[:, 0], corner_points[:, 1], **kwargs)

    if show:
        pyplot.show()


def plot_raster(raster: gdal.Dataset, band: int = 1, axis: pyplot.Axes = None, show: bool = False, **kwargs):
    """
    Plot the given GDAL raster dataset using its georeference information.

    Parameters
    ----------
    raster
        GDAL raster dataset
    band
        raster band (1-indexed)
    axis
        `pyplot` axis to plot to
    show
        whether to show the plot
    """

    if axis is None:
        axis = pyplot.gca()

    raster_band = raster.GetRasterBand(band)
    raster_data = numpy.flip(raster_band.ReadAsArray(), axis=0)
    raster_data[raster_data == raster_band.GetNoDataValue()] = numpy.nan

    geotransform = raster.GetGeoTransform()
    if geotransform[1] < 0:
        raster_data = numpy.flip(raster_data, axis=1)
    if geotransform[5] < 0:
        raster_data = numpy.flip(raster_data, axis=0)

    axis.matshow(raster_data, extent=raster_bounds(raster)[[0, 2, 1, 3]], aspect='auto', **kwargs)

    if show:
        pyplot.show()
