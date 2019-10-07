from re import match
from typing import Union

import numpy
import rasterio
import rasterio.crs
from affine import Affine
from osgeo import gdal, osr, ogr
from rasterio.features import shapes as rasterio_shapes, geometry_mask
from scipy.spatial import Delaunay, cKDTree
from shapely.geometry import Polygon, MultiLineString, MultiPolygon, shape as shapely_shape, mapping
from shapely.ops import polygonize
from shapely.wkt import loads as load_wkt


def gdal_to_xyz(dataset: gdal.Dataset, index: int = 0, nodata: float = None) -> numpy.array:
    """
    Extract XYZ points from a GDAL dataset.

    Parameters
    ----------
    dataset
        GDAL dataset (point cloud or raster)
    index
        zero-based index of layer / band to read from the given dataset
    nodata
        value to exclude from point creation

    Returns
    -------
    numpy.array
        N x 3 array of XYZ points
    """

    is_raster = dataset.GetProjectionRef() != ''

    if index < 0:
        index += dataset.RasterCount if is_raster else dataset.GetLayerCount()

    if is_raster:
        raster_band = dataset.GetRasterBand(index + 1)
        geotransform = dataset.GetGeoTransform()

        if nodata is None:
            nodata = raster_band.GetNoDataValue()

        return geoarray_to_xyz(raster_band.ReadAsArray(), origin=(geotransform[0], geotransform[3]),
                               resolution=(geotransform[1], geotransform[5]), nodata=nodata)
    else:
        point_layer = dataset.GetLayerByIndex(index)
        num_points = point_layer.GetFeatureCount()

        output_points = numpy.empty((num_points, 3))
        for point_index in range(num_points):
            feature = point_layer.GetFeature(point_index)
            output_point = feature.geometry().GetPoint()

            if output_point[2] != nodata:
                output_points[point_index, :] = output_point

        return output_points


def geoarray_to_xyz(data: numpy.array, origin: (float, float), resolution: (float, float), nodata: float = None) -> numpy.array:
    """
    Extract XYZ points from an array of data using the given raster-like georeference (origin  and resolution).

    Parameters
    ----------
    data
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

    if nodata is None:
        nodata = numpy.nan

    x_values, y_values = numpy.meshgrid(numpy.linspace(origin[0], origin[0] + resolution[0] * data.shape[1], data.shape[1]),
                                        numpy.linspace(origin[1], origin[1] + resolution[1] * data.shape[0], data.shape[0]))

    return numpy.stack((x_values[data != nodata], y_values[data != nodata], data[data != nodata]), axis=1)


def xyz_to_gdal(points: numpy.array, spatial_reference: osr.SpatialReference) -> gdal.Dataset:
    """
    Get a GDAL point cloud dataset from XYZ points.

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    utm_zone
        desired UTM zone
    vertical_datum
        name of desired vertical datum

    Returns
    -------
    gdal.Dataset
        GDAL point cloud dataset
    """

    dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
    layer = dataset.CreateLayer('pts', spatial_reference, geom_type=ogr.wkbPoint)

    for point in points:
        geometry = ogr.Geometry(ogr.wkbPoint)
        geometry.AddPoint(*point)
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(geometry)
        layer.CreateFeature(feature)

    return dataset


def raster_edge(data: numpy.array, nodata: float) -> numpy.array:
    """
    Get the cells of the array bordering `nodata`.

    Parameters
    ----------
    data
        array of raster data
    nodata
        value for no data in raster

    Returns
    -------
    numpy.array
        boolean array of edge cells
    """

    elevation_coverage = array_coverage(data, nodata)

    horizontal_difference = numpy.concatenate((numpy.full((elevation_coverage.shape[0], 1), 0),
                                               numpy.diff(numpy.where(elevation_coverage, 1, 0), axis=1)), axis=1)
    vertical_difference = numpy.concatenate((numpy.full((1, elevation_coverage.shape[1]), 0),
                                             numpy.diff(numpy.where(elevation_coverage, 1, 0), axis=0)), axis=0)

    horizontal_edges = (horizontal_difference == 1) | numpy.roll(horizontal_difference == -1, -1, axis=1)
    vertical_edges = (vertical_difference == 1) | numpy.roll(vertical_difference == -1, -1, axis=0)

    return horizontal_edges | vertical_edges


def raster_edge_points(data: numpy.array, origin: (float, float), resolution: (float, float), nodata: float) -> numpy.array:
    """
    Get the points bordering `nodata` in the given array, using the given raster transform.

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

    return geoarray_to_xyz(numpy.where(raster_edge(data, nodata), data, nodata), origin, resolution, nodata)


def shrink_array(data: numpy.array, nodata: float, cells: int = 1) -> numpy.array:
    """
    Shrink the edges of the given array by setting all values bordering `nodata` to `nodata`.

    Parameters
    ----------
    data
        array of data
    nodata
        value for no data
    cells
        cells by which to shrink the data

    Returns
    -------
    numpy.array
        shrunken array of data
    """

    for _ in range(cells):
        data = numpy.where(raster_edge(data, nodata), nodata, data)

    return data


def shrink_gdal_raster(raster: gdal.Dataset, cells: int = 1, band: int = 0) -> gdal.Dataset:
    """
    Shrink the edges of the given GDAL raster by erasing all values bordering `nodata`.

    Parameters
    ----------
    raster
        GDAL raster
    cells
        cells by which to shrink the raster
    band
        zero-based index of raster band

    Returns
    -------
    gdal.Dataset
        shrunken raster
    """

    raster_band = raster.GetRasterBand(band + 1)
    raster_band.WriteArray(shrink_array(raster_band.ReadAsArray(), raster_band.GetNoDataValue(), cells))
    del raster_band

    return raster


def polygon_vertices(polygon: Union[Polygon, MultiPolygon]) -> numpy.array:
    """
    Get the vertices of the given shape.

    Parameters
    ----------
    polygon
        Shapely polygon (or multipolygon)

    Returns
    -------
    numpy.array
        N x 2 array of XY points corresponding to the vertices of the given polygon
    """

    points = []
    if type(polygon) is Polygon:
        points.append(numpy.stack(polygon.exterior.xy, axis=1))
        for interior in polygon.interiors:
            points.append(numpy.stack(interior.xy, axis=1))
    else:
        for geometry in polygon:
            points.append(polygon_vertices(geometry))

    return numpy.concatenate(points, axis=0)


def alpha_hull(points: numpy.array, max_length: float = None, triangles: Delaunay = None) -> MultiPolygon:
    """
    Calculate the alpha shape (concave hull) of the given points.
    inspired by https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    max_length
        maximum length to include in the convex boundary
    triangles
        Delaunay triangulation of points (optional)

    Returns
    ----------
    MultiPolygon
        alpha shape (concave hull) of points
    """

    if points.shape[0] < 4:
        raise ValueError('need at least 4 points to perform triangulation')

    if max_length is None:
        max_length = maximum_nearest_neighbor_distance(points)

    if triangles is None:
        triangles = Delaunay(points)

    indices = numpy.stack(((triangles.simplices[:, 0], triangles.simplices[:, 1]),
                           (triangles.simplices[:, 1], triangles.simplices[:, 2]),
                           (triangles.simplices[:, 2], triangles.simplices[:, 0])), axis=1).T
    edges = points[indices]
    vectors = numpy.squeeze(numpy.diff(edges, axis=2))
    lengths = numpy.hypot(vectors[:, :, 0], vectors[:, :, 1])
    indices = indices[numpy.all(lengths <= max_length, axis=1)]

    if indices.shape[0] > 0:
        indices = numpy.sort(numpy.reshape(indices, (indices.shape[0] * indices.shape[1], indices.shape[2])), axis=1)
        boundary_edge_indices, counts = numpy.unique(indices, axis=0, return_counts=True)
        boundary_edge_indices = boundary_edge_indices[counts == 1]
        return consolidate_disparate_polygons(polygonize(MultiLineString(points[boundary_edge_indices].tolist())))
    else:
        print(f'no edges were found to be shorter than the specified length {max_length}; attempting with {max_length * 2}')
        return alpha_hull(points, max_length * 2, triangles)


def consolidate_disparate_polygons(polygons: [Union[Polygon, MultiPolygon]]) -> MultiPolygon:
    """
    Create a MultiPolygon from the given polygons, assuming interior polygons are holes and that no edges intersect.

    Parameters
    ----------
    polygons
        list of Shapely polygons

    Returns
    -------
    MultiPolygon
        Shapely multipolygon of all polygons with holes included
    """

    if type(polygons) is not list:
        polygons = list(polygons)

    outer_rings = []
    inner_rings = []

    for polygon_1 in polygons:
        for polygon_2 in polygons:
            if polygon_1 is not polygon_2 and polygon_1.intersects(polygon_2):
                if polygon_1.area > polygon_2.area:
                    outer_polygon = polygon_1
                    inner_polygon = polygon_2
                else:
                    outer_polygon = polygon_2
                    inner_polygon = polygon_1

                if outer_polygon not in outer_rings:
                    outer_rings.append(outer_polygon)
                    inner_rings.append([inner_polygon])
                else:
                    for index, outer_ring in enumerate(outer_rings):
                        if outer_ring is outer_polygon:
                            if inner_polygon not in inner_rings[index]:
                                inner_rings[index].append(inner_polygon)

        if polygon_1 not in outer_rings and polygon_1 not in [inner_ring for outer_ring in inner_rings for inner_ring in outer_ring]:
            outer_rings.append(polygon_1)
            inner_rings.append([])

    return MultiPolygon(Polygon(outer_ring.exterior.coords, [inner_ring.exterior.coords for inner_ring in inner_rings[index]])
                        for index, outer_ring in enumerate(outer_rings))


def rasterize_polygon(polygon: Union[Polygon, MultiPolygon], shape: (float, float), origin: (float, float),
                      resolution: (float, float)) -> (numpy.array, Affine):
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
        northwest corner

    Returns
    -------
    numpy.array, Affine
        boolean array mask (True is within the polygon, False outside) and its affine transform
    """

    if type(polygon) is Polygon:
        polygon = MultiPolygon([polygon])

    transform = affine_from_georeference(origin, resolution)
    return ~geometry_mask(polygon, shape, transform), transform


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

    mask, transform = rasterize_polygon(polygon, shape, origin, resolution)

    output_dataset = gdal.GetDriverByName('MEM').Create('', int(shape[1]), int(shape[0]), 1, gdal.GDT_Float32)
    output_dataset.SetProjection(crs_wkt)
    output_dataset.SetGeoTransform(transform.to_gdal())
    output_band = output_dataset.GetRasterBand(1)
    output_band.SetNoDataValue(nodata)
    output_band.WriteArray(numpy.where(mask, 1, nodata))

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


def overwrite_raster(from_raster: gdal.Dataset, onto_raster: gdal.Dataset) -> gdal.Dataset:
    """
    Overwrite the "onto" raster with data values (non-nodata) from the "from" raster.

    Parameters
    ----------
    from_raster
        raster to apply values from
    onto_raster
        raster to write onto

    Returns
    -------
    gdal.Dataset
        the "onto" raster after being overwritten with values from the "from" raster
    """

    if (from_raster.RasterYSize, from_raster.RasterXSize) != (onto_raster.RasterYSize, onto_raster.RasterXSize):
        raise NotImplementedError('overwriting rasters of different shapes is not supported')

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

    shape = numpy.flip((ne_corner - sw_corner) / numpy.abs(resolution)).astype(int)

    return shape, numpy.concatenate((sw_corner, ne_corner))


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
            crs_wkt = crs_wkt_from_epsg(4326)
    elif match('^EPSG:[0-9]+$', crs_wkt):
        crs_wkt = crs_wkt_from_epsg(int(crs_wkt[5:]))

    return crs_wkt


def gdal_crs_from_utm(utm_zone: int, vertical_datum: str = 'NADV88', horizontal_datum: str = 'NAD83') -> osr.SpatialReference:
    """
    Create an OSR NAD83 spatial reference from the given parameters.

    Parameters
    ----------
    utm_zone
        UTM zone
    vertical_datum
        vertical datum
    horizontal_datum
        horizontal reference frame

    Returns
    -------
    osr.SpatialReference
        OSR spatial reference object
    """

    spatial_reference = osr.SpatialReference()
    spatial_reference.SetWellKnownGeogCS(horizontal_datum)
    spatial_reference.SetUTM(abs(utm_zone), 1 if utm_zone > 0 else 0)
    spatial_reference.SetVertCS(vertical_datum, vertical_datum, 2000)
    return spatial_reference


def array_coverage(array: numpy.array, nodata: float = None) -> numpy.array:
    """
    Get a boolean array of where data exists in the given array.

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

    # TODO find reduced generalization of band coverage
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


def vectorize_raster(raster: Union[gdal.Dataset, str], band: int = 0) -> MultiPolygon:
    """
    Vectorize the extent of the given raster where data exists.

    Parameters
    ----------
    raster
        GDAL raster dataset or filename of raster
    band
        zero-based index of raster band

    Returns
    -------
    MultiPolygon
        Shapely multipolygon of coverage extent
    """

    if type(raster) is gdal.Dataset:
        raster_band = raster.GetRasterBand(band + 1)
        raster_array = raster_band.ReadAsArray()
        del raster_band
        transform = Affine.from_gdal(*raster.GetGeoTransform())
    elif type(raster) is str:
        with rasterio.open(raster) as raster:
            raster_array = raster.read()
            transform = raster.transform
    else:
        raise NotImplementedError(f'unsupported input type {type(raster)}')

    return vectorize_geoarray(raster_array, transform)


def affine_from_georeference(origin: (float, float), resolution: (float, float), rotation: (float, float) = None) -> Affine:
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


def maximum_nearest_neighbor_distance(points: numpy.array) -> float:
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

    return numpy.max(cKDTree(points).query(points, k=2)[0][:, 1])


def bounds_from_opposite_corners(corner_1: (float, float), corner_2: (float, float)) -> (float, float, float, float):
    """
    Get bounds from two XY points.

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


def gdal_raster_bounds(raster: gdal.Dataset) -> (float, float, float, float):
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


def gdal_dataset_is_raster(dataset: gdal.Dataset) -> bool:
    """
    Determine whether the given dataset is a raster.

    Parameters
    ----------
    dataset
        GDAL dataset

    Returns
    -------
    bool
        whether the given dataset is a raster
    """

    return dataset.GetProjectionRef() != ''


def geojson_from_geometry_wkt(geometry_wkt: str) -> dict:
    """
    Get GeoJSON of a geometry from its well-known text.

    Parameters
    ----------
    geometry_wkt
        well-known text of geometry

    Returns
    -------
    dict
        GeoJSON of geometry
    """

    return mapping(load_wkt(geometry_wkt))
