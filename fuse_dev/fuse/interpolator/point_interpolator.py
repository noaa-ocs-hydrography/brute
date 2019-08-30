# -*- coding: utf-8 -*-
"""
point_interpolator.py

grice 20180625

V0.0.3 20190211

This is a collection of functions for creating interpolated bathymetry for the
national bathymetry.

The objective is to interpolate both XYZ data and BAGs.


Sources:
    Make ogr dataset from numpy array:
        https://pcjericks.github.io/py-gdalogr-cookbook/geometry.html
    ogr data set to gdal for gridding:
        http://osgeo-org.1560.x6.nabble.com/gdal-dev-DataSource-Dataset-using-
        gdal-Grid-td5322689.html
    GDAL gridding information:
        http://www.gdal.org/grid_tutorial.html#grid_tutorial_interpolation
    gdal_grid:
        http://www.gdal.org/gdal_grid.html

"""

# __version__ = 'point_interpolator 0.0.1'
import datetime
import re
from concurrent import futures
from typing import Union

import matplotlib.cm
import matplotlib.colors
import matplotlib.pyplot
import numpy
import rasterio.crs
import rasterio.io
import rasterio.mask
import rasterio.transform
import scipy
import scipy.interpolate
import scipy.spatial
import scipy.spatial.distance
import shapely.geometry
from osgeo import gdal, osr
from pykrige.ok import OrdinaryKriging
from scipy.spatial import cKDTree
from shapely.ops import unary_union, polygonize

gdal.UseExceptions()


class PointInterpolator:
    """Interpolation methods for creating a raster from points."""

    def __init__(self, window_scalar: float = 1.5, nodata: float = 1000000.0):
        """
        Set some of the precondition, but make them over writable.

        Parameters
        ----------
        window_scalar
            multiplier to use to expand interpolation radius from minimum point spacing
        nodata
            value to use to represent no data
        """

        self.window_scalar = window_scalar
        self.nodata = nodata

    def interpolate(self, points: gdal.Dataset, method: str, output_resolution: float, vector_file_path: str = None,
                    plot: bool = False, layer_index: int = 0) -> gdal.Dataset:
        """
        Interpolate the provided dataset.

        Currently this is assumed to be a gdal dataset.  At some point perhaps
        this should become something more native to python.

        Parameters
        ----------
        points
            gdal point cloud dataset
        method
            method for interpolation
        output_resolution
            resolution of output grid
        vector_file_path
            path to file containing vector data
        plot
            whether to plot the interpolated result next to the original survey points
        layer_index
            index of layer to read

        Returns
        -------
        gdal.Dataset
            raster dataset of interpolated data
        """

        if method == 'invdist_scilin' and vector_file_path is None:
            raise ValueError('Supporting shapefile required')

        self.input_points = gdal_points_to_array(points)
        points_2d = self.input_points[:, 0:2]

        self.input_bounds = _bounds_from_points(points_2d)
        self.output_shape, self.output_bounds = _shape_from_cell_size(output_resolution, self.input_bounds)
        self.crs_wkt = points.GetLayerByIndex(layer_index).GetSpatialRef().ExportToWkt()
        if self.crs_wkt == '':
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(4326)
            self.crs_wkt = srs.ExportToWkt()

        # set the size of the interpolation window to be the maximum of nearest neighbor distances
        self.window_size = self.window_scalar * numpy.max(cKDTree(points_2d).query(points_2d, k=2)[0][:, 1])

        start_time = datetime.datetime.now()

        if method == 'linear':
            # linear interpolation (`gdal.Grid`)
            interpolated_dataset = self.__interp_points_gdal_linear(points)
        elif method == 'invdist':
            # inverse distance interpolation (`gdal.Grid`)
            interpolated_dataset = self.__interp_points_gdal_invdist(points, self.window_size)
        elif method == 'invdist_scilin':
            # inverse distance interpolation (`gdal.Grid`) with linear interpolation (using `scipy.interpolate.griddata`)
            interpolated_dataset = self.__interp_points_gdal_invdist_scipy_linear(points, output_resolution)
        elif method == 'kriging':
            # interpolate using Ordinary Kriging (`pykrige`) by dividing data into smaller chunks
            chunk_size = (40, 40)
            method = f'{method}_{chunk_size[0]}x{chunk_size[1]}'
            interpolated_dataset = self.__interp_points_pykrige_kriging(self.input_points, output_resolution,
                                                                        chunk_size)
        else:
            raise ValueError(f'Interpolation type "{method}" not recognized.')

        # mask the interpolated data using the input points
        mask = self.__get_mask_like(interpolated_dataset)
        output_dataset = _mask_raster(interpolated_dataset, mask)

        print(f'{method} interpolation took {(datetime.datetime.now() - start_time).total_seconds()} s')

        if plot:
            self.__plot(output_dataset, method)

        return output_dataset

    def __get_mask(self, shape: (float, float), resolution: float, nw_corner: (float, float),
                   crs_wkt: str = None, nodata: float = None) -> gdal.Dataset:
        """
        Get a mask of the survey area with the given raster properties.

        Parameters
        ----------
        shape
            shape of mask
        resolution
            cell size of mask
        nw_corner
            northwest corner
        crs_wkt
            spatial reference of mask
        nodata
            value for no data in mask

        Returns
        -------
        gdal.Dataset
            mask
        """

        if type(nw_corner) is not numpy.array:
            nw_corner = numpy.array(nw_corner)

        if crs_wkt is None or crs_wkt == '':
            crs_wkt = self.crs_wkt

        if nodata is None:
            nodata = self.nodata

        concave_hull = _alpha_hull(self.input_points[:, 0:2], self.window_size)

        with rasterio.io.MemoryFile() as rasterio_memory_file:
            with rasterio_memory_file.open(driver='GTiff', width=shape[1], height=shape[0], count=1,
                                           crs=rasterio.crs.CRS.from_wkt(crs_wkt),
                                           transform=rasterio.transform.Affine.translation(*nw_corner) *
                                                     rasterio.transform.Affine.scale(resolution, resolution),
                                           dtype=rasterio.float64, nodata=numpy.array([nodata]).astype(
                        rasterio.float64).item()) as memory_raster:
                memory_raster.write(numpy.full(shape, 1, dtype=rasterio.float64), 1)

            with rasterio_memory_file.open() as memory_raster:
                masked_data, masked_transform = rasterio.mask.mask(memory_raster, [concave_hull])

        output_dataset = gdal.GetDriverByName('MEM').Create('', shape[1], shape[0], 1, gdal.GDT_Float32)
        output_dataset.SetProjection(crs_wkt)
        output_dataset.SetGeoTransform((masked_transform.c, masked_transform.a, masked_transform.b,
                                        masked_transform.f, masked_transform.d, masked_transform.e))
        output_band = output_dataset.GetRasterBand(1)
        output_band.SetNoDataValue(nodata)
        output_band.WriteArray(masked_data[0, :, :])

        return output_dataset

    def __get_mask_like(self, like_raster: gdal.Dataset, band_index: int = 1) -> gdal.Dataset:
        """
        Get a mask of the survey area with the same properties as the given raster.

        Parameters
        ----------
        like_raster
            raster to copy
        band_index
            raster band (1-indexed)

        Returns
        -------
        gdal.Dataset
            mask
        """

        shape = like_raster.RasterYSize, like_raster.RasterXSize
        geotransform = like_raster.GetGeoTransform()
        resolution = geotransform[1]
        nw_corner = geotransform[0], geotransform[3]
        crs_wkt = like_raster.GetProjectionRef()
        if re.match('^EPSG:[0-9]+$', crs_wkt):
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(crs_wkt[-4:]))
            crs_wkt = srs.ExportToWkt()
        nodata = like_raster.GetRasterBand(band_index).GetNoDataValue()
        return self.__get_mask(shape, resolution, nw_corner, crs_wkt, nodata)

    def __interp_points_gdal_linear(self, points: gdal.Dataset, nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating linearly.

        Parameters
        ----------
        points
            GDAL point cloud dataset
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid

        """

        if nodata is None:
            nodata = self.nodata

        algorithm = f'linear:radius=0:nodata={int(nodata)}'
        return gdal.Grid('', points, format='MEM', width=self.output_shape[1], height=self.output_shape[0],
                         outputBounds=self.output_bounds, algorithm=algorithm)

    def __interp_points_gdal_invdist(self, points: gdal.Dataset, radius: float = None,
                                     nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance.

        Parameters
        ----------
        points
            GDAL point cloud dataset
        radius
            size of interpolation window
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if radius is None:
            radius = self.window_size

        if nodata is None:
            nodata = self.nodata

        algorithm = f'invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:nodata={int(nodata)}'
        return gdal.Grid('', points, format='MEM', width=self.output_shape[1], height=self.output_shape[0],
                         outputBounds=self.output_bounds, algorithm=algorithm)

    def __interp_points_gdal_invdist_scipy_linear(self, points: gdal.Dataset, radius: float = None,
                                                  nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance and then linearly (using SciPy).

        Parameters
        ----------
        points
            GDAL point cloud dataset
        radius
            size of interpolation window
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if radius is None:
            radius = self.window_size

        if nodata is None:
            nodata = self.nodata

        # interpolate using inverse distance in GDAL
        interpolated_raster = self.__interp_points_gdal_invdist(points, radius, nodata)
        output_x, output_y = numpy.meshgrid(numpy.arange(interpolated_raster.RasterXSize),
                                            numpy.arange(interpolated_raster.RasterYSize))

        interpolated_data = interpolated_raster.ReadAsArray()
        y_values, x_values = numpy.where(interpolated_data != nodata)
        z_values = interpolated_data[y_values[:], x_values[:]]
        del interpolated_data

        # interpolate linearly in SciPy
        interpolated_data = scipy.interpolate.griddata((x_values, y_values), z_values, (output_x, output_y),
                                                       method='linear', fill_value=nodata)
        interpolated_data[numpy.isnan(interpolated_data)] = nodata

        band_1 = interpolated_raster.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_data)
        del band_1

        return interpolated_raster

    def __interp_points_pykrige_kriging(self, points: numpy.array, output_resolution: float,
                                        chunk_size: (float, float), nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via kriging (using PyKrige).

        Parameters
        ----------
        points
            N x 3 array of XYZ points
        output_resolution
            resolution of output grid
        chunk_size
            size of each chunk with which to divide the interpolation task, in units
            setting this value to less than 40x20 has adverse consequences on memory consumption
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if type(chunk_size) is not numpy.array:
            chunk_size = numpy.array(chunk_size)

        if nodata is None:
            nodata = self.nodata

        assert not numpy.any(chunk_size < numpy.array((40, 20))), 'a small chunk size may freeze your computer'

        chunk_shape = numpy.array(numpy.round(numpy.flip(chunk_size / output_resolution)), numpy.int)
        expand_indices = numpy.round(chunk_shape / 2).astype(numpy.int)

        input_points_sw = numpy.array((self.input_bounds[0], self.input_bounds[1]))

        interpolated_grid_x = numpy.linspace(self.output_bounds[0], self.output_bounds[2], self.output_shape[1])
        interpolated_grid_y = numpy.linspace(self.output_bounds[1], self.output_bounds[3], self.output_shape[0])
        interpolated_grid_shape = numpy.array((len(interpolated_grid_y), len(interpolated_grid_x)), numpy.int)
        interpolated_grid_values = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)
        interpolated_grid_variance = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)

        chunks = numpy.array(numpy.floor(interpolated_grid_shape / chunk_shape), numpy.int)

        chunk_grid_index = numpy.array((0, 0), numpy.int)
        with futures.ProcessPoolExecutor() as concurrency_pool:
            running_futures = {}

            for chunk_index in range(int(numpy.product(chunks))):
                # determine indices of the lower left and upper right bounds of the chunk within the output grid
                grid_start_index = numpy.array(chunk_grid_index * chunk_shape)
                grid_end_index = grid_start_index + chunk_shape
                grid_slice = tuple(slice(grid_start_index[axis], grid_end_index[axis])
                                   for axis in range(len(chunk_shape)))

                # expand the interpolation area past the edges of the chunk to minimize edge formation at the boundary
                interpolation_sw = input_points_sw + numpy.flip((grid_start_index - expand_indices) * output_resolution)
                interpolation_ne = input_points_sw + numpy.flip((grid_end_index + expand_indices) * output_resolution)
                interpolation_points = self.input_points[numpy.where(
                    (self.input_points[:, 0] >= interpolation_sw[0]) & (
                            self.input_points[:, 0] < interpolation_ne[0]) & (
                            self.input_points[:, 1] >= interpolation_sw[1]) & (
                            self.input_points[:, 1] < interpolation_ne[1]))[0]]

                if interpolation_points.shape[0] >= 3:
                    current_future = concurrency_pool.submit(_krige_onto_grid, interpolation_points,
                                                             interpolated_grid_x[grid_slice[1]],
                                                             interpolated_grid_y[grid_slice[0]])
                    running_futures[current_future] = grid_slice

                if chunk_grid_index[0] >= chunks[0]:
                    chunk_grid_index[0] = 0
                    chunk_grid_index[1] += 1
                else:
                    chunk_grid_index[0] += 1

            for _, completed_future in enumerate(futures.as_completed(running_futures)):
                grid_slice = running_futures[completed_future]
                try:
                    chunk_interpolated_values, chunk_interpolated_variance = completed_future.result()
                except Exception as error:
                    # I have absolutely no idea why this error happens, but wrapping it in `try ... except` seems to work without negative consequences
                    print(error)
                interpolated_grid_values[grid_slice] = chunk_interpolated_values
                interpolated_grid_variance[grid_slice] = chunk_interpolated_variance

        interpolated_grid_uncertainty = numpy.sqrt(interpolated_grid_variance) * 2.5
        del interpolated_grid_variance

        output_dataset = gdal.GetDriverByName('MEM').Create('temp', self.output_shape[1], self.output_shape[0], 2,
                                                            gdal.GDT_Float32)

        output_dataset.SetGeoTransform((self.output_bounds[0], output_resolution, 0.0,
                                        self.output_bounds[1], 0.0, output_resolution))
        output_dataset.SetProjection(self.crs_wkt)

        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_grid_values)

        band_2 = output_dataset.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(interpolated_grid_uncertainty)

        del band_1, band_2
        return output_dataset

    def __plot(self, interpolated_raster: Union[numpy.array, gdal.Dataset], interpolation_method: str,
               nodata: float = None, band_index: int = 1):
        """
        Plot preinterpolated points and an interpolated raster side-by-side on synchronized subplots.

        Parameters
        ----------
        interpolated_raster
            array of interpolated data
        interpolation_method
            method of interpolation
        band_index
            raster band (1-indexed)
        """

        if nodata is None:
            nodata = self.nodata

        # if the input raster is a GDAL raster dataset, extract the data from the first band
        if type(interpolated_raster) is gdal.Dataset:
            interpolated_raster = numpy.flip(interpolated_raster.GetRasterBand(band_index).ReadAsArray(), axis=0)

        # replace `nodata` values with NaN
        interpolated_raster[interpolated_raster == nodata] = numpy.nan

        # get minimum and maximum values for all three dimensions
        min_x, min_y, max_x, max_y = self.input_bounds
        z_values = numpy.concatenate((numpy.ravel(interpolated_raster), self.input_points[:, 2]))
        min_z = numpy.nanmin(z_values)
        max_z = numpy.nanmax(z_values)

        # create a new Matplotlib figure window
        figure = matplotlib.pyplot.figure()

        # create subplots
        survey_points_axis = figure.add_subplot(1, 2, 1)
        survey_points_axis.set_title('survey points')
        survey_points_axis.set_xlim(min_x, max_x)
        survey_points_axis.set_ylim(min_y, max_y)
        interpolated_raster_axis = figure.add_subplot(1, 2, 2, sharex=survey_points_axis, sharey=survey_points_axis)
        interpolated_raster_axis.set_title(f'{interpolation_method} interpolation to raster')

        # plot data
        survey_points_axis.scatter(self.input_points[:, 0], self.input_points[:, 1], c=self.input_points[:, 2], s=1,
                                   vmin=min_z, vmax=max_z)

        interpolated_raster_axis.imshow(interpolated_raster, extent=(min_x, max_x, min_y, max_y), aspect='auto',
                                        vmin=min_z, vmax=max_z)

        # create colorbar
        figure.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min_z, vmax=max_z)),
                        ax=(interpolated_raster_axis, survey_points_axis))

        # pause program execution and show the figure
        matplotlib.pyplot.show()

    def __plot_concave_hull(self):
        concave_hull = _alpha_hull(self.input_points[:, 0:2], self.window_size)
        triangles = scipy.spatial.Delaunay(self.input_points[:, 0:2])
        matplotlib.pyplot.plot(*concave_hull.exterior.xy, c='r')
        matplotlib.pyplot.triplot(self.input_points[:, 0], self.input_points[:, 1], triangles=triangles.simplices,
                                  c='g')
        matplotlib.pyplot.scatter(self.input_points[:, 0], self.input_points[:, 1], c='b')
        matplotlib.pyplot.show()


def _mask_raster(raster: gdal.Dataset, mask: gdal.Dataset, mask_band_index: int = 1,
                 mask_value: float = None) -> gdal.Dataset:
    """
    Mask a raster using the nodata values in the given mask.
    Both rasters are assumed to be collocated, as all operations are conducted on the pixel level.

    Parameters
    ----------
    raster
        raster to apply mask to
    mask
        mask to apply to the raster
    mask_value
        value to use as mask
    mask_band_index
        raster band (1-indexed)

    Returns
    -------
    gdal.Dataset
        masked raster dataset
    """

    mask_band = mask.GetRasterBand(mask_band_index)
    mask_array = mask_band.ReadAsArray()

    if mask_value is None:
        mask_value = mask_band.GetNoDataValue()

    for band_index in range(1, raster.RasterCount + 1):
        raster_band = raster.GetRasterBand(band_index)
        raster_array = raster_band.ReadAsArray()
        raster_array[mask_array == mask_value] = raster_band.GetNoDataValue()
        raster_band.WriteArray(raster_array)

    return raster


def _bounds_from_points(points: numpy.array) -> (float, float, float, float):
    """
    Calculate the X and Y bounds of

    Parameters
    ----------
    points
        array of XY or XYZ points

    Returns
    -------
    tuple
        min X, min Y, max X, max Y
    """

    return numpy.min(points[:, 0]), numpy.min(points[:, 1]), numpy.max(points[:, 0]), numpy.max(points[:, 1])


def _shape_from_cell_size(resolution: float, bounds: (float, float, float, float)) -> (
        (int, int), (float, float, float, float)):
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
    tuple, tuple
        shape, bounds
    """

    min_x, min_y, max_x, max_y = bounds
    cols = int(numpy.ceil(float((max_x - min_x) / resolution)))
    cell_remainder_x = min_x % resolution
    rounded_min_x = min_x - cell_remainder_x
    rounded_max_x = rounded_min_x + float(cols * resolution)
    rows = int(numpy.ceil(float((max_y - min_y) / resolution)))
    cell_remainder_y = min_y % resolution
    rounded_min_y = min_y - cell_remainder_y
    rounded_max_y = rounded_min_y + float(rows * resolution)
    return (rows, cols), (rounded_min_x, rounded_min_y, rounded_max_x, rounded_max_y)


def gdal_points_to_array(points: gdal.Dataset, layer_index: int = 0) -> numpy.array:
    """
    Take a gdal vector xyz point cloud and return a numpy array.

    Parameters
    ----------
    dataset
        a GDAL point cloud dataset
    layer_index
        index of layer to read

    Returns
    -------
    numpy array
        XYZ points
    """

    point_layer = points.GetLayerByIndex(layer_index)
    num_points = point_layer.GetFeatureCount()
    output_points = numpy.empty((num_points, 3))

    for point_index in range(num_points):
        feature = point_layer.GetFeature(point_index)
        output_points[point_index, :] = feature.geometry().GetPoint()

    return output_points


def _krige_onto_grid(input_points: numpy.array, output_x: numpy.array, output_y: numpy.array) -> (
        numpy.array, numpy.array):
    """
    Perform kriging interpolation from the given set of points onto a regular grid.

    Parameters
    ----------
    input_points
        N x 3 array of XYZ points
    output_x
        X values of output grid
    output_y
        Y values of output grid
    crop
        list of slices to which to crop output

    Returns
    -------
    interpolated values, variance
    """

    interpolator = OrdinaryKriging(input_points[:, 0], input_points[:, 1], input_points[:, 2], variogram_model='linear',
                                   verbose=False, enable_plotting=False)
    return interpolator.execute('grid', output_x, output_y)


def _boundary_edges(points: numpy.array, radius: float = 1.0):
    """
    Get the edges of the Delaunay triangulation that exist on the boundary.

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    radius
        radius of triangulation
    """

    if points.shape[0] < 4:
        raise ValueError('need at least 4 points to perform triangulation')

    def add_edge(edges, indices):
        """ add an edge between the indices, if it is not already in the edge list """
        if indices not in edges and reversed(indices) not in edges:
            edges.add(indices)
        else:
            # remove the edge if it already exists
            edges.remove(indices)

    triangles = scipy.spatial.Delaunay(points)
    edges = set()
    for index_a, index_b, index_c in triangles.vertices:
        point_a = points[index_a]
        point_b = points[index_b]
        point_c = points[index_c]

        length_a = numpy.hypot(*(point_a - point_b))
        length_b = numpy.hypot(*(point_b - point_c))
        length_c = numpy.hypot(*(point_c - point_a))

        semiperimeter = (length_a + length_b + length_c) / 2.0
        area = numpy.sqrt(
            semiperimeter * (semiperimeter - length_a) * (semiperimeter - length_b) * (semiperimeter - length_c))
        circumcircle_radius = length_a * length_b * length_c / (4.0 * area)

        if circumcircle_radius < radius:
            add_edge(edges, (index_a, index_b))
            add_edge(edges, (index_b, index_c))
            add_edge(edges, (index_c, index_a))

    return points[numpy.array(tuple(edges))].tolist()


def _alpha_hull(points: numpy.array, max_length: float = None):
    """
    Calculate the alpha shape (concave hull) of the given points.
    inspired by https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    max_length
        maximum edge length to include in the convex boundary
    """

    if points.shape[0] < 4:
        raise ValueError('need at least 4 points to perform triangulation')

    if max_length is None:
        max_length = numpy.max(cKDTree(points).query(points, k=2)[0][:, 1])

    triangles = scipy.spatial.Delaunay(points)
    boundary_edges = set()
    for index_a, index_b, index_c in triangles.simplices:
        point_a = points[index_a]
        point_b = points[index_b]
        point_c = points[index_c]

        length_a = numpy.hypot(*(point_a - point_b))
        length_b = numpy.hypot(*(point_b - point_c))
        length_c = numpy.hypot(*(point_c - point_a))

        semiperimeter = (length_a + length_b + length_c) / 2
        area = numpy.sqrt(
            semiperimeter * (semiperimeter - length_a) * (semiperimeter - length_b) * (semiperimeter - length_c))
        circumcircle_diameter = length_a * length_b * length_c / (area * 2)

        if circumcircle_diameter < max_length:
            for indices in ((index_a, index_b), (index_b, index_c), (index_c, index_a)):
                if indices not in boundary_edges and reversed(indices) not in boundary_edges:
                    boundary_edges.add(indices)
                else:
                    # remove the edge if it already exists (if it is shared by another triangle)
                    boundary_edges.remove(indices)

    boundary_edges = points[numpy.array(tuple(boundary_edges))].tolist()
    return unary_union(list(polygonize(shapely.geometry.MultiLineString(boundary_edges))))
