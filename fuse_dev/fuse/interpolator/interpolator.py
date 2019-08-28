# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

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
from osgeo import gdal, osr, ogr
from pykrige.ok import OrdinaryKriging
from scipy.spatial import cKDTree
from shapely.ops import unary_union, polygonize


class PointInterpolator:
    """An abstraction for data interpolation."""

    def __init__(self, dataset: gdal.Dataset, default_nodata: float = 1000000.0, window_scalar: float = 1,
                 layer_index: int = 0):
        """
        Create a new point interpolator object for the given GDAL point cloud dataset.

        Parameters
        ----------
        dataset
            GDAL point cloud dataset
        default_nodata
            value to set for no data in the output interpolation
        window_scalar
            multiplier to use to expand interpolation radius from minimum point spacing
        layer_index
            index of layer to read in given dataset
        """

        self.dataset = dataset
        self.default_nodata = default_nodata

        self.points = gdal_points_to_array(self.dataset)
        points_2d = self.points[:, 0:2]

        # set window size as the minimum distance to all nearest neighbors
        self.window_size = window_scalar * numpy.min(cKDTree(points_2d).query(points_2d, k=2)[0][:, 1])

        # get the concave hull of the survey points
        self.concave_hull = _alpha_hull(points_2d, self.window_size * 3)

        # get the bounds (min X, min Y, max X, max Y)
        self.input_bounds = numpy.min(points_2d[:, 0]), numpy.min(points_2d[:, 1]), numpy.max(
            points_2d[:, 0]), numpy.max(points_2d[:, 1])

        # get the well-known text of the CRS
        self.crs_wkt = _get_gdal_crs_wkt(self.dataset, layer_index)

    def interpolate(self, method: str, output_resolution: float, output_nodata: float = None,
                    plot: bool = True) -> gdal.Dataset:
        """
        Take a gdal dataset and run the interpolation, returning a gdal raster.

        Parameters
        ----------
        method
            interpolation method
        metadata
            dictionary of metadata

        Returns
        -------
            interpolated GDAL dataset and a dictionary of metadata
        """

        output_shape, output_bounds = _shape_from_cell_size(output_resolution, self.input_bounds)

        if output_nodata is None:
            output_nodata = self.default_nodata

        start_time = datetime.datetime.now()

        if method == 'linear':
            # linear interpolation (`gdal.Grid`)
            interpolated_dataset = self.__linear_gdal(output_shape, output_bounds, output_nodata)
        elif method == 'invdist':
            # inverse distance interpolation (`gdal.Grid`)
            interpolated_dataset = self.__invdist_gdal(output_shape, output_bounds, output_nodata)
        elif method == 'invdist_linear':
            # inverse distance interpolation (`gdal.Grid`) with linear interpolation (using `scipy.interpolate.griddata`)
            interpolated_dataset = self.__invdist_gdal_linear_scipy(output_shape, output_bounds, output_nodata)
        elif method == 'kriging':
            # interpolate using Ordinary Kriging (`pykrige`) by dividing data into smaller chunks
            chunk_size = (40, 40)
            method = f'{method}_{chunk_size}m'
            interpolated_dataset = self.__kriging_pykrige(output_shape, output_bounds, chunk_size, output_nodata)
        else:
            raise ValueError(f'Interpolation type "{method}" not recognized.')

        # mask the interpolated data using the input points
        output_dataset = _mask_raster(interpolated_dataset, self.__get_mask_like(interpolated_dataset))

        print(f'{method} interpolation took {(datetime.datetime.now() - start_time).total_seconds()} s')

        if plot:
            self.__plot(output_dataset, method)

        return output_dataset

    def __get_mask(self, shape: (float, float), resolution: float, nw_corner: (float, float), nodata: float,
                   crs_wkt: str = None) -> gdal.Dataset:
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
            nodata = self.default_nodata

        with rasterio.io.MemoryFile() as rasterio_memory_file:
            with rasterio_memory_file.open(driver='GTiff', width=shape[1], height=shape[0], count=1,
                                           crs=rasterio.crs.CRS.from_wkt(crs_wkt),
                                           transform=rasterio.transform.Affine.translation(*nw_corner) *
                                                     rasterio.transform.Affine.scale(resolution, resolution),
                                           dtype=rasterio.float64, nodata=numpy.array([nodata]).astype(
                        rasterio.float64).item()) as memory_raster:
                memory_raster.write(numpy.full(shape, 1, dtype=rasterio.float64), 1)

            with rasterio_memory_file.open() as memory_raster:
                masked_data, masked_transform = rasterio.mask.mask(memory_raster, [self.concave_hull])

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
        nodata = like_raster.GetRasterBand(band_index).GetNoDataValue()
        crs_wkt = like_raster.GetProjectionRef()
        if re.match('^EPSG:[0-9]+$', crs_wkt):
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(int(crs_wkt[-4:]))
            crs_wkt = srs.ExportToWkt()
        return self.__get_mask(shape, resolution, nw_corner, nodata, crs_wkt)

    def __linear_gdal(self, output_shape: (int, int), output_bounds: (float, float, float, float),
                      output_nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating linearly.

        Parameters
        ----------
        output_shape
            shape (rows, cols) of output grid
        output_bounds
            bounds (min X, min Y, max X, max Y) of output grid
        output_nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid

        """

        if output_nodata is None:
            output_nodata = self.default_nodata

        return gdal.Grid('', self.dataset, format='MEM', width=output_shape[1], height=output_shape[0],
                         outputBounds=output_bounds, algorithm=f'linear:radius=0:nodata={int(output_nodata)}')

    def __linear_scipy(self, output_shape: (int, int), output_bounds: (float, float, float, float),
                       output_nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating linearly.

        Parameters
        ----------
        output_shape
            shape (rows, cols) of output grid
        output_bounds
            bounds (min X, min Y, max X, max Y) of output grid
        output_nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid

        """

        if type(output_shape) is not numpy.array:
            output_shape = numpy.array(output_shape)

        if output_nodata is None:
            output_nodata = self.default_nodata

        output_x, output_y = numpy.meshgrid(numpy.arange(output_shape[1]), numpy.arange(output_shape[0]))

        # interpolate linearly in SciPy
        interpolated_data = scipy.interpolate.griddata((self.points[:, 0], self.points[:, 1]), self.points[:, 2],
                                                       (output_x, output_y), method='linear', fill_value=output_nodata)
        interpolated_data[numpy.isnan(interpolated_data)] = output_nodata

        output_resolution = (numpy.array(output_bounds[2:4]) - numpy.array(output_bounds[0:2])) / output_shape

        output_dataset = gdal.GetDriverByName('MEM').Create('temp', output_shape[1], output_shape[0], 2,
                                                            gdal.GDT_Float32)
        output_dataset.SetGeoTransform((output_bounds[0], output_resolution, 0.0,
                                        output_bounds[1], 0.0, output_resolution))
        output_dataset.SetProjection(self.crs_wkt)
        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(output_nodata)
        band_1.WriteArray(interpolated_data)
        del band_1
        return output_dataset

    def __invdist_gdal(self, output_shape: (int, int), output_bounds: (float, float), output_nodata: float = None,
                       radius: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance.

        Parameters
        ----------
        output_shape
            shape (rows, cols) of output grid
        output_bounds
            bounds (min X, min Y, max X, max Y) of output grid
        output_nodata
            value for no data in output grid
        radius
            size of interpolation window

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if radius is None:
            radius = self.window_size

        if output_nodata is None:
            output_nodata = self.default_nodata

        return gdal.Grid('', self.dataset, format='MEM', width=output_shape[1], height=output_shape[0],
                         outputBounds=output_bounds,
                         algorithm=f'invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:nodata={int(output_nodata)}')

    def __invdist_gdal_linear_scipy(self, output_shape: (int, int), output_bounds: (float, float, float, float),
                                    output_nodata: float = None, radius: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance and then linearly (using SciPy).

        Parameters
        ----------
        output_shape
            shape (rows, cols) of output grid
        output_bounds
            bounds (min X, min Y, max X, max Y) of output grid
        output_nodata
            value for no data in output grid
        radius
            size of interpolation window

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if radius is None:
            radius = self.window_size

        if output_nodata is None:
            output_nodata = self.default_nodata

        # interpolate using inverse distance in GDAL
        interpolated_raster = self.__invdist_gdal(output_shape, output_bounds, output_nodata, radius)
        output_x, output_y = numpy.meshgrid(numpy.arange(interpolated_raster.RasterXSize),
                                            numpy.arange(interpolated_raster.RasterYSize))

        interpolated_data = interpolated_raster.ReadAsArray()
        y_values, x_values = numpy.where(interpolated_data != output_nodata)
        z_values = interpolated_data[y_values, x_values]
        del interpolated_data

        # interpolate linearly in SciPy
        interpolated_data = scipy.interpolate.griddata((x_values, y_values), z_values, (output_x, output_y),
                                                       method='linear', fill_value=output_nodata)
        interpolated_data[numpy.isnan(interpolated_data)] = output_nodata

        band_1 = interpolated_raster.GetRasterBand(1)
        band_1.SetNoDataValue(output_nodata)
        band_1.WriteArray(interpolated_data)
        del band_1

        return interpolated_raster

    def __kriging_pykrige(self, output_shape: (int, int), output_bounds: (float, float, float, float),
                          chunk_size: (float, float), output_nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via kriging (using PyKrige).

        Parameters
        ----------
        output_shape
            shape (rows, cols) of output grid
        output_bounds
            bounds (min X, min Y, max X, max Y) of output grid
        output_resolution
            resolution of output grid
        chunk_size
            size of each chunk with which to divide the interpolation task, in units; setting this value to less than 40x20 has adverse consequences on memory consumption
        output_nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if type(output_shape) is not numpy.array:
            output_shape = numpy.array(output_shape)

        if type(chunk_size) is not numpy.array:
            chunk_size = numpy.array(chunk_size)

        if output_nodata is None:
            output_nodata = self.default_nodata

        assert not numpy.any(chunk_size < numpy.array((40, 20))), 'a small chunk size may freeze your computer'

        output_resolution = (numpy.array(output_bounds[2:4]) - numpy.array(output_bounds[0:2])) / output_shape

        # define function to perform kriging (using PyKrige) from given N x 3 points onto the given regular grid
        def grid_krige(input_points: numpy.array, output_x: [float], output_y: [float]) -> (numpy.array, numpy.array):
            return OrdinaryKriging(input_points[:, 0], input_points[:, 1], input_points[:, 2], variogram_model='linear',
                                   verbose=False, enable_plotting=False).execute('grid', output_x, output_y)

        chunk_shape = numpy.array(numpy.round(numpy.flip(chunk_size / output_resolution)), numpy.int)
        expand_indices = numpy.round(chunk_shape / 2).astype(numpy.int)

        input_points_sw = numpy.array((self.input_bounds[0], self.input_bounds[1]))

        interpolated_grid_x = numpy.linspace(output_bounds[0], output_bounds[2], output_shape[1])
        interpolated_grid_y = numpy.linspace(output_bounds[1], output_bounds[3], output_shape[0])
        interpolated_grid_shape = numpy.array((len(interpolated_grid_y), len(interpolated_grid_x)), numpy.int)
        interpolated_grid_values = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)
        interpolated_grid_variance = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)

        chunks = numpy.array(numpy.floor(interpolated_grid_shape / chunk_shape), numpy.int)

        chunk_grid_index = numpy.array((0, 0), numpy.int)
        with futures.ProcessPoolExecutor() as concurrency_pool:
            running_futures = {}

            for chunk_index in range(int(numpy.product(chunks))):
                # determine indices of the lower left and upper right bounds of the chunk within the output grid
                grid_start_index = chunk_grid_index * chunk_shape
                grid_end_index = grid_start_index + chunk_shape
                grid_slice = tuple(slice(grid_start_index[axis], grid_end_index[axis])
                                   for axis in range(len(chunk_shape)))

                # expand the interpolation area past the edges of the chunk to make the chunk contextually aware
                interpolation_sw = input_points_sw + numpy.flip((grid_start_index - expand_indices) * output_resolution)
                interpolation_ne = input_points_sw + numpy.flip((grid_end_index + expand_indices) * output_resolution)

                # get points within the interpolation area
                interpolation_points = self.points[numpy.where((self.points[:, 0] >= interpolation_sw[0]) &
                                                               (self.points[:, 0] <= interpolation_ne[0]) &
                                                               (self.points[:, 1] >= interpolation_sw[1]) &
                                                               (self.points[:, 1] <= interpolation_ne[1]))[0]]

                # only interpolate if there are points
                if interpolation_points.shape[0] >= 3:
                    current_future = concurrency_pool.submit(grid_krige, interpolation_points,
                                                             interpolated_grid_x[grid_slice[1]],
                                                             interpolated_grid_y[grid_slice[0]])
                    running_futures[current_future] = grid_slice

                # go to next chunk
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
                    # TODO figure out why this causes an error
                    print(error)
                interpolated_grid_values[grid_slice] = chunk_interpolated_values
                interpolated_grid_variance[grid_slice] = chunk_interpolated_variance

        interpolated_grid_uncertainty = numpy.sqrt(interpolated_grid_variance) * 2.5
        del interpolated_grid_variance

        output_dataset = gdal.GetDriverByName('MEM').Create('temp', output_shape[1], output_shape[0], 2,
                                                            gdal.GDT_Float32)

        output_dataset.SetGeoTransform((output_bounds[0], output_resolution, 0.0,
                                        output_bounds[1], 0.0, output_resolution))
        output_dataset.SetProjection(self.crs_wkt)

        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(output_nodata)
        band_1.WriteArray(interpolated_grid_values)

        band_2 = output_dataset.GetRasterBand(2)
        band_2.SetNoDataValue(output_nodata)
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
            nodata = self.default_nodata

        # if the input raster is a GDAL raster dataset, extract the data from the first band
        if type(interpolated_raster) is gdal.Dataset:
            interpolated_raster = numpy.flip(interpolated_raster.GetRasterBand(band_index).ReadAsArray(), axis=0)

        # replace `nodata` values with NaN
        interpolated_raster[interpolated_raster == nodata] = numpy.nan

        # get minimum and maximum values for all three dimensions
        min_x, min_y, max_x, max_y = self.input_bounds
        z_values = numpy.concatenate((numpy.ravel(interpolated_raster), self.points[:, 2]))
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
        survey_points_axis.scatter(self.points[:, 0], self.points[:, 1], c=self.points[:, 2], s=1,
                                   vmin=min_z, vmax=max_z)

        interpolated_raster_axis.imshow(interpolated_raster, extent=(min_x, max_x, min_y, max_y), aspect='auto',
                                        vmin=min_z, vmax=max_z)

        # create colorbar
        figure.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=min_z, vmax=max_z)),
                        ax=(interpolated_raster_axis, survey_points_axis))

        # pause program execution and show the figure
        matplotlib.pyplot.show()

    def __plot_concave_hull(self):
        triangles = scipy.spatial.Delaunay(self.points[:, 0:2])
        matplotlib.pyplot.plot(*self.concave_hull.exterior.xy, c='r')
        matplotlib.pyplot.triplot(self.points[:, 0], self.points[:, 1], triangles=triangles.simplices, c='g')
        matplotlib.pyplot.scatter(self.points[:, 0], self.points[:, 1], c='b')
        matplotlib.pyplot.show()


class RasterInterpolator(PointInterpolator):
    def __init__(self, raster: gdal.Dataset, coverage_raster_files: [str], file_size: int, band_index: int = 1,
                 default_nodata: float = None, window_scalar: float = 1):
        """
        Create a new raster interpolator object for the given GDAL raster dataset.

        Parameters
        ----------
        raster
            GDAL raster dataset
        coverage_raster_files
            paths to coverage rasters
        file_size
            size of the file in megabytes
        band_index
            raster band (1-indexed)
        default_nodata
            value to set for no data in the output interpolation
        window_scalar
            multiplier to use to expand interpolation radius from the default (minimum nearest-neighbor distance)
        """

        raster_nodata = raster.GetNoDataValue()

        if default_nodata is None:
            default_nodata = raster_nodata

        raster_shape = raster.RasterYSize, raster.RasterXSize

        geotransform = raster.GetGeotransform()
        raster_x, raster_y = numpy.meshgrid(
            numpy.linspace(geotransform[0], geotransform[0] + geotransform[1] * raster_shape[1], raster_shape[1]),
            numpy.linspace(geotransform[3], geotransform[3] + geotransform[5] * raster_shape[0], raster_shape[0]))
        raster_data = raster.GetRasterBand(band_index).ReadAsArray()

        raster_x[raster_data == raster_nodata] = numpy.nan
        raster_y[raster_data == raster_nodata] = numpy.nan
        raster_data[raster_data == raster_nodata] = numpy.nan

        spatial_reference = osr.SpatialReference(wkt=_get_gdal_crs_wkt(raster))
        point_dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = point_dataset.CreateLayer('pts', spatial_reference, geom_type=ogr.wkbPoint)

        for point_index, point in enumerate(raster_data):
            geometry = ogr.Geometry(ogr.wkbPoint)
            geometry.AddPoint(raster_x[point_index], raster_y[point_index], point)
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(geometry)
            layer.CreateFeature(feature)

        super().__init__(point_dataset, default_nodata, window_scalar, layer_index=0)
        self.coverage_raster_files = coverage_raster_files
        self.file_size = file_size


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


def _alpha_hull(points: numpy.array, radius: float = 1.0):
    """
    Calculate the alpha shape (concave hull) of the given points.
    inspired by https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    radius
        radius of triangulation
    """

    if points.shape[0] < 4:
        raise ValueError('need at least 4 points to perform triangulation')

    triangles = scipy.spatial.Delaunay(points)
    boundary_edges = set()
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
            for indices in ((index_a, index_b), (index_b, index_c), (index_c, index_a)):
                if indices not in boundary_edges and reversed(indices) not in boundary_edges:
                    boundary_edges.add(indices)
                else:
                    # remove the edge if it already exists (that is, if it is shared by another triangle)
                    boundary_edges.remove(indices)

    boundary_edges = points[numpy.array(tuple(boundary_edges))].tolist()
    return unary_union(list(polygonize(shapely.geometry.MultiLineString(boundary_edges))))


def _mask_raster(raster: gdal.Dataset, mask: gdal.Dataset, mask_value: float = None,
                 mask_band_index: int = 1) -> gdal.Dataset:
    """
    Mask a raster using the given mask values in the given mask band.
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
        raster band of mask (1-indexed)

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


def _epsg_to_wkt(epsg: int) -> str:
    """
    Use OSR to get the well-known text of a CRS from its EPSG code.

    Parameters
    ----------
    epsg
        EPSG code of CRS

    Returns
    -------
        well-known text of CRS
    """

    spatial_reference = osr.SpatialReference()
    spatial_reference.ImportFromEPSG(epsg)
    return spatial_reference.ExportToWkt()


def _get_gdal_crs_wkt(dataset: gdal.Dataset, layer_index: int = 0) -> str:
    """
    extract the well-known text of the CRS from the given GDAL dataset

    Parameters
    ----------
    dataset
        GDAL dataset (either raster or vector)
    layer_index
        index of vector layer

    Returns
    -------
        well-known text of CRS
    """

    crs_wkt = dataset.GetProjectionRef()

    if crs_wkt == '':
        vector_layer = dataset.GetLayerByIndex(layer_index)
        if vector_layer is not None:
            crs_wkt = vector_layer.GetSpatialRef().ExportToWKT()
        else:
            crs_wkt = _epsg_to_wkt(4326)
    elif re.match('^EPSG:[0-9]+$', crs_wkt):
        crs_wkt = _epsg_to_wkt(int(crs_wkt[5:]))

    return crs_wkt
