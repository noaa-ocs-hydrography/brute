# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

from concurrent import futures
from datetime import datetime
from re import match
from typing import Union

import fiona
import numpy
import rasterio
from affine import Affine
from matplotlib import pyplot
from matplotlib.cm import ScalarMappable, get_cmap
from matplotlib.colors import Normalize
from osgeo import gdal, osr
from pykrige.ok import OrdinaryKriging
from rasterio.crs import CRS
from rasterio.features import shapes as rasterio_shapes
from rasterio.io import MemoryFile
from rasterio.mask import mask
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import zoom
from scipy.spatial import cKDTree, Delaunay
from shapely.geometry import Polygon, MultiPolygon, MultiLineString, shape as shapely_shape
from shapely.ops import unary_union, polygonize


class Interpolator:
    """An abstraction for data interpolation."""

    def __init__(self, dataset: gdal.Dataset, default_nodata: float = None, window_scalar: float = 1.5, index: int = 0,
                 sidescan_raster_filenames: [str] = None):
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
        index
            index of layer / band to read from the given dataset
        """

        self.dataset = dataset
        self.is_raster = dataset.GetProjectionRef() != ''
        self.index = index

        if self.index == 0 and self.is_raster:
            self.index = 1
        elif self.index < 0:
            self.index += self.dataset.RasterCount + 1 if self.is_raster else self.dataset.GetLayerCount()

        self.default_nodata = default_nodata

        self.crs_wkt = _gdal_crs_wkt(self.dataset, self.index)

        if self.is_raster:
            if sidescan_raster_filenames is None:
                raise ValueError('raster interpolation requires coverage rasters')

            elevation_band = dataset.GetRasterBand(self.index)
            elevation_nodata = elevation_band.GetNoDataValue()
            elevation_data = elevation_band.ReadAsArray()

            raster_shape = dataset.RasterYSize, dataset.RasterXSize
            raster_geotransform = dataset.GetGeoTransform()
            raster_nw_corner = numpy.array((raster_geotransform[0], raster_geotransform[3]))
            raster_resolution = numpy.array((raster_geotransform[1], raster_geotransform[5]))

            self.input_bounds = numpy.array([raster_geotransform[0],
                                             raster_geotransform[3] - raster_shape[0] + abs(raster_geotransform[5]),
                                             raster_geotransform[0] + raster_shape[1] + raster_geotransform[1],
                                             raster_geotransform[3]])

            if self.default_nodata is None:
                self.default_nodata = elevation_nodata

            self.window_size = numpy.mean(window_scalar * raster_resolution)

            sidescan_coverage = _combined_coverage_within_window([raster_filename for raster_filename in sidescan_raster_filenames if
                                                                  '.tif' in raster_filename or '.gpkg' in raster_filename],
                                                                 raster_nw_corner, raster_resolution, raster_shape)
            elevation_coverage = _coverage_mask(elevation_data, elevation_nodata)

            # self.elevation_region = _alpha_hull(
            #     _geoarray_to_points(elevation_coverage & sidescan_coverage, raster_nw_corner, raster_resolution,
            #                         False)[::50, 0:2]).simplify(100)

            # TODO change this from interpolating within the entire region to only interpolating within the intersection between sidescan and elevation coverage (afterwards adding back original elevation outside sidescan)
            self.interpolation_region = _vectorize_geoarray(sidescan_coverage | elevation_coverage, raster_nw_corner,
                                                            raster_resolution, False)

            # if type(self.interpolation_region) in (GeometryCollection, MultiPolygon):
            #     x = []
            #     y = []
            #
            #     for feature in self.interpolation_region:
            #         current_x, current_y = feature.exterior.xy
            #         x.extend(current_x)
            #         y.extend(current_y)
            #
            #     self.interpolation_region = _alpha_hull(numpy.stack((x, y), axis=1))

            self.points = _geoarray_to_points(numpy.where(sidescan_coverage, elevation_data, elevation_nodata),
                                              raster_nw_corner, raster_resolution, elevation_nodata)
        else:
            if self.default_nodata is None:
                self.default_nodata = 1000000.0

            self.points = gdal_to_xyz(self.dataset, self.index, self.default_nodata)
            points_2d = self.points[:, 0:2]

            # get the bounds (min X, min Y, max X, max Y)
            self.input_bounds = numpy.min(points_2d[:, 0]), numpy.min(points_2d[:, 1]), \
                                numpy.max(points_2d[:, 0]), numpy.max(points_2d[:, 1])

            # set the size of the interpolation window size to be the minimum of all nearest neighbor distances
            self.window_size = window_scalar * numpy.max(cKDTree(points_2d).query(points_2d, k=2)[0][:, 1])

            # get the concave hull of the survey points
            self.interpolation_region = _alpha_hull(points_2d, self.window_size)

    def interpolate(self, method: str, output_resolution: (float, float), output_nodata: float = None, plot: bool = False) -> gdal.Dataset:
        """
        Generate a raster from the input data using the given interpolation method.

        Parameters
        ----------
        method
            interpolation method
        output_resolution
            resolution of output grid
        output_nodata
            value to set for no data in interpolated output
        plot
            whether to plot the interpolated result

        Returns
        -------
            interpolated GDAL dataset
        """

        if output_nodata is None:
            output_nodata = self.default_nodata

        output_shape, output_bounds = _shape_from_cell_size(output_resolution, self.input_bounds)

        output_shape = numpy.array(output_shape)
        output_bounds = numpy.array(output_bounds)

        start_time = datetime.now()

        if self.is_raster:
            if method not in ('linear', 'linear_scipy', 'kriging'):
                raise NotImplementedError(f'raster interpolation method "{method}" is not supported')

            if method == 'linear':
                method = f'{method}_scipy'

        if method == 'linear':
            # linear interpolation (`gdal.Grid`)
            interpolated_dataset = self.__linear_gdal(output_shape, output_bounds, output_nodata)
        elif method == 'linear_scipy':
            # linear interpolation (`scipy.interpolate.griddata`)
            interpolated_dataset = self.__linear_scipy(output_shape, output_bounds, output_nodata)
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
            raise NotImplementedError(f'interpolation method "{method}" is not supported')

        print(f'{method} interpolation took {(datetime.now() - start_time).total_seconds()} s')

        # mask the interpolated data to the interpolation region
        interpolation_mask = self.__get_mask_like(self.interpolation_region, interpolated_dataset)
        interpolated_dataset = _mask_raster(interpolated_dataset, interpolation_mask)

        if plot:
            self.__plot(interpolated_dataset, method, show=True)

        return interpolated_dataset

    def __get_mask(self, mask_polygon: Polygon, shape: (float, float), resolution: (float, float),
                   nw_corner: (float, float), nodata: float, crs_wkt: str = None) -> gdal.Dataset:
        """
        Get a mask of the survey area with the given raster properties.

        Parameters
        ----------
        mask_polygon
            Shapely polygon with which to create mask
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

        with MemoryFile() as rasterio_memory_file:
            with rasterio_memory_file.open(driver='GTiff', width=shape[1], height=shape[0], count=1,
                                           crs=CRS.from_wkt(crs_wkt),
                                           transform=_get_affine(nw_corner, resolution),
                                           dtype=rasterio.float64,
                                           nodata=numpy.array([nodata]).astype(rasterio.float64).item()) as memory_raster:
                memory_raster.write(numpy.full(shape, 1, dtype=rasterio.float64), 1)

            with rasterio_memory_file.open() as memory_raster:
                masked_data, masked_transform = mask(memory_raster, [mask_polygon])

        output_dataset = gdal.GetDriverByName('MEM').Create('', shape[1], shape[0], 1, gdal.GDT_Float32)
        output_dataset.SetProjection(crs_wkt)
        output_dataset.SetGeoTransform((masked_transform.c, masked_transform.a, masked_transform.b,
                                        masked_transform.f, masked_transform.d, masked_transform.e))
        output_band = output_dataset.GetRasterBand(1)
        output_band.SetNoDataValue(nodata)
        output_band.WriteArray(masked_data[0, :, :])

        return output_dataset

    def __get_mask_like(self, mask_region: Polygon, like_raster: gdal.Dataset,
                        band_index: int = 1) -> gdal.Dataset:
        """
        Get a mask of the survey area with the same properties as the given raster.

        Parameters
        ----------
        mask_region
            Shapely polygon with which to create mask
        like_raster
            raster to copy
        band_index
            raster band (1-indexed)

        Returns
        -------
        gdal.Dataset
            mask
        """

        geotransform = like_raster.GetGeoTransform()
        raster_band = like_raster.GetRasterBand(band_index)
        return self.__get_mask(mask_polygon=mask_region,
                               shape=(like_raster.RasterYSize, like_raster.RasterXSize),
                               resolution=(geotransform[1], geotransform[5]),
                               nw_corner=(geotransform[0], geotransform[3]),
                               nodata=raster_band.GetNoDataValue(),
                               crs_wkt=_gdal_crs_wkt(like_raster))

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

        if type(output_bounds) is not numpy.array:
            output_bounds = numpy.array(output_bounds)

        if output_nodata is None:
            output_nodata = self.default_nodata

        output_resolution = (output_bounds[2:4] - output_bounds[0:2]) / numpy.flip(output_shape)

        # interpolate using SciPy griddata
        output_x, output_y = numpy.meshgrid(numpy.arange(output_bounds[0], output_bounds[2] + output_resolution[0], output_resolution[0]),
                                            numpy.arange(output_bounds[1], output_bounds[3] + output_resolution[1], output_resolution[1]))
        interpolated_data = griddata((self.points[:, 0], self.points[:, 1]), self.points[:, 2], (output_x, output_y),
                                     method='linear', fill_value=output_nodata)

        output_data = numpy.flip(interpolated_data, axis=0)
        output_data[output_data == output_nodata] = numpy.nan
        pyplot.imshow(output_data, extent=output_bounds[[0, 2, 1, 3]])
        _plot_raster(self.dataset, show=True)

        output_raster = gdal.GetDriverByName('MEM').Create('', int(interpolated_data.shape[1]),
                                                           int(interpolated_data.shape[0]), 1, gdal.GDT_Float32)
        output_raster.SetGeoTransform((output_bounds[0], output_resolution[0], 0.0,
                                       output_bounds[1], 0.0, output_resolution[1]))
        output_raster.SetProjection(self.crs_wkt)
        band_1 = output_raster.GetRasterBand(1)
        band_1.SetNoDataValue(output_nodata)
        band_1.WriteArray(interpolated_data)
        del band_1
        return output_raster

    def __invdist_gdal(self, output_shape: (int, int), output_bounds: (float, float, float, float),
                       output_nodata: float = None, radius: float = None) -> gdal.Dataset:
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
        interpolated_data = griddata((x_values, y_values), z_values, (output_x, output_y), method='linear',
                                     fill_value=output_nodata)
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

        if type(output_bounds) is not numpy.array:
            output_bounds = numpy.array(output_bounds)

        if type(chunk_size) is not numpy.array:
            chunk_size = numpy.array(chunk_size)

        if output_nodata is None:
            output_nodata = self.default_nodata

        assert not numpy.any(chunk_size < numpy.array((40, 20))), 'a small chunk size may freeze your computer'

        output_resolution = (output_bounds[2:4] - output_bounds[0:2]) / numpy.flip(output_shape)

        chunk_shape = numpy.array(numpy.round(chunk_size / output_resolution), numpy.int)
        expand_indices = numpy.round(chunk_shape).astype(numpy.int)

        input_points_sw = numpy.array(self.input_bounds[:2])

        interpolated_grid_x = numpy.linspace(output_bounds[0], output_bounds[2], output_shape[1])
        interpolated_grid_y = numpy.linspace(output_bounds[1], output_bounds[3], output_shape[0])
        interpolated_grid_shape = numpy.array((len(interpolated_grid_y), len(interpolated_grid_x)), numpy.int)
        interpolated_grid_values = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)
        interpolated_grid_variance = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)

        chunk_grid_shape = numpy.array(numpy.floor(interpolated_grid_shape / chunk_shape), numpy.int)
        chunk_grid_index = numpy.array((0, 0), numpy.int)
        with futures.ProcessPoolExecutor() as concurrency_pool:
            running_futures = {}

            for _ in range(numpy.product(chunk_grid_shape)):
                # determine indices of the lower left and upper right bounds of the chunk within the output grid
                grid_start_index = chunk_grid_index * chunk_shape
                grid_end_index = grid_start_index + chunk_shape

                # expand the interpolation area past the edges of the chunk to make the chunk contextually aware
                interpolation_sw = input_points_sw + numpy.flip(grid_start_index - expand_indices) * output_resolution
                interpolation_ne = input_points_sw + numpy.flip(grid_end_index + expand_indices) * output_resolution

                # get points within the interpolation area
                interpolation_points = self.points[numpy.where((self.points[:, 0] >= interpolation_sw[0]) &
                                                               (self.points[:, 0] <= interpolation_ne[0]) &
                                                               (self.points[:, 1] >= interpolation_sw[1]) &
                                                               (self.points[:, 1] <= interpolation_ne[1]))[0]]

                # only interpolate if there are points
                if interpolation_points.shape[0] >= 3:
                    grid_slice = tuple(slice(grid_start_index[dimension], grid_end_index[dimension]) for dimension in
                                       range(len(chunk_shape)))
                    current_future = concurrency_pool.submit(_krige_points_onto_grid, interpolation_points,
                                                             interpolated_grid_x[grid_slice[1]],
                                                             interpolated_grid_y[grid_slice[0]])
                    running_futures[current_future] = grid_slice

                # go to next chunk
                if chunk_grid_index[0] >= chunk_grid_shape[0]:
                    chunk_grid_index[0] = 0
                    chunk_grid_index[1] += 1
                else:
                    chunk_grid_index[0] += 1

            for completed_future in futures.as_completed(running_futures):
                grid_slice = running_futures[completed_future]
                chunk_interpolated_values, chunk_interpolated_variance = completed_future.result()
                interpolated_grid_values[grid_slice] = chunk_interpolated_values
                interpolated_grid_variance[grid_slice] = chunk_interpolated_variance

        interpolated_grid_uncertainty = numpy.sqrt(interpolated_grid_variance) * 2.5
        del interpolated_grid_variance

        output_raster = gdal.GetDriverByName('MEM').Create('', int(output_shape[1]), int(output_shape[0]), 2,
                                                           gdal.GDT_Float32)
        output_raster.SetGeoTransform((output_bounds[0], output_resolution[0], 0.0,
                                       output_bounds[1], 0.0, output_resolution[1]))
        output_raster.SetProjection(self.crs_wkt)

        band_1 = output_raster.GetRasterBand(1)
        band_1.SetNoDataValue(output_nodata)
        band_1.WriteArray(interpolated_grid_values)

        band_2 = output_raster.GetRasterBand(2)
        band_2.SetNoDataValue(output_nodata)
        band_2.WriteArray(interpolated_grid_uncertainty)

        del band_1, band_2
        return output_raster

    def __plot(self, interpolated_raster: Union[numpy.array, gdal.Dataset], interpolation_method: str,
               nodata: float = None, band_index: int = 1, show: bool = False):
        """
        Plot preinterpolated points and an interpolated raster side-by-side on synchronized subplots.

        Parameters
        ----------
        interpolated_raster
            array of interpolated data
        interpolation_method
            method of interpolation
        nodata
            value for no data in given data
        band_index
            raster band (1-indexed)
        """

        # if the input raster is a GDAL raster dataset, extract the data from the first band
        if type(interpolated_raster) is gdal.Dataset:
            geotransform = interpolated_raster.GetGeoTransform()
            min_x = geotransform[0]
            max_x = geotransform[0] + interpolated_raster.RasterXSize + geotransform[1]
            min_y = geotransform[3] - interpolated_raster.RasterYSize + abs(geotransform[5])
            max_y = geotransform[3]
            interpolated_band = interpolated_raster.GetRasterBand(band_index)
            default_nodata = interpolated_band.GetNoDataValue()
            interpolated_raster = numpy.flip(interpolated_band.ReadAsArray(), axis=0)
        else:
            min_x, min_y, max_x, max_y = self.input_bounds
            default_nodata = self.default_nodata

        if nodata is None:
            nodata = default_nodata

        # replace `nodata` values with NaN
        interpolated_raster[interpolated_raster == nodata] = numpy.nan

        # get minimum and maximum values for all three dimensions
        z_values = numpy.concatenate((numpy.ravel(interpolated_raster),
                                      numpy.ravel(self.dataset.ReadAsArray()[band_index - 1]) if self.is_raster else self.points[:, 2]))
        min_z = numpy.nanmin(z_values)
        max_z = numpy.nanmax(z_values)

        # create a new Matplotlib figure window
        figure = pyplot.figure()

        # create subplots
        original_data_axis = figure.add_subplot(1, 2, 1)
        original_data_axis.set_title('survey data')
        original_data_axis.set_xlim(min_x, max_x)
        original_data_axis.set_ylim(min_y, max_y)
        interpolated_data_axis = figure.add_subplot(1, 2, 2, sharex=original_data_axis, sharey=original_data_axis)
        interpolated_data_axis.set_title(f'{interpolation_method} interpolation to raster')

        # plot data
        if self.is_raster:
            _plot_raster(self.dataset, self.index, aspect='auto', vmin=min_z, vmax=max_z)
        else:
            original_data_axis.scatter(self.points[:, 0], self.points[:, 1], c=self.points[:, 2], s=1, vmin=min_z,
                                       vmax=max_z)

        interpolated_data_axis.imshow(interpolated_raster, extent=(min_x, max_x, min_y, max_y), aspect='auto',
                                      vmin=min_z, vmax=max_z)

        # create colorbar
        figure.colorbar(ScalarMappable(norm=Normalize(vmin=min_z, vmax=max_z)),
                        ax=(interpolated_data_axis, original_data_axis))

        # pause program execution and show the figure
        if show:
            pyplot.show()

    def __plot_concave_hull(self, axis: pyplot.Axes = None, show: bool = False, **kwargs):
        if axis is None:
            axis = pyplot.gca()

        axis.plot(*self.interpolation_region.exterior.xy, c='r')
        axis.scatter(self.points[:, 0], self.points[:, 1], c='g')
        triangles = Delaunay(self.points[:, 0:2])
        axis.triplot(self.points[:, 0], self.points[:, 1], triangles=triangles.simplices, c='b', **kwargs)

        if show:
            pyplot.show()


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
            raise ValueError(f'invalid index for raster band "{index}"')

        raster_band = dataset.GetRasterBand(index)
        geotransform = dataset.GetGeoTransform()

        if nodata is None:
            nodata = raster_band.GetNoDataValue()

        return _geoarray_to_points(raster_band.ReadAsArray(), nw_corner=(geotransform[0], geotransform[3]),
                                   resolution=(geotransform[1], geotransform[5]), nodata=nodata)


def _geoarray_to_points(grid: numpy.array, nw_corner: (float, float), resolution: (float, float),
                        nodata: float = None) -> numpy.array:
    """
    Extract XYZ points from an array of data using the given geographic reference.

    Parameters
    ----------
    grid
        array of gridded data
    nw_corner
        X, Y coordinates of northwest corner
    resolution
        cell size
    nodata
        value to exclude from point creation from the input grid

    Returns
    -------
        N x 3 array of XYZ points
    """

    x_values, y_values = numpy.meshgrid(numpy.linspace(nw_corner[0], nw_corner[0] + resolution[0] * grid.shape[1], grid.shape[1]),
                                        numpy.linspace(nw_corner[1], nw_corner[1] + resolution[1] * grid.shape[0], grid.shape[0]))

    if nodata is not None:
        x_values = x_values[grid != nodata]
        y_values = y_values[grid != nodata]
        grid = grid[grid != nodata]

    return numpy.stack((x_values, y_values, grid), axis=1)


def _alpha_hull(points: numpy.array, max_length: float = None) -> Polygon:
    """
    Calculate the alpha shape (concave hull) of the given points.
    inspired by https://sgillies.net/2012/10/13/the-fading-shape-of-alpha.html

    Parameters
    ----------
    points
        N x 3 array of XYZ points
    max_length
        maximum length to include in the convex boundary
    """

    if points.shape[0] < 4:
        raise ValueError('need at least 4 points to perform triangulation')

    if max_length is None:
        max_length = numpy.max(cKDTree(points).query(points, k=2)[0][:, 1])

    triangles = Delaunay(points)
    boundary_edges = set()
    for index_a, index_b, index_c in triangles.simplices:
        point_a = points[index_a]
        point_b = points[index_b]
        point_c = points[index_c]

        vectors = numpy.stack((point_a - point_b, point_b - point_c, point_c - point_a))
        lengths = numpy.hypot(vectors[:, 0], vectors[:, 1])
        semiperimeter = numpy.sum(lengths) / 2

        if numpy.product(lengths) / (2 * numpy.sqrt(semiperimeter * numpy.product(semiperimeter - lengths))) < max_length:
            for indices in ((index_a, index_b), (index_b, index_c), (index_c, index_a)):
                if indices not in boundary_edges and reversed(indices) not in boundary_edges:
                    boundary_edges.add(indices)
                else:
                    # remove the edge if it already exists (if it is shared by another triangle)
                    boundary_edges.remove(indices)

    boundary_edges = points[numpy.array(tuple(boundary_edges))].tolist()
    return unary_union(list(polygonize(MultiLineString(boundary_edges))))


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


def _shape_from_cell_size(resolution: (float, float), bounds: (float, float, float, float)) -> ((int, int), (float, float, float, float)):
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
        shape and bounds (rounded to cell)
    """

    if type(resolution) is not numpy.array:
        resolution = numpy.array(resolution)

    min_x, min_y, max_x, max_y = bounds

    sw_corner = numpy.array(bounds[:2])
    ne_corner = numpy.array(bounds[2:])

    cell_remainder = (ne_corner - sw_corner) % resolution
    ne_corner -= cell_remainder

    shape = numpy.flip((ne_corner - sw_corner) / resolution).astype(int)

    return shape, (*sw_corner, *ne_corner)


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


def _gdal_crs_wkt(dataset: gdal.Dataset, layer_index: int = 0) -> str:
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
            crs_wkt = vector_layer.GetSpatialRef().ExportToWkt()
        else:
            crs_wkt = _epsg_to_wkt(4326)
    elif match('^EPSG:[0-9]+$', crs_wkt):
        crs_wkt = _epsg_to_wkt(int(crs_wkt[5:]))

    return crs_wkt


def _combined_coverage_within_window(coverage_raster_filenames: [str], output_nw_corner: (float, float), output_resolution: (float, float),
                                     output_shape: (int, int)) -> numpy.array:
    """
    Return a boolean mask of the given rasters where data exists within the given window.

    Parameters
    ----------
    coverage_raster_filenames
        list of filenames of coverage rasters
    output_nw_corner
        XY coordinates of northwest corner of output mask
    output_resolution
        XY resolution of output mask
    output_shape
        shape of output mask

    Returns
    -------
        boolean array indicating where data exists
    """

    if type(output_nw_corner) is not numpy.array:
        output_nw_corner = numpy.array(output_nw_corner)

    output_se_corner = output_nw_corner + (output_resolution * numpy.flip(output_shape))

    if type(output_resolution) is not numpy.array:
        output_resolution = numpy.array(output_resolution)

    if type(output_shape) is not numpy.array:
        output_shape = numpy.array(output_shape)

    output_coverage_masks = []

    for coverage_raster_filename in coverage_raster_filenames:
        with rasterio.open(coverage_raster_filename) as coverage_raster:
            coverage_transform = coverage_raster.transform
            coverage_shape = coverage_raster.shape

        coverage_resolution = numpy.array((coverage_transform.a, coverage_transform.e))

        coverage_nw_corner = coverage_transform.c, coverage_transform.f
        coverage_se_corner = coverage_nw_corner + (coverage_resolution * numpy.flip(coverage_shape))

        # ensure the coverage array intersects the given window
        if not (output_nw_corner[0] > coverage_se_corner[0] or output_se_corner[0] < coverage_nw_corner[0] or
                output_se_corner[1] > coverage_nw_corner[1] or output_nw_corner[1] < coverage_se_corner[1]):
            with rasterio.open(coverage_raster_filename) as coverage_raster:
                coverage_data = coverage_raster.read()

            coverage_mask = numpy.where(_coverage_mask(coverage_data), 1, 0)

            # resample coverage data to match BAG resolution
            resolution_ratio = coverage_resolution / output_resolution
            if numpy.any(resolution_ratio != 1):
                coverage_mask = zoom(coverage_mask, zoom=resolution_ratio, order=3, prefilter=False)

            # clip resampled coverage data to the bounds of the BAG
            output_array = numpy.full(output_shape, 0)

            nw_index_delta = numpy.round((output_nw_corner - coverage_nw_corner) / output_resolution).astype(int)
            se_index_delta = numpy.round((output_se_corner - coverage_nw_corner) / output_resolution).astype(int)

            # indices to be written onto the output array
            output_array_index_slices = [slice(0, None), slice(0, None)]

            # BAG leftmost X is to the left of coverage leftmost X
            if nw_index_delta[0] < 0:
                output_array_index_slices[1] = slice(nw_index_delta[0] * -1, output_array_index_slices[1].stop)
                nw_index_delta[0] = 0

            # BAG topmost Y is above coverage topmost Y
            if nw_index_delta[1] < 0:
                output_array_index_slices[0] = slice(nw_index_delta[1] * -1, output_array_index_slices[0].stop)
                nw_index_delta[1] = 0

            # BAG rightmost X is to the right of coverage rightmost X
            if se_index_delta[0] > coverage_mask.shape[1]:
                output_array_index_slices[1] = slice(output_array_index_slices[1].start, coverage_mask.shape[1] - se_index_delta[0])
                se_index_delta[0] = coverage_mask.shape[1]

            # BAG bottommost Y is lower than coverage bottommost Y
            if se_index_delta[1] > coverage_mask.shape[0]:
                output_array_index_slices[0] = slice(output_array_index_slices[0].start, coverage_mask.shape[0] - se_index_delta[1])
                se_index_delta[1] = coverage_mask.shape[0]

            coverage_mask = coverage_mask[nw_index_delta[1]: se_index_delta[1], nw_index_delta[0]: se_index_delta[0]]

            # write the relevant coverage data to a slice of the output array corresponding to the coverage extent
            output_array[output_array_index_slices[0], output_array_index_slices[1]] = coverage_mask
            output_coverage_masks.append(output_array)

    # collapse coverage masks into single raster array
    return numpy.logical_or.reduce(output_coverage_masks)


def _coverage_mask(raster_data: numpy.array, nodata: float = None):
    """
    Return a boolean mask of the given array where data exists.

    Parameters
    ----------
    raster_data
        array of gridded data
    nodata
        value where there is no data in the given array

    Returns
    -------
        boolean array indicating where data exists
    """

    # nodata for sidescan coverage is usually the maximum value
    if nodata is None:
        nodata = numpy.max(raster_data)

    if raster_data.shape[0] == 1:
        raster_data = raster_data[0, :, :]

    # TODO find reduced generalization of multiple bands
    if raster_data.shape[0] == 3:
        raster_mask = (raster_data[0, :, :] != nodata) | (raster_data[1, :, :] != nodata) | (raster_data[2, :, :] != nodata)
    else:
        raster_mask = raster_data != nodata

    return raster_mask


def _vectorize_geoarray(raster_data: numpy.array, nw_corner: (float, float), resolution: (float, float),
                        nodata: float = None) -> MultiPolygon:
    """
    Vectorize the extent of the given raster where data exists.

    Parameters
    ----------

    raster_data
        array of gridded data
    nw_corner
        XY coordinates of the northwest corner
    resolution
        XY cell size
    nodata
        value where there is no data in the given array

    Returns
    -------
        Shapely polygon or multipolygon of coverage extent
    """

    raster_mask = _coverage_mask(raster_data, nodata)
    return MultiPolygon(shapely_shape(shape[0]) for shape in
                        rasterio_shapes(numpy.where(raster_mask, 1, 0), mask=raster_mask, transform=_get_affine(nw_corner, resolution)))


def _vectorize_raster_coverage(raster_filename: str) -> MultiPolygon:
    """
    Vectorize the extent of the given raster where data exists.

    Parameters
    ----------
    raster_filename
        file path to raster
    nodata
        value where there is no data in the given raster file

    Returns
    -------
        Shapely polygon or multipolygon of coverage extent
    """

    with rasterio.open(raster_filename) as raster:
        raster_data = raster.read()
        raster_transform = raster.transform

    nw_corner = raster_transform.c, raster_transform.f
    resolution = raster_transform.a, raster_transform.e

    return _vectorize_geoarray(raster_data, nw_corner, resolution)


def _write_region(output_filename: str, region: MultiPolygon, crs_wkt: str):
    geometry_type = 'Polygon' if type(region) is Polygon else 'MultiPolygon'

    # TODO generalize polygon / multipolygon
    coordinates = [[numpy.stack(polygon.exterior.xy, axis=1).tolist()] for polygon in region]

    with fiona.open(output_filename, 'w', 'GPKG', {'geometry': geometry_type, 'properties': {'name': 'str'}},
                    crs_wkt=crs_wkt) as output_vector_file:
        output_vector_file.write(
            {'geometry': {'type': geometry_type, 'coordinates': coordinates}, 'properties': {'name': 'elevation data'}})


def _get_affine(nw_corner: (float, float), resolution: (float, float)) -> Affine:
    """
    Calculate affine transformation from the given geographic parameters.

    Parameters
    ----------
    nw_corner
        XY coordinates of the northwest corner
    resolution
        XY cell size

    Returns
    -------
        Affine transformation
    """

    if type(resolution) is float:
        resolution = tuple(resolution)

    return Affine.translation(nw_corner[0], nw_corner[1]) * Affine.scale(resolution[0], resolution[-1])


def _krige_points_onto_grid(points: numpy.array, grid_x: [float], grid_y: [float]) -> (numpy.array, numpy.array):
    """
    Perform kriging (using PyKrige) from the given points onto a regular grid.

    Parameters
    ----------
    points
        N x 3 array of points to interpolate
    grid_x
        X coordinates of output grid
    grid_y
        Y coordinates of output grid

    Returns
    -------
        interpolated values and uncertainty
    """

    return OrdinaryKriging(points[:, 0], points[:, 1], points[:, 2], variogram_model='linear', verbose=False,
                           enable_plotting=False).execute('grid', grid_x, grid_y)


def _plot_regions(regions: [Polygon], colors: [str] = None, axis: pyplot.Axes = None, show: bool = False, **kwargs):
    if axis is None:
        axis = pyplot.gca()

    if colors is None:
        cmap = get_cmap('gist_rainbow')
        colors = [cmap(color_index / len(regions)) for color_index in range(len(regions))]

    for color_index, geometry in enumerate(regions):
        color = colors[color_index]
        if type(geometry) is Polygon:
            axis.plot(*geometry.exterior.xy, c=color, **kwargs)
        else:
            for polygon in geometry:
                axis.plot(*polygon.exterior.xy, c=color, **kwargs)

    if show:
        pyplot.show()


def _plot_bounds(sw_corner: (float, float), ne_corner: (float, float), axis: pyplot.Axes = None, show: bool = False, **kwargs):
    if axis is None:
        axis = pyplot.gca()

    corner_points = numpy.array([sw_corner, (ne_corner[0], sw_corner[1]), ne_corner, (sw_corner[0], ne_corner[1]), sw_corner])

    axis.plot(corner_points[:, 0], corner_points[:, 1], **kwargs)

    if show:
        pyplot.show()


def _plot_raster(raster: gdal.Dataset, band_index: int = 1, axis: pyplot.Axes = None, show: bool = False, **kwargs):
    if axis is None:
        axis = pyplot.gca()

    raster_geotransform = raster.GetGeoTransform()
    raster_band = raster.GetRasterBand(band_index)
    raster_data = raster_band.ReadAsArray()
    raster_data[raster_data == raster_band.GetNoDataValue()] = numpy.nan

    axis.imshow(raster_data, extent=[raster_geotransform[0], raster_geotransform[0] + raster.RasterXSize + raster_geotransform[1],
                                     raster_geotransform[3] - raster.RasterYSize + abs(raster_geotransform[5]), raster_geotransform[3]],
                **kwargs)

    if show:
        pyplot.show()


class ExtentError(Exception):
    pass
