# -*- coding: utf-8 -*-
"""
interpolator.py

Created on Mon Feb 11 12:55:51 2019

@author: grice

An abstraction for data interpolation.
"""

from concurrent import futures
from datetime import datetime

import numpy
import rasterio
import rasterio.crs
from fuse.utilities import bounds_from_opposite_corners, gdal_crs_wkt, raster_bounds, array_coverage, georeference_to_affine, alpha_hull, \
    raster_edge_points, vectorize_geoarray, geoarray_to_points, gdal_to_xyz, shape_from_cell_size, raster_mask_like, apply_raster_mask, \
    plot_raster, raster_mask, overwrite_raster
from matplotlib import pyplot
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from osgeo import gdal
from pykrige.ok import OrdinaryKriging
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import zoom
from scipy.spatial import cKDTree

CATZOC = {
    'A1': [.01, .5],
    'A2': [.02, 1],
    'B': [.02, 1],
    'C': [.05, 2]
}


class Interpolator:
    """Interpolator for bathymetric surveys."""

    def __init__(self, dataset: gdal.Dataset, nodata: float = None, window_scalar: float = 1.5, band: int = 0,
                 sidescan_rasters: [str] = None, catzoc: str = 'B'):
        """
        Create a new point interpolator object for the given GDAL point cloud dataset.

        Parameters
        ----------
        dataset
            GDAL point cloud dataset
        nodata
            value to set for no data in the output interpolation
        window_scalar
            multiplier to use to expand interpolation radius from minimum point spacing
        band
            index of layer / band to read from the given dataset
        sidescan_rasters
            filenames of sidescan rasters
        catzoc
            CATZOC score
        """

        self.dataset = dataset
        self.is_raster = dataset.GetProjectionRef() != ''
        self.index = band
        self.catzoc = catzoc.upper()
        assert self.catzoc in CATZOC, f'CATZOC score "{self.catzoc}" not supported'

        if self.is_raster:
            self.index += 1
        elif self.index < 0:
            self.index += self.dataset.RasterCount + 1 if self.is_raster else self.dataset.GetLayerCount()

        self.nodata = nodata

        self.crs_wkt = gdal_crs_wkt(self.dataset, self.index)

        if self.is_raster:
            if sidescan_rasters is None:
                raise ValueError('raster interpolation requires coverage rasters')

            elevation_band = dataset.GetRasterBand(self.index)
            elevation_nodata = elevation_band.GetNoDataValue()
            elevation_data = elevation_band.ReadAsArray()
            del elevation_band

            uncertainty_band = dataset.GetRasterBand(self.index)
            uncertainty_nodata = uncertainty_band.GetNoDataValue()
            uncertainty_data = uncertainty_band.ReadAsArray()
            del uncertainty_band

            raster_shape = dataset.RasterYSize, dataset.RasterXSize
            raster_geotransform = dataset.GetGeoTransform()
            raster_origin = numpy.array((raster_geotransform[0], raster_geotransform[3]))
            raster_resolution = numpy.array((raster_geotransform[1], raster_geotransform[5]))

            self.input_bounds = raster_bounds(self.dataset)

            if self.nodata is None:
                self.nodata = elevation_nodata

            # multiplier cannot be too low or it segments the hull, nor too high lest the hull not respect survey concavity
            # TODO maybe make multiplier dynamic?
            self.window_size = numpy.mean(numpy.abs(raster_resolution)) * window_scalar * 30

            sidescan_coverage = _combined_coverage_within_window((raster_filename for raster_filename in sidescan_rasters if
                                                                  '.tif' in raster_filename or '.gpkg' in raster_filename),
                                                                 raster_origin, raster_resolution, raster_shape)
            elevation_coverage = array_coverage(elevation_data, elevation_nodata)

            transform = georeference_to_affine(raster_origin, raster_resolution)

            elevation_region = alpha_hull(raster_edge_points(elevation_coverage, raster_origin, raster_resolution, False)[:, :2],
                                          self.window_size)
            sidescan_region = vectorize_geoarray(sidescan_coverage, transform, False)

            if not sidescan_region.is_valid:
                sidescan_region = sidescan_region.buffer(0)

            self.interpolation_region = elevation_region.intersection(sidescan_region)

            self.points = geoarray_to_points(numpy.where(sidescan_coverage, elevation_data, elevation_nodata), raster_origin,
                                             raster_resolution, elevation_nodata)
            self.uncertainty_points = geoarray_to_points(numpy.where(sidescan_coverage, uncertainty_data, elevation_nodata), raster_origin,
                                                         raster_resolution, uncertainty_nodata)
        else:
            if self.nodata is None:
                self.nodata = 1000000.0

            self.points = gdal_to_xyz(self.dataset, self.index, self.nodata)
            points_2d = self.points[:, :2]

            # get the bounds (min X, min Y, max X, max Y)
            self.input_bounds = numpy.concatenate((numpy.min(points_2d, axis=0), numpy.max(points_2d, axis=0)))

            # set the size of the interpolation window size to be the minimum of all nearest neighbor distances
            self.window_size = window_scalar * numpy.max(cKDTree(points_2d).query(points_2d, k=2)[0][:, 1])

            # get the concave hull of the survey points
            self.interpolation_region = alpha_hull(points_2d, self.window_size)

    def interpolate(self, method: str, resolution: (float, float), nodata: float = None, plot: bool = False) -> gdal.Dataset:
        """
        Generate a raster from the input data using the given interpolation method.

        Parameters
        ----------
        method
            interpolation method
        resolution
            resolution of output grid
        nodata
            value to set for no data in interpolated output
        plot
            whether to plot the interpolated result

        Returns
        -------
            interpolated GDAL dataset
        """

        if nodata is None:
            nodata = self.nodata

        if self.is_raster:
            output_shape = self.dataset.RasterYSize, self.dataset.RasterXSize
            output_bounds = self.input_bounds
        else:
            output_shape, output_bounds = shape_from_cell_size(resolution, self.input_bounds)

        start_time = datetime.now()

        if self.is_raster:
            if method not in ('linear', 'kriging'):
                raise NotImplementedError(f'raster interpolation method "{method}" is not supported')

        if method == 'linear':
            # linear interpolation (`gdal.Grid`)
            interpolated_dataset = self.__linear(output_shape, output_bounds, nodata)
        elif method == 'invdist':
            # inverse distance interpolation (`gdal.Grid`)
            interpolated_dataset = self.__invdist_gdal(output_shape, output_bounds, nodata)
        elif method == 'kriging':
            # interpolate using Ordinary Kriging (`pykrige`) by dividing data into smaller chunks
            interpolated_dataset = self.__kriging_pykrige(output_shape, output_bounds, (40, 40), nodata, 1, 2)
        else:
            raise NotImplementedError(f'interpolation method "{method}" is not supported')

        print(f'{method} interpolation took {(datetime.now() - start_time).total_seconds()} s')

        # mask the interpolated data to the interpolation region
        interpolation_mask = raster_mask_like(self.interpolation_region, interpolated_dataset)
        interpolated_dataset = apply_raster_mask(interpolated_dataset, interpolation_mask)

        if self.is_raster:
            interpolated_dataset = overwrite_raster(self.dataset, interpolated_dataset)

        if plot:
            self.plot(interpolated_dataset, method, show=True)

        return interpolated_dataset

    def __linear(self, shape: (int, int), bounds: (float, float, float, float), nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating linearly.

        Parameters
        ----------
        shape
            shape (rows, cols) of output grid
        bounds
            bounds (min X, min Y, max X, max Y) of output grid
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid

        """

        if type(shape) is not numpy.array:
            shape = numpy.array(shape)

        if type(bounds) is not numpy.array:
            bounds = numpy.array(bounds)

        if nodata is None:
            nodata = self.nodata

        if self.is_raster:
            # interpolate using SciPy griddata
            output_x, output_y = numpy.meshgrid(numpy.linspace(bounds[0], bounds[2], shape[1]),
                                                numpy.linspace(bounds[1], bounds[3], shape[0]))
            interpolated_values = numpy.flip(griddata((self.points[:, 0], self.points[:, 1]), self.points[:, 2], (output_x, output_y),
                                                      method='linear', fill_value=nodata), axis=0)
            geotransform = self.dataset.GetGeoTransform()
        else:
            interpolated_raster = gdal.Grid('', self.dataset, format='MEM', width=shape[1], height=shape[0], outputBounds=bounds,
                                            algorithm=f'linear:radius=0:nodata={int(nodata)}')
            interpolated_band = interpolated_raster.GetRasterBand(1)
            interpolated_values = interpolated_band.ReadAsArray()
            geotransform = interpolated_raster.GetGeoTransform()

        output_raster = gdal.GetDriverByName('MEM').Create('', int(shape[1]), int(shape[0]), 2, gdal.GDT_Float32)
        output_raster.SetGeoTransform(geotransform)
        output_raster.SetProjection(self.crs_wkt)

        band_1 = output_raster.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_values)
        del band_1

        uncertainty = self.__uncertainty(numpy.where(interpolated_values != nodata, interpolated_values, numpy.nan), self.catzoc)

        if self.is_raster:
            uncertainty_band = self.dataset.GetRasterBand(self.index + 1)
            input_uncertainty = uncertainty_band.ReadAsArray()
            uncertainty_nodata = uncertainty_band.GetNoDataValue()
            del uncertainty_band

            input_uncertainty[input_uncertainty == uncertainty_nodata] = numpy.nan
            input_uncertainty = numpy.flip(input_uncertainty, axis=0)
            uncertainty = numpy.sqrt(input_uncertainty ** 2 + uncertainty ** 2)

        uncertainty[numpy.isnan(uncertainty)] = nodata

        band_2 = output_raster.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(uncertainty)
        del band_2

        return output_raster

    def __invdist_gdal(self, shape: (int, int), bounds: (float, float, float, float), nodata: float = None,
                       radius: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance.

        Parameters
        ----------
        shape
            shape (rows, cols) of output grid
        bounds
            bounds (min X, min Y, max X, max Y) of output grid
        nodata
            value for no data in output grid
        radius
            size of interpolation window

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if type(shape) is not numpy.array:
            shape = numpy.array(shape)

        if type(bounds) is not numpy.array:
            bounds = numpy.array(bounds)

        if radius is None:
            radius = self.window_size

        if nodata is None:
            nodata = self.nodata

        interpolated_raster = gdal.Grid('', self.dataset, format='MEM', width=shape[1], height=shape[0], outputBounds=bounds,
                                        algorithm=f'invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:nodata={int(nodata)}')
        interpolated_band = interpolated_raster.GetRasterBand(1)
        interpolated_values = interpolated_band.ReadAsArray()

        output_resolution = (bounds[2:] - bounds[:2]) / numpy.flip(shape)

        output_raster = gdal.GetDriverByName('MEM').Create('', int(shape[1]), int(shape[0]), 2, gdal.GDT_Float32)
        output_raster.SetGeoTransform((bounds[0], output_resolution[0], 0.0, bounds[1], 0.0, output_resolution[1]))
        output_raster.SetProjection(self.crs_wkt)

        band_1 = output_raster.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_values)
        del band_1

        uncertainty = self.__uncertainty(numpy.where(interpolated_values != nodata, interpolated_values, numpy.nan), self.catzoc)

        if self.is_raster:
            uncertainty_band = self.dataset.GetRasterBand(self.index + 1)
            input_uncertainty = uncertainty_band.ReadAsArray()
            uncertainty_nodata = uncertainty_band.GetNoDataValue()
            del uncertainty_band

            input_uncertainty[input_uncertainty == uncertainty_nodata] = numpy.nan
            input_uncertainty = numpy.flip(input_uncertainty, axis=0)
            uncertainty = numpy.sqrt(input_uncertainty ** 2 + uncertainty ** 2)

        uncertainty[numpy.isnan(uncertainty)] = nodata

        band_2 = output_raster.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(uncertainty)
        del band_2

        return output_raster

    def __kriging_pykrige(self, shape: (int, int), bounds: (float, float, float, float), chunk: (float, float), nodata: float = None,
                          expansion: float = 1, multiplier: int = 1) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via kriging (using PyKrige).

        Parameters
        ----------
        shape
            shape (rows, cols) of output grid
        bounds
            bounds (min X, min Y, max X, max Y) of output grid
        output_resolution
            resolution of output grid
        chunk
            size of each chunk with which to divide the interpolation task, in units;
            setting this value to less than 40x20 has adverse consequences on memory consumption
        nodata
            value for no data in output grid
        expansion
            indices by which to expand the chunk size to its interpolation points catchment area;
            a value of 1 means only the points within the chunk are considered in interpolation
        multiplier
            how many chunks to fit within a row / column in the chunk grid;
            values above 1 will result in overlap

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if type(shape) is not numpy.array:
            shape = numpy.round(shape).astype(int)

        if type(bounds) is not numpy.array:
            bounds = numpy.array(bounds)

        if type(chunk) is not numpy.array:
            chunk = numpy.array(chunk)

        if nodata is None:
            nodata = self.nodata

        assert not numpy.any(chunk < numpy.array((20, 20))), 'a small chunk size may freeze your computer'

        resolution = (bounds[2:] - bounds[:2]) / numpy.flip(shape)

        chunk_shape = numpy.array(numpy.round(chunk / resolution), numpy.int)

        input_points_sw = numpy.array(self.input_bounds[:2])

        output_x = numpy.linspace(bounds[0], bounds[2], shape[1])
        output_y = numpy.linspace(bounds[1], bounds[3], shape[0])

        interpolated_values = numpy.full(shape, fill_value=numpy.nan)
        interpolated_uncertainty = numpy.full(shape, fill_value=0)
        interpolated_variance = numpy.full(shape, fill_value=numpy.nan)

        chunk_grid_shape = numpy.array(numpy.floor(shape / chunk_shape), numpy.int)

        mask = raster_mask(self.interpolation_region, shape, resolution, input_points_sw, nodata, self.crs_wkt)
        mask_band = mask.GetRasterBand(1)
        mask = numpy.flip(mask_band.ReadAsArray(), axis=0) != nodata
        del mask_band

        expansion_indices = (expansion * chunk_shape / 2).astype(int)

        print(f'dividing output grid into {chunk_grid_shape * multiplier} chunks; ' +
              f'{chunk_shape} chunk with {expansion * chunk_shape} catchment area')

        with futures.ProcessPoolExecutor() as concurrency_pool:
            running_futures = {}
            running_uncertainty_futures = {}

            # determine indices of the chunk bounds within the output grid, expanding the interpolation area past the edges of the chunk
            for chunk_row in range((chunk_grid_shape[0] + 1) * multiplier):
                start_row = int(chunk_row * chunk_shape[0] / multiplier)
                end_row = start_row + chunk_shape[0]

                interpolation_south = input_points_sw[1] + (start_row - expansion_indices[0]) * resolution[1]
                interpolation_north = input_points_sw[1] + (end_row + expansion_indices[0]) * resolution[1]

                row_indices = (self.points[:, 1] >= interpolation_south) & (self.points[:, 1] <= interpolation_north)

                row_points = self.points[row_indices]

                if self.is_raster:
                    row_uncertainty_points = self.uncertainty_points[row_indices]

                for chunk_col in range((chunk_grid_shape[1] + 1) * multiplier):
                    start_col = int(chunk_col * chunk_shape[1] / multiplier)
                    end_col = start_col + chunk_shape[1]

                    interpolation_west = input_points_sw[0] + (start_col - expansion_indices[1]) * resolution[0]
                    interpolation_east = input_points_sw[0] + (end_col + expansion_indices[1]) * resolution[0]

                    col_indices = (row_points[:, 0] >= interpolation_west) & (row_points[:, 0] <= interpolation_east)

                    interpolation_points = row_points[col_indices]

                    # only interpolate if there are points
                    if interpolation_points.shape[0] >= 3:
                        grid_slice = slice(start_row, end_row), slice(start_col, end_col)
                        current_future = concurrency_pool.submit(_krige_points_onto_grid, interpolation_points, output_x[grid_slice[1]],
                                                                 output_y[grid_slice[0]], mask=mask[grid_slice])
                        running_futures[current_future] = grid_slice

                        if self.is_raster:
                            uncertainty_points = row_uncertainty_points[col_indices]
                            current_uncertainty_future = concurrency_pool.submit(_krige_points_onto_grid, uncertainty_points,
                                                                                 output_x[grid_slice[1]], output_y[grid_slice[0]],
                                                                                 mask=mask[grid_slice])
                            running_uncertainty_futures[current_uncertainty_future] = grid_slice

            for completed_future in futures.as_completed(running_futures):
                grid_slice = running_futures[completed_future]

                try:
                    chunk_interpolated_values, chunk_interpolated_variance = completed_future.result()
                    if type(chunk_interpolated_values) is numpy.ma.MaskedArray:
                        chunk_interpolated_values = chunk_interpolated_values.filled(nodata)
                    if type(chunk_interpolated_variance) is numpy.ma.MaskedArray:
                        chunk_interpolated_variance = chunk_interpolated_variance.filled(nodata)
                    interpolated_values[grid_slice] = chunk_interpolated_values
                    interpolated_variance[grid_slice] = chunk_interpolated_variance
                except ValueError as error:
                    print(f'malformed slice of {shape}: {grid_slice} ({error})')

            if self.is_raster:
                for completed_future in futures.as_completed(running_uncertainty_futures):
                    grid_slice = running_uncertainty_futures[completed_future]

                    try:
                        chunk_interpolated_uncertainty, _ = completed_future.result()
                        if type(chunk_interpolated_uncertainty) is numpy.ma.MaskedArray:
                            chunk_interpolated_uncertainty = chunk_interpolated_uncertainty.filled(nodata)
                        interpolated_uncertainty[grid_slice] = chunk_interpolated_uncertainty
                    except ValueError as error:
                        print(f'malformed slice of {shape}: {grid_slice} ({error})')

        interpolated_uncertainty = numpy.sqrt(interpolated_variance + interpolated_uncertainty ** 2)
        del interpolated_variance

        interpolated_values[numpy.isnan(interpolated_values)] = nodata
        interpolated_uncertainty[numpy.isnan(interpolated_uncertainty)] = nodata

        output_raster = gdal.GetDriverByName('MEM').Create('', int(shape[1]), int(shape[0]), 2, gdal.GDT_Float32)
        output_raster.SetGeoTransform((bounds[0], resolution[0], 0.0, bounds[1], 0.0, resolution[1]))
        output_raster.SetProjection(self.crs_wkt)

        band_1 = output_raster.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_values)
        del band_1

        band_2 = output_raster.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(interpolated_uncertainty)
        del band_2

        return output_raster

    def __uncertainty(self, values: numpy.array, catzoc: str) -> numpy.array:
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

    def plot(self, raster: gdal.Dataset, method: str, nodata: float = None, band: int = 1, show: bool = False):
        """
        Plot preinterpolated points and an interpolated raster side-by-side on synchronized subplots.

        Parameters
        ----------
        raster
            array (or GDAL raster) of interpolated data
        method
            method of interpolation
        nodata
            value for no data in given data
        band
            raster band (1-indexed)
        """

        raster_band = raster.GetRasterBand(band)

        if nodata is None:
            nodata = raster_band.GetNoDataValue()

        raster_data = raster_band.ReadAsArray()
        input_data = self.dataset.GetRasterBand(band).ReadAsArray() if self.is_raster else self.points[:, 2]

        # replace `nodata` values with NaN
        input_data[input_data == nodata] = numpy.nan
        raster_data[raster_data == nodata] = numpy.nan

        # get minimum and maximum values for all three dimensions
        z_values = numpy.concatenate((numpy.ravel(raster_data), numpy.ravel(input_data)))
        min_z = numpy.nanmin(z_values)
        max_z = numpy.nanmax(z_values)

        # create a new figure window with two subplots
        figure = pyplot.figure()
        left_axis = figure.add_subplot(1, 2, 1)
        left_axis.set_title('survey data')
        right_axis = figure.add_subplot(1, 2, 2, sharex=left_axis, sharey=left_axis)
        right_axis.set_title(f'{method} interpolation to raster')

        # plot data
        if self.is_raster:
            plot_raster(self.dataset, self.index, left_axis, vmin=min_z, vmax=max_z)
        else:
            left_axis.scatter(self.points[:, 0], self.points[:, 1], c=self.points[:, 2], s=1, vmin=min_z, vmax=max_z)

        plot_raster(raster, band, right_axis, vmin=min_z, vmax=max_z)
        right_axis.axes.get_yaxis().set_visible(False)

        # create colorbar
        figure.colorbar(ScalarMappable(norm=Normalize(vmin=min_z, vmax=max_z)), ax=(right_axis, left_axis))

        # pause program execution and show the figure
        if show:
            pyplot.show()


def _combined_coverage_within_window(raster_filenames: [str], origin: (float, float), resolution: (float, float),
                                     shape: (int, int)) -> numpy.array:
    """
    Return a boolean mask of the given rasters where data exists within the given window.

    Parameters
    ----------
    raster_filenames
        list of filenames of coverage rasters
    origin
        XY coordinates of northwest corner of output mask
    resolution
        XY resolution of output mask
    shape
        shape of output mask

    Returns
    -------
        boolean array indicating where data exists
    """

    if type(origin) is not numpy.array:
        origin = numpy.array(origin)

    if type(resolution) is not numpy.array:
        resolution = numpy.array(resolution)

    if type(shape) is not numpy.array:
        shape = numpy.array(shape)

    opposite_origin = origin + (resolution * numpy.flip(shape))
    bounds = bounds_from_opposite_corners(origin, opposite_origin)
    sw_corner = bounds[:2]
    ne_corner = bounds[2:]

    output_coverage_masks = []
    for raster_filename in raster_filenames:
        try:
            with rasterio.open(raster_filename) as raster:
                coverage_transform = raster.transform
                coverage_shape = raster.shape

            coverage_resolution = numpy.array((coverage_transform.a, coverage_transform.e))
            coverage_origin = numpy.array((coverage_transform.c, coverage_transform.f))
            coverage_bounds = bounds_from_opposite_corners(coverage_origin,
                                                           coverage_origin + (coverage_resolution * numpy.flip(coverage_shape)))
            coverage_sw_corner = coverage_bounds[:2]
            coverage_ne_corner = coverage_bounds[2:]

            # ensure the coverage array intersects the given window
            if numpy.all((coverage_sw_corner < ne_corner) & (coverage_ne_corner > sw_corner)):
                with rasterio.open(raster_filename) as raster:
                    coverage_data = raster.read()

                # nodata for sidescan coverage is usually the maximum value
                coverage_mask = numpy.where(array_coverage(coverage_data, numpy.max(coverage_data)), 1, 0)

                # resample coverage data to match BAG resolution
                resolution_ratio = coverage_resolution / resolution
                if numpy.any(resolution_ratio != 1):
                    coverage_mask = zoom(coverage_mask, zoom=resolution_ratio, order=3, prefilter=False)

                origin_index_delta = numpy.round((origin - coverage_origin) / resolution).astype(int)
                opposite_origin_index_delta = numpy.round((opposite_origin - coverage_origin) / resolution).astype(int)

                # indices to be written onto the output array
                output_array_index_slices = [slice(0, None), slice(0, None)]

                # BAG leftmost X is to the left of coverage leftmost X
                if origin_index_delta[0] < 0:
                    output_array_index_slices[1] = slice(origin_index_delta[0] * -1, output_array_index_slices[1].stop)
                    origin_index_delta[0] = 0

                # BAG topmost Y is above coverage topmost Y
                if origin_index_delta[1] < 0:
                    output_array_index_slices[0] = slice(origin_index_delta[1] * -1, output_array_index_slices[0].stop)
                    origin_index_delta[1] = 0

                # BAG rightmost X is to the right of coverage rightmost X
                if opposite_origin_index_delta[0] > coverage_mask.shape[1]:
                    output_array_index_slices[1] = slice(output_array_index_slices[1].start,
                                                         coverage_mask.shape[1] - opposite_origin_index_delta[0])
                    opposite_origin_index_delta[0] = coverage_mask.shape[1]

                # BAG bottommost Y is lower than coverage bottommost Y
                if opposite_origin_index_delta[1] > coverage_mask.shape[0]:
                    output_array_index_slices[0] = slice(output_array_index_slices[0].start,
                                                         coverage_mask.shape[0] - opposite_origin_index_delta[1])
                    opposite_origin_index_delta[1] = coverage_mask.shape[0]

                # write the relevant coverage data to a slice of the output array corresponding to the coverage extent
                output_array = numpy.full(shape, 0)
                output_array[output_array_index_slices[0], output_array_index_slices[1]] = coverage_mask[origin_index_delta[1]:
                                                                                                         opposite_origin_index_delta[1],
                                                                                           origin_index_delta[0]:
                                                                                           opposite_origin_index_delta[0]]
                output_coverage_masks.append(output_array)
        except:
            print(f'error opening file {raster_filename}')

    # collapse coverage masks into single raster array
    return numpy.logical_or.reduce(output_coverage_masks)


def _krige_points_onto_grid(points: numpy.array, x: [float], y: [float], **kwargs) -> (numpy.array, numpy.array):
    """
    Perform kriging (using PyKrige) from the given points onto a regular grid.

    Parameters
    ----------
    points
        N x 3 array of points to interpolate
    x
        X coordinates of output grid
    y
        Y coordinates of output grid

    Returns
    -------
        interpolated values and uncertainty
    """

    interpolator = OrdinaryKriging(points[:, 0], points[:, 1], points[:, 2], variogram_model='linear')
    return interpolator.execute('grid', x, y, **kwargs)
