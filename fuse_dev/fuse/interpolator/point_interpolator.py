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
import os
from concurrent import futures
from typing import Tuple, Union

import matplotlib.pyplot as plt
import numpy
import scipy
import scipy.spatial.distance
from osgeo import gdal, ogr, osr
from pykrige.ok import OrdinaryKriging


class PointInterpolator:
    """Interpolation methods for creating a raster from points."""

    def __init__(self, window_scalar: float = 1.1, nodata: float = 1000000.0):
        """
        Set some of the precondition, but make them over writable.
        """

        self.window_scalar = window_scalar
        self.nodata = nodata

    def interpolate(self, gdal_points: gdal.Dataset, method: str, output_resolution: float,
                    vector_file_path: str = None, shrink: bool = True) -> gdal.Dataset:
        """
        Interpolate the provided dataset.

        Currently this is assumed to be a gdal dataset.  At some point perhaps
        this should become something more native to python.

        Parameters
        ----------
        gdal_points
            gdal point cloud dataset
        method
            method for interpolation
        output_resolution
            resolution of output grid
        vector_file_path
            path to file containing vector data
        shrink
            whether to skrink interpolated data back to original extent

        Returns
        -------
        gdal.Dataset
            raster dataset of interpolated data

        """

        if method == 'invdist_scilin' and vector_file_path is None:
            raise ValueError('Supporting shapefile required')

        input_points = gdal_points_to_array(gdal_points)
        _, _, max_spacing = _point_spacing(input_points)
        window_size = max_spacing * self.window_scalar

        if method == 'linear':
            # linear interpolation using triangulation
            interpolated_dataset = self.__interp_points_gdal_linear(gdal_points, output_resolution)

            # trim interpolation back using the original dataset as a mask
            output_dataset = _mask_raster(interpolated_dataset,
                                          self.__get_mask(gdal_points, output_resolution, window_size))
        elif method == 'invdist':
            # inverse distance interpolation
            interpolated_dataset = self.__interp_points_gdal_invdist(gdal_points, output_resolution, window_size)

            # shrink the coverage back on the edges
            if shrink:
                output_dataset = _shrink_coverage(interpolated_dataset, output_resolution, window_size)
            else:
                output_dataset = interpolated_dataset
        elif method == 'invdist_scilin':
            interpolated_dataset = self.__interp_points_gdal_invdist_scipy_linear(gdal_points, output_resolution,
                                                                                  window_size)
            # shrink the coverage back on the edges
            if shrink:
                output_dataset = _shrink_coverage(interpolated_dataset, output_resolution, window_size)
            else:
                output_dataset = interpolated_dataset

            # mask raster to extent from vector file
            output_dataset = _mask_raster(output_dataset,
                                          rasterize_like(vector_file_path, interpolated_dataset, output_resolution))
        elif method == 'kriging':
            # interpolate using Ordinary Kriging
            interpolated_dataset = self.__interp_points_pykrige_kriging(gdal_points, output_resolution, window_size)

            # shrink the coverage back on the edges
            if shrink:
                output_dataset = _shrink_coverage(interpolated_dataset, output_resolution, window_size)
            else:
                output_dataset = interpolated_dataset
        else:
            raise ValueError(f'Interpolation type "{method}" not recognized.')

        self.__plot(input_points, output_dataset, method, nodata=1000000.0)

        return output_dataset

    def __get_mask(self, dataset: gdal.Dataset, resolution: float, window: float) -> gdal.Dataset:
        """
        Currently using the shrunk invdist method as a mask.

        Parameters
        ----------
        dataset
            TODO write description
        resolution
            TODO write description
        window
            TODO write description

        Returns
        -------
        gdal.Dataset
            mask
        """

        return _shrink_coverage(self.__interp_points_gdal_invdist(dataset, resolution, window), resolution, window)

    def __interp_points_gdal_linear(self, gdal_points: gdal.Dataset, resolution: float,
                                    nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating linearly.

        Parameters
        ----------
        gdal_points
            GDAL point cloud dataset
        resolution
            cell size of output grid
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid

        """

        if nodata is None:
            nodata = self.nodata

        input_points = gdal_points_to_array(gdal_points)
        (rows, cols), bounds = _shape_from_cell_size(resolution, _bounds_from_points(input_points))
        algorithm = f"linear:radius=0:nodata={int(nodata)}"
        return gdal.Grid('', gdal_points, format='MEM', width=cols, height=rows, outputBounds=bounds,
                         algorithm=algorithm)

    def __interp_points_gdal_invdist(self, gdal_points: gdal.Dataset, resolution: float, radius: float,
                                     nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance.

        Parameters
        ----------
        gdal_points
            GDAL point cloud dataset
        resolution
            cell size of output grid
        radius
            size of interpolation window
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if nodata is None:
            nodata = self.nodata

        input_points = gdal_points_to_array(gdal_points)
        (rows, cols), bounds = _shape_from_cell_size(resolution, _bounds_from_points(input_points))
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        return gdal.Grid('', gdal_points, format='MEM', width=cols, height=rows, outputBounds=bounds,
                         algorithm=algorithm)

    def __interp_points_gdal_invdist_scipy_linear(self, gdal_points: gdal.Dataset, resolution: float, radius: float,
                                                  nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance and then linearly (using SciPy).

        Parameters
        ----------
        gdal_points
            GDAL point cloud dataset
        resolution
            cell size of output grid
        radius
            size of interpolation window
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if nodata is None:
            nodata = self.nodata

        # interpolate using inverse distance in GDAL
        interpolated_raster = self.__interp_points_gdal_invdist(gdal_points, resolution, radius, nodata)
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
        band_1.SetNoDataValue(float(nodata))
        band_1.WriteArray(interpolated_data)
        del band_1

        return interpolated_raster

    def __interp_points_pykrige_kriging(self, gdal_points: gdal.Dataset, resolution: float,
                                        nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via kriging (using PyKrige).

        Parameters
        ----------
        gdal_points
            GDAL point cloud dataset
        resolution
            cell size of output grid
        nodata
            value for no data in output grid

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if nodata is None:
            nodata = self.nodata

        input_points = gdal_points_to_array(gdal_points)
        (rows, cols), bounds = _shape_from_cell_size(resolution, _bounds_from_points(input_points))

        output_x = numpy.arange(bounds[0], bounds[2], resolution)
        output_y = numpy.arange(bounds[1], bounds[3], resolution)
        interpolated_array = numpy.empty((len(output_y), len(output_x)), dtype=float)
        variance = numpy.empty((len(output_y), len(output_x)), dtype=float)

        side_length = 5
        total_parts = side_length ** 2

        row, col = 0, 0

        min_x, min_y, max_x, max_y = bounds

        x_range = max_x - min_x
        y_range = max_y - min_y

        start_time = datetime.datetime.now()

        with futures.ProcessPoolExecutor() as concurrency_pool:
            running_futures = {}

            for part_index in range(total_parts):
                print(f'processing chunk {part_index + 1} of {total_parts}')

                grid_x_start = int(col * len(output_x) / side_length)
                grid_x_end = int((col + 1) * len(output_x) / side_length - 1)
                grid_y_start = int(row * len(output_y) / side_length)
                grid_y_end = int((row + 1) * len(output_y) / side_length - 1)

                chunk_x_min = min_x + ((x_range / side_length) * part_index)
                chunk_x_max = min_x + ((x_range / side_length) * (part_index + 1))
                chunk_y_min = min_y + ((y_range / side_length) * part_index)
                chunk_y_max = min_y + ((y_range / side_length) * (part_index + 1))

                chunk_points = input_points[numpy.where(
                    (input_points[:, 0] >= chunk_x_min) & (input_points[:, 0] < chunk_x_max) & (
                            input_points[:, 1] >= chunk_y_min) & (input_points[:, 1] < chunk_y_max))[0]]

                print(f'found {chunk_points.shape[0]} points in chunk')

                if chunk_points.shape[0] > 0:
                    interpolator = OrdinaryKriging(chunk_points[:, 0], chunk_points[:, 1], chunk_points[:, 2],
                                                   variogram_model='linear', verbose=False, enable_plotting=False)

                    current_future = concurrency_pool.submit(interpolator.execute, 'grid',
                                                             output_x[grid_x_start:grid_x_end],
                                                             output_y[grid_y_start:grid_y_end])
                    running_futures[current_future] = [slice(grid_y_start, grid_y_end),
                                                       slice(grid_x_start, grid_x_end)]
                else:
                    interpolated_array[grid_y_start:grid_y_end, grid_x_start:grid_x_end] = numpy.full(
                        (grid_y_end - grid_y_start, grid_x_end - grid_x_start), numpy.nan)
                    variance[grid_y_start:grid_y_end, grid_x_start:grid_x_end] = numpy.full(
                        (grid_y_end - grid_y_start, grid_x_end - grid_x_start), numpy.nan)

                if col >= side_length - 1:
                    col = 0
                    row += 1
                else:
                    col += 1

            for completed_future in futures.as_completed(running_futures):
                chunk_slices = running_futures[completed_future]

                current_interpolated_data, current_variance = completed_future.result()

                interpolated_array[chunk_slices] = current_interpolated_data
                variance[chunk_slices] = current_variance

        duration = (datetime.datetime.now() - start_time)
        print(f'interpolating {total_parts} chunks took {duration.total_seconds()} s')

        uncertainty = numpy.sqrt(variance) * 2.5

        memory_raster_driver = gdal.GetDriverByName('MEM')
        output_dataset = memory_raster_driver.Create('temp', cols, rows, 2, gdal.GDT_Float32)

        input_geotransform = gdal_points.GetGeoTransform()
        output_dataset.SetGeoTransform((input_geotransform[0], resolution, 0.0, 0.0, input_geotransform[3], resolution))
        output_dataset.SetProjection(gdal_points.GetProjection())

        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_array)

        band_2 = output_dataset.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(uncertainty)

        del band_1, band_2
        return output_dataset

    def __plot(self, survey_points: numpy.array, interpolated_raster: Union[numpy.array, gdal.Dataset],
               interpolation_method: str, nodata: float = None):
        """

        Parameters
        ----------
        survey_points
            array of XYZ points from original survey
        interpolated_raster
            array of interpolated data
        interpolation_method
            method of interpolation
        """

        if nodata is None:
            nodata = self.nodata

        if type(interpolated_raster) is gdal.Dataset:
            interpolated_raster = numpy.flip(interpolated_raster.ReadAsArray(), axis=0)
        if len(interpolated_raster.shape) == 3:
            interpolated_raster = interpolated_raster[0, :, :]
        interpolated_raster[interpolated_raster == nodata] = numpy.nan

        min_x, min_y, max_x, max_y = _bounds_from_points(survey_points)

        figure = plt.figure()

        original_data_axis = figure.add_subplot(1, 2, 1)
        original_data_axis.set_title('survey points')
        original_data_axis.scatter(survey_points[:, 0], survey_points[:, 1], c=survey_points[:, 2], s=1)
        original_data_axis.set_xlim(min_x, max_x)
        original_data_axis.set_ylim(min_y, max_y)

        interpolated_data_axis = figure.add_subplot(1, 2, 2)
        interpolated_data_axis.set_title(f'{interpolation_method} interpolation to raster')
        interpolated_data_axis.imshow(interpolated_raster, extent=(min_x, max_x, min_y, max_y), aspect='auto')

        plt.show()


def _shrink_coverage(dataset: gdal.Dataset, resolution: float, radius: float) -> gdal.Dataset:
    """
    Shrink coverage of a dataset by the original coverage radius.

    Parameters
    ----------
    dataset
        interpolated raster
    resolution
        cell size of raster
    radius
        radius to shrink by

    Returns
    -------
    gdal.Dataset
        shrunken dataset
    """

    for band_index in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(1)
        nodata_value = band.GetNoDataValue()
        band_data = band.ReadAsArray()

        # convert nodata values to NaN
        band_data[band_data == nodata_value] = numpy.nan

        # divide the window size by the resolution to get the number of cells by which to shrink the coverage
        for _ in numpy.arange(int(numpy.round(radius / resolution))):
            ew_diffs_nan_indices = numpy.where(numpy.isnan(numpy.diff(band_data, axis=0)))
            ns_diffs_nan_indices = numpy.where(numpy.isnan(numpy.diff(band_data, axis=1)))
            band_data[ew_diffs_nan_indices] = numpy.nan
            band_data[ew_diffs_nan_indices[0] + 1, ew_diffs_nan_indices[1]] = numpy.nan
            band_data[ns_diffs_nan_indices] = numpy.nan
            band_data[ns_diffs_nan_indices[0], ns_diffs_nan_indices[1] + 1] = numpy.nan

        # convert NaN to nodata values
        band_data[numpy.isnan(band_data)] = nodata_value
        band.WriteArray(band_data)

    return dataset


def _mask_raster(raster: gdal.Dataset, raster_mask: gdal.Dataset) -> gdal.Dataset:
    """
    Mask a raster using the nodata values in the given mask.
    Both rasters are assumed to be collocated, as all operations are conducted on the pixel level.

    Parameters
    ----------
    raster
        raster to apply mask to
    raster_mask
        mask to apply to the raster

    Returns
    -------
    gdal.Dataset
        masked raster dataset
    """

    input_data = raster.ReadAsArray()
    raster_band = raster.GetRasterBand(1)
    mask_band = raster_mask.GetRasterBand(1)
    raster_nodata_value = raster_band.GetNoDataValue()
    mask_nodata_value = mask_band.GetNoDataValue()

    input_data[raster_mask.ReadAsArray() == mask_nodata_value] = raster_nodata_value
    raster_band.WriteArray(input_data)
    return raster


def _point_spacing(points: numpy.array) -> Tuple[float, float, float]:
    """
    Take a numpy xyz array and return the min, mean, and max spacing
    between different points in the XY direction for interpolation without
    holes. The returned units will be the same as the provided horizontal
    coordinate system.

    Parameters
    ----------
    points
        array of XY or XYZ points

    Returns
    -------
    Tuple[float, float, float]
        min, mean, and max distances between neighboring points
    """

    pairwise_distances = scipy.spatial.distance.pdist(points, 'euclidean')
    return numpy.min(pairwise_distances), numpy.mean(pairwise_distances), numpy.max(pairwise_distances)


def _bounds_from_points(points: numpy.array) -> Tuple[float, float, float, float]:
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


def _shape_from_cell_size(resolution: float, bounds: Tuple[float, float, float, float]) -> Tuple[
    Tuple[int, int], Tuple[float, float, float, float]]:
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


def _vector_file_to_gdal_dataset(vector_file_path: str, spatial_reference_wkt: str,
                                 geotransform: Tuple[float, float, float, float, float, float], resolution: int,
                                 shape: Tuple[int, int], nodata: float = 0) -> gdal.Dataset:
    """
    Import shapefile

    Parameters
    ----------
    vector_file_path
        path to file containing vector data
    spatial_reference_wkt
        well-known text of a spatial reference system
    geotransform
        tuple of the desired GDAL geotransform of the dataset
    resolution
        cell size of the output dataset
    shape
        shape of the the output dataset
    nodata
        value to be substituted for no data in the output dataset

    Returns
    -------
    gdal.Dataset
        raster dataset

    """

    print('getShpRast', vector_file_path)
    fName = os.path.split(vector_file_path)[-1]
    splits = os.path.splitext(fName)
    name = splits[0]
    # tif = f'{splits[0]}.tif'

    # Open the data source and read in the extent
    source_ds = ogr.Open(vector_file_path)
    source_layer = source_ds.GetLayer()
    source_srs = source_layer.GetSpatialRef()

    for input_feature in source_layer:
        if input_feature is not None:
            input_geometry = input_feature.GetGeometryRef()
            #                print(geom.ExportToWkt())
            output_geometry = ogr.CreateGeometryFromWkt(input_geometry.ExportToWkt())
            #                print(source_srs, to_proj, sep='\n')
            coordTrans = osr.CoordinateTransformation(source_srs, spatial_reference_wkt)
            output_geometry.Transform(coordTrans)
            driver = ogr.GetDriverByName('Memory')
            ds = driver.CreateDataSource('temp')
            layer = ds.CreateLayer(name, spatial_reference_wkt, ogr.wkbMultiPolygon)

            # Add one attribute
            layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
            defn = layer.GetLayerDefn()

            # Create a new feature (attribute and geometry)
            output_feature = ogr.Feature(defn)
            output_feature.SetField('id', 123)

            # Make a geometry, from Shapely object
            output_feature.SetGeometry(output_geometry)

            layer.CreateFeature(output_feature)
            del output_feature, input_geometry
            break

    x_min, x_max, y_min, y_max = output_geometry.GetEnvelope()
    meta = ([x_min, y_max], [x_max, y_min])
    print(meta)

    # Create the destination data source
    x_dim = int((x_max - x_min) / resolution)
    y_dim = int((y_max - y_min) / resolution)
    print(x_dim, y_dim)
    target_ds = gdal.GetDriverByName('MEM').Create('', x_dim, y_dim, gdal.GDT_Byte)
    x_orig, y_orig = geotransform[0], geotransform[3]
    target_gt = (x_orig, resolution, 0, y_orig, 0, resolution)
    target_ds.SetGeoTransform(target_gt)
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])
    newarr = band.ReadAsArray()

    # clip resampled coverage data to the bounds of the BAG
    output_array = numpy.full(shape, nodata)

    cov_ul, cov_lr = numpy.array(x_orig), numpy.array(y_orig)
    bag_ul, bag_lr = numpy.array(x_min), numpy.array(y_min)

    if bag_ul[0] > cov_lr[0] or bag_lr[0] < cov_ul[0] or bag_lr[1] > cov_ul[1] or bag_ul[1] < cov_lr[1]:
        raise ValueError('bag dataset is outside the bounds of coverage dataset')

    ul_index_delta = numpy.round((bag_ul - cov_ul) / numpy.array(resolution)).astype(int)
    lr_index_delta = numpy.round((bag_lr - cov_ul) / numpy.array(resolution)).astype(int)

    # indices to be written onto the output array
    output_array_index_slices = [slice(0, None), slice(0, None)]

    # BAG leftmost X is to the left of coverage leftmost X
    if ul_index_delta[0] < 0:
        output_array_index_slices[1] = slice(ul_index_delta[0] * -1, output_array_index_slices[1].stop)
        ul_index_delta[0] = 0

    # BAG topmost Y is above coverage topmost Y
    if ul_index_delta[1] < 0:
        output_array_index_slices[0] = slice(ul_index_delta[1] * -1, output_array_index_slices[0].stop)
        ul_index_delta[1] = 0

    # BAG rightmost X is to the right of coverage rightmost X
    if lr_index_delta[0] > newarr.shape[1]:
        output_array_index_slices[1] = slice(output_array_index_slices[1].start,
                                             newarr.shape[1] - lr_index_delta[0])
        lr_index_delta[0] = newarr.shape[1]

    # BAG bottommost Y is lower than coverage bottommost Y
    if lr_index_delta[1] > newarr.shape[0]:
        output_array_index_slices[0] = slice(output_array_index_slices[0].start,
                                             newarr.shape[0] - lr_index_delta[1])
        lr_index_delta[1] = newarr.shape[0]

    # write the relevant coverage data to a slice of the output array corresponding to the coverage extent
    output_array[output_array_index_slices[0], output_array_index_slices[1]] = newarr[
                                                                               ul_index_delta[1]:lr_index_delta[1],
                                                                               ul_index_delta[0]:lr_index_delta[0]]

    del newarr, band, source_ds, target_ds

    ax_y, ax_x = output_array.shape

    target_ds = gdal.GetDriverByName('MEM').Create('', ax_x, ax_y, 1, gdal.GDT_Float32)
    x_orig, y_orig = geotransform[0], geotransform[3]
    target_gt = (x_orig, resolution, 0, y_orig, 0, resolution)
    target_ds.SetGeoTransform(target_gt)
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(output_array)

    return target_ds


def rasterize_like(vector_file_path: str, like_raster: gdal.Dataset, resolution: int) -> gdal.Dataset:
    """
    Burn given shapefile onto a raster with the same size and location as the given raster.

    Parameters
    ----------
    vector_file_path
        path to file containing vector data
    like_raster
        GDAL raster dataset
    resolution
        cell size of grid

    Returns
    -------
    gdal.Dataset
        raster with vector data burned in
    """

    geotransform = like_raster.GetGeoTransform()
    spatial_reference_wkt = osr.SpatialReference(wkt=like_raster.GetProjectionRef())
    spatial_reference_wkt.MorphFromESRI()

    return _vector_file_to_gdal_dataset(vector_file_path, spatial_reference_wkt, geotransform, int(resolution),
                                        (like_raster.RasterYSize, like_raster.RasterXSize))


def gdal_points_to_array(gdal_points: gdal.Dataset, layer_index: int = 0) -> numpy.array:
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

    point_layer = gdal_points.GetLayerByIndex(layer_index)
    num_points = point_layer.GetFeatureCount()
    output_points = numpy.empty((num_points, 3))

    for point_index in range(num_points):
        feature = point_layer.GetFeature(point_index)
        output_points[point_index, :] = feature.geometry().GetPoint()

    return output_points
