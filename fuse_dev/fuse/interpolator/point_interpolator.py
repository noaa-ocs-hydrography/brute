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


def plot_interpolation(survey_points: numpy.array, interpolated_raster: Union[numpy.array, gdal.Dataset],
                       interpolation_method: str):
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

    bounds = get_bounds(survey_points)

    figure = plt.figure()
    figure.suptitle(interpolation_method)

    original_data_axis = figure.add_subplot(1, 2, 1)
    interpolated_data_axis = figure.add_subplot(1, 2, 2)

    original_data_axis.set_title('survey points')
    interpolated_data_axis.set_title('interpolated raster')

    original_data_axis.scatter(survey_points[:, 0], survey_points[:, 1], c=survey_points[:, 2], s=0.5)
    interpolated_data_axis.imshow(interpolated_raster, extent=(bounds[0], bounds[2], bounds[1], bounds[3]))

    plt.show()


def point_spacing(points: numpy.array) -> Tuple[float, float, float]:
    """
    Take a numpy xyz array and return the min, mean, and max spacing
    between different points in the XY direction for interpolation without
    holes.  The returned units will be the same as the provided horizontal
    coordinate system.

    Parameters
    ----------
    points
        array of XY or XYZ points

    Returns
    -------
    Tuple[float, float, float]
        minimum, mean, and maximum

    """

    pairwise_distances = scipy.spatial.distance.pdist(points, 'euclidean')
    return numpy.min(pairwise_distances), numpy.mean(pairwise_distances), numpy.max(pairwise_distances)


def get_bounds(points: numpy.array) -> Tuple[float, float, float, float]:
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


def gdal_points_to_array(point_cloud_dataset: gdal.Dataset, layer_index: int = 0) -> numpy.array:
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

    point_layer = point_cloud_dataset.GetLayerByIndex(layer_index)
    num_points = point_layer.GetFeatureCount()
    output_points = numpy.empty((num_points, 3))

    for point_index in range(num_points):
        feature = point_layer.GetFeature(point_index)
        output_points[point_index, :] = feature.geometry().GetPoint()

    return output_points


def _compare_vals(value: float, min_value: float, max_value: float) -> Tuple[float, float]:
    """
    This is a small utility for inspecting values and seeing they
    contribute to a min or max estimate.

    Parameters
    ----------
    value: float :
        TODO write description
    min_value: float :
        TODO write description
    max_value: float :
        TODO write description

    Returns
    -------
    Tuple[float, float]
        new minimum and maximum values

    """

    if numpy.isnan(min_value) or value < min_value:
        min_value = value

    if numpy.isnan(max_value) or value > max_value:
        max_value = value

    return min_value, max_value


def _getShpRast(shapefile_path: str, output_wkt: str, output_geotransform: tuple, output_resolution: int, output_shape,
                nodata: float = 0) -> Tuple[gdal.Dataset, tuple]:
    """
    Import shapefile

    Parameters
    ----------
    shapefile_path: str :
        Shapefile file location
    output_wkt: str:
        WKT string with destination spatial reference system
    output_geotransform: tuple:
        gdal.GeoTransform object of the interpolated dataset
    output_resolution: int :
        Resolution of the input dataset
    output_shape : tuple(int, int)
        Shape of the the input dataset
    nodata: float :
        Nodata value for the ouput gdal.Dataset object (Default value = 0)

    Returns
    -------
    Tuple[gdal.Dataset, tuple]
        dataset and geotransform

    """

    print('getShpRast', shapefile_path)
    fName = os.path.split(shapefile_path)[-1]
    splits = os.path.splitext(fName)
    name = splits[0]
    # tif = f'{splits[0]}.tif'

    # Open the data source and read in the extent
    source_ds = ogr.Open(shapefile_path)
    source_layer = source_ds.GetLayer()
    source_srs = source_layer.GetSpatialRef()

    for input_feature in source_layer:
        if input_feature is not None:
            input_geometry = input_feature.GetGeometryRef()
            #                print(geom.ExportToWkt())
            output_geometry = ogr.CreateGeometryFromWkt(input_geometry.ExportToWkt())
            #                print(source_srs, to_proj, sep='\n')
            coordTrans = osr.CoordinateTransformation(source_srs, output_wkt)
            output_geometry.Transform(coordTrans)
            driver = ogr.GetDriverByName('Memory')
            ds = driver.CreateDataSource('temp')
            layer = ds.CreateLayer(name, output_wkt, ogr.wkbMultiPolygon)

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
    x_dim = int((x_max - x_min) / output_resolution)
    y_dim = int((y_max - y_min) / output_resolution)
    print(x_dim, y_dim)
    target_ds = gdal.GetDriverByName('MEM').Create('', x_dim, y_dim, gdal.GDT_Byte)
    x_orig, y_orig = output_geotransform[0], output_geotransform[3]
    target_gt = (x_orig, output_resolution, 0, y_orig, 0, output_resolution)
    target_ds.SetGeoTransform(target_gt)
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])
    newarr = band.ReadAsArray()

    # clip resampled coverage data to the bounds of the BAG
    output_array = numpy.full(output_shape, nodata)

    cov_ul, cov_lr = numpy.array(x_orig), numpy.array(y_orig)
    bag_ul, bag_lr = numpy.array(x_min), numpy.array(y_min)

    if bag_ul[0] > cov_lr[0] or bag_lr[0] < cov_ul[0] or bag_lr[1] > cov_ul[1] or bag_ul[1] < cov_lr[1]:
        raise ValueError('bag dataset is outside the bounds of coverage dataset')

    ul_index_delta = numpy.round((bag_ul - cov_ul) / numpy.array(output_resolution)).astype(int)
    lr_index_delta = numpy.round((bag_lr - cov_ul) / numpy.array(output_resolution)).astype(int)

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
    x_orig, y_orig = output_geotransform[0], output_geotransform[3]
    target_gt = (x_orig, output_resolution, 0, y_orig, 0, output_resolution)
    target_ds.SetGeoTransform(target_gt)
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(output_array)

    target_gt = target_ds.GetGeoTransform()

    return target_ds, target_gt


def _get_shape_mask(self, grid: gdal.Grid, shapefile: str, resolution: float) -> gdal.Dataset:
    """
    TODO write description

    Parameters
    ----------
    grid: gdal.Grid :
        TODO write description
    shapefile: str :
        TODO write description
    resolution: float :
        TODO write description

    Returns
    -------
    type
        shape mask

    """

    band = grid.GetRasterBand(1)
    shape = (grid.RasterYSize, grid.RasterXSize)
    raster = band.ReadAsArray()
    del band
    grid_ref = grid.GetProjectionRef()
    grid_gt = grid.GetGeoTransform()
    proj = osr.SpatialReference(wkt=grid_ref)
    proj.MorphFromESRI()
    print(grid_ref, grid_gt)
    shape_ds, shape_gt = _getShpRast(shapefile, proj, grid_gt, int(resolution), shape)
    print('transformed', shape_gt)
    return shape_ds


class PointInterpolator:
    """Interpolation methods for creating a raster from points."""

    def __init__(self, window_scalar: float = 1.1):
        """
        Set some of the precondition, but make them over writable.
        """

        self.window_scalar = window_scalar

    def interpolate(self, point_cloud_dataset: gdal.Dataset, method: str, output_resolution: float,
                    shapefile_path: str = None, shrink: bool = True) -> gdal.Dataset:
        """
        Interpolate the provided dataset.

        Currently this is assumed to be a gdal dataset.  At some point perhaps
        this should become something more native to python.

        Parameters
        ----------
        point_cloud_dataset
            gdal point cloud dataset
        method
            method for interpolation
        output_resolution
            resolution of output grid
        shapefile_path
            path to shapefile
        shrink
            whether to skrink interpolated data back to original extent

        Returns
        -------
        gdal.Dataset
            raster dataset of interpolated data

        """

        if method == 'invdist_scilin' and shapefile_path is None:
            raise ValueError('Supporting shapefile required')

        input_points = gdal_points_to_array(point_cloud_dataset)
        _, _, max_spacing = point_spacing(input_points)
        window_size = max_spacing * self.window_scalar

        if method == 'linear':
            # linear interpolation using triangulation
            interpolated_dataset = self._interp_points_gdal_linear(point_cloud_dataset, output_resolution)

            # trim interpolation back using the original dataset as a mask
            output_dataset = self._mask_with_raster(interpolated_dataset,
                                                    self._get_mask(point_cloud_dataset, output_resolution, window_size))
        elif method == 'invlin':
            interpolated_dataset = self._interp_points_gdal_invdist_scipy_linear(point_cloud_dataset, output_resolution,
                                                                                 window_size)

            output_dataset = self._mask_with_raster(
                self._shrink_coverage(interpolated_dataset, output_resolution,
                                      window_size) if shrink else interpolated_dataset,
                _get_shape_mask(interpolated_dataset, shapefile_path, output_resolution))
        elif method == 'invdist':
            # inverse distance interpolation
            interpolated_dataset = self._interp_points_gdal_invdist(point_cloud_dataset, output_resolution, window_size)

            # shrink the coverage back on the edges
            output_dataset = self._shrink_coverage(interpolated_dataset, output_resolution,
                                                   window_size) if shrink else interpolated_dataset
        elif method == 'kriging':
            # interpolate using Ordinary Kriging
            interpolated_dataset = self._interp_points_kriging(point_cloud_dataset, output_resolution, window_size)

            # shrink the coverage back on the edges
            output_dataset = self._shrink_coverage(interpolated_dataset, output_resolution,
                                                   window_size) if shrink else interpolated_dataset
        else:
            raise ValueError(f'Interpolation type "{method}" not recognized.')

        return output_dataset

    def _get_mask(self, dataset: gdal.Dataset, resolution: float, window: float) -> gdal.Dataset:
        """
        Currently using the shrunk invdist method as a mask.

        Parameters
        ----------
        dataset: gdal.Dataset :
            TODO write description
        resolution: float :
            TODO write description
        window: float :
            TODO write description

        Returns
        -------
        gdal.Dataset
            mask

        """

        return self._shrink_coverage(self._interp_points_gdal_invdist(dataset, resolution, window), resolution, window)

    def _interp_points_gdal_linear(self, dataset: gdal.Dataset, resolution: float,
                                   nodata: float = 1000000) -> gdal.Grid:
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.

        Parameters
        ----------
        dataset: gdal.Dataset :
            TODO write description
        resolution: float :
            TODO write description
        nodata: float :
            TODO write description (Default value = 1000000)

        Returns
        -------
        type
            interpolated dataset

        """

        input_points = gdal_points_to_array(dataset)
        rows, cols, bounds = self._get_nodes3(resolution, get_bounds(input_points))

        algorithm = f"linear:radius=0:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=cols, height=rows, outputBounds=bounds,
                                algorithm=algorithm)

        interpolated_data = numpy.flip(interp_data.ReadAsArray(), axis=0)
        interpolated_data[interpolated_data == nodata] = numpy.nan

        plot_interpolation(input_points, interpolated_data, 'linear')

        return interp_data

    def _interp_points_gdal_invdist(self, dataset: gdal.Dataset, resolution: float, radius: float,
                                    nodata: float = 1000000) -> gdal.Dataset:
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.

        Parameters
        ----------
        dataset: gdal.Dataset :
            TODO write description
        resolution: float :
            TODO write description
        radius: float :
            TODO write description
        nodata: float :
            TODO write description(Default value = 1000000)

        Returns
        -------
        type
            interpolated grid

        """

        # Find the bounds of the provided data
        xmin, xmax, ymin, ymax = numpy.nan, numpy.nan, numpy.nan, numpy.nan
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()

        for n in numpy.arange(count):
            f = lyr.GetFeature(n)
            x, y, z = f.geometry().GetPoint()
            xmin, xmax = _compare_vals(x, xmin, xmax)
            ymin, ymax = _compare_vals(y, ymin, ymax)

        numrows, numcolumns, bounds = self._get_nodes3(resolution, (xmin, ymin, xmax, ymax))
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)
        return interp_data

    def _interp_points_gdal_invdist_scipy_linear(self, dataset: gdal.Dataset, resolution: float, radius: float,
                                                 nodata: float = 1000000) -> gdal.Dataset:
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.

        Parameters
        ----------
        dataset: gdal.Dataset :
            TODO write description
        resolution: float :
            TODO write description
        radius: float :
            TODO write description
        nodata: float :
            TODO write description (Default value = 1000000)

        Returns
        -------
        type
            interpolated dataset

        """

        print('_gdal_invdist_scilin_interp_points')
        # Find the bounds of the provided data
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()

        input_points = []

        for n in numpy.arange(count):
            f = lyr.GetFeature(n)
            input_points.append(f.geometry().GetPoint())

        input_points = numpy.array(input_points)
        xmin = numpy.min(input_points[:, 0])
        xmax = numpy.max(input_points[:, 0])
        ymin = numpy.min(input_points[:, 1])
        ymax = numpy.max(input_points[:, 1])

        numrows, numcolumns, bounds = self._get_nodes3(resolution, (xmin, ymin, xmax, ymax))
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)
        arr = interp_data.ReadAsArray()
        ycoord, xcoord = numpy.where(arr != nodata)
        zvals = arr[ycoord[:], xcoord[:]]
        del arr

        print('start')
        xa, ya = numpy.arange(numcolumns), numpy.arange(numrows)
        xi, yi = numpy.meshgrid(xa, ya)
        interp_grid = scipy.interpolate.griddata((xcoord, ycoord), zvals, (xi, yi), method='linear', fill_value=nodata)

        interp_grid[numpy.isnan(interp_grid)] = nodata

        plot_interpolation(input_points, interp_grid, 'IDW + linear')
        print('stop')

        band = interp_data.GetRasterBand(1)
        band.SetNoDataValue(float(nodata))
        band.WriteArray(interp_grid)
        del band

        return interp_data

    def _interp_points_kriging(self, point_cloud_dataset: gdal.Dataset, resolution: float, radius: float,
                               nodata: float = 1000000.0) -> gdal.Dataset:
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.
        """

        input_points = gdal_points_to_array(point_cloud_dataset)
        num_rows, num_columns, bounds = self._get_nodes3(resolution, get_bounds(input_points))

        output_grid_x = numpy.arange(bounds[0], bounds[2], resolution)
        output_grid_y = numpy.arange(bounds[1], bounds[3], resolution)
        interpolated_data = numpy.empty((len(output_grid_y), len(output_grid_x)), dtype=float)
        variance = numpy.empty((len(output_grid_y), len(output_grid_x)), dtype=float)

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

                grid_x_start = int(col * len(output_grid_x) / side_length)
                grid_x_end = int((col + 1) * len(output_grid_x) / side_length - 1)
                grid_y_start = int(row * len(output_grid_y) / side_length)
                grid_y_end = int((row + 1) * len(output_grid_y) / side_length - 1)

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
                                                             output_grid_x[grid_x_start:grid_x_end],
                                                             output_grid_y[grid_y_start:grid_y_end])
                    running_futures[current_future] = [slice(grid_y_start, grid_y_end),
                                                       slice(grid_x_start, grid_x_end)]
                else:
                    interpolated_data[grid_y_start:grid_y_end, grid_x_start:grid_x_end] = numpy.full(
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

                interpolated_data[chunk_slices] = current_interpolated_data
                variance[chunk_slices] = current_variance

        duration = (datetime.datetime.now() - start_time)
        print(f'interpolating {total_parts} chunks took {duration.total_seconds()} s')

        uncertainty = numpy.sqrt(variance) * 2.5

        plot_interpolation(input_points, interpolated_data, 'kriging')

        memory_raster_driver = gdal.GetDriverByName('MEM')
        output_dataset = memory_raster_driver.Create('temp', num_columns, num_rows, 2, gdal.GDT_Float32)

        input_geotransform = point_cloud_dataset.GetGeoTransform()
        output_dataset.SetGeoTransform((input_geotransform[0], resolution, 0.0, 0.0, input_geotransform[3], resolution))
        output_dataset.SetProjection(point_cloud_dataset.GetProjection())

        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_data)

        band_2 = output_dataset.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(uncertainty)

        del band_1, band_2
        return output_dataset

    def _get_nodes(self, resolution: float, bounds: Tuple[float, float, float, float]) -> Tuple[
        int, int, Tuple[float, float, float, float]]:
        """
        Get the bounds and number of rows and columns. The desired resolution
        and the data max and min in X and Y are required.  This algorithm uses
        the average of the data min and max and centers the grid on this
        location.

        Parameters
        ----------
        resolution: float :
            TODO write description
        bounds: Tuple[float, float, float, float] :
            TODO write description

        Returns
        -------
        type
            number of rows, number of columns, and bounds

        """

        xmin, ymin, xmax, ymax = bounds
        numcolumns = int(numpy.ceil((1. * (xmax - xmin) / resolution)))
        xmean = numpy.mean([xmin, xmax])
        xnmin = xmean - numcolumns / 2. * resolution
        xnmax = xmean + numcolumns / 2. * resolution
        numrows = int(numpy.ceil((1. * (ymax - ymin) / resolution)))
        ymean = numpy.mean([ymin, ymax])
        ynmin = ymean - numrows / 2. * resolution
        ynmax = ymean + numrows / 2. * resolution
        return numrows, numcolumns, (xnmin, ynmin, xnmax, ynmax)

    def _get_nodes2(self, resolution: float, bounds: Tuple[float, float, float, float]) -> Tuple[
        int, int, Tuple[float, float, float, float]]:
        """
        Get the bounds and number of rows and columns. The desired resolution
        and the data max and min in X and Y are required.  This algorithm use
        the data min as the anchor for the rows and columns, thus the max data
        values may contain sparse data.

        Parameters
        ----------
        resolution: float :
            TODO write description
        bounds: Tuple[float, float, float, float] :
            TODO write description

        Returns
        -------
        type
            number of rows, number of columns, and bounds

        """

        xmin, ymin, xmax, ymax = bounds
        numcolumns = int(numpy.ceil((1. * (xmax - xmin) / resolution)))
        xnmin = xmin
        xnmax = xmin + 1. * numcolumns * resolution
        numrows = int(numpy.ceil((1. * (ymax - ymin) / resolution)))
        ynmin = ymin
        ynmax = ymin + 1. * numrows * resolution
        return numrows, numcolumns, (xnmin, ynmin, xnmax, ynmax)

    def _get_nodes3(self, resolution: float, bounds: Tuple[float, float, float, float]) -> Tuple[
        int, int, Tuple[float, float, float, float]]:
        """
        Get the bounds and number of rows and columns. The desired resolution
        and the data max and min in X and Y are required.  This algorithm use
        the data min as the anchor for the rows and columns, but floored to the
        nearest multiple of the resolution.

        Parameters
        ----------
        resolution: float :

        bounds: Tuple[float, float, float, float] :


        Returns
        -------
        type
            number of rows, number of columns, and bounds

        """

        xmin, ymin, xmax, ymax = bounds
        numcolumns = int(numpy.ceil((1. * (xmax - xmin) / resolution)))
        xrem = xmin % resolution
        xnmin = xmin - xrem
        xnmax = xnmin + 1. * numcolumns * resolution
        numrows = int(numpy.ceil((1. * (ymax - ymin) / resolution)))
        yrem = ymin % resolution
        ynmin = ymin - yrem
        ynmax = ynmin + 1. * numrows * resolution
        return numrows, numcolumns, (xnmin, ynmin, xnmax, ynmax)

    def _mask_with_raster(self, dataset: gdal.Dataset, maskraster: gdal.Dataset) -> gdal.Dataset:
        """
        Read two rasters and use the nodata values from one to mask the other.
        These files are assumed to be collocated as all operation are conducted
        on the pixel level.

        Parameters
        ----------
        dataset: gdal.Dataset :

        maskraster: gdal.Dataset :


        Returns
        -------
        type
            masked dataset

        """

        data = dataset.ReadAsArray()
        data_rb = dataset.GetRasterBand(1)
        data_nd = data_rb.GetNoDataValue()
        mask = maskraster.ReadAsArray()
        mask_rb = maskraster.GetRasterBand(1)
        mask_nd = mask_rb.GetNoDataValue()
        idx = numpy.nonzero(mask == mask_nd)
        data[idx] = data_nd
        data_rb.WriteArray(data)
        return dataset

    def _shrink_coverage(self, dataset: gdal.Dataset, resolution: float, radius: float) -> gdal.Dataset:
        """
        Shrink coverage of a dataset by the original coverage radius.

        Parameters
        ----------
        dataset: gdal.Dataset :

        resolution: float :

        radius: float :


        Returns
        -------
        type
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
