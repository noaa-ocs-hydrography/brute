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
from concurrent import futures
from typing import Tuple, Union

import matplotlib.pyplot as plt
import numpy
import scipy
import scipy.interpolate
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

    def interpolate(self, points: gdal.Dataset, method: str, output_resolution: float, vector_file_path: str = None,
                    shrink: bool = True, plot: bool = False) -> gdal.Dataset:
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
        shrink
            whether to skrink interpolated data back to original extent
        plot
            whether to plot the interpolated result next to the original survey points

        Returns
        -------
        gdal.Dataset
            raster dataset of interpolated data
        """

        if method == 'invdist_scilin' and vector_file_path is None:
            raise ValueError('Supporting shapefile required')

        self.input_points = gdal_points_to_array(points)
        self.input_shape, self.input_bounds = _shape_from_cell_size(output_resolution,
                                                                    _bounds_from_points(self.input_points))
        self.input_spatial_reference = points.GetProjection()

        if self.input_spatial_reference == '':
            self.input_spatial_reference = osr.SpatialReference()
            self.input_spatial_reference.ImportFromEPSG(3857)

        _, _, max_spacing = _point_spacing(self.input_points)
        window_size = max_spacing * self.window_scalar

        start_time = datetime.datetime.now()

        if method == 'linear':
            # linear interpolation using triangulation
            interpolated_dataset = self.__interp_points_gdal_linear(points)

            # trim interpolation back using the original dataset as a mask
            mask = self.__get_mask((interpolated_dataset.RasterYSize, interpolated_dataset.RasterXSize),
                                   output_resolution, (self.input_bounds[0], self.input_bounds[3]),
                                   self.input_spatial_reference, self.nodata)
            output_dataset = _mask_raster(interpolated_dataset, mask)
        elif method == 'invdist':
            # inverse distance interpolation
            interpolated_dataset = self.__interp_points_gdal_invdist(points, window_size)

            # shrink the coverage back on the edges
            if shrink:
                output_dataset = _shrink_coverage(interpolated_dataset, output_resolution, window_size)
            else:
                output_dataset = interpolated_dataset
        elif method == 'invdist_scilin':
            interpolated_dataset = self.__interp_points_gdal_invdist_scipy_linear(points, output_resolution,
                                                                                  window_size)
            # shrink the coverage back on the edges
            if shrink:
                interpolated_dataset = _shrink_coverage(interpolated_dataset, output_resolution, window_size)

            # mask raster to extent using polygon in vector file
            raster_mask = rasterize_like(vector_file_path, interpolated_dataset, output_resolution)
            output_dataset = _mask_raster(interpolated_dataset, raster_mask)
        elif method == 'kriging':
            chunks_per_side = 128

            # interpolate using Ordinary Kriging
            interpolated_dataset = self.__interp_points_pykrige_kriging(points, window_size,
                                                                        chunks_per_side=chunks_per_side)
            output_dataset = _mask_raster(interpolated_dataset, self.__get_mask())

            method = f'{method}_{chunks_per_side}x{chunks_per_side}'
        else:
            raise ValueError(f'Interpolation type "{method}" not recognized.')

        print(f'{method} interpolation took {(datetime.datetime.now() - start_time).total_seconds()} s')

        if plot:
            self.__plot(output_dataset, method)

        return output_dataset

    def __get_mask(self, shape: Tuple[float, float], output_resolution: float, output_nw: Tuple[float, float],
                   output_spatial_reference: osr.SpatialReference, output_nodata: float = None) -> gdal.Dataset:
        """
        Currently using the shrunk invdist method as a mask.

        Parameters
        ----------
        dataset
            TODO write description
        output_resolution
            TODO write description
        window
            TODO write description

        Returns
        -------
        gdal.Dataset
            mask
        """

        if output_nodata is None:
            output_nodata = self.nodata

        gdal_memory_driver = gdal.GetDriverByName('MEM')
        ogr_memory_driver = ogr.GetDriverByName('Memory')

        # mask raster to extent using polygon in convex hull
        convex_hull_vertices = self.input_points[:, 0:2][
            scipy.spatial.ConvexHull(self.input_points[:, 0:2]).vertices]
        convex_hull_vertices = [tuple(point) for point in convex_hull_vertices.tolist()]

        ring = ogr.Geometry(ogr.wkbLinearRing)
        for vertex in convex_hull_vertices:
            ring.AddPoint(*vertex)

        convex_hull = ogr.Geometry(ogr.wkbPolygon)
        convex_hull.AddGeometry(ring)

        output_dataset = gdal_memory_driver.Create('', shape[0], shape[1], 1, gdal.GDT_Float32)
        output_dataset.SetGeoTransform((output_nw[1], output_resolution, 0, output_nw[0], 0, output_resolution))

        vector_dataset = ogr_memory_driver.CreateDataSource('temp')
        vector_layer = vector_dataset.CreateLayer('temp_layer', output_spatial_reference, ogr.wkbMultiPolygon)
        vector_layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        output_feature = ogr.Feature(vector_layer.GetLayerDefn())
        output_feature.SetField('id', 1)
        output_feature.SetGeometry(convex_hull)
        vector_layer.CreateFeature(output_feature)

        gdal.RasterizeLayer(output_dataset, [1], vector_layer, burn_values=[1])

        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(output_nodata)

        return output_dataset

    def __interp_points_gdal_linear(self, gdal_points: gdal.Dataset, nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating linearly.

        Parameters
        ----------
        gdal_points
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

        algorithm = f"linear:radius=0:nodata={int(nodata)}"
        return gdal.Grid('', gdal_points, format='MEM', width=self.input_shape[1], height=self.input_shape[0],
                         outputBounds=self.input_bounds, algorithm=algorithm)

    def __interp_points_gdal_invdist(self, gdal_points: gdal.Dataset, radius: float,
                                     nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance.

        Parameters
        ----------
        gdal_points
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

        if nodata is None:
            nodata = self.nodata

        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}:angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        return gdal.Grid('', gdal_points, format='MEM', width=self.input_shape[1], height=self.input_shape[0],
                         outputBounds=self.input_bounds, algorithm=algorithm)

    def __interp_points_gdal_invdist_scipy_linear(self, gdal_points: gdal.Dataset, radius: float,
                                                  nodata: float = None) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via inverse distance and then linearly (using SciPy).

        Parameters
        ----------
        gdal_points
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

        if nodata is None:
            nodata = self.nodata

        # interpolate using inverse distance in GDAL
        interpolated_raster = self.__interp_points_gdal_invdist(gdal_points, radius, nodata)
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

    def __interp_points_pykrige_kriging(self, gdal_points: gdal.Dataset, nodata: float = None,
                                        chunks_per_side: int = 128) -> gdal.Dataset:
        """
        Create a regular raster grid from the given points, interpolating via kriging (using PyKrige).

        Parameters
        ----------
        gdal_points
            GDAL point cloud dataset
        nodata
            value for no data in output grid
        chunks_per_side
            number of chunks per side in the square with which to divide the interpolation task
            setting this value to less than 5 has adverse consequences on memory consumption

        Returns
        -------
        gdal.Dataset
            interpolated grid
        """

        if nodata is None:
            nodata = self.nodata

        if chunks_per_side <= 4:
            raise ValueError(f'number of chunks per side should be above 4')

        data_sw = numpy.array((self.input_bounds[0], self.input_bounds[1]))

        interpolated_grid_x = numpy.arange(self.input_bounds[0], self.input_bounds[2], self.resolution)
        interpolated_grid_y = numpy.arange(self.input_bounds[1], self.input_bounds[3], self.resolution)
        interpolated_grid_shape = numpy.array((len(interpolated_grid_y), len(interpolated_grid_x)), numpy.int)
        interpolated_grid_values = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)
        interpolated_grid_variance = numpy.full(interpolated_grid_shape, fill_value=numpy.nan)

        chunk_shape = numpy.array(numpy.round(interpolated_grid_shape / chunks_per_side), numpy.int)
        expand_indices = numpy.round(chunk_shape / 2).astype(numpy.int)

        chunk_grid_index = numpy.array((0, 0), numpy.int)
        with futures.ProcessPoolExecutor() as concurrency_pool:
            running_futures = {}

            chunks = chunks_per_side ** 2
            for chunk_index in range(chunks):
                # print(f'processing chunk {chunk_index + 1} of {chunks}')

                # determine indices of the lower left and upper right bounds of the chunk within the output grid
                grid_start_index = numpy.array(chunk_grid_index * chunk_shape)
                grid_end_index = grid_start_index + chunk_shape
                grid_slice = tuple(
                    slice(grid_start_index[axis], grid_end_index[axis]) for axis in range(len(chunk_shape)))

                # expand the interpolation area past the edges of the chunk to minimize edge formation at the boundary
                interpolation_sw = data_sw + numpy.flip((grid_start_index - expand_indices) * self.resolution)
                interpolation_ne = data_sw + numpy.flip((grid_end_index + expand_indices) * self.resolution)
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

                if chunk_grid_index[1] >= chunks_per_side - 1:
                    chunk_grid_index[1] = 0
                    chunk_grid_index[0] += 1
                else:
                    chunk_grid_index[1] += 1

            for _, completed_future in enumerate(futures.as_completed(running_futures)):
                grid_slice = running_futures[completed_future]
                chunk_interpolated_values, chunk_interpolated_variance = completed_future.result()
                interpolated_grid_values[grid_slice] = chunk_interpolated_values
                interpolated_grid_variance[grid_slice] = chunk_interpolated_variance

        interpolated_grid_uncertainty = numpy.sqrt(interpolated_grid_variance) * 2.5
        del interpolated_grid_variance

        interpolated_grid_values = numpy.flip(interpolated_grid_values, axis=0)
        interpolated_grid_uncertainty = numpy.flip(interpolated_grid_uncertainty, axis=0)

        self.__plot(interpolated_grid_values, f'kriging_{chunks_per_side}x{chunks_per_side}')

        memory_raster_driver = gdal.GetDriverByName('MEM')
        output_dataset = memory_raster_driver.Create('temp', self.input_shape[1], self.input_shape[0], 2,
                                                     gdal.GDT_Float32)

        input_geotransform = gdal_points.GetGeoTransform()
        output_dataset.SetGeoTransform(
            (input_geotransform[0], self.resolution, 0.0, 0.0, input_geotransform[3], self.resolution))
        output_dataset.SetProjection(gdal_points.GetProjection())

        band_1 = output_dataset.GetRasterBand(1)
        band_1.SetNoDataValue(nodata)
        band_1.WriteArray(interpolated_grid_values)

        band_2 = output_dataset.GetRasterBand(2)
        band_2.SetNoDataValue(nodata)
        band_2.WriteArray(interpolated_grid_uncertainty)

        del band_1, band_2
        return output_dataset

    def __plot(self, interpolated_raster: Union[numpy.array, gdal.Dataset], interpolation_method: str,
               nodata: float = None):
        """
        Plot preinterpolated points and an interpolated raster side-by-side on synchronized subplots.

        Parameters
        ----------
        interpolated_raster
            array of interpolated data
        interpolation_method
            method of interpolation
        """

        if nodata is None:
            nodata = self.nodata

        if type(interpolated_raster) is gdal.Dataset:
            interpolated_raster = numpy.flip(interpolated_raster.ReadAsArray(0), axis=0)
        if len(interpolated_raster.shape) == 3:
            interpolated_raster = interpolated_raster[0, :, :]
        interpolated_raster[interpolated_raster == nodata] = numpy.nan

        min_x, min_y, max_x, max_y = self.input_bounds

        figure = plt.figure()

        survey_points_axis = figure.add_subplot(1, 2, 1)
        survey_points_axis.set_title('survey points')
        survey_points_axis.scatter(self.input_points[:, 0], self.input_points[:, 1], c=self.input_points[:, 2], s=1)
        survey_points_axis.set_xlim(min_x, max_x)
        survey_points_axis.set_ylim(min_y, max_y)

        interpolated_raster_axis = figure.add_subplot(1, 2, 2, sharex=survey_points_axis)
        interpolated_raster_axis.set_title(f'{interpolation_method} interpolation to raster')
        interpolated_raster_axis.imshow(interpolated_raster, extent=(min_x, max_x, min_y, max_y), aspect='auto')

        plt.show()


def _shrink_coverage(raster: gdal.Dataset, resolution: float, radius: float) -> gdal.Dataset:
    """
    Shrink coverage of a dataset by the original coverage radius.

    Parameters
    ----------
    raster
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

    for band_index in range(1, raster.RasterCount + 1):
        band = raster.GetRasterBand(1)
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

    return raster


def _mask_raster(raster: gdal.Dataset, mask: gdal.Dataset, mask_value: float = None) -> gdal.Dataset:
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

    Returns
    -------
    gdal.Dataset
        masked raster dataset
    """

    raster_band = raster.GetRasterBand(1)
    mask_band = mask.GetRasterBand(1)

    if mask_value is None:
        mask_value = mask_band.GetNoDataValue()

    raster_data = raster_band.ReadAsArray()
    mask_data = mask_band.ReadAsArray()

    raster_data[mask_data == mask_value] = raster_band.GetNoDataValue()

    raster_band.WriteArray(raster_data)
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


def rasterize_like(vector_file_path: str, like_raster: gdal.Dataset, output_resolution: int) -> gdal.Dataset:
    """
    Burn vector data to a raster with the properties of the given raster.

    Parameters
    ----------
    vector_file_path
        path to file containing vector data
    like_raster
        GDAL raster dataset to emulate
    output_resolution
        cell size of output raster

    Returns
    -------
    gdal.Dataset
        raster with vector data burned in
    """

    input_geotransform = like_raster.GetGeoTransform()
    input_spatial_reference = like_raster.GetProjectionRef()

    if input_spatial_reference == '':
        output_spatial_reference = osr.SpatialReference(
            'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
    else:
        output_spatial_reference = osr.SpatialReference(wkt=input_spatial_reference)
        output_spatial_reference.MorphFromESRI()

    return rasterize(vector_file_path, output_spatial_reference, (input_geotransform[0], input_geotransform[3]),
                     output_resolution)


def rasterize(vector_file_path: str, output_spatial_reference: osr.SpatialReference, output_nw: Tuple[float, float],
              output_resolution: int, output_nodata: float = 0) -> gdal.Dataset:
    """
    Burn vector data to a raster with the given properties.

    Parameters
    ----------
    vector_file_path
        path to file containing vector data
    output_spatial_reference
        output spatial reference system as well-known text
    output_nw
        output northwest corner
    output_resolution
        output cell size
    output_shape
        output shape
    output_nodata
        value to be substituted for no data

    Returns
    -------
    gdal.Dataset
        raster dataset
    """

    if type(output_nw) is not numpy.array:
        output_nw = numpy.array(output_nw)

    # Open the data source and read in the extent
    input_dataset = ogr.Open(vector_file_path)
    input_layer = input_dataset.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()

    input_x_min, input_x_max, input_y_min, input_y_max = None, None, None, None
    temp_layer = None

    ogr_memory_driver = ogr.GetDriverByName('Memory')

    # extract the first feature to a temporary layer
    for input_feature in input_layer:
        if input_feature is not None:
            input_geometry = input_feature.GetGeometryRef()
            output_geometry = ogr.CreateGeometryFromWkt(input_geometry.ExportToWkt())
            output_geometry.Transform(osr.CoordinateTransformation(input_spatial_reference, output_spatial_reference))
            input_x_min, input_x_max, input_y_min, input_y_max = output_geometry.GetEnvelope()

            temp_dataset = ogr_memory_driver.CreateDataSource('temp')
            temp_layer = temp_dataset.CreateLayer('temp_layer', output_spatial_reference, ogr.wkbMultiPolygon)
            temp_layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
            output_feature = ogr.Feature(temp_layer.GetLayerDefn())
            output_feature.SetField('id', 1)
            output_feature.SetGeometry(output_geometry)
            temp_layer.CreateFeature(output_feature)
            break

    del input_dataset

    input_sw = numpy.array((input_x_min, input_y_min))
    input_ne = numpy.array((input_x_max, input_y_max))
    output_shape = tuple(numpy.round((input_ne - input_sw) / output_resolution).astype(numpy.int))

    # burn the single feature to the output dataset
    rasterized_dataset = gdal.GetDriverByName('MEM').Create('', int(output_shape[0]), int(output_shape[1]),
                                                            gdal.GDT_Byte)
    rasterized_dataset.SetGeoTransform((output_nw[1], output_resolution, 0, output_nw[0], 0, output_resolution))
    band_1 = rasterized_dataset.GetRasterBand(1)
    band_1.SetNoDataValue(output_nodata)
    gdal.RasterizeLayer(rasterized_dataset, [1], temp_layer, burn_values=[1])
    rasterized_array = band_1.ReadAsArray()

    # clip resampled output data to the bounds of input
    output_array = numpy.full(output_shape, output_nodata)
    output_se = numpy.array((output_nw[0] + (output_resolution * output_shape[1]),
                             output_nw[1] - (output_resolution * output_shape[0])))
    input_nw = numpy.array((input_x_min, input_y_max))
    input_se = numpy.array((input_x_max, input_y_min))

    if input_nw[0] > output_se[0] or input_se[0] < output_nw[0] or \
            input_se[1] > output_nw[1] or input_nw[1] < output_se[1]:
        raise ValueError('input data is outside output grid bounds')

    index_delta_nw = numpy.round((input_nw - output_nw) / numpy.array(output_resolution)).astype(int)
    index_delta_se = numpy.round((input_se - output_nw) / numpy.array(output_resolution)).astype(int)

    # indices to be written onto the output array
    output_index_slices = [slice(0, None), slice(0, None)]

    # check if western bound of input is further west than that of output
    if index_delta_nw[0] < 0:
        output_index_slices[1] = slice(index_delta_nw[0] * -1, output_index_slices[1].stop)
        index_delta_nw[0] = 0

    # check if northern bound of input is further north than that of output
    if index_delta_nw[1] < 0:
        output_index_slices[0] = slice(index_delta_nw[1] * -1, output_index_slices[0].stop)
        index_delta_nw[1] = 0

    # check if eastern bound of input is further east than that of output
    if index_delta_se[0] > rasterized_array.shape[1]:
        output_index_slices[1] = slice(output_index_slices[1].start, rasterized_array.shape[1] - index_delta_se[0])
        index_delta_se[0] = rasterized_array.shape[1]

    # check if southern bound of input is further south than that of output
    if index_delta_se[1] > rasterized_array.shape[0]:
        output_index_slices[0] = slice(output_index_slices[0].start, rasterized_array.shape[0] - index_delta_se[1])
        index_delta_se[1] = rasterized_array.shape[0]

    # write relevant input data to a slice of the output array corresponding to the input extent
    output_array[output_index_slices[0], output_index_slices[1]] = rasterized_array[index_delta_nw[1]:index_delta_se[1],
                                                                   index_delta_nw[0]:index_delta_se[0]]

    del rasterized_dataset, rasterized_array

    output_dataset = ogr_memory_driver.Create('', output_array.shape[1], output_array.shape[0], 1, gdal.GDT_Float32)
    output_dataset.SetGeoTransform((output_nw[1], output_resolution, 0, output_nw[0], 0, output_resolution))

    band_1 = rasterized_dataset.GetRasterBand(1)
    band_1.SetNoDataValue(output_nodata)
    band_1.WriteArray(output_array)

    return output_dataset


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


def _krige_onto_grid(input_points: numpy.array, output_x: numpy.array, output_y: numpy.array) -> Tuple[
    numpy.array, numpy.array]:
    """
    Perform kriging interpolation from the given set of points onto a regular grid.

    Parameters
    ----------
    input_points
        array of points (N x M)
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
