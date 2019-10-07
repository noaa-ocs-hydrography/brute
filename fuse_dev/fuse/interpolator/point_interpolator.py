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

import os
from typing import Tuple, Any

import matplotlib.pyplot as plt
import numpy as np
import scipy
from osgeo import gdal, ogr, osr
from fuse.utilities import gdal_to_xyz, maximum_nearest_neighbor_distance, shape_from_cell_size


def _compare_vals(val: float, valmin: float, valmax: float) -> Tuple[float, float]:
    """
    This is a small utility for inspecting values and seeing they
    contribute to a min or max estimate.

    Parameters
    ----------
    val: float :
        TODO write description
    valmin: float :
        TODO write description
    valmax: float :
        TODO write description

    Returns
    -------

    """

    if np.isnan(valmin):
        valmin = val
    elif val < valmin:
        valmin = val

    if np.isnan(valmax):
        valmax = val
    elif val > valmax:
        valmax = val

    return valmin, valmax


class PointInterpolator:
    """Interpolation methods for creating a raster from points."""

    def __init__(self, window_scalar: float = 1.1):
        """
        Set some of the precondition, but make them over writable.
        """

        self.window_scale = window_scalar

    def interpolate(self, dataset: gdal.Dataset, interpolation_type: str, resolution: float, shapefile: str = None,
                    shrink: bool = True) -> gdal.Dataset:
        """
        Interpolate the provided dataset.

        Currently this is assumed to be a gdal dataset.  At some point perhaps
        this should become something more native to python.

        Parameters
        ----------
        dataset: gdal.Dataset :
            TODO write description
        interpolation_type: str :
            TODO write description
        resolution: float :
            TODO write description
        shapefile: str :
            TODO write description (Default value = None)
        shrink: bool :
            TODO write description (Default value = True)

        Returns
        -------
        type
            interpolated dataset

        """

        linear = False
        natural = False
        invlin = False
        invdist = False

        if interpolation_type == 'linear':
            linear = True
        elif interpolation_type == 'invdist_scilin':
            if shapefile is None:
                raise ValueError('Supporting shapefile required')
            invlin = True
        elif interpolation_type == 'invdist':
            invdist = True
        else:
            raise ValueError('interpolation type not implemented.')

        data_array = gdal_to_xyz(dataset)[:, :2]

        # get the bounds (min X, min Y, max X, max Y)
        input_bounds = np.concatenate((np.min(data_array, axis=0), np.max(data_array, axis=0)))

        # set the size of the interpolation window size to be the minimum of all nearest neighbor distances
        window = self.window_scale * maximum_nearest_neighbor_distance(data_array)

        if linear:
            # do the triangulation interpolation
            ds2 = self._gdal_linear_interp_points(dataset, resolution, input_bounds)
        elif invlin:
            ds2 = self._gdal_invdist_scilin_interp_points(dataset, resolution,
                                                          input_bounds, window)
            if shrink:
                ds4 = self._shrink_coverage(ds2, resolution, window)
        elif invdist:
            # do the inverse distance interpolation
            ds3 = self._gdal_invdist_interp_points(dataset, resolution, input_bounds, window)
            # shrink the coverage back on the edges and in the holidays on the
            # inv dist
            if shrink:
                ds4 = self._shrink_coverage(ds3, resolution, window)
        else:
            print('No interpolation method recognized')

        if linear:
            # trim the triangulated interpolation back using the inv dist as a
            # mask
            ds3 = self._get_mask(dataset, resolution, input_bounds, window)
            ds5 = self._mask_with_raster(ds2, ds3)
        elif invlin:
            ds3 = self._get_shape_mask(ds2, shapefile, resolution)
            if shrink:
                ds5 = self._mask_with_raster(ds4, ds3)
            else:
                ds5 = self._mask_with_raster(ds2, ds3)

        # write the files out using the above function
        if linear or invlin:
            return ds5
        elif invdist:
            if shrink:
                return ds4
            else:
                return ds3
        else:
            raise ValueError('Interpolation type not understood')

    def _gdal2vector(self, dataset: gdal.Dataset) -> np.array:
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

    def _get_point_spacing(self, dataset: np.array) -> Tuple[float, float, float]:
        """
        Take a numpy xyz array and return the min, mean, and max spacing
        between different points in the XY direction for interpolation without
        holes.  The returned units will be the same as the provided horizontal
        coordinate system.

        Parameters
        ----------
        dataset: np.array :
            TODO write description

        Returns
        -------
        type
            minimum, mean, and maximum

        """

        count = len(dataset)
        min_dist = np.zeros(count) + np.inf

        # roll the array through, comparing all points and saving the minimum dist.
        for n in np.arange(1, count):
            tmp = np.roll(dataset, n, axis=0)
            dist = (np.sqrt(np.square(dataset[:, 0] - tmp[:, 0])
                            + np.square(dataset[:, 1] - tmp[:, 1])))
            idx = np.nonzero(dist < min_dist)[0]
            if len(idx) > 0:
                min_dist[idx] = dist[idx]

        return min_dist.min(), min_dist.mean(), min_dist.max()

    def _get_mask(self, dataset: gdal.Dataset, resolution: float, input_bounds: (float, float, float, float), window: float) -> gdal.Dataset:
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
        type
            mask

        """

        data = self._gdal_invdist_interp_points(dataset, resolution, input_bounds, window)
        mask = self._shrink_coverage(data, resolution, window)
        return mask

    def _getShpRast(self, file: str, to_proj: str, to_gt: tuple, to_res: int, to_shape, nodata: float = 0) -> Tuple[
        gdal.Dataset, Any]:
        """
        Import shapefile

        Parameters
        ----------
        file: str :
            Shapefile file location
        to_proj: str:
            WKT object with destination spatial reference system
        to_gt: tuple:
            gdal.GeoTransform object of the interpolated dataset
        to_res: int :
            Resolution of the input dataset
        shape : tuple(int, int)
            Shape of the the input dataset
        nodata: float :
            Nodata value for the ouput gdal.Dataset object (Default value = 0)

        Returns
        -------
        type
            dataset and geotransform

        """

        print('getShpRast', file)
        fName = os.path.split(file)[-1]
        splits = os.path.splitext(fName)
        name = splits[0]
        # tif = f'{splits[0]}.tif'

        # Open the data source and read in the extent
        source_ds = ogr.Open(file)
        source_layer = source_ds.GetLayer()
        source_srs = source_layer.GetSpatialRef()

        for feature in source_layer:
            if feature is not None:
                geom = feature.GetGeometryRef()
                #                print(geom.ExportToWkt())
                ds_geom = ogr.CreateGeometryFromWkt(geom.ExportToWkt())
                #                print(source_srs, to_proj, sep='\n')
                coordTrans = osr.CoordinateTransformation(source_srs, to_proj)
                ds_geom.Transform(coordTrans)
                driver = ogr.GetDriverByName('Memory')
                ds = driver.CreateDataSource('temp')
                layer = ds.CreateLayer(name, to_proj, ogr.wkbMultiPolygon)

                # Add one attribute
                layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
                defn = layer.GetLayerDefn()

                # Create a new feature (attribute and geometry)
                feat = ogr.Feature(defn)
                feat.SetField('id', 123)

                # Make a geometry, from Shapely object
                feat.SetGeometry(ds_geom)

                layer.CreateFeature(feat)
                del feat, geom  # destroy these
                break

        x_min, x_max, y_min, y_max = ds_geom.GetEnvelope()
        meta = ([x_min, y_max], [x_max, y_min])
        print(meta)

        # Create the destination data source
        x_dim = int((x_max - x_min) / to_res)
        y_dim = int((y_max - y_min) / to_res)
        print(x_dim, y_dim)
        target_ds = gdal.GetDriverByName('MEM').Create('', x_dim, y_dim, gdal.GDT_Byte)
        x_orig, y_orig = to_gt[0], to_gt[3]
        target_gt = (x_orig, to_res, 0, y_orig, 0, to_res)
        target_ds.SetGeoTransform(target_gt)
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)

        # Rasterize
        gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])
        newarr = band.ReadAsArray()

        # clip resampled coverage data to the bounds of the BAG
        output_array = np.full(to_shape, nodata)

        cov_ul, cov_lr = np.array(x_orig), np.array(y_orig)
        bag_ul, bag_lr = np.array(x_min), np.array(y_min)

        if bag_ul[0] > cov_lr[0] or bag_lr[0] < cov_ul[0] or bag_lr[1] > cov_ul[1] or bag_ul[1] < cov_lr[1]:
            raise ValueError('bag dataset is outside the bounds of coverage dataset')

        ul_index_delta = np.round((bag_ul - cov_ul) / np.array(to_res)).astype(int)
        lr_index_delta = np.round((bag_lr - cov_ul) / np.array(to_res)).astype(int)

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
        x_orig, y_orig = to_gt[0], to_gt[3]
        target_gt = (x_orig, to_res, 0, y_orig, 0, to_res)
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
        plt.figure()
        plt.imshow(raster)
        plt.show()
        del band
        grid_ref = grid.GetProjectionRef()
        grid_gt = grid.GetGeoTransform()
        proj = osr.SpatialReference(wkt=grid_ref)
        proj.MorphFromESRI()
        print(grid_ref, grid_gt)
        shape_ds, shape_gt = self._getShpRast(shapefile, proj, grid_gt, int(resolution), shape)
        print('transformed', shape_gt)
        return shape_ds

    def _gdal_linear_interp_points(self, dataset: gdal.Dataset, resolution: float, input_bounds: (float, float, float, float),
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


        (numrows, numcolumns), bounds = shape_from_cell_size(resolution, input_bounds)
        algorithm = f"linear:radius=0:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)

        print(interp_data.GetGeoTransform())
        return interp_data

    def _gdal_invdist_scilin_interp_points(self, dataset: gdal.Dataset, resolution: float, input_bounds: (float, float, float, float), radius: float,
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
        (numrows, numcolumns), bounds = shape_from_cell_size(resolution, input_bounds)
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}" + \
                    f":angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)
        arr = interp_data.ReadAsArray()
        ycoord, xcoord = np.where(arr != nodata)
        zvals = arr[ycoord[:], xcoord[:]]
        del arr

        print('start')
        xa, ya = np.arange(numcolumns), np.arange(numrows)
        xi, yi = np.meshgrid(xa, ya)
        interp_grid = scipy.interpolate.griddata((xcoord, ycoord), zvals, (xi, yi), method='linear', fill_value=nodata)

        interp_grid[np.isnan(interp_grid)] = nodata
        plt.figure()
        plt.imshow(interp_grid)
        plt.show()
        print('stop')

        band = interp_data.GetRasterBand(1)
        band.SetNoDataValue(float(nodata))
        band.WriteArray(interp_grid)
        del band

        return interp_data

    def _gdal_invdist_interp_points(self, dataset: gdal.Dataset, resolution: float, input_bounds: (float, float, float, float), radius: float,
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

        (numrows, numcolumns), bounds = shape_from_cell_size(resolution, input_bounds)
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}" + \
                    f":angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)
        return interp_data

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
        numcolumns = int(np.ceil((1. * (xmax - xmin) / resolution)))
        xmean = np.mean([xmin, xmax])
        xnmin = xmean - numcolumns / 2. * resolution
        xnmax = xmean + numcolumns / 2. * resolution
        numrows = int(np.ceil((1. * (ymax - ymin) / resolution)))
        ymean = np.mean([ymin, ymax])
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
        numcolumns = int(np.ceil((1. * (xmax - xmin) / resolution)))
        xnmin = xmin
        xnmax = xmin + 1. * numcolumns * resolution
        numrows = int(np.ceil((1. * (ymax - ymin) / resolution)))
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
        numcolumns = int(np.ceil((1. * (xmax - xmin) / resolution)))
        xrem = xmin % resolution
        xnmin = xmin - xrem
        xnmax = xnmin + 1. * numcolumns * resolution
        numrows = int(np.ceil((1. * (ymax - ymin) / resolution)))
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
        idx = np.nonzero(mask == mask_nd)
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

        data = dataset.ReadAsArray()
        data_rb = dataset.GetRasterBand(1)
        nodata = data_rb.GetNoDataValue()
        idx = np.nonzero(data == nodata)
        data[idx] = np.nan
        # divide the window size by the resolution to get the number of cells
        rem_cells = int(np.round(radius / resolution))
        # print(f'Shrinking coverage back {rem_cells} cells.')

        for n in np.arange(rem_cells):
            ew = np.diff(data, axis=0)
            idx_ew = np.nonzero(np.isnan(ew))
            ns = np.diff(data, axis=1)
            idx_ns = np.nonzero(np.isnan(ns))
            data[idx_ew] = np.nan
            data[idx_ew[0] + 1, idx_ew[1]] = np.nan
            data[idx_ns] = np.nan
            data[idx_ns[0], idx_ns[1] + 1] = np.nan

        idx = np.nonzero(np.isnan(data))
        data[idx] = nodata
        data_rb.WriteArray(data)
        return dataset
