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
import numpy
import scipy
import sklearn.gaussian_process
from matplotlib.mlab import griddata as mlab_griddata
from osgeo import gdal, ogr, osr


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

    if numpy.isnan(valmin) or val < valmin:
        valmin = val

    if numpy.isnan(valmax) or val > valmax:
        valmax = val

    return valmin, valmax


class point_interpolator:
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

        if interpolation_type in ['natural', 'invdist_scilin'] and shapefile == None:
            raise ValueError('Supporting shapefile required')

        data_array = self._gdal2vector(dataset)
        _, _, maxrad = self._get_point_spacing(data_array)
        window = maxrad * self.window_scale

        if interpolation_type == 'linear':
            # do the triangulation interpolation
            ds2 = self._gdal_linear_interp_points(dataset, resolution)
        elif interpolation_type == 'natural':
            ds2 = self._gdal_mlab_natural_interp_points(dataset, resolution)
            if shrink:
                ds4 = self._shrink_coverage(ds2, resolution, window)
        elif interpolation_type == 'invlin':
            ds2 = self._gdal_invdist_scilin_interp_points(dataset, resolution,
                                                          window)
            if shrink:
                ds4 = self._shrink_coverage(ds2, resolution, window)
        elif interpolation_type == 'invdist':
            # do the inverse distance interpolation
            ds3 = self._gdal_invdist_interp_points(dataset, resolution, window)
            # shrink the coverage back on the edges and in the holidays on the
            if shrink:
                ds4 = self._shrink_coverage(ds3, resolution, window)
                # shrink the coverage back on the edges and in the holidays on the inv dist
                if shrink:
                    ds4 = self._shrink_coverage(ds3, resolution, window)
        elif interpolation_type == 'kriging':
            ds3 = self._kriging_interp_points(dataset, resolution, window)

            # shrink the coverage back on the edges and in the holidays on the inv dist
            if shrink:
                ds4 = self._shrink_coverage(ds3, resolution, window)
        else:
            print('No interpolation method recognized')

        if interpolation_type == 'linear':
            # trim the triangulated interpolation back using the inv dist as a
            # mask
            ds3 = self._get_mask(dataset, resolution, window)
            ds5 = self._mask_with_raster(ds2, ds3)
        elif interpolation_type in ['natural', 'invlin']:
            ds3 = self._get_shape_mask(ds2, shapefile, resolution)
            if shrink:
                ds5 = self._mask_with_raster(ds4, ds3)
            else:
                ds5 = self._mask_with_raster(ds2, ds3)

        # write the files out using the above function
        if interpolation_type in ['linear', 'natural', 'invlin']:
            return ds5
        elif interpolation_type == 'invdist':
            if shrink:
                return ds4
            else:
                return ds3
        else:
            raise ValueError('Interpolation type not understood')

    def _gdal2vector(self, dataset: gdal.Dataset) -> numpy.array:
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
        data = numpy.zeros((count, 3))

        for n in numpy.arange(count):
            f = lyr.GetFeature(n)
            data[n, :] = f.geometry().GetPoint()

        return data

    def _get_point_spacing(self, dataset: numpy.array) -> Tuple[float, float, float]:
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
        min_dist = numpy.zeros(count) + numpy.inf

        # roll the array through, comparing all points and saving the minimum dist.
        for n in numpy.arange(1, count):
            tmp = numpy.roll(dataset, n, axis=0)
            dist = (numpy.sqrt(numpy.square(dataset[:, 0] - tmp[:, 0]) + numpy.square(dataset[:, 1] - tmp[:, 1])))
            idx = numpy.nonzero(dist < min_dist)[0]
            if len(idx) > 0:
                min_dist[idx] = dist[idx]

        return min_dist.min(), min_dist.mean(), min_dist.max()

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
        type
            mask

        """

        data = self._gdal_invdist_interp_points(dataset, resolution, window)
        mask = self._shrink_coverage(data, resolution, window)
        return mask

    def _getShpRast(self, file: str, to_proj: str, to_gt: tuple, to_res: int, to_y, to_x, nodata: float = 0) -> Tuple[
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
        to_y :
            TODO write description
        to_x:
            TODO write description
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
        plt.figure()
        plt.imshow(newarr)
        plt.show()

        # 6
        print(x_orig, y_orig)
        print(x_min, y_min)
        ollx, olly = x_orig, y_orig
        rllx, rlly = x_min, y_min
        dllx, dlly = 0, 0
        if ollx != rllx:
            dllx = rllx - ollx
            print(ollx, rllx, dllx)
        if olly != rlly:
            dlly = olly - rlly
            print(olly, rlly, dlly)
        print(dllx, dlly)

        # 7
        oShape = (to_y, to_x)
        oSy, oSx = oShape
        rSy, rSx = newarr.shape
        print(oSy, rSy)
        print(oSx, rSx)
        expx, expy = 0, 0
        if newarr.shape != oShape:
            print(oSy - rSy, oSx - rSx)
            if rSy < oSy:
                expy = int(numpy.abs(oSy - rSy))
                print('expy', expy)
            if rSx < oSx:
                expx = int(numpy.abs(oSx - rSx))
                print('expx', expx)
        ay = numpy.full((rSy + expy, rSx + expx), nodata)
        print('expz', ay.shape, oShape)
        rollx = int(dllx / to_res)
        rolly = int(dlly / to_res)

        # 8
        up, left = 0, 0
        down, right = 0, 0
        if dlly < 0:
            down = abs(rolly)
            up = 0
        elif dlly > 0:
            up = -int(rolly)
        if dllx < 0:
            left = -int(rollx)
        elif dllx > 0:
            right = abs(rollx)
            left = 0

        if dllx != 0 or dlly != 0:
            print('rollz', up, left, down, right)
            temp = newarr[up:, left:]
            print(temp.shape)
            plt.imshow(temp)
            plt.show()
            ay[down:temp.shape[0] + down, right:temp.shape[1] + right] = temp[:, :]
            del temp
        else:
            ay[:] = newarr[:]

        print('expz', ay.shape)
        ax = numpy.full(oShape, nodata)
        ax[:] = ay[:oSy, :oSx]

        del newarr, ay, band, source_ds, target_ds

        ax_y, ax_x = ax.shape

        target_ds = gdal.GetDriverByName('MEM').Create('', ax_x, ax_y, 1, gdal.GDT_Float32)
        x_orig, y_orig = to_gt[0], to_gt[3]
        target_gt = (x_orig, to_res, 0, y_orig, 0, to_res)
        target_ds.SetGeoTransform(target_gt)
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        band.WriteArray(ax)

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
        bandy, bandx = grid.RasterYSize, grid.RasterXSize
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
        shape_ds, shape_gt = self._getShpRast(shapefile, proj, grid_gt, int(resolution), bandy, bandx)
        print('transformed', shape_gt)
        return shape_ds

    def _gdal_linear_interp_points(self, dataset: gdal.Dataset, resolution: float,
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
        algorithm = f"linear:radius=0:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)
        return interp_data

    def _gdal_mlab_natural_interp_points(self, dataset: gdal.Dataset, resolution: float,
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
        nodata: float :
            TODO write description (Default value = 1000000)

        Returns
        -------
        type
            interpolated dataset

        """

        print('_gdal_mlab_natural_interp_points')
        # Find the bounds of the provided data
        xmin, xmax, ymin, ymax = numpy.nan, numpy.nan, numpy.nan, numpy.nan
        lyr = dataset.GetLayerByIndex(0)
        proj = lyr.GetSpatialRef().ExportToWkt()
        count = lyr.GetFeatureCount()
        xvals, yvals, zvals = [], [], []

        for n in numpy.arange(count):
            f = lyr.GetFeature(n)
            x, y, z = f.geometry().GetPoint()
            xvals.append(x)
            yvals.append(y)
            zvals.append(z)
            xmin, xmax = _compare_vals(x, xmin, xmax)
            ymin, ymax = _compare_vals(y, ymin, ymax)

        numrows, numcolumns, bounds = self._get_nodes3(resolution, (xmin, ymin, xmax, ymax))
        print(bounds)
        xbound, ybound = bounds[0], bounds[1]
        xvals, yvals, zvals = numpy.array(xvals), numpy.array(yvals), numpy.array(zvals)
        xvals, yvals = (xvals - xbound) / resolution, (yvals - ybound) / resolution

        print('start', xvals, yvals, zvals)
        xi, yi = numpy.arange(numcolumns), numpy.arange(numrows)
        interp_obj = mlab_griddata(xvals, yvals, zvals, xi, yi, interp='nn')
        interp_grid, interp_mask = interp_obj.data, interp_obj.mask
        interp_grid[numpy.isnan(interp_grid)] = nodata
        plt.figure()
        plt.imshow(interp_grid)
        plt.show()
        plt.figure()
        plt.imshow(interp_mask)
        plt.show()
        print('stop')

        interp_data = gdal.GetDriverByName('MEM').Create('', numcolumns, numrows, 1, gdal.GDT_Float32)
        interp_gt = (xbound, resolution, 0, ybound, 0, resolution)
        interp_data.SetGeoTransform(interp_gt)
        interp_data.SetProjection(proj)

        band = interp_data.GetRasterBand(1)
        band.SetNoDataValue(float(nodata))
        band.WriteArray(interp_grid)

        del band

        return interp_data

    def _gdal_invdist_scilin_interp_points(self, dataset: gdal.Dataset, resolution: float, radius: float,
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
        xmin, xmax, ymin, ymax = numpy.nan, numpy.nan, numpy.nan, numpy.nan
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()

        for n in numpy.arange(count):
            f = lyr.GetFeature(n)
            x, y, z = f.geometry().GetPoint()
            xmin, xmax = _compare_vals(x, xmin, xmax)
            ymin, ymax = _compare_vals(y, ymin, ymax)

        numrows, numcolumns, bounds = self._get_nodes3(resolution, (xmin, ymin, xmax, ymax))
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}" + \
                    f":angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
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
        plt.figure()
        plt.imshow(interp_grid)
        plt.show()
        print('stop')

        band = interp_data.GetRasterBand(1)
        band.SetNoDataValue(float(nodata))
        band.WriteArray(interp_grid)
        del band

        return interp_data

    def _gdal_invdist_interp_points(self, dataset: gdal.Dataset, resolution: float, radius: float,
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
        algorithm = f"invdist:power=2.0:smoothing=0.0:radius1={radius}:radius2={radius}" + \
                    f":angle=0.0:max_points=0:min_points=1:nodata={int(nodata)}"
        interp_data = gdal.Grid('', dataset, format='MEM', width=numcolumns, height=numrows, outputBounds=bounds,
                                algorithm=algorithm)
        return interp_data

    def _kriging_interp_points(self, dataset: gdal.Dataset, resolution: float, radius: float,
                               nodata: float = 1000000) -> gdal.Dataset:
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.
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

        # do kriging on points
        gp = sklearn.gaussian_process.GaussianProcessRegressor(kernel=sklearn.gaussian_process.kernels.ConstantKernel)

        # TODO actually make this work
        interpolated_data = gp.fit(x, y)

        # TODO convert back to GDAL dataset
        return interpolated_data

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

        data = dataset.ReadAsArray()
        data_rb = dataset.GetRasterBand(1)
        nodata = data_rb.GetNoDataValue()
        idx = numpy.nonzero(data == nodata)
        data[idx] = numpy.nan
        # divide the window size by the resolution to get the number of cells
        rem_cells = int(numpy.round(radius / resolution))
        # print(f'Shrinking coverage back {rem_cells} cells.')

        for n in numpy.arange(rem_cells):
            ew = numpy.diff(data, axis=0)
            idx_ew = numpy.nonzero(numpy.isnan(ew))
            ns = numpy.diff(data, axis=1)
            idx_ns = numpy.nonzero(numpy.isnan(ns))
            data[idx_ew] = numpy.nan
            data[idx_ew[0] + 1, idx_ew[1]] = numpy.nan
            data[idx_ns] = numpy.nan
            data[idx_ns[0], idx_ns[1] + 1] = numpy.nan

        idx = numpy.nonzero(numpy.isnan(data))
        data[idx] = nodata
        data_rb.WriteArray(data)
        return dataset
