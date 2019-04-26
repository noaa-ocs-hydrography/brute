# -*- coding: utf-8 -*-
"""
point_interpolator.py

grice 20180625

V0.0.3 20190211

This is a collection of functions for creating interpolated bathymetry for the
national bathymetry.

The objective is to interpolate both XYZ data and BAGs.


Sources:
    Make ogr dataset from numpy array: https://pcjericks.github.io/py-gdalogr-cookbook/geometry.html
    ogr data set to gdal for gridding: http://osgeo-org.1560.x6.nabble.com/gdal-dev-DataSource-Dataset-using-gdal-Grid-td5322689.html
    GDAL gridding information: http://www.gdal.org/grid_tutorial.html#grid_tutorial_interpolation
    gdal_grid: http://www.gdal.org/gdal_grid.html
    
"""

# __version__ = 'point_interpolator 0.0.1'

from osgeo import gdal
import numpy as np
import matplotlib as mpl
from mpl_toolkits import natgrid as _nn

class point_interpolator:
    """
    Interpolation methods for creating a raster from points.
    """
    def __init__(self, window_scalar = 1.1):
        """
        Set some of the precondition, but make them over writable.
        """
        self.win_scal = window_scalar
    
    def interpolate(self, dataset, interpolation_type, resolution):
        """
        Interpolate the provided dataset.
        
        Currently this is assumed to be a gdal dataset.  At some point perhaps
        this should become something more native to python.
        """
        linear = False
        invdist = False
        natural = False
        if interpolation_type == 'linear':
            linear = True
        elif interpolation_type == 'invdist':
            invdist = True
        elif interpolation_type == 'natural':
            natural = True
        else:
            raise ValueError('interpolation type not implemented.')
        data_array = self._gdal2vector(dataset)
        minrad, meanrad, maxrad = self._get_point_spacing(data_array)
        window = maxrad * self.window_scal
        if linear:
            # do the triangulation interpolation
            ds2 = self._gdal_linear_interp_points(dataset, resolution)
        elif natural:
            ds2 = self._get_natural_interp(dataset)
        elif invdist:
            # do the inverse distance interpolation
            ds3 = self._gdal_invdist_interp_points(dataset, resolution, window)
            # shrink the coverage back on the edges and in the holidays on the inv dist
            ds4 = self._shrink_coverage(ds3, resolution, window)
        if linear or natural:
            # trim the triangulated interpolation back using the inv dist as a mask
            ds3 = self._get_mask(dataset, resolution, window)
            ds4 = self._shrink_coverage(ds3, resolution, window)
            ds5 = self._mask_with_raster(ds2, ds4)
        # write the files out using the above function
        if linear or natural:
            return ds5
        elif invdist:
            return ds4
        else:
            raise ValueError('Interpolation type not understood')
    
    def _gdal2vector(self, dataset):
        """
        Take a gdal vector xyz point cloud and return a numpy array.
        """
        # get the data out of the gdal data structure
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()
        data = np.zeros((count,3))
        for n in np.arange(count):
            f = lyr.GetFeature(n)
            data[n,:] = f.geometry().GetPoint()
        return data
            
    def _get_point_spacing(self, dataset):
        """
        Take a numpy xyz array and return the min, mean, and max spacing
        between different points in the XY direction for interpolation without 
        holes.  The returned units will be the same as the provided horizontal 
        coordinate system.
        """
        count = len(dataset)
        min_dist = np.zeros(count) + np.inf
        # roll the array through, comparing all points and saving the minimum dist.
        for n in np.arange(1,count):
            tmp = np.roll(dataset,n,axis = 0)
            dist = np.sqrt(np.square(dataset[:,0] - tmp[:,0]) + np.square(dataset[:,1] - tmp[:,1]))
            idx = np.nonzero(dist < min_dist)[0]
            if len(idx) > 0:
                min_dist[idx] = dist[idx]
        return min_dist.min(), min_dist.mean(), min_dist.max()
    
    def _get_mask(self, dataset, resolution, window):
        """
        This is a hack to compute the mask.  Casiano will make this better some
        day.
        
        Return a gdal raster that has nodes where data should be populated with
        1, and all other nodes populated with the "no data" value.
        """
        # Find the bounds of the provided data
        xmin,xmax,ymin,ymax = np.nan,np.nan,np.nan,np.nan
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()
        for n in np.arange(count):
            f = lyr.GetFeature(n)
            x,y,z = f.geometry().GetPoint()
            xmin, xmax = self._compare_vals(x,xmin,xmax)
            ymin, ymax = self._compare_vals(y,ymin,ymax)
        numrows, numcolumns, bounds = self._get_nodes3(resolution, [xmin,ymin,xmax,ymax])
        # casiano's process for figuring out which nodes are to be populated
        # make a grid with the no data value
        # populate the nodes that should have data with one
        # turn the array into a gdal dataset
        mask = self._gdal_invdist_interp_points(dataset, resolution, window)
        return mask
    
    def _gdal_linear_interp_points(self, dataset, resolution, nodata = 1000000):
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.
        """
        # Find the bounds of the provided data
        xmin,xmax,ymin,ymax = np.nan,np.nan,np.nan,np.nan
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()
        for n in np.arange(count):
            f = lyr.GetFeature(n)
            x,y,z = f.geometry().GetPoint()
            xmin, xmax = self._compare_vals(x,xmin,xmax)
            ymin, ymax = self._compare_vals(y,ymin,ymax)
        numrows, numcolumns, bounds = self._get_nodes3(resolution, [xmin,ymin,xmax,ymax])
        algorithm = "linear:radius=0:nodata=" + str(int(nodata))
        interp_data = gdal.Grid('', dataset, format='MEM',
                                width = numcolumns,
                                height = numrows,
                                outputBounds = bounds,
                                algorithm=algorithm)
        return interp_data
    
    def _gdal_invdist_interp_points(self, dataset, resolution, radius, nodata = 1000000):
        """
        Interpolate the provided gdal vector points and return the interpolated
        data.
        """
        # Find the bounds of the provided data
        xmin,xmax,ymin,ymax = np.nan,np.nan,np.nan,np.nan
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()
        for n in np.arange(count):
            f = lyr.GetFeature(n)
            x,y,z = f.geometry().GetPoint()
            xmin, xmax = self._compare_vals(x,xmin,xmax)
            ymin, ymax = self._compare_vals(y,ymin,ymax)
        numrows, numcolumns, bounds = self._get_nodes3(resolution, [xmin,ymin,xmax,ymax])
        algorithm = "invdist:power=2.0:smoothing=0.0:radius1=" + str(radius)+":radius2="+str(radius)+":angle=0.0:max_points=0:min_points=1:nodata=" + str(int(nodata))
        interp_data = gdal.Grid('', dataset, format='MEM', 
                                width = numcolumns,
                                height = numrows,
                                outputBounds = bounds,
                                algorithm=algorithm)
        return interp_data
    
    def _get_natural_interp(self, dataset):
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()
        data = np.zeros((count,3))
        print (lyr, count, data)
    
    def _compare_vals(self, val, valmin, valmax):
        """
        This is a small utility for inspecting values and seeing they contribute
        to a min or max estimate.
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
    
    def _get_nodes(self, resolution, bounds):
        """
        Get the bounds and number of rows and columns. The desired resolution and
        the data max and min in X and Y are required.  This algorithm uses the 
        average of the data min and max and centers the grid on this location.
        """
        xmin, ymin, xmax, ymax = bounds
        numcolumns = int(np.ceil((1. * (xmax - xmin) / resolution)))
        xmean = np.mean([xmin,xmax])
        xnmin = xmean - numcolumns / 2. * resolution
        xnmax = xmean + numcolumns / 2. * resolution
        numrows = int(np.ceil((1. * (ymax - ymin) / resolution)))
        ymean = np.mean([ymin,ymax])
        ynmin = ymean - numrows/ 2. * resolution
        ynmax = ymean + numrows/ 2. * resolution
        return numrows, numcolumns, [xnmin,ynmin,xnmax,ynmax]
    
    def _get_nodes2(self, resolution, bounds):
        """
        Get the bounds and number of rows and columns. The desired resolution and
        the data max and min in X and Y are required.  This algorithm uses the 
        data min as the anchor for the rows and columns, thus the max data values
        may contain sparse data.
        """
        xmin, ymin, xmax, ymax = bounds
        numcolumns = int(np.ceil((1. * (xmax - xmin) / resolution)))
        xnmin = xmin
        xnmax = xmin + 1. * numcolumns * resolution
        numrows = int(np.ceil((1. * (ymax - ymin) / resolution)))
        ynmin = ymin
        ynmax = ymin + 1. * numrows * resolution
        return numrows, numcolumns, [xnmin,ynmin,xnmax,ynmax]
    
    def _get_nodes3(self, resolution, bounds):
        """
        Get the bounds and number of rows and columns. The desired resolution and
        the data max and min in X and Y are required.  This algorithm uses the 
        data min as the anchor for the rows and columns, but floored to the nearest
        multiple of the resolution.
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
        return numrows, numcolumns, [xnmin,ynmin,xnmax,ynmax]
    
    def _mask_with_raster(self, dataset, maskraster):
        """
        Read two rasters and use the nodata values from one to mask the other.
        These files are assumed to be collocated as all operation are conducted on
        the pixel level.
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
    
    def _shrink_coverage(self, dataset, resolution, radius):
        """
        Shrink coverage of a dataset by the original coverage radius.
        """
        data = dataset.ReadAsArray()
        data_rb = dataset.GetRasterBand(1)
        nodata = data_rb.GetNoDataValue()
        idx = np.nonzero(data == nodata)
        data[idx] = np.nan
        # divide the window size by the resolution to get the number of cells
        rem_cells = int(np.round(radius / resolution)) 
        # print ('Shrinking coverage back ' + str(rem_cells) + ' cells.')
        for n in np.arange(rem_cells):
            ew = np.diff(data, axis = 0)
            idx_ew = np.nonzero(np.isnan(ew))
            ns = np.diff(data, axis = 1)
            idx_ns = np.nonzero(np.isnan(ns))
            data[idx_ew] = np.nan
            data[idx_ew[0] + 1, idx_ew[1]] = np.nan
            data[idx_ns] = np.nan
            data[idx_ns[0], idx_ns[1] + 1] = np.nan
        idx = np.nonzero(np.isnan(data))
        data[idx] = nodata
        data_rb.WriteArray(data)
        return dataset
