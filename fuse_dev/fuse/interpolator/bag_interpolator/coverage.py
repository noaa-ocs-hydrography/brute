# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:51 2019

@author: Casiano.Koprowski
"""

import os as _os
import numpy as _np
from osgeo import gdal as _gdal
from osgeo import ogr as _ogr
from osgeo import osr as _osr
from osgeo import gdalconst as _gdalconst

def _maxValue(arr):
    """Returns the most used value in the array as an integer

    Takes an input array and finds the most used value in the array, this
    value is used by the program to assume the array's nodata value

    Parameters
    ----------
    arr : numpy.array
        An input array

    Returns
    -------
    maxVal : int
        Returns the most used value in the array as an integer

    """
    nums, counts = _np.unique(arr, return_counts =True)
    index = _np.where(counts==_np.amax(counts))
    return int(nums[index])


class geotiff:
    """This class serves as the main container for geotiff data.

    In order to read and assign the data relevant to the interpolation process,
    this class heavily depends on OSGEO's GDAL, OGR, and ORS libraries to open
    and populate the data needed

    Parameters
    ----------
    filepath : str
        The complete file path of the input coverage file

    Attributes
    ----------
    nodata : float
        1000000.0 is the nodata value associated with the BAG format
    elevation : numpy.array
        The elevation layer of the BAG
    uncertainty : numpy.array
        The uncertainty layer of the BAG
    shape : tuple
        The (y, x) dimensions of the layer data
    bounds : touple
        The geographic bounds of the data in order (nw, se)=([sx,ny],[nx,sy])
    wkt : str
        The WKT representation of the data CRS

    """
    def __init__(self, filename):
        _ds = _gdal.Open(filename)
        self.bounds, self.resolution = self._getBounds(_ds)
        self.array, self.shape, self.nodata = self._getArrayData(_ds)
        _fName = _os.path.split(filename)[-1]
        _splits = _fName.split('_')
        self.name = '_'.join([x for x in _splits[:2]])
        _ds=None

    def _getBounds(self, gdal_obj):
        ulx, xres, xskew, uly, yskew, yres  = gdal_obj.GetGeoTransform()
        lrx = ulx + (gdal_obj.RasterXSize * xres)
        lry = uly + (gdal_obj.RasterYSize * yres)
        if lrx<ulx:
            s = ulx
            ulx = lrx
            lrx = s
        if uly<lry:
            s = lry
            lry = uly
            uly = s
        meta = ([ulx, uly], [lrx, lry])
        return meta, (xres, yres)

    def _getArrayData(self, gdal_obj):
        band = gdal_obj.GetRasterBand(1)
        arr = band.ReadAsArray()
        band=None
        maxVal = _maxValue(arr)
        shape = arr.shape
        return arr, shape, maxVal

class geopackage:
    """This class serves as the main container for geopackage data.

    In order to read and assign the data relevant to the interpolation process,
    this class heavily depends on OSGEO's GDAL, OGR, and ORS libraries to open
    and populate the data needed

    A geopackage coverage requres a ``to_crs`` WKT spatial reference system in
    order to ensure that the contents of the geopackage can be properly
    rasterized.

    Parameters
    ----------
    filepath : str
        The complete file path of the input coverage file
    to_crs : str
        WKT object with destination spatial reference system

    Attributes
    ----------
    nodata : float
        1000000.0 is the nodata value associated with the BAG format
    elevation : numpy.array
        The elevation layer of the BAG
    uncertainty : numpy.array
        The uncertainty layer of the BAG
    shape : tuple
        The (y, x) dimensions of the layer data
    bounds : touple
        The geographic bounds of the data in order (nw, se)=([sx,ny],[nx,sy])
    wkt : str
        The WKT representation of the data CRS

    """
    def __init__(self, filename, to_crs, pixel_size=1, nodata=255):
        self.wkt = to_crs
        self.resolution = pixel_size
        self.nodata = 255
        self.array, self.shape, self.bounds = self._readAndRasterize(filename,
                                                    to_crs, pixel_size, nodata)

    def _readAndRasterize(self, filename, to_crs, pixel_size, nodata):
        fName = _os.path.split(filename)[-1]
        splits = _os.path.splitext(fName)
        name = splits[0]

        # Create destination Spatial Ref from wkt
        to_srs = _osr.SpatialReference(wkt=to_crs)

        # Open the data source and read in the extent
        source_ds = _ogr.Open(filename)
        source_layer = source_ds.GetLayer()
        source_srs = source_layer.GetSpatialRef()

        for feature in source_layer:
            if feature != None:
                geom = feature.GetGeometryRef()
#                print (geom.ExportToWkt())
                ds_geom = _ogr.CreateGeometryFromWkt(geom.ExportToWkt())
#                print (source_srs, to_proj, sep='\n')
                coordTrans = _osr.CoordinateTransformation(source_srs, to_srs)
                ds_geom.Transform(coordTrans)
                driver = _ogr.GetDriverByName('Memory')
                ds = driver.CreateDataSource('temp')
                layer = ds.CreateLayer(name, to_srs, _ogr.wkbMultiPolygon)

                # Add one attribute
                layer.CreateField(_ogr.FieldDefn('Survey', _ogr.OFTString))
                defn = layer.GetLayerDefn()

                # Create a new feature (attribute and geometry)
                feat = _ogr.Feature(defn)
                feat.SetField('Survey', name)

                # Make a geometry, from Shapely object
                feat.SetGeometry(ds_geom)

                layer.CreateFeature(feat)
                feat = geom = None  # destroy these
                break


        x_min, x_max, y_min, y_max = ds_geom.GetEnvelope()
        bounds = ([x_min,y_max],[x_max,y_min])

        # Create the destination data source
        x_dim = int((x_max - x_min) / pixel_size)
        y_dim = int((y_max - y_min) / pixel_size)
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_dim, y_dim,
                                                        _gdal.GDT_Byte)
        gt = (x_min, pixel_size, 0, y_max, 0, -pixel_size)
        target_ds.SetGeoTransform(gt)
        target_ds.SetProjection(to_srs)
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)

        # Rasterize
        _gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])

        array = band.ReadAsArray()
        shape = array.shape

        band = None
        source_ds = None
        target_ds = None

        return array, shape, bounds