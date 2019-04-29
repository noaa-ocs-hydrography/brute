# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 16:52:49 2018

@author: Casiano.Koprowski

Made in part with code from:
    - hyo.bag Library Pydro by HydrOffice
    - "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
        (https://gis.stackexchange.com/a/278965)]
    - "How to get raster corner..." on on GIS Stack Exchange [Answer by 'James'
        (https://gis.stackexchange.com/a/201320)]

Icon from https://www.ngdc.noaa.gov/mgg/image/2minrelief.html
    via https://commons.wikimedia.org/wiki/File:Atlantic_bathymetry.jpg
"""

import os as _os
import re as _re
import numpy as _np
import tables as _tb
import time as _time
import scipy as _scipy
import shutil as _shutil
import lxml.etree as _et
import matplotlib.pyplot as _plt
import astropy.convolution as _apc

from datetime import datetime as _dt
from string import ascii_lowercase as _al
from scipy.ndimage.interpolation import zoom as _zoom

from osgeo import gdal as _gdal
from osgeo import ogr as _ogr
from osgeo import osr as _osr
from osgeo import gdalconst as _gdalconst

#import autointerp_ui

# this allows GDAL to throw Python Exceptions
_gdal.UseExceptions()

progLoc = _os.getcwd()
print (progLoc)

ns = {
    'bag': 'http://www.opennavsurf.org/schema/bag',
    'gco': 'http://www.isotc211.org/2005/gco',
    'gmd': 'http://www.isotc211.org/2005/gmd',
    'gmi': 'http://www.isotc211.org/2005/gmi',
    'gml': 'http://www.opengis.net/gml/3.2',
    'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
}

ns2 = {
    'gml': 'http://www.opengis.net/gml',
    'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
    'smXML': 'http://metadata.dgiwg.org/smXML',
}

catZones = {
        'A1':[.01,.5],
        'A2/B':[.02,1],
        'C':[.05,2]}

combo = _re.compile(r'COMBINED', _re.IGNORECASE)
interp = _re.compile(r'INTERP', _re.IGNORECASE)

def getTifElev(files):
    '''This function takes a list of GeoTiff filepaths and opens each in order
    to access their values. Returned as a list of tif objects are the order a
    file was found, the path of the file, and the file's values as
    [x, file, extent, arr] and a list of file names

    Upper Left Corner and Lower Right Corner Bounds Directly from:
    "How to get raster corner..." on on GIS Stack Exchange [Answer by 'James'
    (https://gis.stackexchange.com/a/201320)]

    Parameters
    ----------
    files : list
        File paths of the input GeoTiff files

    Returns
    -------
    tifFiles : list
        List of GeoTiff file objects
        Position loaded, file name, raster extents, numpy.array
    names : list
        list of names of imported GeoTiffs

    '''
    print ('getTifElev')
    tifFiles = []
    names = []
    maxVal = 0
    y = 0
    for file in files:
        print (file)
        ds = _gdal.Open(file)
#        print (_gdal.Info(ds))
        ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
        print (xres, yres)
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)
        print (ulx, lrx)
        print (uly, lry)
        if lrx<ulx:
            s = ulx
            ulx = lrx
            lrx = s
        if uly<lry:
            s = lry
            lry = uly
            uly = s
        meta = ([ulx, uly], [lrx, lry])

        fName = _os.path.split(file)[-1]
        splits = fName.split('_')
        name = '_'.join([x for x in splits[:2]])
        names.append(name)
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        print (arr)
        maxVal = maxValue(arr)
        print (maxVal)
        [cols, rows] = arr.shape
        print(arr.shape)
        print (_np.amax(arr), _np.amin(arr))
        print (meta)
        band=None
        ds=None
        tifFiles.append ([y, file, meta, arr])
        y += 1
    return tifFiles, names


def getBagLyrs(fileObj):
    """This function takes a list of BAG filepaths and opens each in order
    to access their values. Returned is a list of bag objects are the order a
    file was found, the path of the file, the file's elevation values, and the
    file's uncertainty values as [x, file, extent, elev] and a list of file
    names

    Parameters
    ----------
    files : string
        File path of the input BAG file

    Returns
    -------
    bag : list
        BAG file object
        Position loaded, file name, grid extents, file size, grid dimensions,
        uncertainty numpy.array, bathymetry numpy.array
    maxVal : int
        Maximum occurance value AKA nodata value

    """
    print ('getBagLyrs')
#    bagFiles = []
    names = []
    y = 0
    print (fileObj)
    fName = _os.path.split(fileObj)[-1]
    splits = fName.split('_')
    name = '_'.join([x for x in splits[:2]])
    names.append(name)
    print (splits)
    size = int(_np.round(_os.path.getsize(fileObj)/1000))
    print (size)
    with _tb.open_file(fileObj, mode = 'r') as bagfile:
        meta_read = [str(x, 'utf-8', 'ignore') for x in bagfile.root.BAG_root.metadata.read()]
#        print (meta_read)
        meta_xml = ''.join(meta_read)
#        print (meta_xml)
        encodeVal = 0
        for x in meta_xml:
            if meta_xml[encodeVal] == '>':
                meta_xml = meta_xml[encodeVal:]
                break
            else:
                encodeVal += 1
        startVal = 0
        for x in meta_xml:
            if meta_xml[startVal] == '<':
                meta_xml = meta_xml[startVal:]
                break
            else:
                startVal += 1
        xml_tree = _et.XML(meta_xml)
        a = 0
        try:
            ret = xml_tree.xpath('//*/gmd:spatialRepresentationInfo/gmd:MD_Georectified/'
                                      'gmd:cornerPoints/gml:Point/gml:coordinates',
                                      namespaces=ns)[0].text.split()
            a = 1
        except (_et.Error, IndexError) as e:
            try:
                ret = xml_tree.xpath('//*/spatialRepresentationInfo/smXML:MD_Georectified/'
                                          'cornerPoints/gml:Point/gml:coordinates',
                                          namespaces=ns2)[0].text.split()
                a = 2
            except (_et.Error, IndexError) as e:
                print("unable to read corners SW and NE: %s" % e)
                return

        try:
            sw = [float(c) for c in ret[0].split(',')]
            ne = [float(c) for c in ret[1].split(',')]
        except (ValueError, IndexError) as e:
            print("unable to read corners SW and NE: %s" % e)
            return
        print (a, ret, sw, ne)
        sx,sy = sw
        nx,ny = ne
        meta = [[sx,ny],[nx,sy]]
        elev = _np.flipud(bagfile.root.BAG_root.elevation.read())
        uncr = _np.flipud(bagfile.root.BAG_root.uncertainty.read())
        shape = elev.shape
#        maxVal = _np.amax(elev)
        maxVal = 1000000.0
        bag = [y, fileObj, meta, size, shape, uncr, elev]
#        bag = [y, fileObj, meta, size, shape]
        print (bag, maxVal)
        bagfile.close()
    return bag, maxVal, bagfile

def getShpRast(file):
    """Import shapefile

    Parameters
    ----------
    file : string
        Shapefile file location

    Todo
    ----
    Inital working idea for importing shapefiles as rasters

    """
    print ('getShpRast')
    file = _ogr.Open(file)
    shape = file.GetLayer()
    fnum = shape.GetFeatureCount()
    x = _np.arrange(fnum)
    polys = []
    for n, x in enumerate(_np.arange(fnum)):
        polys.append(shape.GetFeature(n))
    print (polys)

def write_raster(raster_array, gt, data_obj, outputpath, dtype=_gdal.GDT_UInt32,
                 options=0, color_table=0, nbands=1, nodata=False):
    """Directly From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
    (https://gis.stackexchange.com/a/278965)]

    Parameters
    ----------
    raster_array : numpy.array
        Array to be written to a GeoTiff file
    gt : tuple, gdal.GeoTransform
        Norhtern extent, resolution, 0.0, Western extent, 0.0, -resolution)
    data_obj : gdal.RasterBand
        gdal.RasterBand
    outputpath : string
        Folder to save the GeoTiff raster

    """
    print('write_raster')

    height, width = raster_array.shape

    # Prepare destination file
    driver = _gdal.GetDriverByName("GTiff")
    if options != 0:
        dest = driver.Create(outputpath, width, height, nbands, dtype, options)
    else:
        dest = driver.Create(outputpath, width, height, nbands, dtype)

    # Write output raster
    if color_table != 0:
        dest.GetRasterBand(1).SetColorTable(color_table)

    dest.GetRasterBand(1).WriteArray(raster_array)

    if nodata is not False:
        dest.GetRasterBand(1).SetNoDataValue(nodata)

    # Set transform and projection
    dest.SetGeoTransform(gt)
    wkt = data_obj.GetProjection()
    srs = _osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dest.SetProjection(srs.ExportToWkt())

    # Close output raster dataset
    dest = None

def maxValue(arr):
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
    print ('maxValue')
    nums, counts = _np.unique(arr, return_counts =True)
    index = _np.where(counts==_np.amax(counts))
    print (index, nums[index])
    return int(nums[index])


def tupleGrid(grid, maxVal):
    """Takes an input matrix and an assumed nodata value. The function iterates
    through the matrix and compiles a list of 'edge' points [[x, y, z], ...]
    where:

    1. the current value is not a nodata value and previous value was a nodata value.
        - sets io True, indicating that the next value to compare against should be a nodata value.
    2. the current value is a nodata value and the previous value was not a nodata value.
        - sets io False, indicating that the next value to compare against should not be a nodata value.
    3. returns a list of the points found.

    Parameters
    ----------
    grid : numpy.array
        An input array
    maxVal : int
        The array's nodata value

    Returns
    -------
    np.array
        Array of indecies where nodata values meet data values
        in order x, y, z

    """
    print ('tupleGrid')
    points = []
    a = 0
    for x in range(grid.shape[1]):
        io = False
        for y in range(grid.shape[0]):
            if grid[y,x] == maxVal:
                if grid[y-1,x] != maxVal:
                    val = grid[y-1,x]
                    point = [x, y-1, val]
                    if a == 1:
                        print (point, val)
                        a += 1
                    points.append(point)
                    io = False
                else:
                    pass
            else:
                if io == False:
                    val = grid[y,x]
                    point = [x, y, val]

                    if a == 0:
                        print (point, val)
                        a += 1
                    points.append(point)
                    io = True
    return _np.array(points)

def concatGrid(grids, maxVal, shape):
    """Takes an input of an array of grid objects and the assumed nodata value
    Passes the assumed nodata value and the arrays held within each of the
    listed grid objects to tupleGrid(grid, maxVal) for a return of an array of
    edge points for each grid [[x, y, z], ...]

    Takes the results of both tupleGrid calls and combines them into a single
    array of edge points.  Then uses .view().sort() to sort the combined
    products [[x, y, z], ...] based on ascending x, then column y.

    The results of the sorted point array are then split into a list of points
    [[x, y], ...] and values [[z], ...] and returned

    Parameters
    ----------
    grids : list
        The BAG and GeoTiff objects
    maxVal : int
        The BAG data's nodata value
    shape : tuple
        Dimensions of the input BAG data (y, x)

    Returns
    -------
    grid : numpy.array
        Array of indecies where nodata values meet data values
        in order x, y
    vals : numpy.array
        Array of values where nodata values meet data values
        in order z

    """
    print ('concatGrid')
    print ('tpts', _dt.now())
    tpts = tupleGrid(grids[0], maxVal)
    print ('bpts', _dt.now())
    bpts = tupleGrid(grids[1], maxVal)
    print ('done', _dt.now())
    print ('combo1', _dt.now())
    if len(tpts) == 0 and len(bpts) != 0:
        comb = bpts
        comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
        print (comb)
        grid = _np.hsplit(comb, [2, 4])
        vals = grid[1].squeeze()
        grid = grid[0]
    elif len(bpts) == 0 and len(tpts) != 0:
        comb = tpts
        comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
        print (comb)
        grid = _np.hsplit(comb, [2, 4])
        vals = grid[1].squeeze()
        grid = grid[0]
    elif len(tpts) == 0 and len(bpts) == 0:
        grid, vals = [], []
    else:
        comb = _np.concatenate([tpts, bpts])
        comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
        print (comb)
        grid = _np.hsplit(comb, [2, 4])
        vals = grid[1].squeeze()
        grid = grid[0]
    print (grid, vals)
    return grid, vals

def alignTifs(tifs):
    """Takes an input of an array of tiff objects. The goal of this function
    is to fit the provided tifs to the largest combined area of the inputs.


    1. Find NW and SE corners of the first input array and the number of rows and columns in the input.
    2. Find NW and SE corners of subsiquent input arrays and record differences as the NW differences (nxd, nyd) SE differences (sxd, syd).

        - also modify the original NW and/or SE corner values if applicable.
        - also modify the original number or rows and columns if applicable.

    3. Resize the area of each input array to the new total combined size of all input arrays based on the stored number of rows and columns.
    4. Reposition the data within each input array to maintain original geographic placement based on the difference of the input array's original NW corner and the recorded 'largest' NW boundry of the combined input areas.
    5. Return tiff objects with the newly resized/repositioned tifs and the new NW and SE extents of the new combined area.


    Parameters
    ----------
    tifs : list
        List of GeoTiff objects created by :func:`getTifElev`

    Returns
    -------
    sizedTifs : list
        GeoTiff file objects with numpy.array raster data replaced with data
        that is alligned within the maximum extents of all input GeoTiffs
    tifs : list
        Replaces sizedTifs if there is no resizing needed
    ext : list
        The maximum extents of all input GeoTiffs

    """
    print ('alignTifs')
    x = 0
#    maxVal = 0
    nw, se = [0,0], [0,0]
    rows, cols = 0, 0
    sxd, nyd = 0, 0
    nxd, syd = 0, 0
    for tif in tifs:
        bounds = tif[-2]
        ul = bounds[0]
        lr = bounds[-1]
        if x == 0:
            nw = ul
            se = lr
            print (nw, se)
            x += 1
        else:
            ulx, uly = ul
            lrx, lry = lr
            nwx, nwy = nw
            scx, scy = se
            if ulx <= nwx:
                nxd = nwx - ulx
                nwx = ulx
            elif ulx >= nwx:
                nxd = ulx - nwx
#                nwx = nwx
            if uly >= nwy:
                nyd = uly - nwy
                nwy = uly
            elif uly <= nwy:
                nyd = nwy - uly
#                nwy = nwy
            if lrx >= scx:
                sxd = lrx - scx
                scx = lrx
            elif lrx <= scx:
                sxd = scx - lrx
#                scx = scx
            if lry <= scy:
                syd = lry - scy
                scy = lry
            elif lry >= scy:
                syd = scy - lry
#                scy = scy
            print (nxd, nyd, sxd, syd)
            print (sxd, nyd)
            nw = [nwx, nwy]
            se = [scx, scy]
            cols, rows = int(_np.round(nwy - scy)), int(_np.round(scx - nwx))
            print('cols: ' , nwy - scy, '\nrows: ', scx - nwx)
    print (nw, se)
    sizedTifs = []
    print ('resize?')
    if nxd != 0 or nyd != 0 or sxd != 0 or syd != 0:
        print ('yes', _dt.now())
        ref = _np.full((cols, rows), 0)
        print (ref.shape)
        for tif in tifs:
            maxVal = maxValue(tif[-1])
            sizedTif = tif[:-1]
            print (tif)
            bndx, bndy = tif[2][0]
            nwx, nwy = nw
            print ('old:', bndx, bndy)
            print ('new', nwx, nwy)
            arr = tif[-1]
            bef = arr.shape
            y, x = arr.shape
            print (y, rows)
            print (x, cols)
            if x != cols:
                exp = int(abs(rows - x))
                print (exp)
                add = _np.full((y, exp), maxVal)
                arr = _np.column_stack([arr, add])
                y, x = arr.shape
            if y != rows:
                exp = int(abs(cols - y))
                print (exp)
                add = _np.full((exp, x), maxVal)
                arr = _np.vstack([arr, add])
            print (bef, arr.shape)
            if nwx != bndx:
                rollx = bndx - nwx
                print (rollx)
                arr = _np.roll(arr, int(rollx), axis=1)
            if nwy != bndy:
                rolly = nwy - bndy
                print (rolly)
                arr = _np.roll(arr, int(rolly), axis=0)
            sizedTif.append(arr)
            sizedTifs.append(sizedTif)
#            sizedTifs.append(arr)
        ext = [nw, se]
        print (ext)
        print ('done', _dt.now())
        return sizedTifs, ext
    else:
        print ('same')
        ext = [nw, se]
        return tifs, ext

def polyTifVals(tifs, path, names, extent):
    """Heavy Influence From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
    (https://gis.stackexchange.com/a/278965)]

    This function takes an array input of tif objects, a destination path (for
    ouput saving), and a list the names of the input tif objects.

    Takes the arrays of each tif object and combines them along the z axis of
    each::

        [[za,zb],[za,zb],[za,zb],
         [za,zb],[za,zb],[za,zb],
         [za,zb],[za,zb],[za,zb]]

    Then takes the mean of the values along the z axis. The reslult is then
    modified to a binary raster by making all values that are not the nodata
    value 1 and the values that are nodata 0.

    This binary raster is saved at the input path and also returned with the
    new full path for the output

    Parameters
    ----------
    tifs : list
        List of GeoTiff objects
    path : string
        Folder path for the output files
    names : list
        List of names of imported GeoTiffs
    extent : list
        The maximum extents of all input GeoTiffs

    Returns
    -------
    meanTiff : numpy.array
        Binary raster only showing data and nodata areas for all of the input
        GeoTiff data
    outputtiff : string
        File name and path generated for saving meanTiff data in Tiff format
    outputhdf5 : string
        File name and path generated for saving meanTiff data in HDF5 format
    targs : list
        Additional arguments passed back for use in :func:`alignGrids`

    """
    print ('polyTifVals')
    letters = []
    for letter in _al:
        letters.append(letter)
    maxRes = [0,0]
    maxVal = 0
    tifList = []
    x = 0
    for tif in tifs:
        print (tif[-1].shape)
        maxVal = maxValue(tif[-1])
        print (maxVal)
        [rows, cols] = tif[-1].shape
        if x == 0:
            print ('original', maxRes)
            if rows >= maxRes[0]:
                maxRes[0] = rows
            if cols >= maxRes[1]:
                maxRes[1] = cols
            print ('original', maxRes)
            tifList.append([tif[0], tif[-1]])
            x += 1
        elif x != 0:
            if [tif[-1].shape[0], tif[-1].shape[1]] != maxRes:
                print ('nope')
                pass
            else:
                tifList.append([tif[0], tif[-1]])
                print ('yup')
        tifDict = dict(tifList)
        for i, tif in tifDict.items():
            array = _np.expand_dims(tif,2)
            if i == 0:
                allarrays = array
            else:
                allarrays = _np.concatenate((allarrays, array), axis=2)

        meanTiff = _np.nanmean(allarrays, axis=2)
        if maxVal != 0:
            meanTiff = (meanTiff < maxVal).astype(_np.int)
        else:
            meanTiff = (meanTiff > maxVal).astype(_np.int)

    outputname = path + '\\' + names[0] + '_COMBINEDPOLY'
    outputtiff = outputname + '.tif'
    outputhdf5 = outputname + '.h5'

    print(outputname)
    while True:
        if _os.path.exists(outputtiff):
            _os.remove(outputtiff)
        elif not _os.path.exists(outputtiff):
            break
        if _os.path.exists(outputhdf5):
            _os.remove(outputhdf5)
        elif not _os.path.exists(outputhdf5):
            break
    ny, nx = extent[0]
    reso = 1
    gtran = (ny, reso, 0.0, nx, 0.0, -(reso))
#    write_raster(meanTiff, gtran, gd_obj, outputtiff, options = ['COMPRESS=LZW'])
#    gd_obj = None

    targs = [gtran, tifs[-1][1], outputtiff]
    print ('done?')
    return meanTiff, outputtiff,  outputhdf5, targs

def shift5(arr, y, x, fill_value=_np.nan):
    """https://stackoverflow.com/questions/30399534/shift-elements-in-a-numpy-array"""
    print ('shift5', y)
    result = _np.empty_like(arr)
    if y > 0:
        result[:y,:] = fill_value
        result[y:,:] = arr[:-y,:]
    elif y < 0:
        result[y:,:] = fill_value
        result[:y,:] = arr[-y:,:]
    else:
        result = arr
    print ('shift6', x)
    if x > 0:
        result[:,:x] = fill_value
        result[:,x:] = arr[:,:-x]
    elif x < 0:
        result[:,x:] = fill_value
        result[:,:x] = arr[:,-x:]
    else:
        result = arr
    return result

def alignGrids(bag, tif, maxVal, targs):
    """Takes an input of two arrays representing bag and tif data. These arrays
    hold information like extent, data, and more. The goal of this function is
    to fit and/or shift the provided tif to the area of the bag.

    Steps:

    1. Determine/confirm the resolution of the tif using a mix of file naming, extent, and dimensions.
    2. Determine the resolution of the bag using file naming conventions.
    3. Calculate the tif/bag resolution to determine the value needed to resize the tif.
    4. Resize the tif so that the resolution matches the bag using.
    5. Convert tif values to either numpy.nan or maxVal for use later.
    6. Store the extents of the bag and tif and calculate the difference.
    7. Create an empty array with dimensions equal to the original extents of the tiff data plus the additional space needed to match the total extents of the BAG data.
    8. Align the content of the tiff as needed in order to match the bag and apply the shifted content to the empty array created in Step 7.
    9. Enforce the size of the bag onto the shifted tiff data stored in the array created and modified in Step(s) 7 and 8.


    Parameters
    ----------
    bag : list
        BAG file object
        Position loaded, file name, grid extents, file size, grid dimensions,
        uncertainty numpy.array, bathymetry numpy.array
    tif : list
        GeoTiff file objects with numpy.array raster data replaced with binary
        data only showing data and nodata areas and that is alligned within
        the maximum extents of all input GeoTiffs
    maxVal : int
        Maximum occurance value AKA nodata value
    targs : list
        Additional arguments generated in :func:`polyTifVals` for use in saving
        the BAG aligned Tiff data

    Returns
    -------
    bagRes : float
        Resolution of the BAG data
    grids : list
        The BAG and GeoTiff objects, GeoTiff object's numpy.array raster data
        replaced with BAG aligned Tiff data
    bagBounds : list
        Geographic extent of the BAG data
    bShape : tuple
        Dimensions of the input BAG data (y, x)

    """
    print ('alignGrids')
    print (bag[-1])

    ##1
    tu, tl = tif[-2]
    ty, tx = tif[-1].shape
    tex, tey = int(tl[0] - tu[0]), int(tu[1] - tl[1])
    print (tx, ty)
    print (tex, tey)
    if tex-ty == 0 and tey-tx == 0:
        tifRes = 1
        print ('same tif res',tifRes)
    else:
        tifRes = _np.round(_np.mean([tex/tx, tey/ty]))
        print ('diff tif res',tifRes)

    ##2
    splits = _os.path.split(bag[1])[-1]
    print (splits)
    splits = splits.split('_')[2]
    print (splits)
    bagRes = _re.split('(\d+)',splits)[-2:]
    print (bagRes)

    ## 3
    if bagRes[-1] == 'cm':
        bagRes = int(bagRes[0]) / 100
        zres = tifRes/bagRes
    elif bagRes[-1] =='m':
        bagRes = int(bagRes[0])
        zres = tifRes/bagRes
    print (bagRes, zres)
    _plt.imshow(tif[-1][::100,::100])
    _plt.show()
    ## 4
    print (tif[-1])
    if zres == 1:
        newarr = tif[-1]
    else:
        print ('_zoom', _dt.now())
        newarr = _zoom(tif[-1], zoom=[zres, zres], order=3, prefilter=False)
        print ('zoomed', _dt.now())
    _plt.imshow(newarr[::100,::100])
    _plt.show()

    ## 5
    newarr = newarr.astype('float64')
    newarr[newarr > 0] = _np.nan
    newarr[newarr < 1] = float(maxVal)
    _plt.imshow(newarr[::100,::100])
    _plt.show()
    print (newarr)
    print (tif[-1].shape, newarr.shape)

    ## 6
    bagBounds = bag[2]
    tifBounds = tif[2]
    print (bagBounds)
    print (tifBounds)
    bulx, buly = bagBounds[0]
    blrx, blry = bagBounds[-1]
    tulx, tuly = tifBounds[0]
    tlrx, tlry = tifBounds[-1]
    dulx, duly, dlrx, dlry = 0, 0, 0, 0
    if bulx != tulx:
        dulx = tulx - bulx
        print (bulx, tulx, dulx)
    if buly != tuly:
        duly = buly - tuly
        print (buly, tuly, duly)
    if blrx != tlrx:
        dlrx = blrx - tlrx
        print (blrx, tlrx, dlrx)
    if blry != tlry:
        dlry = tlry - blry
        print (blry, tlry, dlry)
    print (dulx, duly)
    print (dlrx, dlry)

    ## 7
    bShape = bag[4]
    bSy, bSx = bShape
    tSy, tSx = newarr.shape
    print (bSy, tSy)
    print (bSx, tSx)
    expx, expy = 0, 0
    if newarr.shape != bShape:
        print (bSy - tSy, bSx - tSx)
        if tSy < bSy:
            expy = int(_np.abs(bSy - tSy))
            print ('expy', expy)
        if tSx < bSx:
            expx = int(_np.abs(bSx - tSx))
            print ('expx', expx)
    ay = _np.full((tSy+expy,tSx+expx), maxVal)
    print ('expz', ay.shape, bShape)
    rollx = int(dulx*zres)
    rolly = int(duly*zres)

    ## 8
    up, left = 0, 0
    down, right = 0, 0
    if duly < 0:
        up = -int(rolly)
    elif duly > 0:
        down = rolly
        up = 0
    if dulx < 0:
        left = -int(rollx)
    elif dulx > 0:
        right = rollx
        left = 0

    if dulx != 0 or duly != 0:
        print ('rollz', up, left, down, right)
        temp = newarr[up:,left:]
        print (temp.shape)
        _plt.imshow(temp[::100,::100])
        _plt.show()
        ay[down:temp.shape[0]+down,right:temp.shape[1]+right] = temp[:,:]
        temp = None
    else:
        ay[:] = newarr[:]
    print ('expz', ay.shape)
    ax = _np.full(bShape, maxVal)
    ax[:] = ay[:bSy,:bSx]
    newarr = None
    ay = None

    ## 9
    gd_obj = _gdal.Open(targs[1])
    tif.pop()
    print (tif)
    tif.append(ax)
    print (tif, tif[-1].shape)
#    print ('h5')
#    outputhdf5 = tif[1]
#    f = _tb.open_file(outputhdf5, 'w')
#    atom = _tb.Atom.from_dtype(ax.dtype)
#    ds = f.create_carray(f.root, 'data', atom, bShape)
#    ds[:] = ax
#    f.close()
    print ('tiff')
    temp = targs[0]
    gt = (bulx, bagRes, temp[2], buly, temp[4], -(bagRes))
    write_raster(ax, gt, gd_obj, targs[2], options=['COMPRESS=LZW'])
    ax = None

    grids = [tif, bag]
    bShape = bag[-1].shape
    return bagRes, grids, bagBounds, bShape

def comboGrid(grids):
    """Prepares the tif and bag data to be split in concatGrid and returns
    the results

    Parameters
    ----------
    grids : list
        The BAG and GeoTiff objects

    Returns
    -------
    combo : numpy.array
        Array of indecies where nodata values meet data values
        in order x, y
    vals : numpy.array
        Array of values where nodata values meet data values
        in order z


    """
    print ('comboGrid')
    maxVal = maxValue(grids[1][-1])
    shape = grids[1][-1].shape
    tif = grids[0][-1]
    bag = grids[1][-1]
    arrs = [tif, bag]
    combo, vals = concatGrid(arrs, maxVal, shape)
    return combo, vals

def rePrint(grids, ugrids, maxVal, ioVal, debug=0):
    """Uses a mix of interpolated and original bag and tif data in order to
    determine where new interpolated data should be applied. No interpolated
    data is used where original bag data exists.


    bag : numpy.array
        orignal bag data
    interp : numpy.array
        interpolated bag data (smoothed)
    pbag : numpy.array
        interpolated bag data (before smooth)
    poly : numpy.array
        binary raster of tif coverage
    uncr : numpy.array
        original uncertainty
    iuncrt : numpy.array
        interpolated uncertainty
    maxVal : numpy.array
        bag no-data value
    ioVal : numpy.array
        user option to include all data or interpolated data only


    Steps:

    1. Create a binary grid of tif dat, tpoly.
    2. Create a binary grid of bag data, bpoly.
    3. Create a combined grid of the "complete coverage" of data; cpoly, numpy.logical_or of bpoly, tpoly.
    4. Use bpoly to remove the bag coverage from the combined "complete coverage"; dpoly, numpy.logical_xor of bpoly, cpoly.
    5. Use dpoly to determine where interpolated data can be applied and apply original bag data everywhere else; ibag, numpy.where of dpoly, pbag, bag.
    6. Create a binary grid of ibag, npoly.
    7. Use npoly and dpoly to create a binary grid where only data is present; fpoly, numpy.logical_and of dpoly, nopoly.
    8. Finalize results by using fpoly to apply interpolated data where appropriate and original bag data everywhere else.


    Parameters
    ----------
    grids : list
        List of original data objects for BAG and GeoTiff data
    ugrids : list
        List of interpolated data objects for BAG depth and uncertainty data
    ioVal : int, bool
        User input. Determines whether origninal and interpolated or only
        interpolated data is output
    debug : int, bool
        Whether or not polyList includes all or only the last step in the
        evaluation process

    Returns
    -------
    nbag : numpy.array
        Interpolated BAG bathymetry
    nunc : numpy.array
        Interpolated BAG uncertainty
    polyList : list
        List of numpy.array objects containting all or only the last step in
        the evaluation process

    """
    print ('rePrint', _dt.now())
    print (maxVal)
    poly = grids[0][-1]
    bag = grids[-1][-1]
    uncr = grids[-1][-2]
    interp = ugrids[0]
    iuncrt = ugrids[1]
    pbag = ugrids[2]
    rows, cols = bag.shape
    ## 1
    tpoly = _np.nan_to_num(poly)
    tpoly = (tpoly < maxVal).astype(_np.int)
    ## 2
    bpoly = (bag < maxVal).astype(_np.int)
    ## 3
    cpoly = _np.logical_or(bpoly, tpoly)
    ## 4
    dpoly = _np.logical_xor(bpoly, cpoly)
    ## 5
    ibag = _np.where(dpoly, pbag, bag)
    ## 6
    npoly = (ibag < maxVal).astype(_np.int)
    ## 7
    fpoly = _np.logical_and(dpoly, npoly)
    ## 8
    if ioVal == False:
        nbag = _np.where(fpoly, interp, bag)
        nunc = _np.where(fpoly, iuncrt, uncr)
    elif ioVal == True:
        nbag = _np.where(fpoly, interp, maxVal)
        nunc = _np.where(fpoly, iuncrt, maxVal)
    print ('done', _dt.now())
    if debug == 0:
        polyList = [fpoly]
        return nbag, nunc, polyList
    elif debug == 1:
        polyList = [tpoly,bpoly,cpoly,dpoly,npoly,fpoly,ibag]
        return nbag, nunc, polyList

def triangulateSurfaces(grids, combo, vals, uval):
    """Interpolates input data and convolves the ouput of the interpolation

    Takes input BAG and GeoTiff objects; tile or complete data.  This data is
    used to inform the shape/size of the resulting output of the interpolation
    function. If combo and vals are empty, the interpolation process is skipped
    and the tile or complete data is passed back from the function. Otherwise,
    combo, vals and the shape of the BAG data input are passed to
    :func:`scipy.interpolate.griddata` and output to the variable "grid_pre".
    The output of this function is also saved to a seperate variable "grid"
    while the original output is left untouched. The seperate output is then
    passed to :func:`astropy.convolution.convolve` to smooth the output.

    The uncertainty layer is calculated using uval and grid::

        >>> m, b = uval
        >>> uncr = (grid*m)+b

    Parameters
    ----------
    grids : list
        The BAG and GeoTiff objects
    combo : numpy.array
        Array of indecies where nodata values meet data values
        in order x, y
    vals : numpy.array
        Array of values where nodata values meet data values
        in order z
    uval : tuple
        Values for uncertainty calculation

    Returns
    -------
    grid : numpy.array
        The interpolated and convolved bathemetry
    uncr : numpy.array
        The calculated uncertainty from the interpolated and convolved
        bathemetry
    grid_pre : numpy.array
        The interpolated bathemetry before convolution is performed

    """
    print ('triangulateSurfaces')
    bagObj = grids[-1]
    bag = bagObj[-1]
    uncr = bagObj[-2]
    tifObj = grids[0]
    poly = tifObj[-1]
    maxVal = maxValue(bag)
    m, b = uval
    x, y = _np.arange(bag.shape[1]), _np.arange(bag.shape[0])
    xi, yi = _np.meshgrid(x, y)
    if len(combo) != 0 or len(vals) != 0:
        print ('try tri', _dt.now())
        grid_pre = _scipy.interpolate.griddata(combo, vals, (xi, yi),
                                            method='linear', fill_value=maxVal)
        print ('done', _dt.now())
        grid = grid_pre
        grid = _np.asarray(grid, dtype='float64')
        print (grid.shape, poly.shape)
        grid[grid>0] = _np.nan
        print ('interp')
        kernel = _apc.Gaussian2DKernel(3)
        grid = _apc.convolve(grid,kernel)
        print ('interp done')
        grid[_np.isnan(grid)]=maxVal
        grid[grid>=0] = maxVal
        uncr = (grid*m)+b
    else:
        grid = bag
        uncr = uncr
        grid_pre = bag
    return grid, uncr, grid_pre

def bagSave(bag, new, tifs, res, ext, path, newu, polyList, ioVal):
    """Primary function for saving final products of the tool.

    Extended description of function.

    Parameters
    ----------
    bag : list
        BAG file object
    new : numpy.array
        Origninal and interpolated BAG data or only interpolated data
    tifs : list
        GeoTiff file object
    res : float
        Resolution of input BAG data
    ext : list
        Extent of the input BAG data [NW, SE] in [x, y]
    path : str
        Output folder path
    newu : numpy.array
        Origninal and interpolated Uncertainty data or only interpolated data
    polyList : list
        List of binary grids generated by :func:`rePrint` for
        output for use in validating or debugging products
    ioVal : int, bool
        User input. Determines whether origninal and interpolated or only
        interpolated data is output

    """
    for tif in tifs:
        gd_obj = _gdal.Open(tif[1])
        break
    print ('bagSave')
    oextX, oextY = bag[2]
    ny, nx = ext[0]
    sy, sx = ext[-1]
    reso = float(res)
    print (sx, sy, nx, ny)
    gtran = (ny, reso, 0.0, nx, 0.0, -(reso))
    print (gtran)
    fName = bag[1].split('\\')[-1]
    print (fName)
    split = fName.split('_')[:2]
    if res < 1:
        res = str(int(res*100)) + 'cm'
    else:
        res = str(res) + 'm'
    if ioVal == 1:
        bagName = '_'.join([x for x in split]) + '_' + res + '_INTERP_ONLY'
    else:
        bagName = '_'.join([x for x in split]) + '_' + res + '_INTERP_FULL'
    outputpath2 = path + '\\' + bagName + '.bag'
#    print(outputpath)
    while True:
#        if _os.path.exists(outputpath):
#            _os.remove(outputpath)
#        elif not _os.path.exists(outputpath):
#            break
        if _os.path.exists(outputpath2):
            _os.remove(outputpath2)
        elif not _os.path.exists(outputpath2):
            break
    for num in range(len(polyList)):
        outputpath = path + '\\' + bagName + '_' + str(num) +'.tif'
        print(outputpath)
        write_raster(polyList[num], gtran, gd_obj, outputpath,
                     dtype=_gdal.GDT_Float64, nodata=0, options = [ 'COMPRESS=LZW' ])
    polyList = None
    _shutil.copy2(bag[1], outputpath2)
    with _tb.open_file(outputpath2, mode = 'a') as bagfile:
        new = _np.flipud(new)
        bagfile.root.BAG_root.elevation[:,:] = new
        newu = _np.flipud(newu)
        bagfile.root.BAG_root.uncertainty[:,:] = newu
        if ioVal == 1:
            bagfile.root.BAG_root.tracking_list.remove_rows(0,bagfile.root.BAG_root.tracking_list.nrows)
        bagfile.flush()
    bagfile.close()
    gd_obj = None
    print ('done')

def sliceFinder(size, res, shape, var=5000):
    """Uses the file size of the bag to determine if the grid should be tiled.
    If the file is less than 100Mb, the file will not be tiled.  If the file is
    large enough to tile, the number of tiles and index size of each tile will
    be calculated based on the ratio of the total size of each array.

    yChunk = 5000\*sqrt(height/width)
    xChunk = 5000\*sqrt(width/height)

    ny = _np.ceil(height/yChunk)
    nx = _np.ceil(width/xChunk)

    tiles = nx\*ny

    chunkgird is both the arrangement of tiles in relation to the grid and the
    order in which the tiles are processed.

    Parameters
    ----------
    value : int
        Size of the input BAG file
    res : float
        Resolution of the input BAG data
    shape : tuple
        Dimensions of the input BAG data (y, x)
    var : int, optional
        Arbitrary value for determining chunk size

    Returns
    -------
    tiles : int
        Number of tiles used to process the data
    chunkGrid : numpy.array, None
        The arrangement of tiles in relation to the grid and the
        order in which the tiles are proceced
    sliceInfo : list, None
        The value for the tile buffer and dimensions of each tile

    Example
    -------
    >>> chunkGrid = _np.arrange(tiles).reshape((ny, nx))
    array([[ 0,  1,  2,  3,  4,  5],
           [ 6,  7,  8,  9, 10, 11],
           [12, 13, 14, 15, 16, 17],
           [18, 19, 20, 21, 22, 23],
           [24, 25, 26, 27, 28, 29],
           [30, 31, 32, 33, 34, 35]])
    # 36 Total Tiles
    # Tile 3 is index [0,3] and has a value of 2

    """
    print ('sliceFinder')
    if res < 1:
        b = 25/res
    elif res > 1:
        b = 25
    if size <= 100000:
        tiles = 0
        return tiles, None, None
    elif size > 100000:
        yChunk = int(_np.round(var*_np.sqrt(shape[0]/shape[1]))) #y
        xChunk = int(_np.round(var*_np.sqrt(shape[1]/shape[0]))) #x
        ny = int(_np.ceil(shape[0]/yChunk)) #ny
        nx = int(_np.ceil(shape[1]/xChunk)) #nx
        print (ny, nx)
        tiles = ny*nx
        chunckGrid = _np.arange(tiles).reshape((ny, nx))
        sliceInfo = [b, yChunk, xChunk]
        print (tiles, chunckGrid, sliceInfo)
        return tiles, chunckGrid, sliceInfo


def interp(grids, size, res, shape, uval, ioVal):
    """Summary of function

    Extended description of function.

    Parameters
    ----------
    grids : list
        The BAG and GeoTiff objects
    size : int
        Size of the input BAG file
    res : float
        Resolution of the input BAG data
    shape : tuple
        Dimensions of the input BAG data (y, x)
    uval : tuple
        Values for uncertainty calculation
    ioVal : int, bool
        User input. Determines whether origninal and interpolated or only
        interpolated data is output

    Returns
    -------
    ugrids : list
        Complete interpolated data as unitedBag, unitedUncr, and unitedPre

    """
    z, chunkGrid, sliceInfo = sliceFinder(size, res, shape)
    if z > 1:
        tifObjras = grids[0][-1]
        bagObjras = grids[-1][-1]
        uncObjras = grids[-1][-2]
        unitedBag = _np.empty_like(bagObjras)
        unitedUncr = _np.empty_like(bagObjras)
        unitedPre = _np.empty_like(bagObjras)
        bagShape = bagObjras.shape
        _plt.imshow(tifObjras[::100,::100])
        _plt.show()
        for ySlice in range(chunkGrid.shape[0]):
            for xSlice in range(chunkGrid.shape[1]):
                ts = _dt.now()
                index = ySlice, xSlice
                print ('\nTile', chunkGrid[index]+1, 'of', z, '-', ts)
                tile = chunk(sliceInfo, index, bagShape)
                tiffTile = tifObjras[tile.yMin:tile.yMax,tile.xMin:tile.xMax]
#                tiffTile[tiffTile < 0] = chunkGrid[index]
#                tiffTile[tiffTile > 0] = -chunkGrid[index]
                bathTile = bagObjras[tile.yMin:tile.yMax,tile.xMin:tile.xMax]
#                bathTile[bathTile < 1] = chunkGrid[index]
#                bathTile[tiffTile > 0] = -chunkGrid[index]
                uncrTile = uncObjras[tile.yMin:tile.yMax,tile.xMin:tile.xMax]
                gridSplit = [grids[0][:-1], grids[-1][:-2]]
                gridSplit[0].append(tiffTile)
                gridSplit[-1].append(uncrTile)
                gridSplit[-1].append(bathTile)
                tiffTile = None
                bathTile = None
                uncrTile = None
                combo, vals = comboGrid(gridSplit)
                print ('interp is next')
                newBag, newUncr, preBag = triangulateSurfaces(gridSplit, combo, vals, uval)
                unitedBag[tile.yBMin:tile.yBMax,tile.xBMin:tile.xBMax] = newBag[tile.yIMin:tile.yIMax,tile.xIMin:tile.xIMax]
                unitedUncr[tile.yBMin:tile.yBMax,tile.xBMin:tile.xBMax] = newUncr[tile.yIMin:tile.yIMax,tile.xIMin:tile.xIMax]
                unitedPre[tile.yBMin:tile.yBMax,tile.xBMin:tile.xBMax] = preBag[tile.yIMin:tile.yIMax,tile.xIMin:tile.xIMax]
#                unitedBag[tile.yBMin:tile.yBMax,tile.xBMin:tile.xBMax] = bathTile[tile.yIMin:tile.yIMax,tile.xIMin:tile.xIMax]
#                unitedPre[tile.yBMin:tile.yBMax,tile.xBMin:tile.xBMax] = tiffTile[tile.yIMin:tile.yIMax,tile.xIMin:tile.xIMax]
                newBag = None
                newUncr = None
                preBag = None
                td = _dt.now()
                tdelt = td - ts
                print ('Tile complete -', td, '| Tile took:', tdelt)
        tifObjras = None
        bagObjras = None
        uncObjras = None
    else:
        ts = _dt.now()
        print ('\nTile 1 of 1 -', ts)
        combo, vals = comboGrid(grids)
        print ('interp is next')
        unitedBag, unitedUncr, unitedPre = triangulateSurfaces(grids, combo, vals, uval)
        td = _dt.now()
        tdelt = td - ts
        print ('Tile complete -', td, '| Tile took:', tdelt)
#    _plt.imshow(unitedBag[::100,::100])
#    _plt.show()
#    _plt.imshow(unitedPre[::100,::100])
#    _plt.show()
    ugrids = [unitedBag, unitedUncr, unitedPre]
    return ugrids


def main(bagPath, tifPath, desPath, catzoc, ioVal):
    """main function of the interpolation process.  This function handles
    the flow of data in input, tiling, interpolation, and reintigration of the
    final grid.

    Parameters
    ----------
    bagPath : string
        File path of the input BAG file
    tifPath : list
        File paths of the input GeoTiff files
    desPath : string
        Folder path for the output files
    catzoc : string
        Selection of catzoc grade used to calculate interpolated uncertainty
    ioVal : int, bool
        User input. Determines whether origninal and interpolated or only
        interpolated data is output

    Returns
    -------
    msg : string
        String with the total time it took to complete the process

    """
    start = _dt.now()
    print ('start', start)
    tifFiles, names = getTifElev(tifPath)
    if len(tifFiles) > 1:
        tifGrids, extent = alignTifs(tifFiles)
        for tif in tifGrids:
            print(tif, '\n')
    else:
        tifGrids = tifFiles
        print (tifGrids)
        extent = tifGrids[0][2]
    comboTif, name, nameh5, targs = polyTifVals(tifGrids, desPath, names, extent)
    _plt.imshow(comboTif[::100,::100])
    _plt.show()
#    tifGrids = None
    comboArr = [0, nameh5, extent, comboTif]
    print ('\ndone with tif section\n')

    bag, maxVal, bagfile = getBagLyrs(bagPath)
    print(bag, '\n')
    size = bag[3]
    res, grids, ext, shape = alignGrids(bag, comboArr, maxVal, targs)
#    bagfile.close

    ##  Tiling Begins
    uval = catZones.get(catzoc)
    ugrids = interp(grids, size, res, shape, uval, ioVal)

    ## Tiling Ends
    print ('\nsaving is next\n')
    saveBag, saveUnc, polyList = rePrint(grids, ugrids, maxVal, ioVal)
    grids = None
    ugrids = None
    bagSave(bag, saveBag, tifGrids, res, ext, desPath, saveUnc, polyList, ioVal)
    done = _dt.now()
    print ('acually done', done)
    delta = done - start
    msg = 'Done! Took: ' + str(delta)
    print (msg)
    return msg

class chunk():
    """chunk() serves as the data container for individual tile data. It's
    inputs inlcude sliceInfo=[buffer, height, width], chunkSlice=tile[y,x]
    (tile position), and the total height and width of the BAG grid.  This
    information is then used to calculate the indices of each tile (with and
    without a buffer) and the indices of the data within the buffered tiles.

    """

    yMin, yMax, xMin, xMax = 0, 0, 0, 0
    yBMin, yBMax, xBMin, xBMax = 0, 0, 0, 0
    yIMin, yIMax, xIMin, xIMax = 0, 0, 0, 0

    def __init__(self, sliceInfo, chunkSlice, shape):
        '''Breaks the input sliceInfo into individual variables and starts the
        calculation function calcSlice()

        '''
        buffer = sliceInfo[0]
        yChunk = sliceInfo[1]
        xChunk = sliceInfo[2]
        self.calcSlice(buffer, yChunk, xChunk, chunkSlice, shape)

    def calcSlice(self, buffer, yChunk, xChunk, chunkSlice, shape):
        """Takes the complete information of a give individual tile and
        calculates the indices of each tile (with and without a buffer) and the
        indices of the data within the buffered tiles.


        The calculations are the same, regardless of axis:

        1. Determination of the data (with a buffer) for interpolation.::

            Min = 0 #(for tiles on the upper or right borders)
            #or
            Min = (chunk \* index) - buffer
            Max = chunk \* (inex + 1) + buffer
            #or
            Max = axis total

        2. Determination of the data (without a buffer) for use in re-applying data back to the original shape of the data.::

            BMin = 0 #(for tiles on the upper or right borders)
            #or
            BMin = (chunk \* index)
            BMax = chunk \* (inex + 1)
            #or
            BMax = axis total

        3. Determination of where the data (without a buffer) lies within an interpolated tile (with a buffer).::

            IMin = 0
            #or
            IMin = BMin - Min
            ma = Max - BMin
            #or
            ma = axis total
            if ma == 0:
                IMax = chunk + buffer
            if ma != 0:
                IMax = -(yma)


        Example
        -------
        >>> tiles, chunkGrid, sliceInfo = sliceFinder(20, .5, [25710,35010])
        >>> print(tiles)
        36 # 36 Total Tiles
        >>> print(chunkGrid)
        array([[ 0,  1,  2,  3,  4,  5],
               [ 6,  7,  8,  9, 10, 11],
               [12, 13, 14, 15, 16, 17],
               [18, 19, 20, 21, 22, 23],
               [24, 25, 26, 27, 28, 29],
               [30, 31, 32, 33, 34, 35]])
        # Tile 3 is index [0,3] and has a value of 2
        >>> print(sliceInfo)
        [40.0, 4285, 5835]

        For a shape of::

            [25710,35010]
            (4285, 5835)

        Tile 1 of 36::

            chunkSlice = [0,0]
            [0, 4325, 0, 5875]
            [0, 4285, 0, 5835]
            [0, -40, 0, -40]

        Tile 8 of 36::

            chunkSlice = [1,2]
            [4245, 8610, 5795, 11710]
            [4285, 8570, 5835, 11670]
            [40, -40, 40, -40]


        """
        ## 1
        self.yMin = int(max(0, (yChunk * chunkSlice[0]) - buffer))
        self.yMax = int(min(yChunk * (chunkSlice[0]+1) + buffer, shape[0]))
        self.xMin = int(max(0, (xChunk * chunkSlice[1]) - buffer))
        self.xMax = int(min(xChunk * (chunkSlice[1]+1) + buffer, shape[1]))
        slices = [self.yMin, self.yMax, self.xMin, self.xMax]

        ## 2
        self.yBMin = int(max(0, (yChunk * chunkSlice[0])))
        self.yBMax = int(min(yChunk * (chunkSlice[0]+1), shape[0]))
        self.xBMin = int(max(0, (xChunk * chunkSlice[1])))
        self.xBMax = int(min(xChunk * (chunkSlice[1]+1), shape[1]))
        tiles = [self.yBMin, self.yBMax, self.xBMin, self.xBMax]

        ## 3
        self.yIMin = int(max(0,(self.yBMin-self.yMin)))
        yma = int(min((self.yMax-self.yBMax), shape[0]))
        if yma == 0:
            self.yIMax = int(yChunk + buffer)
        else:
            self.yIMax = -int(yma)
        self.xIMin = int(max(0,(self.xBMin-self.xMin)))
        xma = int(min((self.xMax-self.xBMax), shape[1]))
        if xma == 0:
            self.xIMax = int(xChunk + buffer)
        else:
            self.xIMax = -int(xma)
        borders = [self.yIMin, self.yIMax, self.xIMin, self.xIMax]

        print (slices, tiles, borders, sep='\n')
