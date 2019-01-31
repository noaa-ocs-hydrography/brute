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
#from __future__ import absolute_import, division, print_function

import os
import re
import wx
import ast
import sys
import cv2
import scipy
import shutil
import datetime
import numpy as np
import tables as tb

from lxml import etree
from string import ascii_lowercase
from scipy.ndimage.interpolation import zoom
from scipy.spatial import cKDTree
from osgeo import gdal, ogr, osr, gdalconst
import matplotlib.pyplot as plt

import autointerp_ui

# this allows GDAL to throw Python Exceptions
gdal.UseExceptions()

progLoc = os.getcwd()
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

combo = re.compile(r'COMBINED', re.IGNORECASE)
interp = re.compile(r'INTERP', re.IGNORECASE)

def getTifElev(files):
    '''This function takes a list of GeoTiff filepaths and opens each in order
    to access their values. Returned as a list of tif objects are the order a
    file was found, the path of the file, and the file's values as
    [x, file, extent, arr] and a list of file names

    Upper Left Corner and Lower Right Corner Bounds Directly from:
    "How to get raster corner..." on on GIS Stack Exchange [Answer by 'James'
    (https://gis.stackexchange.com/a/201320)]
    '''
    print ('getTifElev')
    tifFiles = []
    names = []
    maxVal = 0
    y = 0
    for file in files:
        print(file)
        ds = gdal.Open(file)
        
        ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
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

        fName = os.path.split(file)[-1]
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
        print (np.amax(arr), np.amin(arr))
        print (meta)
        band=None
        ds=None
        tifFiles.append ([y, file, meta, arr])
        y += 1
    return tifFiles, names


def getBagLyrs(fileObj):
    '''This function takes a list of BAG filepaths and opens each in order
    to access their values. Returned is a list of bag objects are the order a
    file was found, the path of the file, the file's elevation values, and the
    file's uncertainty values as [x, file, extent, elev] and a list of file 
    names
    '''
    print ('getBagLyrs')
#    bagFiles = []
    names = []
    y = 0
    print (fileObj)
    fName = os.path.split(fileObj)[-1]
    splits = fName.split('_')
    name = '_'.join([x for x in splits[:2]])
    names.append(name)
    print (splits)
    with tb.open_file(fileObj, mode = 'r') as bagfile:
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
        xml_tree = etree.XML(meta_xml)
        try:
            ret = xml_tree.xpath('//*/gmd:spatialRepresentationInfo/gmd:MD_Georectified/'
                                      'gmd:cornerPoints/gml:Point/gml:coordinates',
                                      namespaces=ns)[0].text.split()
        except (etree.Error, IndexError) as e:
            try:
                ret = xml_tree.xpath('//*/spatialRepresentationInfo/smXML:MD_Georectified/'
                                          'cornerPoints/gml:Point/gml:coordinates',
                                          namespaces=ns2)[0].text.split()
            except (etree.Error, IndexError) as e:
                print("unable to read corners SW and NE: %s" % e)
                return

        try:
            sw = [float(c) for c in ret[0].split(',')]
            ne = [float(c) for c in ret[1].split(',')]
        except (ValueError, IndexError) as e:
            print("unable to read corners SW and NE: %s" % e)
            return
        print (ret, sw, ne)
        sx,sy = sw
        nx,ny = ne
        meta = [[sx,ny], [nx,sy]]
        elev = np.flipud(bagfile.root.BAG_root.elevation.read())
        uncr = np.flipud(bagfile.root.BAG_root.uncertainty.read())
        print (np.amax(elev), np.amin(elev))
        print (np.amax(uncr), np.amin(uncr))
        bag = [y, fileObj, meta, uncr, elev]
        print (bag)
        bagfile.close()
    return bag

def getShpRast(file):
    print ('getShpRast')
    file = ogr.Open(file)
    shape = file.GetLayer()
    fnum = shape.GetFeatureCount()
    x = np.arrange(fnum)
    polys = []
    for n, x in enumerate(np.arange(fnum)):
        polys.append(shape.GetFeature(n))
    print (polys)

def write_raster(raster_array, gt, data_obj, outputpath, dtype=gdal.GDT_UInt32,
                 options=0, color_table=0, nbands=1, nodata=False):
    '''Directly From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
    (https://gis.stackexchange.com/a/278965)]
    '''
    print('write_raster')

    height, width = raster_array.shape

    # Prepare destination file
    driver = gdal.GetDriverByName("GTiff")
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
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dest.SetProjection(srs.ExportToWkt())

    # Close output raster dataset
    dest = None

def maxValue(arr):
    '''Takes an input array and finds the most used value in the array, this 
    value is used by the program to assume the array's nodata value
    
    returns the most used value in the array as an integer
    '''
    print ('maxValue')
    nums, counts = np.unique(arr, return_counts =True)
    index = np.where(counts==np.amax(counts))
    print (index, nums[index])
    return int(nums[index])


def tupleGrid(grid, maxVal, add='no'):
    '''Takes an input matrix and an assumed nodata value. The function iterates 
    through the matrix and compiles a list of 'edge' points [[x, y, z], ...] where 
        1) the current value is not a nodata value and previous value was a 
        nodata value
            - sets bool(io) True, indicating that the next value to compare 
            against should be a nodata value
        2) the current value is a nodata value and the previous value was not 
        a nodata value
            - sets bool(io) False, indicating that the next value to compare
            against should not be a nodata value
    returns a list of the points found
    
        TODO: figure out a supported way to differentiate from a add input of
        None or input array.  Numpy gives a warning against using typewise
        or elementwise comparasons (ex. if add == None [type]
                                     or if add == 'no' [element])
    '''
    print ('tupleGrid')
    print (add, type(add))
    if add == 'no':
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
        return points
    else:
        points = []
        points2 = []
        a = 0
        for x in range(grid.shape[1]):
            io = False
            for y in range(grid.shape[0]):
                if grid[y,x] == maxVal:
                    if grid[y-1,x] != maxVal:
                        val = grid[y-1,x]
                        val2 = add[y-1,x]
                        point = [x, y-1, val]
                        point2 = [x, y-1, val2]
                        if a == 1:
                            print (point, val)
                            print (point2, val2)
                            a += 1
                        points.append(point)
                        points2.append(point2)
                        io = False
                    else:
                        pass
                else:
                    if io == False:
                        val = grid[y,x]
                        val2 = add[y,x]
                        point = [x, y, val]
                        point2 = [x, y, val2]
                        if a == 0:
                            print (point, val)
                            print (point2, val2)
                            a += 1
                        points.append(point)
                        points2.append(point2)
                        io = True
        return points, points2

def tupleGrid2(grid, maxVal):
    ''' Goal of tupleGrid2:
            - Determine better or more efficient way to handle 'edge' detection
            for the data.
        Why this is hard:
            - Most out of the box edge detectors either interpolate the edges
            for a better look or arrive at the same result due to a more general
            algorithm.
                - This causes the lines of the edges to cross into areas of no-
                data accross the grids
                -This is not desireable because the goal of the algorithm is to
                locate all data that sits on the edge of no-data boundries as a
                form of data reduction during triangulaiton/interpolation.
                    - To populate no-data values into these functions is not 
                    only undesireable but also detrimental to the expected 
                    result
        Current solution:
            - Continue using original 'edge' detection method and limit the size
            of input files in order to prevent program memory overload
    '''
    print ('tupleGrid2')
#    grid = np.nan_to_num(grid)

    return [[],[]]

def concatGrid(grids, maxVal, shape):
    '''Takes an input of an array of grid objects and the assumed nodata value
    Passes the assumed nodata value and the arrays held within each of the 
    listed grid objects to tupleGrid(grid, maxVal) for a return of an array of
    edge points for each grid [[x, y, z], ...]
    
    Takes the results of both tupleGrid calls and combines them into a single 
    array of edge points.  Then uses .view().sort() to sort the combined 
    products [[x, y, z], ...] based on ascending x, then column y.  
    
    The results of the sorted point array are then split into a list of points 
    [[x, y], ...] and values [[z], ...] and returned
    '''
    print ('concatGrid')
    print ('tpts', datetime.datetime.now())
#    print (grids[0])
    tpts = tupleGrid(grids[0], maxVal)
    print ('bpts', datetime.datetime.now())
    bpts = tupleGrid(grids[1], maxVal)#, add=grids[2])
    print ('done', datetime.datetime.now())
    print ('combo1', datetime.datetime.now())
    comb = np.concatenate([tpts, bpts])
    comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
    print (comb)
    grid = np.hsplit(comb, [2, 4])
    vals = grid[1].squeeze()
    grid = grid[0]
    print (grid, vals)
#    print ('combo2', datetime.datetime.now())
#    comb2 = np.concatenate([tpts, upts])
#    comb2.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
#    print (comb2)
#    grid2 = np.hsplit(comb2, [2, 4])
#    vals2 = grid2[1].squeeze()
#    grid2 = grid2[0]
#    print (grid2, vals2)
    return grid, vals#, vals2

def alignTifs(tifs):
    '''Takes an input of an array of tiff objects. The goal of this function 
    is to fit the provided tifs to the largest combined area of the inputs.
    
    Steps:
        1) Find NW and SE corners of the first input array and the number of
        rows and columns in the input
        2) Find NW and SE corners of subsiquent input arrays and record 
        differences as the NW differences (nxd, nyd) SE differences (sxd, syd)
            - also modify the original NW and/or SE corner values if applicable
            - also modify the original number or rows and columns if applicable
        3) Resize the area of each input array to the new total combined size
        of all input arrays based on the stored number of rows and columns
        4) Reposition the data within each input array to maintain original 
        geographic placement based on the difference of the input 
        array's original NW corner and the recorded 'largest' NW boundry of the
        combined input areas
        5) return tiff objects with the newly resized/repositioned tifs and the 
        new NW and SE extents of the new combined area
    '''
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
            rows, cols = tif[-1].shape
            print (nw, se)
            x += 1
        else:
            ulx, uly = ul
            lrx, lry = lr
            nwx, nwy = nw
            scx, scy = se
            if ulx <= nwx:
                nxd = ulx - nwx
                nwx = ulx
            elif ulx >= nwx:
                nxd = ulx - nwx
                nwx = ulx
            if uly >= nwy:
                nyd = uly - nwy
                nwy = uly
            elif uly <= nwy:
                nyd = uly - nwy
                nwy = uly
            if lrx >= scx:
                sxd = lrx - scx
                scx = lrx
            elif lrx <= scx:
                sxd = scx - lrx
                scx = lrx
            if lry <= scy:
                syd = lry - scy
                scy = lry
            elif lry >= scy:
                syd = lry - scy
                scy = lry
            r, c = tif[-1].shape[0], tif[-1].shape[1]
            print (rows, cols)
            print (r, c)
            if rows <= r:
                rows = r
            if cols <= c:
                cols = c
            print (nxd, nyd, sxd, syd)
            print (sxd, nyd)
            nw = [nwx, nwy]
            se = [scx, scy]
            print('cols: ' , nwy - scy, '\nrows: ', scx - nwx)
    print (nw, se)
    sizedTifs = []
    print ('resize?')
    if sxd != 0 or nyd != 0:
        print ('yes', datetime.datetime.now())
        for tif in tifs:
            maxVal = maxValue(tif[-1])
            sizedTif = tif[:-1]
            print (tif)
            bndx, bndy = tif[2][0]
            nwx, nwy = nw
            print (bndx, bndy)
            print (nwx, nwy)
            arr = tif[-1]
            bef = arr.shape
            x, y = arr.shape
            print (x, rows)
            print (y, cols)
            if y != cols:
                exp = cols - y
                print (exp)
#                for a in range(exp):
#                    arr = np.column_stack([arr, [maxVal] * x])
                add = np.full((x, exp), maxVal)
                arr = np.column_stack([arr, add])
                x, y = arr.shape
            if x != rows:
                exp = rows - x
                print (exp)
#                for a in range(exp):
#                    arr = np.vstack([arr, [maxVal] * y])
                add = np.full((exp, y), maxVal)
                arr = np.vstack([arr, add])
            print (bef, arr.shape)
            if nwx != bndx:
                rollx = bndx - nwx
                print (rollx)
                arr = np.roll(arr, int(rollx), axis=1)
            if nwy != bndy:
                rolly = nwy - bndy
                print (rolly)
                arr = np.roll(arr, int(rolly), axis=0)
            sizedTif.append(arr)
            sizedTifs.append(sizedTif)
#            sizedTifs.append(arr)
        ext = [nw, se]
        print ('done', datetime.datetime.now())
        return sizedTifs, ext
    else:
        print ('same')
        ext = [nw, se]
        return tifs, ext

def polyTifVals(tifs, path, names):
    '''Heavy Influence From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
    (https://gis.stackexchange.com/a/278965)]

    This function takes an array input of tif objects, a destination path (for 
    ouput saving), and a list the names of the input tif objects.
    
    Takes the arrays of each tif object and combines them along the z axis of 
    each:
        [[za,zb],[za,zb],[za,zb],
         [za,zb],[za,zb],[za,zb],
         [za,zb],[za,zb],[za,zb]]
    Then takes the mean of the values along the z axis. The reslult is then 
    modified to a binary raster by making all values that are not the nodata
    value 1 and the values that are nodata 0.
    
    This binary raster is saved at the input path and also returned with the 
    new full path for the output
    '''
    print ('polyTifVals')
    letters = []
    for letter in ascii_lowercase:
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
            print (maxRes)
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
        gd_obj = gdal.Open(tifs[i][1])
        array = np.expand_dims(tif,2)
        if i == 0:
            allarrays = array
        else:
            allarrays = np.concatenate((allarrays, array), axis=2)
#    meanTiff = (allarrays < maxVal).astype(np.int)
    meanTiff = np.nanmean(allarrays, axis=2)
#    meanTiff = np.nanmean(meanTiff, axis=2)
    if maxVal != 0:
        meanTiff = (meanTiff < maxVal).astype(np.int)
    else:
        meanTiff = (meanTiff > maxVal).astype(np.int)

    outputpath = path + '/' + names[0] + '_COMBINEDPOLY.tif'

    print(outputpath)
    while True:
        if os.path.exists(outputpath):
            os.remove(outputpath)
        elif not os.path.exists(outputpath):
            break
    write_raster(meanTiff, gd_obj.GetGeoTransform(), gd_obj, outputpath)
    gd_obj = None

#    pointTiff = tupleGrid(meanTiff, 0)

    print ('done?')
    return meanTiff, outputpath

def alignGrids(bag, tif):
    print ('alignGrids')
    print (bag[-1])
    maxVal = maxValue(bag[-1])
    maxTif = maxValue(tif[-1])
    tu, tl = tif[-2]
    tx, ty = tif[-1].shape
    tex, tey = tu[1] - tl[1], tl[0] - tu[0]
    print (tx, ty)
    print (tex, tey)
    if tex-tx == 0 and tey-ty == 0:
        tifRes = 1
        print (tifRes)
    else:
        tifRes = np.round(np.mean([tex/tx, tey/ty]))
        print(tifRes)
#    tifRes = 1
    splits = bag[1].split('/')[-1]
    splits = splits.split('_')[2]
    print (splits)
    bagRes = re.split('(\d+)',splits)[-2:]
    print (bagRes)
    if bagRes[-1] == 'cm':
        bagRes = int(bagRes[0]) / 100
        zres = tifRes/bagRes
    elif bagRes[-1] =='m':
        bagRes = int(bagRes[0])
        zres = tifRes/bagRes
    print (bagRes, zres)
    tmax = 0
    bmax = maxValue(bag[-1])

    print (tif[-1])
    if zres == 1:
        newarr = tif[-1]
    else:
        print ('zoom', datetime.datetime.now())
        newarr = zoom(tif[-1], zoom=[zres, zres], order=3, prefilter=False)        
        print ('zoomed', datetime.datetime.now())

    newarr = newarr.astype('float64')
    newarr[newarr > 0] = np.nan
    newarr[newarr < 1] = float(maxVal)

    print (newarr)

#    maxTif = maxValue(newarr)

    print (tif[-1].shape, newarr.shape)

    tif.pop()
    print (tif)
    tif.append(newarr)
    print (tif, tif[-1].shape)

    bSx, bSy = bag[-1].shape
    tSx, tSy = tif[-1].shape
    bagBounds = bag[2]
    tifBounds = tif[2]
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
    newarr2 = tif[-1]
    print (newarr2.shape)
    expx, expy = 0, 0
    if newarr2.shape != bag[-1].shape:
        print (bSx - tSx, bSy - tSy)
        if tSx < bSx:
            expx = np.abs(bSx - tSx)
            print (expx)
            add = np.full((expx, tSy), maxVal)
            newarr2 = np.vstack([newarr2, add])
            tSx, tSy = newarr2.shape
        if tSy < bSy:
            expy = np.abs(bSy - tSy)
            print (expy)
            add = np.full((tSx, expy), maxVal)
            newarr2 = np.column_stack([newarr2, add])
    rollx = dulx*zres
    rolly = duly*zres
    print (dulx, rollx, duly, rolly)
    if dulx != 0:
        print('x')
        newarr2 = np.roll(newarr2, int(rollx), axis=1)
    if dlry != 0:
        print('y')
        newarr2 = np.roll(newarr2, int(rolly), axis=0)
    newarr2 = newarr2[0:bSx,0:bSy]
    print (newarr.shape, newarr2.shape, (bSx, bSy))
    tif.pop()
    print (tif)
    tif.append(newarr2)
    print (tif, tif[-1].shape)

    grids = [tif, bag]

    return bagRes, grids, bagBounds

def comboGrid(grids):
    print ('comboGrid')
    maxVal = maxValue(grids[1][-1])
    shape = grids[1][-1].shape
    tif = grids[0][-1]
    bag = grids[1][-1]
#    uncr = grids[1][-2]
    arrs = [tif, bag]#, uncr]
    combo, vals = concatGrid(arrs, maxVal, shape)
    return combo, vals

#def comboGrid2(grids):
#    print ('comboGrid2')
#    maxVal = maxValue(grids[1][-1])
#    shape = grids[1][-1].shape
#    arrs = []
#    for grid in grids:
#        arrs.append(grid[-1])
#    combo, vals = concatGrid2(arrs, maxVal, shape)
#    return combo, vals

def rePrint(bag, interp, poly, maxVal, uncr, uval, ioVal):
    print ('rePrint', datetime.datetime.now())
    m, b = uval
    perVal = maxVal*m
    print (maxVal, perVal)
    rows, cols = bag.shape
    tpoly = np.nan_to_num(poly)
    tpoly = (tpoly < maxVal).astype(np.int)
    bpoly = (bag < maxVal).astype(np.int)
    cpoly = np.logical_or(bpoly, tpoly)
    dpoly = np.logical_xor(bpoly, cpoly)
    if ioVal == False:
        nbag = np.where(dpoly, interp, bag)
        nunc = np.where(dpoly, (interp*m)+b, uncr)
    elif ioVal == True:
        nbag = np.where(dpoly, interp, maxVal)
        nunc = np.where(dpoly, (interp*m)+b, maxVal)
    nunc[nunc>perVal] = maxVal
    print ('done', datetime.datetime.now())
    return nbag, nunc, dpoly

def triangulateSurfaces(grids, combo, vals, uval, ioVal):
    print ('triangulateSurfaces')
    bagObj = grids[-1]
    bag = bagObj[-1]
    uncr = bagObj[-2]
    tifObj = grids[0]
    poly = tifObj[-1]
    maxVal = maxValue(bag)
    x, y = np.arange(bag.shape[1]), np.arange(bag.shape[0])
    xi, yi = np.meshgrid(x, y)
    print ('try tri', datetime.datetime.now())
    values = scipy.interpolate.griddata(combo, vals, (xi, yi), 
                                        method='linear', fill_value=maxVal)
    print ('done', datetime.datetime.now())
#    print ('try uncr', datetime.datetime.now())
#    values2 = scipy.interpolate.griddata(combo, uval, (xi, yi), 
#                                        method='linear', fill_value=maxVal)
#    print ('done', datetime.datetime.now())
#    maxVal = maxValue(values)
#    print (values)
    values = np.asarray(values, dtype='float64')
    values[np.isnan(values)]=maxVal
    print (values.shape, poly.shape)
#    values2 = np.asarray(values2, dtype='float64')
#    values2[np.isnan(values2)]=maxVal
    grid, uncr, dpoly = rePrint(bag,values,poly,maxVal,uncr,uval, ioVal)
    return grid, uncr, dpoly

def bagSave(bag, new, tifs, res, ext, path, newu, dpoly):
    for tif in tifs:
        gd_obj = gdal.Open(tif[1])
        break
    print ('bagSave')
    oextX, oextY = bag[2]
    nx, ny = ext[0]
    sx, sy = ext[-1]
    reso = float(res)
    print (nx, ny, sx, sy)
    gtran = (nx, reso, 0.0, ny, 0.0, -(reso))
    print (gtran)
    fName = bag[1].split('\\')[-1]
    print (fName)
    split = fName.split('_')[:2]
    if res < 1:
        res = '50cm'
    else:
        res = str(res) + 'm'
    bagName = '_'.join([x for x in split]) + '_' + res + '_INTERP'
    outputpath = path + '/' + bagName +'.tif'
    outputpath2 = path + '/' + bagName +'.bag'
    print(outputpath)
    while True:
        if os.path.exists(outputpath):
            os.remove(outputpath)
        elif not os.path.exists(outputpath):
            break
        if os.path.exists(outputpath2):
            os.remove(outputpath2)
        elif not os.path.exists(outputpath2):
            break
    write_raster(dpoly, gtran, gd_obj, outputpath, 
                 dtype=gdal.GDT_Float64, nodata=0)
    shutil.copy2(bag[1], outputpath2)
    with tb.open_file(outputpath2, mode = 'a') as bagfile:
        new = np.flipud(new)
        bagfile.root.BAG_root.elevation[:,:] = new
        newu = np.flipud(newu) 
        bagfile.root.BAG_root.uncertainty[:,:] = newu
        bagfile.flush()
    bagfile.close()
    gd_obj = None
    print ('done')

def interp(bagPath, tifPath, desPath, catzoc, ioVal):
    start = datetime.datetime.now()
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
    comboTif, name = polyTifVals(tifGrids, desPath, names)
    comboArr = [0, name, extent, comboTif]

    bag = getBagLyrs(bagPath)
    print(bag, '\n')
    res, grids, ext = alignGrids(bag, comboArr)
    combo, vals = comboGrid(grids)
    tifObj = grids[0]
    poly = tifObj[-1]
    print (combo.shape, poly.shape)
    uval = catZones.get(catzoc)
#    print (uval)
    newBag, newUncr, dpoly = triangulateSurfaces(grids, combo, vals, uval, ioVal)
    bagSave(bag, newBag, tifGrids, res, ext, desPath, newUncr, dpoly)
    done = datetime.datetime.now()
    delta = done - start
    print ('done', done, '\ntook', delta)
    return True

class Form(autointerp_ui.Form):
    '''Load ui and: define tif storage columns, overwrite ui defined fucntions 
    with desired function behaviour
    '''
    def __init__(self, parent):
        autointerp_ui.Form.__init__(self, parent)
        self.insInd = 0
        #Instantiation of 'GeoTIFF File List' box Columns:
        self.list_tif.InsertColumn(0, 'File', width=200)
        self.list_tif.InsertColumn(1, 'Path', width=500)
        
    def programQuit(self, event):
        '''Closes GUI, ends program.
        Maps to Cancel button, File->Quit, and CTRL+Q
        '''
        self.Close()
        
    def itemInsert(self, event):
        '''Adds files selected from the 'Add GeoTIFF File' to the 'GeoTIFF File
        List' box. 'File' holds the name of the file and 'Path' holds the 
        complete file path
        '''
        print (self.picker_tif.GetPath())
        tif = self.picker_tif.GetPath()
        self.gettifList()
        if tif not in self.tifList:  
            name = os.path.split(self.picker_tif.GetPath())
            self.list_tif.InsertItem(self.insInd, name[1])
            self.list_tif.SetItem(self.insInd, 1, tif)
            self.insInd += 1
        
    def itemRemove(self, event):
        '''Removes selected files from the 'GeoTIFF File List' box when the 
        'Remove' button is clicked'''
        selected = self.list_tif.SelectedItemCount
        for x in range(0, selected):
            sel = self.list_tif.GetFirstSelected()
            self.list_tif.DeleteItem(sel)
            if self.insInd > 0:
                self.insInd -= 1
        self.gettifList()
            
    def programProg(self, event):
        '''Collects the GUI field values for use in running the tools main
        function 'interp(bagPath, tifPath, desPath)'
        '''
        bagPath = self.picker_bag.GetPath()
        self.gettifList()
        tifPath = self.tifList
        desPath = self.picker_des.GetPath()
        if desPath == '':
            desPath = os.path.split(bagPath)[0]
        catzoc = self.choice_catzoc.GetString(self.choice_catzoc.GetCurrentSelection())
        ioOut = self.radio_data.GetSelection()
        interp(bagPath, tifPath, desPath, catzoc, ioOut)
            
    def gettifList(self):
        tifCount = self.list_tif.GetItemCount()
        self.tifList = []
        for x in range(0, tifCount):
            self.tifList.append(self.list_tif.GetItemText(x, col=1))
        print (self.tifList)
            
class Done(autointerp_ui.Done):
    def __init__(self, parent):
        autointerp_ui.Done.__init__(self, parent)
        
        
app = wx.App()
frame = Form(None)
icon = wx.Icon()
icon.CopyFromBitmap(wx.Bitmap("autointerp.ico", wx.BITMAP_TYPE_ANY))
frame.SetIcon(icon)
frame.Show()
app.MainLoop()