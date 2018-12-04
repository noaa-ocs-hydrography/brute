# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 16:52:49 2018

@author: Casiano.Koprowski

Made in part with code from:
    - Extract Bag Outlines on Pydro by Dayong.Shen
    - QC Tools 2: Flier Finder on Pydro by HydrOffice
    - QC Tools 2: Bag Metadata Editor on Pydro by HydrOffice
    - "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon' 
        (https://gis.stackexchange.com/a/278965)]
    - "How to get raster corner..." on on GIS Stack Exchange [Answer by 'James' 
        (https://gis.stackexchange.com/a/201320)]
"""
#from __future__ import absolute_import, division, print_function

import os
import re
import sys
import cv2
import scipy
import shutil
import datetime
import numpy as np
import tables as tb
import matplotlib.pyplot as plt

from string import ascii_lowercase
from scipy.ndimage.interpolation import zoom
from osgeo import gdal, ogr, osr, gdalconst

from hyo.bag import BAGFile
from hyo.bag import BAGError
from hyo.bag.helper import Helper
from hyo.bag.meta import Meta

import logging

logger = logging.getLogger()
logger.setLevel(logging.NOTSET)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)  # change to WARNING to reduce verbosity, DEBUG for high verbosity
ch_formatter = logging.Formatter('%(levelname)-9s %(name)s.%(funcName)s:%(lineno)d > %(message)s')
ch.setFormatter(ch_formatter)
logger.addHandler(ch)

# this allows GDAL to throw Python Exceptions
gdal.UseExceptions()

os.chdir('R:/Scripts/svn-nbs/scripts/baginterp')
from autoInterp_ui import Ui_MainWindow
#os.chdir('R:/Scripts/svn-nbs/scripts/baginterp')
#qtCreatorFile = 'autoInterp.ui'
#Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

combo = re.compile(r'COMBINED', re.IGNORECASE)
interp = re.compile(r'INTERP', re.IGNORECASE)

def tifInput(path):
    '''This function takes a file/folder path and identifies and returns a list
    of all the paths of GeoTiff files found.
    '''
    print ('tifInput')
    inFiles = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.tif') or file.lower().endswith('.tiff'):
                print (interp.match(file))
                print (interp.search(file))
                if interp.search(file) != None or combo.search(file) != None:
                    pass
                else:
                    inFiles.append(root + '/' + file)
    return inFiles

def bagInput(path):
    '''This function takes a file/folder path and identifies and returns a list
    of all the paths of BAG files found.
    '''
    print ('bagInput')
    inFiles = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.bag'):
                if interp.search(file) != None or combo.search(file) != None:
                    pass
                else:
                    inFiles.append(root + '/' + file)
    return inFiles

def getTifElev(files):   
    '''This function takes a list of GeoTiff filepaths and opens each in order
    to access their values. Returned as a list of tif objects are the order a
    file was found, the path of the file, and the file's values as 
    [x, file, arr] and a list of file names
    
    Upper Left Corner and Lower Right Corner Bounds Directly from:
    "How to get raster corner..." on on GIS Stack Exchange [Answer by 'James' 
    (https://gis.stackexchange.com/a/201320)]
    '''
    print ('getTifElev')
    tifFiles = []
    names = []
    maxVal = 0
    x = 0
    for file in files:
        print(file)
        ds = gdal.Open(file)
#        print (gdal.Info(ds))
        
        ulx, xres, xskew, uly, yskew, yres  = ds.GetGeoTransform()
        lrx = ulx + (ds.RasterXSize * xres)
        lry = uly + (ds.RasterYSize * yres)
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
        tifFiles.append ([x, file, meta, arr])
        x += 1
    return tifFiles, names


def getBagLyrs(files):
    '''This function takes a list of BAG filepaths and opens each in order
    to access their values. Returned is a list of bag objects are the order a
    file was found, the path of the file, the file's elevation values, and the 
    file's uncertainty values as [x, file, elev, uncr] and a list of file names
    '''
    print ('getBagLyrs')
    bagFiles = []
    names = []
    x = 0
    for file in files:
        print (file)
        bag_0 = BAGFile(file)   
#        print (Meta(bag_0.metadata()))
        meta = Meta(bag_0.metadata()).sw, Meta(bag_0.metadata()).ne
        print (meta)
        bag_0.close()
        fName = os.path.split(file)[-1]
        splits = fName.split('_')
        name = '_'.join([x for x in splits[:2]])
        names.append(name)
        print (splits)
        with tb.open_file(file, mode = 'r') as bagfile:
            elev = np.flipud(bagfile.root.BAG_root.elevation.read())
#            uncr = np.flipud(bagfile.root.BAG_root.uncertainty.read())
            print (np.amax(elev), np.amin(elev))
#            print (np.amax(uncr), np.amin(uncr))

            bagFiles.append([x, file, meta, elev])
            print(elev.shape)
            bagfile.close()
        x += 1
    return bagFiles, names

def write_raster(raster_array, gt, data_obj, outputpath, dtype=gdal.GDT_UInt32,
                 options=0, color_table=0, nbands=1, nodata=False):
    '''
    Directly From:
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
    print ('maxValue')
    nums, counts = np.unique(arr, return_counts =True)
    index = np.where(counts==np.amax(counts))
    print (index, nums[index])
    return int(nums[index])
    
    
def tupleGrid(grid, maxVal):
    print ('tupleGrid')
    points = []
#    zvals = []
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
#                    zvals.append(val)
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
#                    zvals.append(val)
                    points.append(point)
                    io = True
    return points#, zvals

def tupleGrid2(grid, maxVal):
    print ('tupleGrid2')
    edge_horizont = scipy.ndimage.sobel(grid, 0)
    edge_vertical = scipy.ndimage.sobel(grid, 1)
    magnitude = np.hypot(edge_horizont, edge_vertical)
#    print (magnitude, maxValue(magnitude))
#    magnitude = np.where(0<magnitude<.5)
#    print (magnitude, maxValue(magnitude))
#    points = np.take(grid, magnitude)
    
#    return [magnitude, points]
    return [[],[]]

def concatGrid(grids, maxVal, shape):
    print ('concatGrid')
    print ('tpts', datetime.datetime.now())
#    print (grids[0])
    tpts = tupleGrid(grids[0], maxVal)
    print ('bpts', datetime.datetime.now())
    bpts = tupleGrid(grids[1], maxVal)
    print ('done', datetime.datetime.now())
    comb = np.concatenate([tpts, bpts])
    comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0) #<-- returns None
    print (comb)
    grid = np.hsplit(comb, [2, 4])
    vals = grid[1].squeeze()
    grid = grid[0]
    print (grid, vals)
    return grid, vals

def concatGrid2(grids, maxVal, shape):
    print ('concatGrid2')
    print ('tpts', datetime.datetime.now())
    tpts = tupleGrid2(grids[0], maxVal)
    tcom = np.concatenate([tpts])
    print (tcom.shape, tcom)
    print ('bpts', datetime.datetime.now())
    bpts = tupleGrid2(grids[1], maxVal)
    bcom = np.concatenate([bpts])
    print (bcom.shape, bcom)
    print ('done', datetime.datetime.now())
    comb = np.concatenate([tcom, bcom])
    comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
    grid = np.hsplit(comb, [2, 4])
    grid = grid[0]
    vals = grid[1].squeeze()
    print (grid, vals)
    return grid, vals
         
def alignTifs(tifs):
    print ('alignTifs')
    x = 0
#    maxVal = 0
    nw, se = [0,0], [0,0]
    rows, cols = 0, 0
    sxd, nyd = 0, 0
    nxd, syd = 0, 0
    for tif in tifs:
        maxVal = maxValue(tif[-1])
        print (maxVal)
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
            if arr.shape != bef:
                arr = np.roll(arr, int(sxd), axis=0)
                arr = np.roll(arr, int(nyd), axis=1)
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
    '''
    Heavy Influence From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon' 
    (https://gis.stackexchange.com/a/278965)]
    
    This function takes an input of tif
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
        tifRes = np.mean([tex/tx, tey/ty])
        print(tifRes)
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
    
    print ('zoom', datetime.datetime.now())
    newarr = zoom(tif[-1], zoom=[zres, zres], order=3, prefilter=False)
    newarr = newarr.astype('float64')
    
    newarr[newarr > 0] = np.nan
    newarr[newarr < 1] = maxVal
    
    print ('zoomed', datetime.datetime.now())
    maxTif = maxValue(newarr)

    print (tif[-1].shape, newarr.shape)

    tif.pop()
    tif.append(newarr)
    
    nw, se = [0,0], [0,0]
    rows, cols = 0, 0
    sxd, nyd = 0, 0
    nxd, syd = 0, 0
    bSx, bSy = bag[-1].shape
    bagBounds = bag[2]
    tifBounds = tif[2]
    bulx, buly = bagBounds[0]
    blrx, blry = bagBounds[-1]
    tulx, tuly = tifBounds[0]
    tlrx, tlry = tifBounds[-1]
    dulx, duly, dlrx, dlry = 0, 0, 0, 0
    if bulx != tulx:
        dulx = bulx - tulx
    if buly != tuly:
        print (buly, tuly)
        duly = buly - tuly
    if blrx != tlrx:
        dlrx = blrx - tlrx
    if blry != tlry:
        dlry = blry - tlry
    newarr2 = tif[-1]
    tSx, tSy = tif[-1].shape
    if newarr2.shape != bag[-1].shape:
        if tSx < bSx:
            exp = bSx - tSx
            add = np.full((np.abs(exp), np.abs(tSy)), maxVal)
            newarr2 = np.vstack([newarr2, add])
            tSx, tSy = newarr2.shape
        if tSy < bSy:
            exp = bSy - tSy
            add = np.full((np.abs(tSx), np.abs(exp)), maxVal)
            newarr2 = np.column_stack([newarr2, add])
    print (dulx, duly)
    if dulx != 0:
        tif[-1] = np.roll(newarr2, int(dulx), axis=0)
    if duly != 0:
        tif[-1] = np.roll(newarr2, int(duly), axis=1)
    newarr2 = newarr2[0:bSx,0:bSy]
    print (newarr.shape, newarr2.shape, (bSx, bSy))
    tif.pop()
    tif.append(newarr2)
    
    grids = [tif, bag]
    
    return bagRes, grids, bagBounds

def comboGrid(grids):
    print ('comboGrid')
    maxVal = maxValue(grids[1][-1])
    shape = grids[1][-1].shape
    arrs = []
    for grid in grids:
        arrs.append(grid[-1])
    combo, vals = concatGrid(arrs, maxVal, shape)
    return combo, vals

def comboGrid2(grids):
    print ('comboGrid2')
    maxVal = maxValue(grids[1][-1])
    shape = grids[1][-1].shape
    arrs = []
    for grid in grids:
        arrs.append(grid[-1])
    combo, vals = concatGrid2(arrs, maxVal, shape)
    return combo, vals
    

def rePrint(bag, interp, poly, maxVal):
    print ('rePrint', datetime.datetime.now())
    rows, cols = bag.shape
    arr = np.nan * np.empty((rows,cols))
    z = 0
    print (maxValue(poly))
    for a in range(bag.shape[1]):
        for b in range(bag.shape[0]):
            g = bag[b,a]
            i = interp[b,a]
            p = poly[b,a]
            if g != maxVal:
                bag[b,a] = g
            elif p != maxVal:
                if g == maxVal or i != np.nan:
                    bag[b,a] = i
    print ('done', datetime.datetime.now())
    return bag

def rePrint2(bag, interp, poly, maxVal):
    print ('rePrint', datetime.datetime.now())
    rows, cols = bag.shape
    tpoly = np.nan_to_num(poly)
    tpoly = (tpoly < maxVal).astype(np.int)
    bpoly = (bag < maxVal).astype(np.int)
    ipoly = (interp <maxVal).astype(np.int)
    cpoly = np.logical_or(bpoly, tpoly)
    rpoly = np.logical_not(bpoly)
    dpoly = np.logical_xor(bpoly, cpoly)
    interp = np.where(dpoly, interp, bag)
    print ('done', datetime.datetime.now())
    return interp

def triangulateSurfaces(grids, combo, vals):
    print ('triangulateSurfaces') 
    bagObj = grids[-1]
    bag = bagObj[-1]
    tifObj = grids[0]
    poly = tifObj[-1]
    maxVal = maxValue(bag)   
    print ('try tri', datetime.datetime.now())
    x, y = np.arange(bag.shape[1]), np.arange(bag.shape[0])
    xi, yi = np.meshgrid(x, y)
    values = scipy.interpolate.griddata(combo, vals, (xi, yi), method='linear', fill_value=maxVal)
    print ('done', datetime.datetime.now())
    maxVal = maxValue(values)
    print (values)
    ret = 0
    values = np.asarray(values, dtype='float64')
    values[np.isnan(values)]=maxVal
    print (values.shape, poly.shape)
    grid = rePrint(bag,values,poly,maxVal)
    return grid

def triangulateSurfaces2(grids, combo, vals):
    print ('triangulateSurfaces2') 
    bagObj = grids[-1]
    bag = bagObj[-1]
    tifObj = grids[0]
    poly = tifObj[-1]
    maxVal = maxValue(bag)   
    print ('try tri', datetime.datetime.now())
    x, y = np.arange(bag.shape[1]), np.arange(bag.shape[0])
    xi, yi = np.meshgrid(x, y)
    values = scipy.interpolate.griddata(combo, vals, (xi, yi), method='linear', fill_value=maxVal)
    print ('done', datetime.datetime.now())
    maxVal = maxValue(values)
    print (values)
    ret = 0
    values = np.asarray(values, dtype='float64')
    values[np.isnan(values)]=maxVal
    print (values.shape, poly.shape)
    grid = rePrint2(bag,values,poly,maxVal)
    return grid

def bagSplice(bag, combo):
    bagObj = bag[-1]
    tifObj = combo[-1]
    bags, combos = [], []
    return

def bagSave(bag, new, tifs, res, ext):
    for tif in tifs:
        gd_obj = gdal.Open(tif[1])
        break
    print ('bagSave')
    oextX, oextY = bag[2]
    nx, ny = ext[0]
    sx, sy = ext[-1]
    reso = float(res)
    print (nx, ny, sx, sy)
    gtran = (nx, reso, 0.0, sy, 0.0, -(reso))
    print (gtran)
    fName = bag[1].split('/')[-1]
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
    write_raster(new, gtran, gd_obj, outputpath, dtype=gdal.GDT_Float64, nodata=1000000.0)
    shutil.copy2(bag[1], outputpath2)
    with tb.open_file(outputpath2, mode = 'a') as bagfile:
        new = np.flipud(new)
        bagfile.root.BAG_root.elevation[:,:] = new
        bagfile.flush()
    bagfile.close()
    gd_obj = None
    print ('done')

path = 'C:/Users/Casiano.Koprowski/Desktop/Testing Files/BAGs and SSS Mosaics for Interpolation/H12600'
fileList = tifInput(path)
print (fileList)
tifFiles, names = getTifElev(fileList)
tifGrids, extent = alignTifs(tifFiles)
for tif in tifGrids:
    print(tif, '\n')
comboTif, name = polyTifVals(tifGrids, path, names)
comboArr = [0, name, extent, comboTif]

bagList = bagInput(path)
print (bagList)
bagFiles, names = getBagLyrs(bagList)
for bag in bagFiles:
    print(bag, '\n')    
    res, grids, ext = alignGrids(bag, comboArr)
    combo, vals = comboGrid(grids)
#    break
    tifObj = grids[0]
    poly = tifObj[-1]
    print (combo.shape, poly.shape)
    newBag = triangulateSurfaces2(grids, combo, vals)
    bagSave(bag, newBag, tifGrids, res, ext)
    break

