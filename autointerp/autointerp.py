# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 16:52:49 2018

@author: Casiano.Koprowski

Made in part with code from:
    - Extract Bag Outlines on Pydro by Dayong.Shen
    - QC Tools 2: Flier Finder on Pydro by
    - QC Tools 2: Bag Metadata Editor on Pydro by
"""
#from __future__ import absolute_import, division, print_function

import os
import re
import datetime
import scipy
import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import tables as tb
import shapely as sp

from shapely.geometry import MultiPoint
from shapely.ops import cascaded_union
from osgeo import gdal, ogr, osr, gdalconst
from string import ascii_lowercase

# this allows GDAL to throw Python Exceptions
gdal.UseExceptions()

combo = re.compile(r'COMBINED', re.IGNORECASE)

def tifInput(path):
    inFiles = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.tif') or file.lower().endswith('.tiff'):
                inFiles.append(root + '/' + file)
    return inFiles

def bagInput(path):
    inFiles = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.lower().endswith('.bag'):
                inFiles.append(root + '/' + file)
    return inFiles

def getTifElev(files):   
    tifFiles = []
    names = []
    x = 0
    for file in files:
        print(file)
        if not combo.search(file):
            ds = gdal.Open(file)
            fName = os.path.split(file)[-1]
            splits = fName.split('_')
            name = '_'.join([x for x in splits[:2]])
            print (fName, splits, name)
            names.append(name)
            band = ds.GetRasterBand(1)
            arr = band.ReadAsArray()
            [cols, rows] = arr.shape
            print(arr.shape)
            print (np.amax(arr), np.amin(arr))
            band=None
            ds=None  
            tifFiles.append ([x, file, arr])
            x += 1
    return tifFiles, names


def getBagLyrs(files):
    bagFiles = []
    names = []
    x = 0
    for file in files:
        fName = os.path.split(file)[-1]
        splits = fName.split('_')
        name = '_'.join([x for x in splits[:2]])
        print (fName, splits, name)
        names.append(name)
        with tb.open_file(file, mode = 'r') as bagfile:
            elev = bagfile.root.BAG_root.elevation.read()
            uncr = bagfile.root.BAG_root.uncertainty.read()
            print (np.amax(elev), np.amin(elev))
            print (np.amax(uncr), np.amin(uncr))
            bagFiles.append([x, file, elev, uncr])
            print(elev.shape)
            bagfile.close()
        x += 1
    return bagFiles, names

def write_raster(raster_array, gt, data_obj, outputpath, dtype=gdal.GDT_UInt16, options=0, color_table=0, nbands=1, nodata=False):
    '''
    Directly From:
    https://gis.stackexchange.com/questions/278946/what-is-the-simplest-way-to-get-the-average-of-multiple-raster-images-while-igno
    '''
    
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
         
def polyTifVals(tifs, path, names):
    '''
    Heavy Influence From:
    https://gis.stackexchange.com/questions/278946/what-is-the-simplest-way-to-get-the-average-of-multiple-raster-images-while-igno
    
    Needs: interpolation/smoothing to final merged product
    '''
    letters = []
    for letter in ascii_lowercase:
        letters.append(letter)
    maxRes = [0,0]
    maxVal = 0
    tifList = []
    for tif in tifs:
        if np.amax(tif[-1]) >= maxVal:
            maxVal = np.amax(tif[-1])
        [rows, cols] = tif[-1].shape
        if rows >= maxRes[0]:
            maxRes[0] = rows
        if cols >= maxRes[1]:
            maxRes[1] = cols
        tifList.append([tif[0], tif[-1]])
    tifDict = dict(tifList)
    
    for i, tif in tifDict.items(): 
        gd_obj = gdal.Open(tifs[i][1])
        array = np.expand_dims(tif,2)
        if i == 0:
            allarrays = array
        else:
            allarrays = np.concatenate((allarrays, array), axis=2)
    
    polyTiff = (allarrays < maxVal).astype(np.int)
    polyTiff = np.nanmean(polyTiff, axis=2)
    polyTiff = (polyTiff > 0).astype(np.int)
          
    outputpath = path + '/' + names[0] + '_COMBINEDPOLY.tif'
    
    print(outputpath)
    write_raster(polyTiff, gd_obj.GetGeoTransform(), gd_obj, outputpath)
    gd_obj = None
    
    print ('done?')
    return outputpath

def tupleGrid(grid, maxVal):
    print ('tupleGrid')
    points = []
    for x in range(grid.shape[1]):
        for y in range(grid.shape[0]):
            if grid[y,x] == maxVal:
                pass
            else:
                point = (x, y, grid[y,x])
                points.append(point)
    print ('done')
    return points

def triangulateSurfaces(bags):
    print ('triangulateSurfaces')
    grids = []
    for bag in bags:
        maxVal = np.amax(bag[-2])
#        print (maxVal)
        elev = tupleGrid(bag[-2], maxVal)
        print ('MultiPoint Conversion')
        elev = MultiPoint(elev)
        mask = (bag[-2] < maxVal).astype(np.int)
        print ('done')
        grids.append([elev, mask])
        break
    
    triangles = []
    for x in grids:
        print ('try tri', datetime.datetime.now())
        tri = scipy.spatial.Delaunay(x)
        triangles.append(tri)
        print ('done', datetime.datetime.now())
    return triangles


path = 'C:/Users/Casiano.Koprowski/Desktop/Testing Files/BAGs and SSS Mosaics for Interpolation/H12525'
fileList = tifInput(path)
print (fileList)
tifFiles, names = getTifElev(fileList)
for item in tifFiles:
    print(item, '\n')
comboTif = polyTifVals(tifFiles, path, names)

bagList = bagInput(path)
print (bagList)
bagFiles, names = getBagLyrs(bagList)
triBagAreas = triangulateSurfaces(bagFiles)
