# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:43:55 2019

@author: Casiano.Koprowski
"""

import os

from osgeo import gdal
from proc_io import proc_io

progLoc = os.getcwd()
fileLoc = 'R:\\Scripts\\vlab-nbs\\fuse_dev\\fuse'

# bSize = [5000,5000]
# nodata = 100000.0
# prj = '''PROJCS["unnamed",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'''
# origin = ([605260.5, 4494884.5], [620521.5, 4483923.5])
# bag_grid = np.random.randint(-20, high=0, size=bSize,dtype='int32')
# print (bag_grid, bag_grid.shape)
# bag_data = ['test_bag', 4, origin, bSize, prj, nodata, bag_grid]

fileName = 'R:\\Scripts\\Testing Files\\BAGs and SSS Mosaics for Interpolation\H12525\H12525_SSS_100.tif'
'''
PROJCS["unnamed",
GEOGCS["NAD83",
       DATUM["North_American_Datum_1983",
       SPHEROID["GRS 1980",6378137,298.2572221010042,
AUTHORITY["EPSG","7019"]],
AUTHORITY["EPSG","6269"]],
    PRIMEM["Greenwich",0],
UNIT["degree",0.0174532925199433],
AUTHORITY["EPSG","4269"]],
PROJECTION["Transverse_Mercator"],
PARAMETER["latitude_of_origin",0],
PARAMETER["central_meridian",-75],
PARAMETER["scale_factor",0.9996],
PARAMETER["false_easting",500000],
PARAMETER["false_northing",0],
UNIT["metre",1,AUTHORITY["EPSG","9001"]]]
'''
# fileName = 'R:\\Scripts\\Testing Files\\BAGs and SSS Mosaics for Interpolation\\H12600\\H12600_SSS_1m_100__A.tif'
'''
PROJCS["NAD83 / UTM zone 18N",
GEOGCS["NAD83",
       DATUM["North_American_Datum_1983",
       SPHEROID["GRS 1980",6378137,298.2572221010042,
AUTHORITY["EPSG","7019"]],
AUTHORITY["EPSG","6269"]],
    PRIMEM["Greenwich",0],
UNIT["degree",0.0174532925199433],
AUTHORITY["EPSG","4269"]],
PROJECTION["Transverse_Mercator"],
PARAMETER["latitude_of_origin",0],
PARAMETER["central_meridian",-75],
PARAMETER["scale_factor",0.9996],
PARAMETER["false_easting",500000],
PARAMETER["false_northing",0],
UNIT["metre",1,
AUTHORITY["EPSG","9001"]],
AUTHORITY["EPSG","26918"]]
'''

bag_data = gdal.Open(fileName)

# print (gdal.Info(bag_data))

bagWrite = proc_io('gdal', 'csar', fileLoc + '\\rasterdata')
# print (bagWrite)
x = 0
while True:
    outname = f'write_test_{x}.csar'
    filePath = os.path.join(fileLoc, outname)
    print(outname, filePath)
    if os.path.exists(filePath):
        x += 1
    else:
        bagWrite.write(bag_data, outname)  # , z_up=False)
        break
print('done')
