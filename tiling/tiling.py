# -*- coding: utf-8 -*-
"""
tiling.py

Created on Tue Jul 16 04:34:58 2019

@author: grice

A set of utilities for building global quadralateral tesselation schemes.
These schemes are global with differing dimensions in the X (easting) and Y
(northing) directions.  The tile sizes are defined in degrees as

dimension = 180 / (2**R), R = [0, 1, 2, ...35]

The tesselation set starts with the lower left corner of the first cell at
longitude -180 deg and latitude -90 deg.  Tiles are named from west to east in rows
with the naming scheme continuing with the next row north once the last cell in
a row reaches +180 degrees on the right hand boundry.
"""

import os
import string
import numpy as np
import osgeo.ogr as ogr
import osgeo.osr as osr

_basedigits = string.digits + string.ascii_uppercase

def get_tile_name(tile_number: int, name_len: int = 7) -> str:
    """
    Provide the name for a tile number.  The tile number is a flattened 
    series of the 2 dimensional tile set.
    """
    base = len(_basedigits)
    name = []
    while tile_number > 0:
        tile_number, res = divmod(tile_number, base)
        name.append(_basedigits[res])
    name.reverse()
    name = ''.join(name)
    name = name.zfill(name_len)
    return name

def generate_tile_labels(number_of_tiles: int):
    """
    Create a generator for the tile names within a certain tle number range.
    """
    for n in range(number_of_tiles):
        name = get_tile_name(n)
        yield name
            
def get_dimension(r: int) -> float:
    """
    Return the dimension in degrees as defined by
    
    dimension = 180 / 2**r
    """
    return 180. / (2**r)

def get_num_cols(r: int) -> int:
    """
    Return the number of columns for a given "R" value assuming coverage from
    -180 to +180.
    """
    return 2**(r+1)

def get_num_rows(r: int) -> int:
    """
    Return the number of rows for a given "R" value and assuming coverage from
    -90 to + 90.
    """
    return 2**(r)

def dig2num(d: str) -> int:
    """
    Return the number associated with a digit.
    """
    return _basedigits.index(d)

def num2dig(n: int) -> str:
    """
    Return the digit associated with a number.
    """
    return _basedigits[n]

def get_tile_bounds(xy: str):
    """
    Return the tile set with names for a given two character resolution name.
    """
    try:
        assert len(xy) == 2
    except AssertionError:
        raise ValueError('The provided name must be two characters')
    xyu = xy.upper()
    x,y = xyu
    xr = dig2num(x)
    yr = dig2num(y)
    xn = get_num_cols(xr)
    yn = get_num_rows(yr)
    n = xn * yn
    names = generate_tile_labels(n)
    xb = np.linspace(-180., 180., xn + 1)
    yb = np.linspace(-90., 90., yn + 1)
    return names, xb, yb
    
def build_resolution_scheme(xy_res: str, path = '.'):
    """
    Create a geopackage with the tile set.
    """
    field_name = 'TileID'
    name,xb,yb = get_tile_bounds(xy_res)
    driver = ogr.GetDriverByName("GPKG")
    outfilename = os.path.join(path, f'{xy_res}_tesselation.gpkg')
    ds = driver.CreateDataSource(outfilename)
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    lyr = ds.CreateLayer('Tessellation', srs, ogr.wkbPolygon)
    field = ogr.FieldDefn(field_name, ogr.OFTString)
    lyr.CreateField(field)
    fd = lyr.GetLayerDefn()
    for m in range(len(yb) - 1):
        for n in range(len(xb) - 1):
            tile = _create_tile(xb[n], yb[m], xb[n+1], yb[m+1])
            f = ogr.Feature(fd)
            f.SetGeometry(tile)
            f.SetField(field_name, next(name))
            lyr.CreateFeature(f)
            f = None
    ds = None
    
def _create_tile(x_min: float, y_min: float, x_max: float, y_max: float):
    """
    Return an ogr polygon representing a single tile.
    """
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(x_min, y_min)
    ring.AddPoint(x_min, y_max)
    ring.AddPoint(x_max, y_max)
    ring.AddPoint(x_max, y_min)
    ring.AddPoint(x_min, y_min)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly
