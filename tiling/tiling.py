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
from shapely.geometry import Polygon
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

def get_tile_set_bounds(xy: str):
    """
    Return the tile bounds for a given two character resolution name.
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
    xb = np.linspace(-180., 180., xn + 1)
    yb = np.linspace(-90., 90., yn + 1)
    return xb, yb

def get_tile_from_point(xy: str, point: list):
    """
    Get the bounds of a tile that includes a point.
    """
    xb,yb = get_tile_set_bounds(xy)
    x,y = point
    pxb, idx = _get_closest_bounds(xb, x)
    pyb, idy = _get_closest_bounds(yb, y)
    bounds = [pxb.min(), pyb.min(), pxb.max(), pyb.max()]
    lin_id = _get_subn(len(xb), [x], [y])[0,0]
    name = get_tile_name(lin_id)
    return bounds, name

def build_gpkg(xy_res: str, bbox: list, path: str = '.'):
    """
    Create a geopackage with the tile set.
    """
    # get the bounds subset
    xb,yb = get_tile_set_bounds(xy_res)
    x_min, y_min, x_max, y_max = bbox
    subx, idx = _get_subbounds(xb, x_min, x_max)
    suby, idy = _get_subbounds(yb, y_min, y_max)
    # get the linear indicies for the tiles to get the names
    rx = dig2num(xy_res[0])
    ncols = get_num_cols(rx)
    idn = _get_subn(ncols, idx, idy)
    c = 0
    # setup the gdal object
    id_field = 'TileID'
    name_field = 'TileName'
    driver = ogr.GetDriverByName("MEMORY")
    ds = driver.CreateDataSource("tmp")
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    lyr = ds.CreateLayer('Tessellation', srs, ogr.wkbPolygon)
    fid = ogr.FieldDefn(id_field, ogr.OFTString)
    fname = ogr.FieldDefn(name_field, ogr.OFTString)
    lyr.CreateField(fid)
    lyr.CreateField(fname)
    fd = lyr.GetLayerDefn()
    for m in range(len(suby) - 1):
        for n in range(len(subx) - 1):
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(subx[n], suby[m])
            ring.AddPoint(subx[n], suby[m+1])
            ring.AddPoint(subx[n+1], suby[m+1])
            ring.AddPoint(subx[n+1], suby[m])
            ring.AddPoint(subx[n], suby[m])
            tile = ogr.Geometry(ogr.wkbPolygon)
            tile.AddGeometry(ring)
            f = ogr.Feature(fd)
            f.SetGeometry(tile)
            f.SetField(id_field, str(idn[m,n]))
            f.SetField(name_field, get_tile_name(idn[m,n]))
            c += 1
            lyr.CreateFeature(f)
            f = None
            tile = None
    outfilename = os.path.join(path, f'{xy_res}_tesselation.gpkg')
    ogr.GetDriverByName("GPKG").CopyDataSource(ds, outfilename)
    ds = None

def get_shapely(xy_res: str, bbox: list):
    """
    Return a list of shapely polygons that intersect or are contained by the
    provided bounding box.
    
    Things that go over the date line will probably be goofy...
    """
    # get the bounds subset
    xb,yb = get_tile_set_bounds(xy_res)
    x_min, y_min, x_max, y_max = bbox
    subx, idx = _get_subbounds(xb, x_min, x_max)
    suby, idy = _get_subbounds(yb, y_min, y_max)
    # get the linear indicies for the tiles to get the names
    rx = dig2num(xy_res[0])
    ncols = get_num_cols(rx)
    idn = _get_subn(ncols, idx, idy)
    c = 0
    collect = []
    names = []
    for m in range(len(suby)-1):
        for n in range(len(subx)-1):
            p = Polygon([(subx[n], suby[m]), 
                         (subx[n], suby[m + 1]),
                         (subx[n + 1], suby[m + 1]),
                         (subx[n + 1], suby[m])])
            collect.append(p)
            names.append(get_tile_name(idn[m,n]))
            c += 1
    return collect, names
    
def _get_subbounds(bounds, a_min: float, a_max: float):
    """
    Return the subset that intersect or are contained by the bounds.
    """
    idx = np.nonzero((bounds > a_min) & (bounds < a_max))[0]
    idx_min = idx.min() - 1
    idx_max = idx.max() + 1
    idx = np.insert(idx, 0, [idx_min])
    idx = np.append(idx, [idx_max])
    subset = bounds[idx]
    return subset, idx

def _get_subn(dimx: int, idx, idy):
    """
    Get the global indicies for the subtiles.
    """
    idx = np.array(idx)
    idy = np.array(idy)
    n_subx = len(idx)
    n_suby = len(idy)
    idn = np.zeros(n_subx * n_suby, dtype = np.int)
    for n,m in enumerate(idy):
        idn[n * n_subx:(n+1) * n_subx] = m * dimx + idx
    idn.shape = (n_suby, n_subx)
    return idn

def _get_closest_bounds(bounds, a: float):
    """
    Get the bounds around a location from the entire bounds set.
    """
    b = np.zeros(2)
    # get the first point
    rel = np.abs(bounds - a)
    id1 = rel.argmin()
    b[0] = bounds[id1]
    # get the second point
    rel[id1] = np.nan
    id2 = np.nanargmin(rel)
    b[1] = bounds[id2]
    ida = min(id1, id2)
    return b, ida
    
    