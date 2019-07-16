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

import string
import numpy as np

_basedigits = string.digits + string.ascii_uppercase

def get_tile_name(tile_number: int, name_len: int = 7) -> str:
    """
    Provide the name for a tile number.  The tile number is a flattened 
    series of the 2 dimensional tile set.
    """
    base = len(_basedigits)
    name = []
    while tile_number > 0:
        res = tile_number % base
        name.append(_basedigits[res])
        tile_number -= res
        tile_number = int(tile_number / base)
    name.reverse()
    name = ''.join(name)
    name = name.zfill(name_len)
    return name

def generate_tile_labels(number_of_tiles: int):
    """
    Create a generator for the tile names within a certain tle number range.
    """
    n = 0
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

def get_tile_set(xy: str):
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
    xb = np.linspace(-180., 180., xn)
    yb = np.linspace(-90., 90., yn)