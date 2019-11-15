# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:31:40 2019

@author: Casiano.Koprowski
"""

from datetime import datetime as _dt
from typing import Union

import astropy.convolution as _apc
import numpy as _np
import scipy as _scipy

#from fuse.utilities import array_edge_points


def tupleGrid(grid: _np.array, nodata: int):
    """
    Takes an input matrix and an assumed nodata value. The function iterates
    through the matrix and compiles a list of 'edge' points [[x, y, z], ...]
    where:

    1. the current value is not a nodata value and previous value was a nodata value.
        - sets io True, indicating that the next value to compare against should be a nodata value.
    2. the current value is a nodata value and the previous value was not a nodata value.
        - sets io False, indicating that the next value to compare against should not be a nodata value.
    3. returns a list of the points found.

    Parameters
    ----------
    grid
        An input array
    nodata
        The array's nodata value

    Returns
    -------
    numpy.array
        Array of indecies where nodata values meet data values in order x, y, z
    """

    points = []
    a = 0
    for x in range(grid.shape[1]):
        io = False
        for y in range(grid.shape[0]):
            if grid[y, x] == nodata:
                if grid[y - 1, x] != nodata:
                    val = grid[y - 1, x]
                    point = [x, y - 1, val]
                    if a == 1:
                        a += 1
                    points.append(point)
                    io = False
                else:
                    pass
            else:
                if not io:
                    val = grid[y, x]
                    point = [x, y, val]

                    if a == 0:
                        a += 1
                    points.append(point)
                    io = True
    return _np.array(points)


def regulararray_to_xyz(data: _np.array, nodata: float = None) -> _np.array:
    """
    Extract XYZ points from an array of data

    Parameters
    ----------
    data
        2D array of gridded data
    nodata
        value to exclude from point creation from the input grid

    Returns
    -------
    numpy.array
        N x 3 array of XYZ points
    """

    if nodata is None:
        nodata = _np.nan

    x_values, y_values = _np.meshgrid(_np.arange(data.shape[1]), _np.arange(data.shape[0]))

    return _np.stack((x_values[data != nodata], y_values[data != nodata], data[data != nodata]), axis=1)


def array_coverage(array: _np.array, nodata: float = None) -> _np.array:
    """
    Get a boolean array of where data exists in the given array.

    Parameters
    ----------
    array
        array of gridded data with dimensions (Z)YX
    nodata
        value where there is no data in the given array

    Returns
    -------
    numpy.array
        array of booleans indicating where data exists
    """

    if len(array.shape) > 2:
        array = _np.squeeze(array)

    if nodata is None:
        nodata = _np.nan

    # TODO find reduced generalization of band coverage
    if array.shape[0] == 3:
        coverage = (array[0, :, :] != nodata) | (array[1, :, :] != nodata) | (array[2, :, :] != nodata)
    else:
        coverage = array != nodata

    return coverage


def raster_edge(data: _np.array, nodata: float) -> _np.array:
    """
    Get the cells of the array bordering `nodata`.

    Parameters
    ----------
    data
        array of raster data
    nodata
        value for no data in raster

    Returns
    -------
    numpy.array
        boolean array of edge cells
    """

    elevation_coverage = array_coverage(data, nodata)

    horizontal_difference = _np.concatenate((_np.full((elevation_coverage.shape[0], 1), 0),
                                               _np.diff(_np.where(elevation_coverage, 1, 0), axis=1)), axis=1)
    vertical_difference = _np.concatenate((_np.full((1, elevation_coverage.shape[1]), 0),
                                             _np.diff(_np.where(elevation_coverage, 1, 0), axis=0)), axis=0)

    horizontal_edges = (horizontal_difference == 1) | _np.roll(horizontal_difference == -1, -1, axis=1)
    vertical_edges = (vertical_difference == 1) | _np.roll(vertical_difference == -1, -1, axis=0)

    return horizontal_edges | vertical_edges


def array_edge_points(data: _np.array, nodata: float) -> _np.array:
    """
    Get the points bordering `nodata` in the given array.

    Parameters
    ----------
    data
        array of raster data
    nodata
        value for no data in raster

    Returns
    -------
    numpy.array
        N x 3 array of points
    """

    return regulararray_to_xyz(_np.where(raster_edge(data, nodata), data, nodata), nodata)


def concatGrid(arr_1, arr_2, nodata: int, no_nan: bool = False, split: bool = True):
    """
    Takes an input of an array of grid objects and the assumed nodata value
    Passes the assumed nodata value and the arrays held within each of the
    listed grid objects to :func:`tupleGrid` for a return of an array of
    edge points for each grid [[x, y, z], ...]

    Takes the results of both tupleGrid calls and combines them into a single
    array of edge points.  Then uses .view().sort() to sort the combined
    products [[x, y, z], ...] based on ascending x, then column y.

    The results of the sorted point array are then split into a list of points
    [[x, y], ...] and values [[z], ...] and returned

    Parameters
    ----------
    arr_1
        param arr_2:
    nodata
        The BAG data's nodata value

    Returns
    -------

    """
    if arr_1[arr_1 != nodata].size != 0 and arr_2[arr_2 != nodata].size != 0:
        if not no_nan:
            points_1 = array_edge_points(arr_1, nodata)
            points_2 = array_edge_points(arr_2, nodata)
            comb = _np.concatenate([points_1, points_2])
        else:
            comb = array_edge_points(arr_1, nodata)
        comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
        grid = _np.hsplit(comb, [2, 4])
        if split:
            z = grid[1].squeeze()
            xy = grid[0]
            return xy, z
        else:
            return grid
    else:
        if split:
            return [], []
        else:
            return []


def rePrint(bag_elev: _np.array, bag_uncr: _np.array, cov_array: _np.array, ugrids: list, maxVal: _np.array,
            ioVal: Union[int, bool], debug: Union[int, bool] = False):
    """
    Uses a mix of interpolated and original bag and tif data in order to
    determine where new interpolated data should be applied. No interpolated
    data is used where original bag data exists.

    Steps:

    1. Create a binary grid of cov data, tpoly.
    2. Create a binary grid of bag data, bpoly.
    3. Create a combined grid of the "complete coverage" of data; cpoly, numpy.logical_or of bpoly, tpoly.
    4. Use bpoly to remove the bag coverage from the combined "complete coverage"; dpoly, numpy.logical_xor of bpoly, cpoly.
    5. Use dpoly to determine where interpolated data can be applied and apply original bag data everywhere else; ibag, numpy.where of dpoly, pbag, bag.
    6. Create a binary grid of ibag, npoly.
    7. Use npoly and dpoly to create a binary grid where only data is present; fpoly, numpy.logical_and of dpoly, nopoly.
    8. Finalize results by using fpoly to apply interpolated data where appropriate and original bag data everywhere else.

    Parameters
    ----------
    bag_elev :
        orignal bag data
    bag_uncr :
        bag uncertainty layer
    cov_array :
        array of non-bathymetry coverage
    ugrids :
        List of interpolated data objects for BAG depth and uncertainty data
    maxVal :
        bag no-data value
    ioVal :
        User input. Determines whether origninal and interpolated or only interpolated data is output
    debug :
        Whether or not polyList includes all or only the last step in the evaluation process (Default value = False)

    Returns
    -------

    """

    print('rePrint', _dt.now())
    print(maxVal)
    poly = cov_array
    bag = bag_elev
    uncr = bag_uncr
    interp = ugrids[0]
    iuncrt = ugrids[1]
    pbag = ugrids[2]
    rows, cols = bag.shape
    # 1
    tpoly = _np.nan_to_num(poly)
    del poly
    tpoly = (tpoly < maxVal).astype(_np.int)
    # 2
    bpoly = (bag < maxVal).astype(_np.int)
    # 3
    cpoly = _np.logical_or(bpoly, tpoly)
    del tpoly
    # 4
    dpoly = _np.logical_xor(bpoly, cpoly)
    del cpoly, bpoly
    # 5
    ibag = _np.where(dpoly, pbag, bag)
    # 6
    npoly = (ibag < maxVal).astype(_np.int)
    del ibag
    # 7
    fpoly = _np.logical_and(dpoly, npoly)
    del dpoly, npoly
    # 8
    if not ioVal:
        nbag = _np.where(fpoly, interp, bag)
        nunc = _np.where(fpoly, iuncrt, uncr)
    elif ioVal:
        nbag = _np.where(fpoly, interp, maxVal)
        nunc = _np.where(fpoly, iuncrt, maxVal)
    del fpoly
    print('done', _dt.now())
    return nbag, nunc


class Interpolate:
    """
    Interpolates input data and convolves the ouput of the interpolation, if
    applicable.

    Attributes
    ----------
    bathy : numpy.array
        The resultant bathemetry data
    uncrt : numpy.arry
        The resultant uncertainty data
    unint : numpy.array
        The pre-gaussian resultant bathemetry data
    """

    def __init__(self, method: str, bathy: _np.array, uncrt: _np.array, covrg: _np.array, catzoc: tuple = None,
                 nodata: float = 1000000.0):
        """
        Takes input bathy and coverage arrays (tile or complete data) as well as
        the uncertainty array.  This data is used to inform the shape/size of the
        resulting output of the interpolation function.

        The input bathy and covrg arrays are passed to :func:`concatGrid` in order
        to gather the edges of the each of the arrays and combine them into a
        single list of combined xy and z egde values. If xy and z are empty, the
        interpolation process is skipped and the original tile or complete data is
        passed back from the function.

        Parameters
        ----------
        method : str
            The interpolation method
        bathy : numpy.array
            The input bathemetry data
        uncrt : numpy.arry
            The input uncertainty data
        covrg : numpy.array
            The input coverage data
        catzoc : tuple, optional
            The input values for uncertainty calculation
        nodata : float, optional
            The default value is 1000000.0, the nodata value associated with
            the BAG format
        """

        xi, yi = _np.meshgrid(_np.arange(bathy.shape[1]), _np.arange(bathy.shape[0]))
        if method == 'linear':
            xy, z = concatGrid(bathy, covrg, nodata)
            print(xy, z)
            if len(xy) != 0:
                self.bathy, self.uncrt, self.unint = self._linear(xy, z, xi, yi, catzoc, nodata)
            else:
                self.bathy, self.uncrt, self.unint = bathy, uncrt, bathy
        elif method == 'kriging':
            xyz = concatGrid(bathy, covrg, nodata, no_nan=True, split=False)
            if len(xyz) != 0:
                self.bathy, self.uncrt, self.unint = self._kriging(xyz, xi, yi, catzoc, nodata)
            else:
                self.bathy, self.uncrt, self.unint = bathy, uncrt, bathy
        else:
            raise ValueError(f'Interpolation type "{method}" not recognized.')

    def _linear(self, xy, z, xi, yi, uval, nodata):
        """
        Linear interpolation of input points

        The lists xy, z, and the shape of
        the BAG data input are passed to :func:`scipy.interpolate.griddata` an
        output to the variable "grid_pre". The output of this function is also
        saved to a seperate variable "grid" while the original output is left
        untouched. The seperate output is then passed to
        :func:`astropy.convolution.convolve` to smooth the output.

        The uncertainty layer is calculated using catzoc and grid:
        >>> m, b = catzoc
        >>> uncrt = (grid*m)+b

        Returns
        -------

        """

        m, b = uval
        grid_pre = _scipy.interpolate.griddata(xy, z, (xi, yi), method='linear', fill_value=nodata)
        grid = grid_pre
        grid = _np.asarray(grid, dtype='float64')
        grid[grid > 0] = _np.nan
        kernel = _apc.Gaussian2DKernel(3)
        grid = _apc.convolve(grid, kernel)
        grid[_np.isnan(grid)] = nodata
        grid[grid >= 0] = nodata
        uncr = (grid * m) + b
        return grid, uncr, grid_pre

    def _kriging(self, xyz, xi, yi, nodata):
        raise NotImplementedError('The kriging method has not been implemented yet!')


def sliceFinder(size: int, shape: (int, int), res: float, var: int = 5000):
    """
    Uses the file size of the bag to determine if the grid should be tiled.
    If the file is less than 100Mb, the file will not be tiled.  If the file is
    large enough to tile, the number of tiles and index size of each tile will
    be calculated based on the ratio of the total size of each array.

    yChunk = var*sqrt(height/width)
    xChunk = var*sqrt(width/height)

    ny = _np.ceil(height/yChunk)
    nx = _np.ceil(width/xChunk)

    tiles = nx*ny

    chunkgird is both the arrangement of tiles in relation to the grid and the
    order in which the tiles are processed.

    Parameters
    ----------
    size :
        Size of the input BAG file
    shape :
        Dimensions of the input BAG data (y, x)
    res :
        Resolution of the input BAG data
    var :
        Arbitrary value for determining chunk size (Default value = 5000)

    Returns
    -------

    >>> chunkGrid = _np.arrange(tiles).reshape((ny, nx))
    array([[ 0,  1,  2,  3,  4,  5],
           [ 6,  7,  8,  9, 10, 11],
           [12, 13, 14, 15, 16, 17],
           [18, 19, 20, 21, 22, 23],
           [24, 25, 26, 27, 28, 29],
           [30, 31, 32, 33, 34, 35]])
    # 36 Total Tiles
    # Tile 3 is index [0, 2] and has a value of 2
    """

    print('sliceFinder')
    if res < 1 and size <= 100000:
        size /= res if size > 50000 else res/2
    if res == 1 and size <= 100000:
        size += res
    elif res > 1 and res < 4 and size <= 100000:
        size *= res

    if res < 1:
        b = int(25 / res)
    elif res >= 1 and res <= 4:
        b = int(25 * res)
    elif res >= 4:
        b = 25

    if size <= 100000:
        tiles = 0
        return tiles, None, None
    elif size > 100000:
        yChunk = int(_np.round(var * _np.sqrt(shape[0] / shape[1])))  # y
        xChunk = int(_np.round(var * _np.sqrt(shape[1] / shape[0])))  # x
        ny = int(_np.ceil(shape[0] / yChunk))  # ny
        nx = int(_np.ceil(shape[1] / xChunk))  # nx
        print(ny, nx)
        tiles = ny * nx
        chunckGrid = _np.arange(tiles).reshape((ny, nx))
        sliceInfo = [b, yChunk, xChunk]
        print(tiles, chunckGrid, sliceInfo)
        return tiles, chunckGrid, sliceInfo


def chunk(arr, tile, mode=None, copy=None):
    """
    TODO write description

    Parameters
    ----------
    arr :
        param tile:
    mode :
        Default value = None)
    copy :
        Default value = None)
    tile :


    Returns
    -------

    """

    if mode == 'a':
        arr = arr[tile.yMin:tile.yMax, tile.xMin:tile.xMax]
        return arr
    elif mode == 'b':
        arr = arr[tile.yIMin:tile.yIMax, tile.xIMin:tile.xIMax]
        return arr
    elif mode == 'c':
        copy[tile.yBMin:tile.yBMax, tile.xBMin:tile.xBMax] = arr[tile.yIMin:tile.yIMax, tile.xIMin:tile.xIMax]
        return copy
    elif mode is None:
        raise ValueError("Mode value required.")


class BagTile:
    """
    tiles() serves as the data container for individual tile data. It's
    inputs inlcude sliceInfo=[buffer, height, width], chunkSlice=tile[y,x]
    (tile position), and the total height and width of the BAG grid.  This
    information is then used to calculate the indices of each tile (with and
    without a buffer) and the indices of the data within the buffered tiles.

    Parameters
    ----------

    Returns
    -------

    """

    def __init__(self, sliceInfo, chunkSlice, shape):
        """
        Takes the complete information of a give individual tile and
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
        # Tile 3 is index [0, 2] and has a value of 2
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

        :param sliceInfo:
        :param chunkSlice:
        :param shape:
        """

        self.yMin, self.yMax, self.xMin, self.xMax = 0, 0, 0, 0
        self.yBMin, self.yBMax, self.xBMin, self.xBMax = 0, 0, 0, 0
        self.yIMin, self.yIMax, self.xIMin, self.xIMax = 0, 0, 0, 0
        buffer = sliceInfo[0]
        yChunk = sliceInfo[1]
        xChunk = sliceInfo[2]
        self._calcSlice(buffer, yChunk, xChunk, chunkSlice, shape)

    def _calcSlice(self, buffer, yChunk, xChunk, chunkSlice, shape):
        """
        Takes the complete information of a give individual tile and
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

        Parameters
        ----------
        buffer :
            param yChunk:
        xChunk :
            param chunkSlice:
        shape :

        yChunk :

        chunkSlice :


        Returns
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
        """

        # 1
        self.yMin = int(max(0, (yChunk * chunkSlice[0]) - buffer))
        self.yMax = int(min(yChunk * (chunkSlice[0] + 1) + buffer, shape[0]))
        self.xMin = int(max(0, (xChunk * chunkSlice[1]) - buffer))
        self.xMax = int(min(xChunk * (chunkSlice[1] + 1) + buffer, shape[1]))
        slices = [self.yMin, self.yMax, self.xMin, self.xMax]

        # 2
        self.yBMin = int(max(0, (yChunk * chunkSlice[0])))
        self.yBMax = int(min(yChunk * (chunkSlice[0] + 1), shape[0]))
        self.xBMin = int(max(0, (xChunk * chunkSlice[1])))
        self.xBMax = int(min(xChunk * (chunkSlice[1] + 1), shape[1]))
        tiles = [self.yBMin, self.yBMax, self.xBMin, self.xBMax]

        # 3
        self.yIMin = int(max(0, (self.yBMin - self.yMin)))
        yma = int(min((self.yMax - self.yBMax), shape[0]))
        if yma == 0:
            self.yIMax = int(yChunk + buffer)
        else:
            self.yIMax = -int(yma)
        self.xIMin = int(max(0, (self.xBMin - self.xMin)))
        xma = int(min((self.xMax - self.xBMax), shape[1]))
        if xma == 0:
            self.xIMax = int(xChunk + buffer)
        else:
            self.xIMax = -int(xma)
        borders = [self.yIMin, self.yIMax, self.xIMin, self.xIMax]

        print(slices, tiles, borders, sep='\n')
