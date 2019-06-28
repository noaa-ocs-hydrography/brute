# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:31:40 2019

@author: Casiano.Koprowski
"""
import numpy as _np
import scipy as _scipy
import astropy.convolution as _apc
from datetime import datetime as _dt

import matplotlib.pyplot as plt


def tupleGrid(grid, nodata):
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
    nodata : float
        The array's nodata value

    Returns
    -------
    numpy.array
        Array of indecies where nodata values meet data values
        in order x, y, z

    """
    points = []
    a = 0
    for x in range(grid.shape[1]):
        io = False
        for y in range(grid.shape[0]):
            if grid[y, x] == nodata:
                if grid[y-1, x] != nodata:
                    val = grid[y-1, x]
                    point = [x, y-1, val]
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


def concatGrid(arr_1, arr_2, nodata):
    """Takes an input of an array of grid objects and the assumed nodata value
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
    grids : list
        The BAG and GeoTiff objects
    maxVal : int
        The BAG data's nodata value

    Returns
    -------
    xy : numpy.array
        Array of indecies where nodata values meet data values
        in order x, y
    z : numpy.array
        Array of values where nodata values meet data values
        in order z

    """
    points_1 = tupleGrid(arr_1, nodata)
    points_2 = tupleGrid(arr_2, nodata)
    if len(points_1) == 0 or len(points_2) == 0:
        xy, z = [], []
    else:
        comb = _np.concatenate([points_1, points_2])
        comb.view('i8,i8,i8').sort(order=['f0', 'f1'], axis=0)
        grid = _np.hsplit(comb, [2, 4])
        z = grid[1].squeeze()
        xy = grid[0]
    return xy, z


def rePrint(bag_elev, bag_uncr, cov_array, ugrids, maxVal, ioVal, debug=False):
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
    tpoly = (tpoly < maxVal).astype(_np.int)
    # 2
    bpoly = (bag < maxVal).astype(_np.int)
    # 3
    cpoly = _np.logical_or(bpoly, tpoly)
    # 4
    dpoly = _np.logical_xor(bpoly, cpoly)
    # 5
    ibag = _np.where(dpoly, pbag, bag)
    # 6
    npoly = (ibag < maxVal).astype(_np.int)
    # 7
    fpoly = _np.logical_and(dpoly, npoly)
    # 8
    if not ioVal:
        nbag = _np.where(fpoly, interp, bag)
        nunc = _np.where(fpoly, iuncrt, uncr)
    elif ioVal:
        nbag = _np.where(fpoly, interp, maxVal)
        nunc = _np.where(fpoly, iuncrt, maxVal)
    print('done', _dt.now())
    polyList = [tpoly, bpoly, cpoly, dpoly, npoly, fpoly, ibag]
    plt.figure()
    for rast in polyList:
        plt.imshow(rast)
        plt.show()
    if not debug:
        # polyList = [fpoly, cpoly]
        return nbag, nunc, cpoly.astype(_np.int)
    elif debug:
        return nbag, nunc, polyList


class linear:
    """Interpolates input data and convolves the ouput of the interpolation, if
    applicable.

    Takes input bathy and coverage arrays (tile or complete data) as well as
    the uncertainty array.  This data is used to inform the shape/size of the
    resulting output of the interpolation function.

    The input bathy and covrg arrays are passed to :func:`concatGrid` in order
    to gather the edges of the each of the arrays and combine them into a
    single list of combined xy and z egde values. If xy and z are empty, the
    interpolation process is skipped and the original tile or complete data is
    passed back from the function. Otherwise, the lists xy, z, and the shape of
    the BAG data input are passed to :func:`scipy.interpolate.griddata` an
    output to the variable "grid_pre". The output of this function is also
    saved to a seperate variable "grid" while the original output is left
    untouched. The seperate output is then passed to
    :func:`astropy.convolution.convolve` to smooth the output.

    The uncertainty layer is calculated using catzoc and grid::

        >>> m, b = catzoc
        >>> uncrt = (grid*m)+b

    Parameters
    ----------
    bathy : numpy.array
        The input bathemetry data
    uncrt : numpy.array
        The input uncertainty data
    covrg : numpy.array
        The input coverage data
    catzoc : tuple
        The input values for uncertainty calculation
    nodata : float, optional
        The default value is 1000000.0, the nodata value associated with the
        BAG format

    Returns
    -------
    bathy : numpy.array
        The interpolated and convolved bathemetry
    uncrt : numpy.array
        The calculated uncertainty from the interpolated and convolved
        bathemetry
    unint : numpy.array
        The interpolated bathemetry before convolution is performed

    """
    def __init__(self, bathy, uncrt, covrg, catzoc, nodata=1000000.0):
        x, y = _np.arange(bathy.shape[1]), _np.arange(bathy.shape[0])
        xi, yi = _np.meshgrid(x, y)
        xy, z = concatGrid(bathy, covrg, nodata)
        print(xy, z)
        if len(xy) != 0:
            self.bathy, self.uncrt, self.unint = self._interpolate(xy, z, xi,
                                                                   yi, catzoc,
                                                                   nodata)
        else:
            self.bathy, self.uncrt, self.unint = bathy, uncrt, bathy

    def _interpolate(self, xy, z, xi, yi, uval, nodata):
        m, b = uval
        grid_pre = _scipy.interpolate.griddata(xy, z, (xi, yi),
                                               method='linear',
                                               fill_value=nodata)
        grid = grid_pre
        grid = _np.asarray(grid, dtype='float64')
        grid[grid > 0] = _np.nan
        kernel = _apc.Gaussian2DKernel(3)
        grid = _apc.convolve(grid, kernel)
        grid[_np.isnan(grid)] = nodata
        grid[grid >= 0] = nodata
        uncr = (grid*m)+b
        return grid, uncr, grid_pre


def sliceFinder(size, shape, res, var=5000):
    """Uses the file size of the bag to determine if the grid should be tiled.
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
    print('sliceFinder')
    if res < 1:
        b = 25/res
    elif res >= 1:
        b = 25
    if size <= 100000:
        tiles = 0
        return tiles, None, None
    elif size > 100000:
        yChunk = int(_np.round(var*_np.sqrt(shape[0]/shape[1])))  # y
        xChunk = int(_np.round(var*_np.sqrt(shape[1]/shape[0])))  # x
        ny = int(_np.ceil(shape[0]/yChunk))  # ny
        nx = int(_np.ceil(shape[1]/xChunk))  # nx
        print(ny, nx)
        tiles = ny*nx
        chunckGrid = _np.arange(tiles).reshape((ny, nx))
        sliceInfo = [b, yChunk, xChunk]
        print(tiles, chunckGrid, sliceInfo)
        return tiles, chunckGrid, sliceInfo


def chunk(arr, tile, mode=None, copy=None):
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


class tile:
    """tiles() serves as the data container for individual tile data. It's
    inputs inlcude sliceInfo=[buffer, height, width], chunkSlice=tile[y,x]
    (tile position), and the total height and width of the BAG grid.  This
    information is then used to calculate the indices of each tile (with and
    without a buffer) and the indices of the data within the buffered tiles.

    """
    def __init__(self, sliceInfo, chunkSlice, shape):
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
        self.yMin, self.yMax, self.xMin, self.xMax = 0, 0, 0, 0
        self.yBMin, self.yBMax, self.xBMin, self.xBMax = 0, 0, 0, 0
        self.yIMin, self.yIMax, self.xIMin, self.xIMax = 0, 0, 0, 0
        buffer = sliceInfo[0]
        yChunk = sliceInfo[1]
        xChunk = sliceInfo[2]
        self._calcSlice(buffer, yChunk, xChunk, chunkSlice, shape)

    def _calcSlice(self, buffer, yChunk, xChunk, chunkSlice, shape):
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
        # 1
        self.yMin = int(max(0, (yChunk * chunkSlice[0]) - buffer))
        self.yMax = int(min(yChunk * (chunkSlice[0]+1) + buffer, shape[0]))
        self.xMin = int(max(0, (xChunk * chunkSlice[1]) - buffer))
        self.xMax = int(min(xChunk * (chunkSlice[1]+1) + buffer, shape[1]))
        slices = [self.yMin, self.yMax, self.xMin, self.xMax]

        # 2
        self.yBMin = int(max(0, (yChunk * chunkSlice[0])))
        self.yBMax = int(min(yChunk * (chunkSlice[0]+1), shape[0]))
        self.xBMin = int(max(0, (xChunk * chunkSlice[1])))
        self.xBMax = int(min(xChunk * (chunkSlice[1]+1), shape[1]))
        tiles = [self.yBMin, self.yBMax, self.xBMin, self.xBMax]

        # 3
        self.yIMin = int(max(0, (self.yBMin-self.yMin)))
        yma = int(min((self.yMax-self.yBMax), shape[0]))
        if yma == 0:
            self.yIMax = int(yChunk + buffer)
        else:
            self.yIMax = -int(yma)
        self.xIMin = int(max(0, (self.xBMin-self.xMin)))
        xma = int(min((self.xMax-self.xBMax), shape[1]))
        if xma == 0:
            self.xIMax = int(xChunk + buffer)
        else:
            self.xIMax = -int(xma)
        borders = [self.yIMin, self.yIMax, self.xIMin, self.xIMax]

        print(slices, tiles, borders, sep='\n')
