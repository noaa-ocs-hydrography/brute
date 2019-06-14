# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:30:21 2019

@author: Casiano.Koprowski
"""

import numpy as _np
from datetime import datetime as _dt

from . import bag as _bag
from . import coverage as _cvg
from . import interpolator as _itp
from fuse.proc_io.proc_io import proc_io

_catZones = {
        'A1':(.01,.5),
        'A2/B':(.02,1),
        'C':(.05,2)}

class intitialize:
    def __init__(self, outlocation, mode, catzoc, io):
        self._outlocation = outlocation
        self._uval = _catZones.get(catzoc)
        self._io = io
        if mode == 'linear':
            pass
        else:
            raise ValueError('Interpolation type not implemented.')

    def linear(self, filepath, coverage_list):
        bag = _bag.bag_file()
        bag.open_file(filepath, 'hyo')
        bag.generate_name(self._outlocation, self._io)
        coverage = _cvg.unified_coverage(coverage_list, bag.wkt, bag.name)
        coverage = _cvg.align2grid(coverage, bag.bounds, bag.shape,
                                   bag.resolution, bag.nodata)

        z, tiles, tile_info = _itp.sliceFinder(bag.size, bag.shape,
                                               bag.resolution[0])
        if z > 1:
            unitedBag = _np.empty_like(bag.elevation)
            unitedUnc = _np.empty_like(bag.elevation)
            unitedPre = _np.empty_like(bag.elevation)
            bagShape = bag.shape
            for ySlice in range(tiles.shape[0]):
                for xSlice in range(tiles.shape[1]):
                    ts = _dt.now()
                    index = ySlice, xSlice
                    print ('\nTile', tiles[index]+1, 'of', z, '-', ts)
                    tile = _itp.tile(tile_info, index, bagShape)
                    covgTile = _itp.chunk(coverage.array, tile, mode='a')
                    bathTile = _itp.chunk(bag.elevation, tile, mode='a')
                    uncrTile = _itp.chunk(bag.uncertainty, tile, mode='a')
                    print ('interp is next')
                    interp = _itp.linear(bathTile,uncrTile,covgTile,self._uval)
                    covgTile = bathTile = uncrTile = None
                    unitedBag = _itp.chunk(interp.bathy, tile, mode='c',
                                           copy=unitedBag)
                    unitedUnc = _itp.chunk(interp.uncrt, tile, mode='c',
                                            copy=unitedUnc)
                    unitedPre = _itp.chunk(interp.unint, tile, mode='c',
                                           copy=unitedPre)
                    interp = None
                    td = _dt.now()
                    tdelt = td - ts
                    print ('Tile complete -', td, '| Tile took:', tdelt)
            ugrids = [unitedBag, unitedUnc, unitedPre]
        else:
            ts = _dt.now()
            print ('\nTile 1 of 1 -', ts)
            print ('interp is next')
            interp = _itp.linear(bag.elevation, bag.uncertainty,
                                 coverage.array, self._uval)
            ugrids = [interp.bathy, interp.uncrt, interp.unint]
            td = _dt.now()
            tdelt = td - ts
            print ('Tile complete -', td, '| Tile took:', tdelt)
        bag.elevation, bag.uncertainty, coverage.array = _itp.rePrint(bag.elevation, bag.uncertainty,
                                                           coverage.array, ugrids, bag.nodata, self._io)
        print (coverage.array)

        save = _bag.gdal_create('MLLW')
#        save.components2gdal([bag.elevation, bag.uncertainty], bag.shape,
#                             bag.bounds, bag.resolution, bag.wkt, bag.nodata)
        save.bag2gdal(bag)

        writer = proc_io('gdal', 'bag')
        writer.write(save.dataset, bag.outfilename)

        _cvg.write_vector(coverage, self._outlocation)

        coverage = None
        bag = None
        save = None


