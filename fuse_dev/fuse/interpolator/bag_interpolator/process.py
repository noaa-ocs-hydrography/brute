# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:30:21 2019

@author: Casiano.Koprowski
"""

from datetime import datetime as _dt

import numpy as _np
from fuse.proc_io.proc_io import ProcIO

from . import bag as _bag
from . import coverage as _cvg
from . import interpolator as _itp

_catZones = {
    'A1': (.01, .5),
    'A2/B': (.02, 1),
    'C': (.05, 2)
}


class RasterInterpolator:
    """TODO write description"""

    def __init__(self):
        ...

    def interpolate(self, dataset, interpolation_type: str, resolution: float,
                    support_files: list, catzoc: str = 'A2/B',
                    io: bool = False):
        if interpolation_type == 'linear':
            self._linear(dataset, support_files, catzoc, io)
        else:
            raise ValueError('Interpolation type not implemented:'
                             + f'{interpolation_type}')

    def _linear(dataset, support_files: list, catzoc: str, io: bool):
        """
        Linear interpolation of a raster using supporting files as a mask

        Parameters
        ----------
        dataset :
            A gdal Dataset
        support_files

        Returns
        -------
        dataset
            A gdal Dataset

        """
        uval = _catZones.get(catzoc)

        bag = _bag.BagFile()
        bag.from_gdal(dataset)

        coverage = _cvg.UnifiedCoverage(support_files, bag.wkt)
        coverage = _cvg.align2grid(coverage, bag.bounds, bag.shape, bag.resolution, bag.nodata)

        z, tiles, tile_info = _itp.sliceFinder(bag.size, bag.shape, bag.resolution)

        if z > 1:
            unitedBag = _np.empty_like(bag.elevation)
            unitedUnc = _np.empty_like(bag.elevation)
            unitedPre = _np.empty_like(bag.elevation)
            bagShape = bag.shape

            for ySlice in range(tiles.shape[0]):
                for xSlice in range(tiles.shape[1]):
                    ts = _dt.now()
                    index = ySlice, xSlice
                    print(f'\nTile {tiles[index] + 1} of {z} - {ts}')
                    tile = _itp.BagTile(tile_info, index, bagShape)
                    covgTile = _itp.chunk(coverage.array, tile, mode='a')
                    bathTile = _itp.chunk(bag.elevation, tile, mode='a')
                    uncrTile = _itp.chunk(bag.uncertainty, tile, mode='a')
                    print('interp is next')
                    interp = _itp.LinearInterpolator(bathTile, uncrTile, covgTile, uval)
                    del covgTile, bathTile, uncrTile
                    unitedBag = _itp.chunk(interp.bathy, tile, mode='c',
                                           copy=unitedBag)
                    unitedUnc = _itp.chunk(interp.uncrt, tile, mode='c',
                                           copy=unitedUnc)
                    unitedPre = _itp.chunk(interp.unint, tile, mode='c',
                                           copy=unitedPre)
                    del interp
                    td = _dt.now()
                    tdelt = td - ts
                    print('Tile complete -', td, '| Tile took:', tdelt)

            ugrids = [unitedBag, unitedUnc, unitedPre]
            del unitedBag, unitedUnc, unitedPre
        else:
            ts = _dt.now()
            print('\nTile 1 of 1 -', ts)
            print('interp is next')
            interp = _itp.LinearInterpolator(bag.elevation, bag.uncertainty,
                                             coverage.array, uval)
            ugrids = [interp.bathy, interp.uncrt, interp.unint]
            del interp
            td = _dt.now()
            tdelt = td - ts
            print('Tile complete -', td, '| Tile took:', tdelt)

        bag.elevation, bag.uncertainty, coverage.array = _itp.rePrint(bag.elevation, bag.uncertainty,
                                                                      coverage.array, ugrids, bag.nodata, io)
        print(coverage.array)

        save = _bag.BagToGDALConverter('MLLW')
        save.bag2gdal(bag)

        del coverage, bag, ugrids

        return save.dataset

class Intitializor:
    """TODO write description"""

    def __init__(self, outlocation: str, mode: str, catzoc: str, io: bool):
        """

        Parameters
        ----------
        outlocation
        mode
        catzoc
        io
        """

        self._outlocation = outlocation
        self._uval = _catZones.get(catzoc)
        self._io = io

        if mode == 'linear':
            pass
        else:
            raise ValueError('Interpolation type not implemented.')

    def linear(self, filepath: str, coverage_list: list):
        """
        TODO write description

        Parameters
        ----------
        filepath: str :
            TODO write description
        coverage_list: list :
            TODO write description

        Returns
        -------

        """

        bag = _bag.BagFile()
        bag.open_file(filepath, 'hack')
        bag.generate_name(self._outlocation, self._io)
        coverage = _cvg.UnifiedCoverage(coverage_list, bag.wkt, bag.name)
        coverage = _cvg.align2grid(coverage, bag.bounds, bag.shape, bag.resolution, bag.nodata)

        z, tiles, tile_info = _itp.sliceFinder(bag.size, bag.shape, bag.resolution[0])

        if z > 1:
            unitedBag = _np.empty_like(bag.elevation)
            unitedUnc = _np.empty_like(bag.elevation)
            unitedPre = _np.empty_like(bag.elevation)
            bagShape = bag.shape

            for ySlice in range(tiles.shape[0]):
                for xSlice in range(tiles.shape[1]):
                    ts = _dt.now()
                    index = ySlice, xSlice
                    print(f'\nTile {tiles[index] + 1} of {z} - {ts}')
                    tile = _itp.BagTile(tile_info, index, bagShape)
                    covgTile = _itp.chunk(coverage.array, tile, mode='a')
                    bathTile = _itp.chunk(bag.elevation, tile, mode='a')
                    uncrTile = _itp.chunk(bag.uncertainty, tile, mode='a')
                    print('interp is next')
                    interp = _itp.LinearInterpolator(bathTile, uncrTile, covgTile, self._uval)
                    del covgTile, bathTile, uncrTile
                    unitedBag = _itp.chunk(interp.bathy, tile, mode='c',
                                           copy=unitedBag)
                    unitedUnc = _itp.chunk(interp.uncrt, tile, mode='c',
                                           copy=unitedUnc)
                    unitedPre = _itp.chunk(interp.unint, tile, mode='c',
                                           copy=unitedPre)
                    del interp
                    td = _dt.now()
                    tdelt = td - ts
                    print('Tile complete -', td, '| Tile took:', tdelt)

            ugrids = [unitedBag, unitedUnc, unitedPre]
            del unitedBag, unitedUnc, unitedPre
        else:
            ts = _dt.now()
            print('\nTile 1 of 1 -', ts)
            print('interp is next')
            interp = _itp.LinearInterpolator(bag.elevation, bag.uncertainty,
                                             coverage.array, self._uval)
            ugrids = [interp.bathy, interp.uncrt, interp.unint]
            del interp
            td = _dt.now()
            tdelt = td - ts
            print('Tile complete -', td, '| Tile took:', tdelt)

        bag.elevation, bag.uncertainty, coverage.array = _itp.rePrint(bag.elevation, bag.uncertainty,
                                                                      coverage.array, ugrids, bag.nodata, self._io)
        print(coverage.array)

        save = _bag.BagToGDALConverter('MLLW')
        #        save.components2gdal([bag.elevation, bag.uncertainty], bag.shape,
        #                             bag.bounds, bag.resolution, bag.wkt, bag.nodata)
        save.bag2gdal(bag)

        writer = ProcIO('gdal', 'bag')
        print(save.dataset.GetGeoTransform())
        writer.write(save.dataset, bag.outfilename)

        _cvg.write_vector(coverage, self._outlocation)

        del coverage, bag, save, ugrids
