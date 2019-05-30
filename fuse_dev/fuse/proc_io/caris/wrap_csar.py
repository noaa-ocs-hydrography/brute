# -*- coding: utf-8 -*-
"""
wrap_csar.py

Created on Thu Feb 14 15:11:39 2019

@author: grice
"""
import os
import sys
import pickle
import numpy as np
import logging
import caris.coverage as cc

def write_csar(dataset, m):
    """
    Convert a gdal dataset into a csar.
    http://www.teledynecaris.com/en/support/caris-python-api/5-1/coverage/raster/intro.html#creating-a-raster
    """
    print ('write_csar')
    dataset = np.array(dataset)
    print (m, dataset, dataset.shape)
    z_dir = cc.Direction.HEIGHT
    layerName = "Elevation"
    if not m['z_up']:
        z_dir = cc.Direction.DEPTH
        layerName = "Depth"
    band_info = cc.BandInfo(name=layerName,
                     type = cc.DataType.FLOAT32,
                     tuple_length = 1,
                     direction = z_dir,
                     units = 'm',
                     category = cc.Category.SCALAR,
                     level_policy = cc.LevelPolicy.BICUBIC)
    resolution = [m['resy'], m['resx']]
    origin = [m['originx'],m['originy']]
#    origin = [0,0]
    dimensions = [m['dimx'], m['dimy']]
    crs = m['crs']
    name = m['outfilename']
    bands = [band_info]

    if os.path.exists(name):
        os.remove(name)
        os.remove(name+'0')
    raster = cc.create_raster(name, crs, origin, resolution, dimensions, bands)

    idx = (dataset < m['nodata']).astype(np.int)
#    idx = np.nonzero(dataset == m['nodata'])
    if not m['z_up']:
        dataset = np.where(idx, dataset, raster.band_info['Depth'].ndv)
    #    # write the data into the csar container
        band_dtype = raster.band_info['Depth'].numpy_dtype
        area = ((0,0),(dimensions[0],dimensions[1]))
        raster.write("Depth", area, dataset.astype(band_dtype))
    else:
        dataset = np.where(idx, dataset, raster.band_info['Elevation'].ndv)
    #    # write the data into the csar container
        band_dtype = raster.band_info['Elevation'].numpy_dtype
        area = ((0,0),(dimensions[0],dimensions[1]))
        raster.write("Elevation", area, dataset.astype(band_dtype))

    raster = None

def check_metadata(meta):
    """
    Check to make sure the required metadata keys are available in the
    provided metadata dictionary.
    """
    req_attrib = {'resx',
                  'resy',
                  'originx',
                  'originy',
                  'dimx',
                  'dimy',
                  'crs',
                  'nodata',
                  'outfilename',
                  'z_up',
                  }
    mkeys = ''
    for key in req_attrib:
        if key not in meta:
            mkeys = mkeys + key + ', '
    if len(mkeys) > 0:
        raise ValueError('Metadata missing to write csar %s' % mkeys)

def main():
    """
    Parse the arguments and send them to the write method.
    """
    # check to make sure the file exists
    data = np.load(sys.argv[1])
    # check to make sure the metadata file exists
    with open(sys.argv[2], 'rb') as metafile:
        metadata = pickle.load(metafile)
    # read the metadata into variables and send to the write method
    check_metadata(metadata)
    write_csar(data, metadata)


if __name__ == '__main__':
    main()
