# -*- coding: utf-8 -*-
"""
wrap_csar.py

Created on Thu Feb 14 15:11:39 2019

@author: grice
"""

import os
import pickle
import sys

import caris.coverage as cc
import numpy as np

def write_raster(dataset, m: dict):
    """
    Convert a gdal dataset into a csar.
    http://www.teledynecaris.com/en/support/caris-python-api/5-1/coverage/raster/intro.html#creating-a-raster

    Parameters
    ----------
    dataset :
        param m:
    dataset: gdal.Dataset :

    m: dict :


    Returns
    -------

    """
    # set the names of the layers in an assumed order
    band_names = ['Elevation', 'Uncertainty','NBSa','NBSb','NBSc']
    # make sure the file does not exist before we write because it corrupts the file.
    csarname = m['outfilename']
    while os.path.exists(csarname):
        try:
            os.remove(csarname)
        except:
            raise RuntimeError('Unable to remove csar file: {}'.format(csarname))
    while os.path.exists(csarname + '0'):
        try:
            os.remove(csarname + '0')
        except:
            raise RuntimeError('Unable to remove csar0 file: {}'.format(csarname))
    # make sure the data has the assumed shape and get the number of layers
    if len(dataset.shape) == 3:
        numlayers,dimy,dimx = dataset.shape
    elif len(dataset.shape) == 2:
        print('Single Band Grid')
        dimy, dimx = dataset.shape
        numlayers = 1
        dataset.shape = (1,shape[0],shape[1])
    else:
        mgs = 'Array with shape {} provided to csar write process'.format(dataset.shape)
        raise ValueError(msg)
    if numlayers > len(band_names):
        raise ValueError('Too many layers provided: {}'.format(numlayers))
    # get the metadata for building the layers
    resolution = [m['resx'], m['resy']]
    origin = [m['originx'], m['originy']]
    dimensions = [dimx, dimy]
    area = ((0, 0), (dimx, dimy))
    crs = m['crs']
    ndv = int(m['nodata'])

    bands = []
    for n in range(numlayers):
        # get the min / max, and change the vertical direction if needed.
        dataset_band = np.flipud(dataset[n,:,:])
        ndv_idx = np.nonzero(dataset_band == ndv)
        dataset_band[ndv_idx] = np.nan
        if n == 0:
            z_type = cc.Direction.HEIGHT
            d_min = np.nanmax(dataset_band) 
            d_max = np.nanmin(dataset_band)
            if not m['z_up']:
                print('Reversing vertical direction of data')
                dataset_band *= -1
        else:
            z_type = cc.Direction.NAP
            d_min = np.nanmin(dataset_band) 
            d_max = np.nanmax(dataset_band) 
       
        # build the band and create the layer
        
        band_info = cc.BandInfo(name=band_names[n], type=cc.DataType.FLOAT32, tuple_length=1, direction=z_type, units='m', category=cc.Category.SCALAR, minimum=d_min, maximum=d_max, level_policy=cc.LevelPolicy.MAX)
        
        bands.append(band_info)
        
    raster = cc.create_raster(csarname, crs, origin, resolution, dimensions, bands)
    for n in range(numlayers):
        dataset_band = np.flipud(dataset[n,:,:])
        # set no data value and write the data to the layer
        dataset_band[ndv_idx] = raster.band_info[band_names[n]].ndv
        band_dtype = band_info.numpy_dtype
        raster.write(band_names[n], area, dataset_band.astype(band_dtype))
    del raster


def write_points(dataset, m: dict):
    """
    Convert a set of GDAL points to a CSAR point cloud.  The provided data is
    assumed to be a depth (positive down) and is assigned to a height
    (positive up).

    Parameters
    ----------
    dataset :
        param m:
    dataset: gdal.Dataset :

    m: dict :


    Returns
    -------

    """

    print('write_cloud')
    print(m)
    outfilename = m['outfilename']

    while os.path.exists(outfilename) and os.path.exists(outfilename + '0'):
        try:
            os.remove(outfilename)
        except:
            pass
        try:
            os.remove(outfilename + '0')
        except:
            pass

    crs = m['crs']

    print(m['z_up'], type(m['z_up']))

    d_min = np.nanmax(dataset[:, 2])
    d_max = np.nanmin(dataset[:, 2])
    print(d_min, d_max)
    # build CSAR bands
    bandInfo = {}  # Define our bands below
    z_dir = cc.Direction.HEIGHT
    layerName = "Elevation"
    print(m['z_up'], layerName, z_dir)
    bandInfo[layerName] = cc.BandInfo(type=cc.DataType.FLOAT64, tuple_length=1, name=layerName, direction=z_dir,
                                      units='m', category=cc.Category.SCALAR, ndv=-1.0, minimum=d_min, maximum=d_max)
    bandInfo['Position'] = cc.BandInfo(type=cc.DataType.FLOAT64, tuple_length=3, name='Position',
                                       direction=cc.Direction.NAP, units='', category=cc.Category.SCALAR,
                                       ndv=(-1.0, -1.0, 0.0))

    # set up the CSAR
    opts = cc.Options()
    opts.open_type = cc.OpenType.WRITE
    opts.position_band_name = 'Position'
    opts.band_info = bandInfo
    opts.extents = ((dataset[:, 0].min(), dataset[:, 1].min(), dataset[:, 2].min()),
                    (dataset[:, 0].max(), dataset[:, 1].max(), dataset[:, 2].max()))
    opts.wkt_cosys = crs

    # Create data for iterator
    if not m['z_up']:
        blocks = [{layerName: list(-dataset[:, 2]), 'Position': list(dataset)}]
    else:
        blocks = [{layerName: list(dataset[:, 2]), 'Position': list(dataset)}]

    opts.iterator = lambda: iter(blocks)
    # print(outfilename)
    pc = cc.Cloud(filename=outfilename, options=opts)

    del pc


def check_metadata(meta: dict, meta_type: str):
    """
    Check to make sure the required metadata keys are available in the
    provided metadata dictionary.

    Parameters
    ----------
    meta :
        param meta_type:
    meta: dict :

    meta_type: str :


    Returns
    -------

    """

    if meta_type == 'gdal':
        req_attrib = {'resx', 'resy', 'originx', 'originy', 'dimx', 'dimy', 'crs', 'nodata', 'outfilename', 'z_up'}
        mkeys = ''

        for key in req_attrib:
            if key not in meta:
                mkeys = mkeys + key + ', '

        if len(mkeys) > 0:
            raise ValueError('Metadata missing to write csar {}'.format(mkeys))
    elif meta_type == 'point':
        req_attrib = {'crs', 'outfilename', 'z_up'}
        mkeys = ''

        for key in req_attrib:
            if key not in meta:
                mkeys = mkeys + key + ', '

        if len(mkeys) > 0:
            raise ValueError('Metadata missing to write csar {}'.format(mkeys))
    else:
        raise ValueError('Unknown metadata typ: {}'.format(meta_type))


def main():
    """Parse the arguments and send them to the write method."""

    # check to make sure the file exists
    data = np.load(sys.argv[1])

    # check to make sure the metadata file exists
    with open(sys.argv[2], 'rb') as metafile:
        metadata = pickle.load(metafile)

    # write type argument
    outfile_type = sys.argv[3]

    # read the metadata into variables and send to the write method
    check_metadata(metadata, outfile_type)

    if outfile_type == 'gdal':
        write_raster(data, metadata)
    elif outfile_type == 'point':
        write_points(data, metadata)


if __name__ == '__main__':
    main()
