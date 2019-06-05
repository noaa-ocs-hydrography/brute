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

    print (name+'0')

    while os.path.exists(name) and os.path.exists(name+'0'):
        try:
            os.remove(name)
        except:
            pass
        try:
            os.remove(name+'0')
        except:
            pass
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

def write_cloud(dataset, m):
    """
    Convert a set of GDAL points to a CSAR point cloud.  The provided data is
    assumed to be a depth (positive down) and is assigned to a height
    (positive up).

    Parameters
    ----------
    dataset : numpy.array
        An array of [x,y,z] data points
    m : dict
        a dictionary containing the 'crs' info, 'outfilename', and 'z_up' to
        determine directionality of the data

    """
    print ('write_cloud')
    print (m)
    outfilename = m['outfilename']
    print (outfilename+'0')
    while os.path.exists(outfilename) and os.path.exists(outfilename+'0'):
        try:
            os.remove(outfilename)
        except:
            pass
        try:
            os.remove(outfilename+'0')
        except:
            pass

    crs = m['crs']

    print(m['z_up'], type(m['z_up']))

    # build CSAR bands
    bandInfo = {} # Define our bands below
    z_dir = cc.Direction.HEIGHT
    layerName = "Elevation"
    if not m['z_up']:
        z_dir = cc.Direction.DEPTH
        layerName = "Depth"
    print(m['z_up'], layerName, z_dir)
    bandInfo[layerName] = cc.BandInfo(type = cc.DataType.FLOAT64,
                                     tuple_length = 1,
                                     name = layerName,
                                     direction = z_dir,
                                     units = 'm',
                                     category = cc.Category.SCALAR,
                                     ndv = -1.0)
    bandInfo['Position'] = cc.BandInfo(type = cc.DataType.FLOAT64,
                                       tuple_length = 3,
                                       name = 'Position',
                                       direction = cc.Direction.NAP,
                                       units = '',
                                       category = cc.Category.SCALAR,
                                       ndv = (-1.0, -1.0, 0.0))
    # set up the CSAR
    opts = cc.Options();
    opts.open_type = cc.OpenType.WRITE
    opts.position_band_name = 'Position'
    opts.band_info = bandInfo
    opts.extents = ((dataset[:,0].min(),dataset[:,1].min(),dataset[:,2].min()),
                    (dataset[:,0].max(),dataset[:,1].max(),dataset[:,2].max()))
    opts.wkt_cosys = crs
    # Create data for iterator
    if not m['z_up']:
        blocks = [ { layerName: list(dataset[:,2]), 'Position': list(dataset) } ]
    else:
        x = dataset[:,0]
        y = dataset[:,1]
        z = -dataset[:,2]
        blocks = [ { layerName: list(z), 'Position': list(zip(x,y,[0]*len(z)))} ]


    opts.iterator = lambda: iter(blocks)
    # print(outfilename)
    pc = cc.Cloud(filename = outfilename, options=opts)

    pc = None

def check_metadata(meta, meta_type):
    """
    Check to make sure the required metadata keys are available in the
    provided metadata dictionary.
    """
    if meta_type == 'gdal':
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
    elif meta_type == 'point':
        req_attrib = {'crs',
                      'outfilename',
                      'z_up',
                      }
        mkeys = ''
        for key in req_attrib:
            if key not in meta:
                mkeys = mkeys + key + ', '
        if len(mkeys) > 0:
            raise ValueError('Metadata missing to write csar %s' % mkeys)
    else:
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
    # write type argument
    outfile_type = sys.argv[3]
    # read the metadata into variables and send to the write method
    check_metadata(metadata, outfile_type)
    if outfile_type == 'gdal':
        write_csar(data, metadata)
    elif outfile_type == 'point':
        write_cloud(data, metadata)


if __name__ == '__main__':
    main()
