# -*- coding: utf-8 -*-
"""
proc_io.py

Created on Thu Feb 14 15:13:52 2019

@author: grice
"""
import sys, os
import pickle
import subprocess
import numpy as np
import tables as tb
from osgeo import gdal
gdal.UseExceptions()


class proc_io:
    """
    A class to abstract the reading and writing of bathymetry.
    """
    def __init__(self, in_data_type, out_data_type, work_dir):
        """
        Initialize with the data type to be worked.
        """
        self._in_data_type = in_data_type
        self._out_data_type = out_data_type
        self._work_dir = work_dir
    
    def write(self, dataset, outfilename, z_up = True):
        """
        Write the provided data to the predefined data type.
        """
        if self._in_data_type =='gdal':
            data, metadata = self._convert_gdal(dataset)
        elif self._in_data_type =='bag':
            data, metadata = self._read_bag(dataset)
        else:
            raise ValueError('input data type unknown: ' + 
                             str(self._in_data_type))
        metadata['outfilename'] = outfilename
        metadata['z_up'] = z_up
        if self._out_data_type == 'csar':
            self._write_csar(dataset, metadata)
        else:
            raise ValueError('writer type unknown: ' + 
                             str(self._out_data_type))
            
    def _convert_gdal(self, dataset):
        """
        Convert the gdal dataset into a numpy array and a dictionary and
        return.
        """
        meta = {}
        # get the logisitics for converting the gdal dataset to csar
        gt = dataset.GetGeoTransform()
        meta['resx'] = gt[1]
        meta['resy'] = gt[5]
        meta['originx'] = gt[0]
        meta['originy'] = gt[3]
        meta['dimx'] = dataset.RasterXSize
        meta['dimy'] = dataset.RasterYSize
        meta['crs'] = dataset.GetProjection()
        rb = dataset.GetRasterBand(1) # should this be hardcoded for 1?
        meta['nodata'] = rb.GetNoDataValue()
        # get the gdal data raster
        data = dataset.ReadAsArray()
        return data, meta
    
    def _read_bag(self, dataset):
        meta = {}
        # get the logisitics for converting the bag dataset to csar
#        drv = ogr.GetDriverByName('BAG')
        gt = dataset.GetGeoTransform()
        meta['resx'] = gt[1]
        meta['resy'] = gt[5]
        meta['originx'] = gt[0]
        meta['originy'] = gt[3]
        meta['dimx'] = dataset.RasterXSize
        meta['dimy'] = dataset.RasterYSize
        meta['crs'] = dataset.GetProjection()
        rb = dataset.GetRasterBand(1) # should this be hardcoded for 1?
        meta['nodata'] = rb.GetNoDataValue()
        # get the gdal data raster
        data = dataset.ReadAsArray()
        return data, meta
            
    def _write_csar(self, dataset, metadata):
        """
        Convert the provided numpy array into a csar file using the provided
        metadata.
        
        The data and metadata are saved out to a file and then loaded into the
        wrapper around the csar writer.
        """
        # save the provided dataset and metadata to a file
        datafilename = os.path.join(self._work_dir, 'rasterdata')
        np.save(datafilename, dataset)
        metafilename = os.path.join(self._work_dir, 'metadata')
        with open(metafilename, 'wb') as metafile:
            pickle.dump(metadata, metafile)
        # set the locations for running the wrap_csar script
        start = os.path.realpath(os.path.dirname(__file__))
        write_csar = os.path.join(start, 'wrap_csar.py')
        activate_file = _retrieve_activate_batch()

        if os.path.exists(write_csar):
            args = ["cmd.exe", "/C", "set pythonpath=", "&&",  # run shell (/K: leave open (debugging), /C close the shell)
                    activate_file, "NBS35", "&&",  # activate the Caris 3.5 virtual environment
                    'python', write_csar,  # call the script
                    '"' + datafilename.replace("&", "^&") + '"',  # surface path
                    '"' + metafilename.replace("&", "^&") + '"',  # metadata path
                    ]

            subprocess.Popen(' '.join(args), creationflags=subprocess.CREATE_NEW_CONSOLE)
        else:
            print("Unable to create %s" % metadata['outfilename'])

# helper function to retrieve the path to the "Scripts" folder in PydroXL
def _retrieve_scripts_folder():
    install_prefix = sys.exec_prefix
    folder_path = os.path.realpath(os.path.join(install_prefix, os.pardir, os.pardir, "Scripts"))
    if not os.path.exists(folder_path):
        raise RuntimeError("The Scripts folder does not exist at: %s" % folder_path)
    return folder_path

# helper function to retrieve the path to the "activate.bat" batch file in PydroXL
def _retrieve_activate_batch():

    scripts_prefix = _retrieve_scripts_folder()
    file_path = os.path.realpath(os.path.join(scripts_prefix, "activate.bat"))
    if not os.path.exists(file_path):
        raise RuntimeError("The activate file does not exist at: %s" % file_path)
    return file_path