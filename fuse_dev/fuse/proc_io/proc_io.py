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
from osgeo import gdal, osr
gdal.UseExceptions()

__version__ = 'Test'

def maxValue(arr):
    '''Takes an input array and finds the most used value in the array, this 
    value is used by the program to assume the array's nodata value
    
    returns the most used value in the array as an integer
    '''
    print ('maxValue')
    nums, counts = np.unique(arr, return_counts =True)
    index = np.where(counts==np.amax(counts))
    print (index, nums[index])
    return int(nums[index])

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
            data, metadata = self._convert_gdal(dataset)
        else:
            raise ValueError('input data type unknown: ' + 
                             str(self._in_data_type))
        metadata['outfilename'] = outfilename
        metadata['z_up'] = z_up
        if self._out_data_type == 'csar':
            self._write_csar(data, metadata)
        elif self._out_data_type == 'bag':
            self._write_bag(data, metadata)
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
        print (gt)
        meta['resx'] = gt[1]
        meta['resy'] = gt[1]
        meta['originx'] = gt[0]
        meta['originy'] = gt[3] - dataset.RasterYSize
        meta['dimx'] = dataset.RasterXSize
        meta['dimy'] = dataset.RasterYSize
        print (meta)
        meta['crs'] = dataset.GetProjection()
        rb = dataset.GetRasterBand(1) # should this be hardcoded for 1?
#        meta['nodata'] = rb.GetNoDataValue()
        # get the gdal data raster
        data = np.flipud(rb.ReadAsArray())
        maxVal = maxValue(data)
        meta['nodata'] = maxVal
        return data, meta
            
    def _write_csar(self, data, metadata):
        """
        Convert the provided numpy array into a csar file using the provided
        metadata.
        
        The data and metadata are saved out to a file and then loaded into the
        wrapper around the csar writer.
        """
        # save the provided dataset and metadata to a file
        if os.path.exists(self._work_dir):
            pass
        else:
            os.mkdir(self._work_dir)
        datafilename = os.path.join(self._work_dir, 'rasterdata.npy')
        np.save(datafilename, data)
        metafilename = os.path.join(self._work_dir, 'metadata')
        with open(metafilename, 'wb') as metafile:
            pickle.dump(metadata, metafile)
        # set the locations for running the wrap_csar script
        start = os.path.realpath(os.path.dirname(__file__))
        write_csar = os.path.join(start, 'wrap_csar.py')
        activate_file = _retrieve_activate_batch()

        if os.path.exists(write_csar):
            args = ["cmd.exe", "/K", "set pythonpath=", "&&",  # run shell (/K: leave open (debugging), /C close the shell)
                    activate_file, "NBS35", "&&",  # activate the Caris 3.5 virtual environment
                    'python', write_csar,  # call the script
                    '"' + datafilename.replace("&", "^&") + '"',  # surface path
                    '"' + metafilename.replace("&", "^&") + '"',  # metadata path
                    ]

            subprocess.Popen(' '.join(args), creationflags=subprocess.CREATE_NEW_CONSOLE)
        else:
            print("Unable to create %s" % metadata['outfilename'])
            
    def _write_bag(data, metadata, dtype=gdal.GDT_UInt32,
                 options=0, color_table=0, nbands=1, nodata=False):
        pass
#        """Directly From:
#        "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
#        (https://gis.stackexchange.com/a/278965)]
#    
#        Parameters
#        ----------
#        raster_array : numpy.array
#            Array to be written to a GeoTiff file
#        gt : tuple, gdal.GeoTransform
#            Norhtern extent, resolution, 0.0, Western extent, 0.0, -resolution)
#        data_obj : gdal.RasterBand
#            gdal.RasterBand
#        outputpath : string
#            Folder to save the GeoTiff raster
#    
#        """
#        print('write_raster')
#    
#        height, width = data.shape
#    
#        # Prepare destination file
#        driver = gdal.GetDriverByName("BAG")
#        if options != 0:
#            dest = driver.Create(metadata['outfilename'], width, height, nbands, dtype, options)
#        else:
#            dest = driver.Create(metadata['outfilename'], width, height, nbands, dtype)
#    
#        # Write output raster
#        if color_table != 0:
#            dest.GetRasterBand(1).SetColorTable(color_table)
#    
#        dest.GetRasterBand(1).WriteArray(data)
#    
#        if nodata is not False:
#            dest.GetRasterBand(1).SetNoDataValue(nodata)
#    
#        # Set transform and projection
#        dest.SetGeoTransform(gt)
#        wkt = data_obj.GetProjection()
#        srs = osr.SpatialReference()
#        srs.ImportFromWkt(wkt)
#        dest.SetProjection(srs.ExportToWkt())
#    
#        # Close output raster dataset
#        dest = None

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