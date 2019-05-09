# -*- coding: utf-8 -*-
"""
proc_io.py

Created on Thu Feb 14 15:13:52 2019

@author: grice
"""
import sys, os
import pickle
import subprocess
from tempfile import TemporaryDirectory as tempdir
import logging
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

__version__ = 'Test'

class proc_io:
    """
    A class to abstract the reading and writing of bathymetry.
    """
    def __init__(self, in_data_type, out_data_type, work_dir = None,
                 z_up = True, nodata=100000.0, caris_env_name = 'NBS35',
                 overwrite = True):
        """
        Initialize with the data type to be worked.
        """
        self._in_data_type = in_data_type
        self._out_data_type = out_data_type
        self._z_up = z_up
        self._write_nodata = nodata
        self._caris_environment_name = caris_env_name
        self.overwrite = overwrite
        if work_dir == None:
            self._work_dir = tempdir()
            self._work_dir_name = self._work_dir.name
        else:
            self._work_dir_name = work_dir
        self._logger = logging.getLogger('fuse')

    def write(self, dataset, outfilename, metadata = None):
        """
        Write the provided data to the predefined data type.
        """
        self._logger.log(logging.DEBUG, 
                         'Begin {} write'.format(self._out_data_type))
        if os.path.exists(outfilename) and self.overwrite:
            self._logger.log(logging.DEBUG, 'Overwriting ' + outfilename)
            os.remove(outfilename)
            if self._out_data_type == 'bag':
                caris_xml = outfilename +'.aux.xml'
                if os.path.exists(caris_xml):
                    os.remove(caris_xml)
        if self._out_data_type == 'csar':
            self._write_csar(dataset, outfilename)
        elif self._out_data_type == 'bag':
            self._write_bag(dataset, outfilename, metadata)
        else:
            raise ValueError('writer type unknown: ' +
                             str(self._out_data_type))
            
    def _write_csar(self, dataset, outfilename):
        """
        Convert the provided gdal dataset into a csar file.

        The data and metadata are saved out to a file and then loaded into the
        wrapper around the csar writer.
        """
        conda_env_name = self.caris_environment_name
        # put the provided data into the right form for the csar conversion.
        if self._in_data_type =='gdal':
            dataset = self._set_gdalndv(dataset)
            data, metadata = self._gdal2array(dataset)
        else:
            raise ValueError('input data type unknown: ' +
                             str(self._in_data_type))
        metadata['outfilename'] = outfilename
        metadata['z_up'] = self._z_up
        conda_env_path = _retrieve_env_path(conda_env_name)
        python_path = os.path.join(conda_env_path, 'python')
        # save the provided dataset and metadata to a file
        datafilename = os.path.join(self._work_dir_name, 'rasterdata.npy')
        np.save(datafilename, data)
        metafilename = os.path.join(self._work_dir_name, 'metadata')
        with open(metafilename, 'wb') as metafile:
            pickle.dump(metadata, metafile)
        # set the locations for running the wrap_csar script
        start = os.path.realpath(os.path.dirname(__file__))
        write_csar = os.path.join(start, 'wrap_csar.py')
        activate_file = _retrieve_activate_batch()
        logfilename = self._get_logfilename()
        if os.path.exists(write_csar):
            args = ["cmd.exe", "/C", "set pythonpath= &&", # setup the commandline
                    activate_file, conda_env_name, "&&",  # activate the Caris 3.5 virtual environment
                    python_path, write_csar,  # call the script
                    '"' + datafilename.replace("&", "^&") + '"',  # surface path
                    '"' + metafilename.replace("&", "^&") + '"',  # metadata path
                    '"' + logfilename.replace("&", "^&") + '"',
                    ]
            args = ' '.join(args)
            self._logger.log(logging.DEBUG, args)
            self._stop_logfile()
            try:
                proc = subprocess.Popen(args, creationflags=subprocess.CREATE_NEW_CONSOLE)
            except:
                err = 'Error executing: {}'.foramt(args)
                print(err)
                self._logger.log(logging.DEBUG, err)
            self._start_logfile(logfilename)
            try:
                stdout, stderr = proc.communicate()
                self._logger.log(logging.DEBUG, stdout)
                self._logger.log(logging.DEBUG, stderr)
            except:
                err = 'Error in handling error output'
                print(err)
                self._logger.log(logging.DEBUG, err)
            if not os.path.exists(metadata['outfilename']):
                err = "Unable to create {}".format(metadata['outfilename'])
                self._logger.log(logging.DEBUG, err)
                raise RuntimeError(err)
        else:
            err = "Unable to overwrite {}".format(metadata['outfilename'])
            self._logger.log(logging.DEBUG, err)
            raise RuntimeError(err)

    def _write_bag(self, dataset, outfilename, metadata = None):
        """
        Parameters
        ----------
        dataset : A georeferenced gdal raster object.
        
        outfilename : A string defining the path and filename of the file to be
            written.
        
        Kwarg
        -----
        metadata : A dictionary containing the metadata to be written to to
            the standard gdal bag xml.  Keys to the written must corrispond to
            the standard xml tags.  Defaults to None which writes nothing.
        """
        if self._in_data_type =='gdal':
            dataset = self._set_gdalndv(dataset)
        else:
            raise ValueError('input data type unknown: ' +
                             str(self._in_data_type))
        if metadata is not None:
            raise NotImplementedError('bag xml metadata write has not been implemented')
        # Prepare destination file
        driver = gdal.GetDriverByName("BAG")
        # write and close output raster dataset
        dest = driver.CreateCopy(outfilename, dataset)
        dest.GetRasterBand(1).SetColorTable(0) # what is this?
        dest = None
        self._logger.log(logging.DEBUG, 'BAG file created')

    def _gdal2array(self, dataset):
        """
        Convert the gdal dataset into a numpy array and a dictionary of
        metadata of the geotransform information and return.
        
        The gdal dataset should have he no data value set appropriately.
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
        meta['nodata'] = rb.GetNoDataValue()
        # get the gdal data raster
        data = np.flipud(rb.ReadAsArray())
        return data, meta
    
    def _set_gdalndv(self, dataset):
        """
        Update the gdal raster object no data value and the raster no data
        values in to corrispond with the object no data value.
        """
        # check the no data value
        rb = dataset.GetRasterBand(1) # should this be hardcoded for 1?
        ndv = rb.GetNoDataValue()
        if self._write_nodata != ndv:
            data = rb.ReadAsArray()
            data = np.where(data==ndv, self._write_nodata, data)
            dataset.GetRasterBand(1).WriteArray(data)
        return dataset

    def _get_logfilename(self):
        """
        Return the log filename.
        """
        if len(self._logger.handlers) > 1:
            raise ValueError('Not sure which hanlder to use for logging csar work. Using first')
        elif len(self._logger.handlers) <1:
            handlefilename = 'casiano.log'
        else:
            h = self._logger.handlers[0]
            handlefilename = h.baseFilename
        return handlefilename

    def _stop_logfile(self):
        """
        Get the logger filename, stop logging to it, and return the filename.
        """
        # remove handlers that might have existed from previous files
        if len(self._logger.handlers) <1:
            pass
        else:
            h = self._logger.handlers[0]
            self._logger.removeHandler(h)

    def _start_logfile(self, handlefilename):
        """
        Add a handler to the logger at the provided filename.
        """
        # create file handler for this filename
        fh = logging.FileHandler(handlefilename)
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        self._logger.addHandler(fh)

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

def _retrieve_env_path(env_name):
    """
    Given a conda environement name, find the environment.
    """
    current_env_loc = os.environ['conda_prefix']
    desired_env_loc = os.path.join(current_env_loc, os.pardir, env_name)
    if os.path.exists(desired_env_loc):
        return desired_env_loc
    else:
        raise ValueError('{} environment does not exist in current conda installation'.format(env_name))
