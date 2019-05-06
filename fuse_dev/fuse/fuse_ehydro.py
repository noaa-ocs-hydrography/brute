# -*- coding: utf-8 -*-
"""
fuse_ehydro.py

Created on Thu Jan 31 10:30:11 2019

@author: grice
"""

import os as _os
import logging as _logging
import fuse.fuse_base_class as _fbc
import fuse.meta_review.meta_review_ehydro as _mre
import fuse.raw_read.usace as _usace
import fuse.datum_transform.transform as _trans
import fuse.interpolator.interpolator as _interp
from fuse.proc_io.proc_io import proc_io

class fuse_ehydro(_fbc.fuse_base_class):
    _cols = ['from_filename',
            'from_path',
            'to_filename',
            'start_date',
            'end_date',
            'from_fips',
            'from_horiz_datum',
            'from_horiz_units',
            'from_horiz_unc',
            'to_horiz_datum',
            'from_vert_datum',
            'from_vert_key',
            'from_vert_units',
            'from_vert_unc',
            'to_vert_datum',
            'to_vert_units',
            'agency',
            'source_indicator',
            'source_type',
            'complete_coverage',
            'complete_bathymetry',
            'vert_uncert_fixed',
            'vert_uncert_vari',
            'horiz_uncert',
            'feat_size',
            'feat_detect',
            'feat_least_depth',
            'interpolated',
            'script_version',
            ]
    
    def __init__(self, config_filename):
        super().__init__(config_filename)
        self._meta_obj = _mre.meta_review_ehydro(self._config['metapath'],
                                                fuse_ehydro._cols)
        self._set_data_reader()
        self._set_data_transform()
        self._set_data_interpolator()
        self._set_data_writer()
        self._meta = {} # initialize the metadata holder
        self.logger = _logging.getLogger('fuse')
        self.logger.setLevel(_logging.DEBUG)
        
    def _set_data_reader(self):
        """
        Use information from the config file to set the reader to use for
        converting the raw data to usable metadata and bathymetry.
        """
        try:
            reader_type = self._config['raw_reader_type'].casefold()
            if reader_type == 'cenan':
                self._reader = _usace.cenan.read_raw()
            elif reader_type == 'cemvn':
                self._reader = _usace.cemvn.read_raw()
            elif reader_type == 'cesaj':
                self._reader = _usace.cesaj.read_raw()
            elif reader_type == 'cesam':
                self._reader = _usace.cesam.read_raw()
            elif reader_type == 'ceswg':
                self._reader = _usace.ceswg.read_raw()
            else:
                raise ValueError('reader type not implemented')
        except:
            raise ValueError("No reader type found in the configuration file.")
            
    def _set_data_transform(self):
        """
        Set up the datum transformation engine. 
        """
        self._transform = _trans.transform(self._config, self._reader)
        
    def _set_data_interpolator(self):
        """
        Set up the interpolator engine.
        """
        engine = self._config['interpolation_engine']
        res = float(self._config['to_resolution'])
        method = self._config['interpolation_method']
        self._interpolator = _interp.interpolator(engine, method, res)
        
    def _set_data_writer(self):
        """
        Set up the location and method to write tranformed and interpolated
        data.
        """
        ext = self._config['bathymetry_intermediate_file']
        self._writer = proc_io('gdal', ext)
        
    def read(self, infilename):
        """
        Extract metadata from the provided eHydro file path and write the metadata
        to the specified metadata file.  The bathymetry will be interpolated and
        writen to a CSAR file in the specificed csarpath.
        """
        self._meta = {}
        self._set_log(infilename)
        # get the metadata
        meta = self._reader.read_metadata(infilename)
        meta['to_horiz_datum'] = self._config['to_horiz_datum']
        meta['to_vert_datum'] = self._config['to_vert_datum']
        meta['to_vert_units'] = 'metres'
        meta['interpolated'] = 'True'
        #meta['script_version'] = meta['script_version'] #+ ',' + __version__ + i2c.__version__
        self._meta.update(meta)
        # write the metadata
        self._meta_obj.write_meta_record(meta)
        if 'from_fips' in self._meta:
            pass
            #self.process(infilename)
        # reading the bathymetry is not required > goes directly to datum trans
                
    def process(self, infilename):
        """
        Do the datum transformtion and interpolation.
        """
        self._get_stored_meta(infilename)
        self._set_log(infilename)
        if 'from_fips' in self._meta:
            # convert the bathy
            outpath = self._config['outpath']
            infilepath, infilebase = _os.path.split(infilename)
            infileroot, ext = _os.path.splitext(infilebase)
            outfilename = _os.path.join(outpath, infileroot)
            new_ext = self._config['bathymetry_intermediate_file']
            outfilename = outfilename + '.' + new_ext 
            # oddly _transform becomes the bathymetry reader here...
            # return a gdal dataset in the right datums for combine
            dataset = self._transform.translate(infilename, self._meta)
            # take a gdal dataset for interpolation and return a gdal dataset
            dataset = self._interpolator.interpolate(dataset)
            self._writer.write(dataset, outfilename)
            self._meta['to_filename'] = outfilename
        else:
            self.logger.log(_logging.DEBUG, 'No fips code found')
            
    def _set_log(self, infilename):
        """
        Set the object logging object and file.
        """
        metapath, metafile = _os.path.split(self._config['metapath'])
        filepath, filename = _os.path.split(infilename)
        fname, ext = _os.path.splitext(filename)
        logname = _os.path.join(metapath, fname + '.log')
        self._meta['logfilename'] = logname
        # remove handlers that might have existed from previous files
        for h in self.logger.handlers:
            self.logger.removeHandler(h)
        # create file handler for this filename
        fh = _logging.FileHandler(logname)
        fh.setLevel(_logging.DEBUG)
        self.logger.addHandler(fh)
        
            
    def _get_stored_meta(self, infilename):
        """
        Get the metadata in a local dictionary so that it can be used within
        the instantiated object.
        """
        # file name is the key rather than the path
        path, f = _os.path.split(infilename)
        if 'from_filename' not in self._meta:
            self._meta = self._meta_obj.read_meta_record(f)
        elif self._meta['from_filename'] is not infilename:
            self._meta = self._meta_obj.read_meta_record(f)
        # need to catch if this file is not in the metadata record yet here.
