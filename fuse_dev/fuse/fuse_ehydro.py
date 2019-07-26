# -*- coding: utf-8 -*-
"""
fuse_ehydro.py

Created on Thu Jan 31 10:30:11 2019

@author: grice
"""

import logging as _logging
import os as _os

import fuse.datum_transform.transform as _trans
import fuse.fuse_base_class as _fbc
import fuse.interpolator.interpolator as _interp
import fuse.meta_review.meta_review_ehydro as _mre
import fuse.raw_read.usace as _usace
from fuse.proc_io.proc_io import proc_io


class fuse_ehydro(_fbc.fuse_base_class):
    """TODO write description"""
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
        self._db = None
        self._meta = {}  # initialize the metadata holder
        self._pickle_meta = {}  # initialize the survey pickle object
        self.logger = _logging.getLogger('fuse')
        self.logger.setLevel(_logging.DEBUG)

    def _set_data_reader(self):
        """
        Use information from the config file to set the reader to use for
        converting the raw data to usable metadata and bathymetry.

        Parameters
        ----------

        Returns
        -------

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
            elif reader_type == 'cespl':
                self._reader = _usace.cespl.read_raw()
            else:
                raise ValueError('reader type not implemented')
        except:
            raise ValueError("No reader type found in the configuration file.")

    def _set_data_transform(self):
        """Set up the datum transformation engine."""

        self._transform = _trans.transform(self._config, self._reader)

    def _set_data_interpolator(self):
        """Set up the interpolator engine."""

        engine = self._config['interpolation_engine']
        res = float(self._config['to_resolution'])
        method = self._config['interpolation_method']
        self._interpolator = _interp.interpolator(engine, method, res)

    def _set_data_writer(self):
        """
        Set up the location and method to write tranformed and interpolated
        data.

        Parameters
        ----------

        Returns
        -------

        """

        ext = self._config['bathymetry_intermediate_file']
        if ext == 'bag':
            ext2 = 'gpkg'
        else:
            ext2 = ext
        self._writer = proc_io('gdal', ext)
        self._points = proc_io('point', 'csar')

    def _read_pickle(self, infilename: str):
        """


        Parameters
        ----------
        infilename :


        Returns
        -------

        """

        pickle = _usace.parse_usace_pickle.pickle_file(infilename)
        return pickle.pickle_meta

    def read(self, infilename: str):
        """
        Extract metadata from the provided eHydro file path and write the metadata
        to the specified metadata file.  The bathymetry will be interpolated and
        writen to a CSAR file in the specificed csarpath.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        self._meta = {}
        self._set_log(infilename)
        print('calling pickle reader at fuse_ehydro level to pass polygon')
        self._pickle_meta = self._read_pickle('infilename')
        # get the metadata
        meta = self._reader.read_metadata(infilename)
        meta['to_horiz_datum'] = self._config['to_horiz_datum']
        meta['to_vert_datum'] = self._config['to_vert_datum']
        meta['to_vert_units'] = 'metres'
        meta['interpolated'] = 'False'
        # meta['script_version'] = f'{meta['script_version']}' #,{__version__}{i2c.__version__}'
        self._meta.update(meta)
        # write the metadata
        self._meta_obj.write_meta_record(meta)

    def process(self, infilename: str, interpolate=True):
        """
        Do the datum transformtion and interpolation.

        Given the generic need to interpolate USACE data the 'interpolate'
        kwarg is set to True as a hack.  This information should be drawn from
        the data reader since there will be cases where we get full res data
        from the reader and interlation is not necessary.

        Parameters
        ----------
        infilename
        interpolate

        Returns
        -------

        """

        self._get_stored_meta(infilename)
        self._set_log(infilename)
        if 'from_fips' in self._meta:
            # convert the bathy for the original data
            outpath = self._config['outpath']
            infilepath, infilebase = _os.path.split(infilename)
            infileroot, ext = _os.path.splitext(infilebase)
            outfilebase = _os.path.join(outpath, infileroot)
            new_ext = self._config['bathymetry_intermediate_file']
            outfilename = f'{outfilebase}.{new_ext}'
            # oddly _transform becomes the bathymetry reader here...
            # return a gdal dataset in the right datums for combine
            dataset = self._transform.translate(infilename, self._meta)
            resolution = self._config['to_resolution']
            self._points.write(dataset, outfilename)
            self._meta['to_filename'] = outfilename
            self._meta_obj.write_meta_record(self._meta)
            # take a gdal dataset for interpolation and return a gdal dataset
            interpfilename = f'{outfilebase}_{resolution}m_interp.{new_ext}'
            interpkeyfilename = f'{infilebase}.interpolated'
            self._meta_interp = self._meta.copy()
            self._meta_interp['interpolated'] = True
            self._meta_interp['from_filename'] = interpkeyfilename
            self._meta_interp['to_filename'] = interpfilename
            if 'poly_name' in self._pickle_meta:
                shapename = self._pickle_meta['poly_name']
                shapepath = _os.path.join(infilepath, shapename)
                dataset = self._interpolator.interpolate(dataset, shapepath)
            else:
                dataset = self._interpolator.interpolate(dataset)
            self._writer.write(dataset, interpfilename)
            self._meta_obj.write_meta_record(self._meta_interp)
        else:
            self.logger.log(_logging.DEBUG, 'No fips code found')

    def post(self, infilename):
        """
        Make the data available for amalgamation.
        """
        self._set_log(infilename)
        self._get_s57_stored_meta(infilename)
        if len(self._meta) > 0:
            if self._db is None:
                self._connect_to_db()
            procfile = self._meta['to_filename']
            print(self._s57_meta)
            self._db.write(procfile, 'new', self._s57_meta)

    def _connect_to_db(self):
        """
        Connect to the database defined in the configuration dictionary.
        """
        if 'database_location' in self._config:
            db_loc = self._config['database_location']
        else:
            raise ValueError('No database location defined in the configuration file.')
        if 'database_location' in self._config:
            db_name = self._config['database_name']
        else:
            raise ValueError('No database name defined in the configuration file.')
        intype = self._config['bathymetry_intermediate_file']
        self._db = proc_io(intype, 'carisbdb51', db_loc=db_loc, db_name=db_name)

    def disconnect(self):
        """
        Asks proc_io to close the database connection

        """
        if self._db:
            self._db.close_connection()
        else:
            if 'database_location' in self._config:
                db_loc = self._config['database_location']
                raise RuntimeError(f'No connection to {db_loc} to close')
            else:
                raise ValueError('No database location defined in the configuration file.')

    def _set_log(self, infilename: str):
        """
        Set the object logging object and file.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        metapath, metafile = _os.path.split(self._config['metapath'])
        filepath, filename = _os.path.split(infilename)
        fname, ext = _os.path.splitext(filename)
        logname = _os.path.join(metapath, f'{fname}.log')
        self._meta['logfilename'] = logname
        # remove handlers that might have existed from previous files
        for h in self.logger.handlers:
            self.logger.removeHandler(h)
        # create file handler for this filename
        fh = _logging.FileHandler(logname)
        fh.setLevel(_logging.DEBUG)
        formatter = _logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

    def _get_stored_meta(self, infilename: str):
        """
        Get the metadata in a local dictionary so that it can be used within
        the instantiated object.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        # file name is the key rather than the path
        path, f = _os.path.split(infilename)
        if 'from_filename' not in self._meta:
            self._meta = self._meta_obj.read_meta_record(f)
        elif self._meta['from_filename'] is not infilename:
            self._meta = self._meta_obj.read_meta_record(f)
        # need to catch if this file is not in the metadata record yet here.

    def _get_s57_stored_meta(self, infilename: str):
        """
        Get the metadata in a local dictionary so that it can be used within
        the instantiated object.  In this case the metadata is converted to
        an s57 version of the metadata.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        # file name is the key rather than the path
        path, f = _os.path.split(infilename)
        if 'from_filename' not in self._meta:
            self._meta = self._meta_obj.read_meta_record(f)
        elif self._meta['from_filename'] is not infilename:
            self._meta = self._meta_obj.read_meta_record(f)
        if len(self._meta) > 0:
            self._s57_meta = self._meta_obj.row2s57(self._meta)
        # need to catch if this file is not in the metadata record yet here.
