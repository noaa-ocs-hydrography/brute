# -*- coding: utf-8 -*-
"""
fuse_ehydro.py

Created on Thu Jan 31 10:30:11 2019

@author: grice
"""

import logging as _logging
import os as _os

import fuse.datum_transform.transform as _trans
import fuse.fuse_processor as _fbc
import fuse.interpolator.interpolator as _interp
import fuse.meta_review.meta_review_ehydro as _mre
import fuse.raw_read.usace as _usace
from fuse.proc_io.proc_io import ProcIO
from fuse import score


class FuseProcessor_eHydro(_fbc.FuseProcessor):
    """TODO write description"""
    _default_quality_metrics = {'complete_coverage': False,
                                'bathymetry': True,
                                'vert_uncert_fixed': 0.5,
                                'vert_uncert_vari': 0.01,
                                'horiz_uncert_fixed': 5.0,
                                'horiz_uncert_vari': 0.05,
                                'feat_detect': False,
                                }

    _datums = ['from_fips',
              'from_horiz_datum',
              'from_horiz_units',
              'to_horiz_datum',
              'from_vert_datum',
              'from_vert_key',
              'from_vert_units',
              'to_vert_datum',
              'to_vert_units',
              ]

    _quality_metrics = ['from_horiz_unc',
                        'from_vert_unc',
                        'complete_coverage',
                        'bathymetry',
                        'vert_uncert_fixed',
                        'vert_uncert_vari',
                        'horiz_uncert_fixed',
                        'horiz_uncert_vari',
                        'feat_size',
                        'feat_detect',
                        'feat_least_depth',
                        'interpolate',
                        ]

    _paths = ['from_filename',
             'from_path',
             'to_filename',
             ]

    _dates = ['start_date',
              'end_date',
              ]

    _source_info = ['agency',
                    'source_indicator',
                    'source_type',
                    'interpolated',
                    'posted',
                    'license',
                    ]

    _processing_info = ['logfilename',
                        'version_reference',
                        ]

    _scores = ['catzoc',
               'supersession_score',
               ]

    def __init__(self, config_filename):
        super().__init__(config_filename)
        cols = FuseProcessor_eHydro._paths \
            + FuseProcessor_eHydro._dates \
            + FuseProcessor_eHydro._datums \
            + FuseProcessor_eHydro._quality_metrics \
            + FuseProcessor_eHydro._scores \
            + FuseProcessor_eHydro._source_info
        self._meta_obj = _mre.MetaReviewer_eHydro(self._config['metapath'], cols)
        self._set_data_reader()
        self._set_data_transform()
        self._set_data_interpolator()
        self._set_data_writer()
        self._db = None
        self._meta = {}  # initialize the metadata holder
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
                self._reader = _usace.cenan.CENANRawReader()
            elif reader_type == 'cemvn':
                self._reader = _usace.cemvn.CEMVNRawReader()
            elif reader_type == 'cesaj':
                self._reader = _usace.cesaj.CESAJRawReader()
            elif reader_type == 'cesam':
                self._reader = _usace.cesam.CESAMRawReader()
            elif reader_type == 'ceswg':
                self._reader = _usace.ceswg.CESWGRawReader()
            elif reader_type == 'cespl':
                self._reader = _usace.cespl.CESPLRawReader()
            elif reader_type == 'cenae':
                self._reader = _usace.cenae.CENAERawReader()
            else:
                raise ValueError('reader type not implemented')
        except:
            raise ValueError("No reader type found in the configuration file.")

    def _set_data_transform(self):
        """Set up the datum transformation engine."""

        self._transform = _trans.DatumTransformer(self._config, self._reader)

    def _set_data_interpolator(self):
        """Set up the interpolator engine."""

        engine = self._config['interpolation_engine']
        res = float(self._config['to_resolution'])
        method = self._config['interpolation_method']
        self._interpolator = _interp.Interpolator(engine, method, res)

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
        self._writer = ProcIO('gdal', ext)
        self._points = ProcIO('point', ext2)

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

        self._set_log(infilename)
        # get the metadata
        meta = self._reader.read_metadata(infilename)
        meta['to_horiz_datum'] = self._config['to_horiz_datum']
        meta['to_vert_datum'] = self._config['to_vert_datum']
        meta['to_vert_units'] = 'metres'
        meta['interpolated'] = 'False'
        meta['posted'] = False
        if not self._quality_metadata_ready(meta):
            default = FuseProcessor_eHydro._default_quality_metrics
            msg = f'Not all quality metadata was found.  Using default values: {default}'
            self.logger.log(_logging.DEBUG, msg)
            meta = {**default, **meta}
        # write the metadata
        self._meta_obj.write_meta_record(meta)
        self._close_log()

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

        TODO: need to add checks to make sure the metadata is ready.
            Perhaps this should be added to the metadata object?
        """

        metadata = self._get_stored_meta(infilename)
        self._set_log(infilename)
        if 'from_fips' in metadata:
            # convert the bathy for the original data
            outpath = self._config['outpath']
            infilepath, infilebase = _os.path.split(infilename)
            infileroot, ext = _os.path.splitext(infilebase)
            outfilebase = _os.path.join(outpath, infileroot)
            new_ext = self._config['bathymetry_intermediate_file']
            outfilename = f'{outfilebase}.{new_ext}'
            # oddly _transform becomes the bathymetry reader here...
            # return a gdal dataset in the right datums for combine
            dataset = self._transform.translate(infilename, metadata)
            resolution = self._config['to_resolution']
            self._points.write(dataset, outfilename)
            metadata['to_filename'] = outfilename
            self._meta_obj.write_meta_record(metadata)
            if 'interpolate' in metadata:
                interpolate = metadata['interpolate']
                if interpolate == 'True':
                    # take a gdal dataset for interpolation and return a gdal dataset
                    interpfilename = f'{outfilebase}_{resolution}m_interp.{new_ext}'
                    interpkeyfilename = f'{infilebase}.interpolated'
                    meta_interp = metadata.copy()
                    meta_interp['interpolated'] = True
                    meta_interp['from_filename'] = interpkeyfilename
                    meta_interp['to_filename'] = interpfilename
                    if 'poly_name' in meta_interp:
                        shapename = meta_interp['poly_name']
                        shapepath = _os.path.join(infilepath, shapename)
                        dataset = self._interpolator.interpolate(dataset, shapepath)
                    else:
                        dataset = self._interpolator.interpolate(dataset)
                    self._writer.write(dataset, interpfilename)
                    self._meta_obj.write_meta_record(meta_interp)
                elif interpolate == 'False':
                    print(f'{infileroot} - No interpolation required')
            else:
                raise ValueError('metadata has no >interpolate< value')
        else:
            self.logger.log(_logging.DEBUG, 'No fips code found')
        self._close_log()

    def post(self, infilename):
        """
        Make the data available for amalgamation.

        TODO: need to add checks to make sure the metadata is ready.
            Perhaps this should be added to the metadata object?
        """
        metadata = self._get_stored_meta(infilename)
        self._set_log(infilename)
        if self._quality_metadata_ready(metadata):
            if not self._score_metadata_ready(metadata):
                metadata['catzoc'] = score.catzoc(metadata)
                metadata['supersession_score'] = score.supersession(metadata)
            if self._db is None:
                self._connect_to_db()
            procfile = metadata['to_filename']
            metadata['posted'] = True
            s57_meta = self._get_meta_as_s57(metadata)
            self._db.write(procfile, 'new', s57_meta)
            # need to check for proper insertion...
            self._meta_obj.write_meta_record(metadata)
        self._close_log()

    def score(self, infilename, date):
        """
        Provided a date, get the decayed quality metric and insert in the
        database, making the information available for amalgamation.
        """
        metadata = self._get_stored_meta(infilename)
        self._set_log(infilename)
        if metadata['posted'].upper() == 'TRUE':
            dscore = score.decay(metadata, date)
            if self._db == None:
                self._connect_to_db()
            procfile = metadata['to_filename']
            s57_meta = self._get_meta_as_s57(metadata)
            s57_meta['dcyscr'] = dscore
            self._db.write(procfile, 'metadata', s57_meta)
            log = f'Posting new decay score of {dscore} to database.'

        else:
            log = 'Insertion of decay score failed.'
        self.logger.log(_logging.DEBUG, log)
        self._close_log()


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
        self._db = ProcIO(intype, 'carisbdb51', db_loc=db_loc, db_name=db_name)

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
        root, ext = _os.path.splitext(infilename)
        if ext == '.interpolated':
            infilename = root
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

    def _close_log(self):
        """
        Close the object logging file.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """
        # remove handlers
        for h in self.logger.handlers:
            self.logger.removeHandler(h)

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
        return self._meta

    def _get_meta_as_s57(self, metadata):
        """
        The metadata is converted to an s57 version of the metadata.

        Parameters
        ----------
        metadata


        Returns
        -------

        """
        s57_meta = self._meta_obj.row2s57(metadata)
        return s57_meta

    def _quality_metadata_ready(self, metadata):
        """
        Check the metadata to see if the required fields are populated.
        """
        # check the feature metadata
        if 'feat_detect' in metadata:
            if metadata['feat_detect'] == True:
                if 'feat_least_depth' in metadata and 'feat_size' in metadata:
                    feature_ready = True
                else:
                    feature_ready = False
            else:
                feature_ready = True
        else:
            feature_ready = False
        if not feature_ready:
            msg = 'Quality metadata for features is not yet available.'
            self.logger.log(_logging.DEBUG, msg)
        # check the uncertainty metadata
        if 'vert_uncert_fixed' in metadata and 'vert_uncert_vari' in metadata:
            vert_uncert_ready = True
        else:
            vert_uncert_ready = False
        if 'horiz_uncert_fixed' in metadata and 'horiz_uncert_vari' in metadata:
            horiz_uncert_ready = True
        else:
            horiz_uncert_ready = False
        if not vert_uncert_ready or not horiz_uncert_ready:
            msg = 'Quality metadata for uncertainty is not yet available.'
            self.logger.log(_logging.DEBUG, msg)
        # check the coverage
        if 'complete_coverage' in metadata and 'bathymetry' in metadata:
            coverage_ready = True
        else:
            coverage_ready = False
        ready = feature_ready and vert_uncert_ready and horiz_uncert_ready and coverage_ready
        return ready

    def _date_metadata_ready(self, metadata):
        """
        Check the metadata to see if the required fields are populated.
        """
        if 'end_date' not in metadata or 'start_date' not in metadata:
            ready = False
        else:
            ready = True
        return ready

    def _score_metadata_ready(self, metadata):
        if 'catzoc' in metadata and 'supersession_score' in metadata:
            return True
        else:
            return False
