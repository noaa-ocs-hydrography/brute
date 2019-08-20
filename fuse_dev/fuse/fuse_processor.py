# -*- coding: utf-8 -*-
"""
fuse_abstract_obj

Created on Thu Jan 31 10:03:30 2019

@author: grice
"""

import configparser as _cp
import os as _os
import logging as _logging
import fuse.raw_read.usace as _usace
import fuse.raw_read.noaa as _noaa
import fuse.datum_transform.transform as _trans
import fuse.interpolator.interpolator as _interp
import fuse.meta_review as _mr
from fuse.proc_io.proc_io import ProcIO
from fuse import score

_ehydro_quality_metrics = {'complete_coverage': False,
                            'bathymetry': True,
                            'vert_uncert_fixed': 0.5,
                            'vert_uncert_vari': 0.1,
                            'horiz_uncert_fixed': 5.0,
                            'horiz_uncert_vari': 0.05,
                            'feat_detect': False,
                            }

class FuseProcessor:
    """The fuse object."""

    _datums = ['from_horiz_datum',
              'from_horiz_frame',
              'from_horiz_type',
              'from_horiz_units',
              'from_horiz_key',
              'from_vert_datum',
              'from_vert_key',
              'from_vert_units',
              'from_vert_direction',
              'to_horiz_frame',
              'to_horiz_type',
              'to_horiz_units',
              'to_horiz_key',
              'to_vert_datum',
              'to_vert_key',
              'to_vert_units',
              'to_vert_direction',
              ]

    _quality_metrics = ['from_horiz_unc',
                        'from_horiz_resolution',
                        'from_vert_unc',
                        'complete_coverage',
                        'bathymetry',
                        'vert_uncert_fixed',
                        'vert_uncert_vari',
                        'horiz_uncert_fixed',
                        'horiz_uncert_vari',
                        'to_horiz_resolution',
                        'feat_size',
                        'feat_detect',
                        'feat_least_depth',
                        ]

    _paths = ['from_filename',
              'from_path',
              'to_filename',
              'support_files',
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
                        'interpolate',
                        'file_size',
                        ]

    _scores = ['catzoc',
               'supersession_score',
               ]

    def __init__(self, configfilename: str = 'generic.config'):
        """
        Initialize with the metadata file to use and the horizontal and
        vertical datums of the workflow.

        Parameters
        ----------
        configfilename
            test
        """
        self._configfilename = configfilename
        self._config = self._read_configfile(configfilename)
        self.rawdata_path = self._config['rawpaths']
        self.procdata_path = self._config['outpath']
        self._cols = self._paths \
            + self._dates \
            + self._datums \
            + self._quality_metrics \
            + self._scores \
            + self._source_info \
            + self._processing_info
        self._meta_obj = _mr.MetaReviewer(self._config['metapath'], self._cols)
        self._set_data_reader()
        self._set_data_transform()
        self._set_data_interpolator()
        self._set_data_writer()
        self._db = None
        self._meta = {}  # initialize the metadata holder
        self.logger = _logging.getLogger('fuse')
        self.logger.setLevel(_logging.DEBUG)

    def _read_configfile(self, confile: str):
        """
        Read, parse, and return the configuration information in the provided
        file.  The actual format of this file is ....

        rawpaths
        outpath
        to_horiz_datum
        to_vert_datum
        metapath

        Parameters
        ----------
        confile :

        confile: str :


        Returns
        -------

        """

        config = {}
        config_file = _cp.ConfigParser()
        config_file.read(confile)
        sections = config_file.sections()
        for section in sections:
            for key in config_file[section]:
                if key == 'rawpaths':
                    rawpaths = []
                    raw = config_file[section][key].split(';')
                    for r in raw:
                        if _os.path.isdir(r):
                            rawpaths.append(r)
                        else:
                            raise ValueError(f'Invalid input path: {r}')
                    config[key] = rawpaths
                elif key == 'outpath':
                    if _os.path.isdir(config_file[section][key]):
                        config[key] = config_file[section][key]
                    else:
                        raise ValueError(f'Invalid input path: {config_file[section][key]}')
                else:
                    config[key] = config_file[section][key]

        if len(config) == 0:
            raise ValueError('Failed to read configuration file.')
        else:
            self._check_config(config)
        return config

    def _check_config(self, config_dict: dict):
        """
        Check to ensure at least the basic elements of the configuration are
        available.

        Parameters
        ----------
        config_dict :

        config_dict: dict :


        Returns
        -------

        """
        options = {'rawpaths': 'path to raw data',
                   'outpath': 'path to output data',
                   'to_horiz_datum' : 'output horizontal datum description',
                   'to_horiz_frame' : 'output horizontal datum frame',
                   'to_horiz_type' : 'output horizontal datum type',
                   'to_horiz_units' : 'output horizontal datum units',
                   'to_vert_key' : 'output vertical datum key',
                   'to_vert_units' : 'output vertical datum units',
                   'to_vert_direction' : 'output vertical datum direction',
                   'to_vert_datum': 'output vertical datum description',
                   'metapath': 'metadata output',
                   }
        for key in options.keys():
            if key not in config_dict:
                raise ValueError(f'No {options[key]} found in configuration file.')

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
                self._read_type = 'ehydro'
            elif reader_type == 'cemvn':
                self._reader = _usace.cemvn.CEMVNRawReader()
                self._read_type = 'ehydro'
            elif reader_type == 'cesaj':
                self._reader = _usace.cesaj.CESAJRawReader()
                self._read_type = 'ehydro'
            elif reader_type == 'cesam':
                self._reader = _usace.cesam.CESAMRawReader()
                self._read_type = 'ehydro'
            elif reader_type == 'ceswg':
                self._reader = _usace.ceswg.CESWGRawReader()
                self._read_type = 'ehydro'
            elif reader_type == 'cespl':
                self._reader = _usace.cespl.CESPLRawReader()
                self._read_type = 'ehydro'
            elif reader_type == 'cenae':
                self._reader = _usace.cenae.CENAERawReader()
                self._read_type = 'ehydro'
            elif reader_type == 'bag':
                self._reader = _noaa.bag.BAGRawReader()
                self._read_type = 'bag'
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
        The file name to use to read the bathymetry and metadata into useable
        forms.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        if self._read_type == 'ehydro':
            self._read_ehydro(infilename)
        elif self._read_type == 'bag':
            self._read_noaa_bag(infilename)
        else:
            raise ValueError('Reader type not implemented')

    def _read_ehydro(self, infilename: str):
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
        raw_meta = self._reader.read_metadata(infilename)
        meta = raw_meta.copy()
        meta['read_type'] = 'ehydro'

        # translate from the reader to common metadata keys for datum transformations
        if 'from_fips' in meta:
            meta['from_horiz_key'] = meta['from_fips']
        if 'from_horiz_units' in meta:
            if meta['from_horiz_units'].upper() == 'US SURVEY FOOT':
                meta['from_horiz_units'] = 'us_ft'
            else:
                raise ValueError(f'Input datum units are unknown: {meta["from_horiz_units"]}')
        if 'from_vert_key' in meta:
            meta['from_vert_key'] = meta['from_vert_key'].lower()
        if 'from_vert_units' in meta:
            if meta['from_vert_units'].upper() == 'US SURVEY FOOT':
                meta['from_vert_units'] = 'us_ft'
            else:
                raise ValueError(f'Input datum units are unknown: {meta["from_vert_units"]}')
        # insert a few default values for datum stuff if it isn't there already
        if 'from_vert_direction' not in meta:
            meta['from_vert_direction'] = 'height'
        if 'from_horiz_frame' not in meta:
            meta['from_horiz_frame'] = 'NAD83'
        if 'from_horiz_type' not in meta:
            meta['from_horiz_type'] = 'spc'
        # get the rest from the config file
        meta['to_horiz_frame'] = self._config['to_horiz_frame']
        meta['to_horiz_type'] = self._config['to_horiz_type']
        meta['to_horiz_units'] = self._config['to_horiz_units']
        if 'to_horiz_key' in self._config:
            meta['to_horiz_key'] = self._config['to_horiz_key']
        meta['to_vert_key'] = self._config['to_vert_key']
        meta['to_vert_units'] = self._config['to_vert_units']
        meta['to_vert_direction'] = self._config['to_vert_direction']
        meta['to_vert_datum'] = self._config['to_vert_datum']
        meta['interpolated'] = 'False'
        meta['posted'] = False
        if not self._quality_metadata_ready(meta):
            default = _ehydro_quality_metrics
            msg = f'Not all quality metadata was found.  Using default values: {default}'
            self.logger.log(_logging.DEBUG, msg)
            meta = {**default, **meta}
        # write the metadata
        self._meta_obj.write_meta_record(meta)
        self._close_log()

    def _read_noaa_bag(self, infilename: str):
        """
        Extract metadata from the provided bag file path and write the metadata
        to the specified metadata file.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        self._set_log(infilename)
        # get the metadata
        raw_meta = self._reader.read_metadata(infilename)
        meta = raw_meta.copy()
        meta['read_type'] = 'bag'
        # translate from the reader to common metadata keys for datum transformations
#        if 'from_fips' in meta:
#            meta['from_horiz_key'] = meta['from_fips']
#        if 'from_horiz_units' in meta:
#            if meta['from_horiz_units'].upper() == 'US SURVEY FOOT':
#                meta['from_horiz_units'] = 'us_ft'
#            else:
#                raise ValueError(f'Input datum units are unknown: {meta["from_horiz_units"]}')
#        if 'from_vert_key' in meta:
#            meta['from_vert_key'] = meta['from_vert_key'].lower()
#        if 'from_vert_units' in meta:
#            if meta['from_vert_units'].upper() == 'US SURVEY FOOT':
#                meta['from_vert_units'] = 'us_ft'
#            else:
#                raise ValueError(f'Input datum units are unknown: {meta["from_vert_units"]}')
#        # insert a few default values for datum stuff if it isn't there already
#        if 'from_vert_direction' not in meta:
#            meta['from_vert_direction'] = 'height'
#        if 'from_horiz_frame' not in meta:
#            meta['from_horiz_frame'] = 'NAD83'
#        if 'from_horiz_type' not in meta:
#            meta['from_horiz_type'] = 'spc'
        # get the rest from the config file
        meta['to_horiz_frame'] = self._config['to_horiz_frame']
        meta['to_horiz_type'] = self._config['to_horiz_type']
        meta['to_horiz_units'] = self._config['to_horiz_units']
        if 'to_horiz_key' in self._config:
            meta['to_horiz_key'] = self._config['to_horiz_key']
        meta['to_vert_key'] = self._config['to_vert_key']
        meta['to_vert_units'] = self._config['to_vert_units']
        meta['to_vert_direction'] = self._config['to_vert_direction']
        meta['to_vert_datum'] = self._config['to_vert_datum']
        meta['interpolated'] = 'False'
        meta['posted'] = False
        if not self._quality_metadata_ready(meta):
            default = _ehydro_quality_metrics
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
        metadata['read_type'] = self._read_type
        self._set_log(infilename)
        if 'from_fips' in metadata or 'from_horiz_frame' in metadata:
            # convert the bathy for the original data
            outpath = self._config['outpath']
            infilepath, infilebase = _os.path.split(infilename)
            infileroot, ext = _os.path.splitext(infilebase)
            outfilebase = _os.path.join(outpath, infileroot)
            metadata['new_ext'] = self._config['bathymetry_intermediate_file']
            # oddly _transform becomes the bathymetry reader here...
            # return a gdal dataset in the right datums for combine
            dataset, transformed = self._transform.translate(infilename, metadata)
            if self._read_type == 'ehydro':
                outfilename = f"{outfilebase}.{metadata['new_ext']}"
                self._points.write(dataset, outfilename)
                metadata['to_filename'] = outfilename
            if self._read_type == 'bag':
                metadata['to_filename'] = infilename
            self._meta_obj.write_meta_record(metadata)
            if 'interpolate' in metadata:
                interpolate = metadata['interpolate']
                if interpolate == 'True':
                    meta_interp = metadata.copy()
                    dataset, meta_interp = self._interpolator.interpolate(dataset, meta_interp)
                    self._writer.write(dataset, meta_interp['to_filename'])
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

    def _datum_metadata_ready(self, metadata):
        """
        Check the metadata to see if the required fields are populated.
        """
        tmp = FuseProcessor._datums.copy()
        tmp.remove('to_horiz_key')
        ready = True
        for key in tmp:
            if key not in metadata:
                ready = False
                break
        return ready

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
