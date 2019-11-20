# -*- coding: utf-8 -*-
"""
fuse_abstract_obj

Created on Thu Jan 31 10:03:30 2019

@author: grice
"""

import configparser as _cp
import logging as _logging
import os as _os

import fuse.datum_transform.transform as _trans
import fuse.interpolator.interpolator as _interp
import fuse.meta_review as _mr
import fuse.raw_read.noaa as _noaa
import fuse.raw_read.usace as _usace
from fuse import score
from fuse.proc_io.proc_io import ProcIO


class FuseProcessor:
    """Bathymetric survey object."""

    _datums = [
        'from_horiz_datum',
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

    _quality_metrics = [
        'from_horiz_unc',
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

    _paths = [
        'from_filename',
        'from_path',
        'to_filename',
        'support_files',
    ]

    _dates = [
        'start_date',
        'end_date',
    ]

    _source_info = [
        'agency',
        'source_indicator',
        'source_type',
        'interpolated',
        'posted',
        'license',
    ]

    _processing_info = [
        'logfilename',
        'version_reference',
        'interpolate',
        'file_size'
    ]

    _scores = [
        'catzoc',
        'supersession_score',
    ]

    def __init__(self, config_filename: str = 'generic.config'):
        """
        Initialize with the metadata file to use and the horizontal and
        vertical datums of the workflow.

        Parameters
        ----------
        config_filename
            path to file with configuration
        """

        self._config_filename = config_filename
        self._config = self._read_configfile(config_filename)
        self.rawdata_path = self._config['rawpaths']
        self.procdata_path = self._config['outpath']
        self._cols = FuseProcessor._paths + FuseProcessor._dates + FuseProcessor._datums + FuseProcessor._quality_metrics + FuseProcessor._scores + FuseProcessor._source_info + FuseProcessor._processing_info

        if 'metatable' in self._config:
            hostname, database, table = self._config['metatable'].split('/')
            self._meta_obj = _mr.MetadataDatabase(hostname, database, table, self._cols)
        elif 'metapath' in self._config:
            self._meta_obj = _mr.MetadataFile(self._config['metapath'], self._cols)
        else:
            raise ConfigParseError('configuration does not specify metadata table or file')

        self.logger = self._set_log('fuse')
        self.logger.info(f'configuration file: {self._config_filename}')
        self.logger.info(f'data input: {self.rawdata_path}')
        self.logger.info(f'data output: {self.procdata_path}')
        self.logger.info(f'metadata output: {self._config["metatable" if "metatable" in self._config else "metapath"]}')

        self._set_data_reader()
        self._set_data_transform()
        self._set_data_interpolator()
        self._set_data_writer()
        self._db = None
        self._meta = {}  # initialize the metadata holder

        self.logger.info('ready to start processing')

    def _read_configfile(self, configuration_file: str):
        """
        Read, parse, and return the configuration information in the provided
        file.  The actual format of this file is ....

        rawpaths
        outpath
        to_horiz_datum
        to_vert_datum
        metapath or metatable

        Parameters
        ----------
        configuration_file
            path to file with configuration

        Returns
        -------
            dictionary of metadata
        """
        if not _os.path.isfile(configuration_file):
            raise ValueError(f'file not found: {configuration_file}')
        config = {}
        config_file = _cp.ConfigParser()
        config_file.read(configuration_file)
        sections = config_file.sections()
        for section in sections:
            config_file_section = config_file[section]
            for key in config_file_section:
                value = config_file_section[key]
                if key == 'rawpaths':
                    rawpaths = []
                    raw = value.split(';')
                    for r in raw:
                        r = r.strip()
                        rawpaths.append(r)
                    config[key] = rawpaths
                else:
                    config[key] = value
        # add the root path
        if 'rootpath' in config:
            root = config['rootpath']
            # raw paths first
            config['rawpaths'] = [_os.path.join(root, path) for path in config['rawpaths']]
            # output paths
            config['outpath'] = _os.path.join(root, config['outpath'])
            # metapath
            if 'metapath' in config:
                config['metapath'] = _os.path.join(root, config['metapath'])
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
        config_dict
            dictionary of metadata keys
        """

        # dictionary of the basic configuration, with their descriptions
        required_config_keys = {
            'rawpaths': 'file path to raw data',
            'outpath': 'file path to output data',
            'to_horiz_datum': 'output horizontal datum',
            'to_horiz_frame': 'frame of output horizontal datum',
            'to_horiz_type': 'type (PCS / GCS) of output horizontal datum',
            'to_horiz_units': 'units of output horizontal datum',
            'to_vert_key': 'EPSG of output vertical datum',
            'to_vert_units': 'units of output vertical datum units',
            'to_vert_direction': 'vertical direction of output vertical datum',
            'to_vert_datum': 'output vertical datum'
        }
        for required_config_key in required_config_keys:
            if required_config_key not in config_dict:
                raise ConfigParseError(
                    f'configuration does not specify a {required_config_keys[required_config_key]} ("{required_config_key}")')

        if 'metapath' not in config_dict and 'metatable' not in config_dict:
            raise ConfigParseError('configuration does not specify a metadata location ("metapath" or "metatable")')

        # check the paths
        for p in config_dict['rawpaths']:
            if not _os.path.isdir(p):
                raise ConfigParseError(f'Invalid input path: {p}')
        if not _os.path.isdir(config_dict['outpath']):
            raise ConfigParseError(f'Invalid output data path: {config_dict["outpath"]}')

    def _set_data_reader(self):
        """
        Use information from the config file to set the reader to use for
        converting the raw data to usable metadata and bathymetry.
        """

        try:
            reader_type = self._config['raw_reader_type'].casefold()
            if reader_type == 'cenan':
                self._reader = _usace.cenan.CENANRawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'cemvn':
                self._reader = _usace.cemvn.CEMVNRawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'cesaj':
                self._reader = _usace.cesaj.CESAJRawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'cesam':
                self._reader = _usace.cesam.CESAMRawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'ceswg':
                self._reader = _usace.ceswg.CESWGRawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'cespl':
                self._reader = _usace.cespl.CESPLRawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'cenae':
                self._reader = _usace.cenae.CENAERawReader(self.logger)
                self._read_type = 'ehydro'
            elif reader_type == 'bag':
                self._reader = _noaa.bag.BAGSurvey(self._config['outpath'], self.logger)
                self._read_type = 'bag'
            elif reader_type == 'bps':
                self._reader = _noaa.bps.BPSRawReader(self.logger)
                self._read_type = 'bps'
            else:
                raise ValueError('reader type not implemented')
        except ValueError:
            raise ValueError("No reader type found in the configuration file.")

    def _set_data_transform(self):
        """Set up the datum transformation engine."""

        self._transform = _trans.DatumTransformer(self._config['vdatum_path'], self._config['java_path'], self._reader)

    def _set_data_interpolator(self):
        """Set up the interpolator engine."""
        engine = self._config['interpolation_engine']
        res = float(self._config['to_resolution'])
        method = self._config['interpolation_method']
        self._interpolator = _interp.Interpolator(engine, method, res)

    def _set_data_writer(self):
        """Set up the location and method to write tranformed and interpolated data."""

        self._raster_extension = self._config['bathymetry_intermediate_file']
        if self._raster_extension == 'bag':
            self._point_extension = 'csar'
        else:
            self._point_extension = self._raster_extension
        self._raster_writer = ProcIO('gdal', self._raster_extension, logger=self.logger)
        self._point_writer = ProcIO('point', self._point_extension, logger=self.logger)

    def read(self, dataid: str) -> [str]:
        """
        Read survey bathymetry and metadata into useable forms.

        Parameters
        ----------
        dataid
            survey name

        Returns
        ----------
        [str]
            input survey path
        """

        self.logger.info(f'reading {dataid}')

        logger = self._set_log(dataid)

        # get the metadata
        try:
            metadata = self._reader.read_metadata(dataid).copy()
            if type(metadata) == dict:
                metadata = [metadata]
        except (RuntimeError, TypeError) as error:
            logger.warning(f'{error.__class__.__name__} {error}')
            return []
        from_id = []
        for m in metadata:
            # get the config file information
            m = self._add_config_metadata(m)
            # check to see if the quality metadata is available.
            if not self._quality_metadata_ready(m, logger):
                logger.warning('Not all quality metadata was found during read.')
            else:
                logger.info('All quality metadata was found during read.')
            # write out the metadata and close the log
            self._meta_obj.insert_records(m)
            from_id.append(m['from_filename'])

        logger.debug(from_id)
        self._close_log(logger)
        return from_id

    def _add_config_metadata(self, metadata):
        """
        Add the metadata contained in the config file to the dictionary.

        Parameters
        ----------
        metadata
            metadata dictionary from the reader.

        Returns
        -------
        dict
            the provided metadata dictionary with the addition of the config
            file metadata.
        """
        metadata['to_horiz_frame'] = self._config['to_horiz_frame']
        metadata['to_horiz_type'] = self._config['to_horiz_type']
        metadata['to_horiz_units'] = self._config['to_horiz_units']
        if 'to_horiz_key' in self._config:
            metadata['to_horiz_key'] = self._config['to_horiz_key']
        metadata['to_vert_key'] = self._config['to_vert_key']
        metadata['to_vert_units'] = self._config['to_vert_units']
        metadata['to_vert_direction'] = self._config['to_vert_direction']
        metadata['to_vert_datum'] = self._config['to_vert_datum']
        return metadata

    def process(self, dataid: str) -> str:
        """
        Do the datum transformtion and interpolation (if required).

        Given the generic need to interpolate USACE data the 'interpolate'
        kwarg is set to True as a hack.  This information should be drawn from
        the data reader since there will be cases where we get full res data
        from the reader and interlation is not necessary.

        TODO: need to add checks to make sure the metadata is ready. Perhaps this should be added to the metadata object?

        Parameters
        ----------
        filename
            filename to process

        Returns
        ----------
        str
            output filename
        """

        self.logger.info(f'processing {dataid}')

        output_filename = ''

        metadata = self._get_stored_meta(dataid)
        logger = self._set_log(dataid)
        logger.info(f'processing {dataid}')
        metadata['read_type'] = self._read_type

        if self._datum_metadata_ready(metadata):
            # convert the bathy for the original data
            frompath = metadata['from_path']
            input_directory = _os.path.splitext(_os.path.basename(frompath))[0]
            metadata['outpath'] = _os.path.join(self._config['outpath'], input_directory)
            metadata['new_ext'] = self._point_extension

            dataset = None
            try:
                dataset, metadata, transformed = self._transform.reproject(frompath, metadata)
            except (ValueError, RuntimeError, IndexError) as error:
                logger.warning(f'transformation error: {error.__class__.__name__} - {error}')
                metadata['interpolate'] = False

            self._meta_obj.insert_records(metadata)

            if 'interpolate' in metadata:
                if metadata['interpolate'] and self._read_type == 'bag':
                    if ('support_files' not in metadata or len(metadata['support_files']) < 1):
                        metadata['interpolate'] = False
                        logger.warning("No coverage files provided; no interpolation can occur")

                if metadata['interpolate']:
                    metadata = self._transform.reproject_support_files(metadata, self._config['outpath'])

                    try:
                        dataset, metadata = self._interpolator.interpolate(dataset, metadata)
                        metadata['interpolated'] = True
                        self._raster_writer.write(dataset, metadata['to_filename'])
                    except (ValueError, RuntimeError, IndexError) as error:
                        logger.warning(f'interpolation error: {error.__class__.__name__} - {error}')

                else:
                    if self._read_type in ['ehydro', 'bps']:
                        metadata['to_filename'] = f'{metadata["outpath"]}.{metadata["new_ext"]}'
                        self._point_writer.write(dataset, metadata['to_filename'])
                    elif self._read_type == 'bag':
                        # only write out the bag if the file was transformed
                        if 'to_filename' not in metadata:
                            metadata['to_filename'] = f"{metadata['outpath']}.{self._raster_extension}"
                            self._raster_writer.write(dataset, metadata['to_filename'])

                    logger.info('No interpolation required')
            else:
                del dataset
                raise ValueError('metadata has no "interpolate" value')

            self._meta_obj.insert_records(metadata)
            output_filename = metadata['to_filename']
        else:
            logger.warning('metadata is missing required datum transformation entries')

        logger.info(f'processed -> {output_filename}')
        self._close_log(logger)
        return output_filename

    def post(self, filename):
        """
        Make the data available for amalgamation.

        TODO: need to add checks to make sure the metadata is ready.
            Perhaps this should be added to the metadata object?
        """

        self.logger.info(f'posting {filename}')

        metadata = self._get_stored_meta(filename)
        logger = self._set_log(filename)
        logger.info(f'posting {filename}')
        if self._quality_metadata_ready(metadata):
            if not self._score_metadata_ready(metadata):
                metadata['catzoc'] = score.catzoc(metadata)
                metadata['supersession_score'] = score.supersession(metadata)
            if self._db is None:
                self._connect_to_db()
            procfile = metadata['to_filename']
            metadata['posted'] = True
            s57_meta = self._metadata_to_s57(metadata)
            self._db.write(procfile, 'new', s57_meta)
            # need to check for proper insertion...
            self._meta_obj.insert_records(metadata)

        logger.info('posted')
        self._close_log(logger)

    def score(self, filename, date):
        """
        Provided a date, get the decayed quality metric and insert in the
        database, making the information available for amalgamation.
        """

        self.logger.info(f'scoring {filename}')

        metadata = self._get_stored_meta(filename)
        logger = self._set_log(filename)
        logger.info(f'scoring {filename}')

        if metadata['posted'].upper() == 'TRUE':
            dscore = score.decay(metadata, date)
            if self._db == None:
                self._connect_to_db()
            procfile = metadata['to_filename']
            s57_meta = self._metadata_to_s57(metadata)
            s57_meta['dcyscr'] = dscore
            self._db.write(procfile, 'metadata', s57_meta)
            logger.info(f'Posting new decay score of {dscore} to database.')
        else:
            logger.warning('Insertion of decay score failed.')

        logger.info('scored')
        self._close_log(logger)

    def _connect_to_db(self):
        """ Connect to the database defined in the configuration dictionary. """

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
        """ Asks proc_io to close the database connection """

        if self._db:
            self._db.close_connection()
        else:
            if 'database_location' in self._config:
                db_loc = self._config['database_location']
                raise RuntimeError(f'No connection to {db_loc} to close')
            else:
                raise ValueError('No database location defined in the configuration file.')

    def _set_log(self, name: str, file_level: str = _logging.DEBUG, console_level: str = _logging.INFO) -> _logging.Logger:
        """
        Set the global logger to the given filename.

        Parameters
        ----------
        name
            name of logger / log file
        file_level
            logging level of logfile
        console_level
            logging level of console

        Returns
        ----------
        logging.Logger
            logging object
        """

        name = _os.path.splitext(_os.path.basename(name))[0]
        log_directory = _os.path.dirname(self._config['metapath']) if 'metapath' in self._config else self._config['outpath']
        log_filename = _os.path.join(log_directory, f'{name}.log')
        self._config['logfilename'] = log_filename

        logger = _logging.Logger(name)

        # remove handlers that might have existed from previous files
        # self._close_log(logger)

        log_format = '%(asctime)s %(name)-30s %(levelname)-8s: %(message)s'

        # create handlers for this filename
        log_file = _logging.FileHandler(log_filename)
        log_file.setLevel(file_level)
        log_file.setFormatter(_logging.Formatter(log_format))
        logger.addHandler(log_file)

        console = _logging.StreamHandler()
        console.setLevel(console_level)
        console.setFormatter(_logging.Formatter(log_format))
        logger.addHandler(console)

        return logger

    def _close_log(self, logger: _logging.Logger):
        """ Close the object logging file. """

        # remove handlers
        for handler in logger.handlers:
            logger.removeHandler(handler)

    def _get_stored_meta(self, filename: str) -> dict:
        """
        Get the metadata in a local dictionary so that it can be used within
        the instantiated object.

        Parameters
        ----------
        filename
            filename of metadata file

        Returns
        -------
            dictionary of metadata
        """
        try:
            # file name is the key rather than the path
            f = _os.path.basename(filename)
            if 'from_filename' not in self._meta or self._meta['from_filename'] != filename:
                self._meta = self._meta_obj[f]
            else:
                # need to catch if this file is not in the metadata record yet here.
                raise KeyError(f'File not referenced in stored metadata: {f}')
            return self._meta
        except KeyError:
            return {}

    def _metadata_to_s57(self, metadata) -> dict:
        """
        The metadata is converted to an s57 version of the metadata.

        Parameters
        ----------
        metadata
            dictionary of metadata

        Returns
        -------
            dictionary of metadata
        """

        return _mr.csv_to_s57(metadata)

    def _datum_metadata_ready(self, metadata) -> bool:
        """
        Check the metadata to see if the required fields are populated.

        Parameters
        ----------
        metadata
            dictionary of metadata

        Returns
        ----------
            whether metadata has all datum fields
        """

        datum_keys = FuseProcessor._datums.copy()

        # if geographic input remove the need for a zone key
        if 'from_horiz_type' in metadata and metadata['from_horiz_type'] == 'geo' and 'from_horiz_key' in datum_keys:
            datum_keys.remove('from_horiz_key')

        return all(key in metadata for key in datum_keys if key != 'to_horiz_key')

    def _quality_metadata_ready(self, metadata, logger: _logging.Logger = None):
        """
        Check the metadata to see if the required fields are populated.

        Parameters
        ----------
        metadata
            dictionary of metadata
        logger
            logging object

        Returns
        ----------
            whether metadata has all quality fields
        """

        if logger is None:
            logger = self._set_log(metadata['from_filename'])

        # check the feature metadata
        if 'feat_detect' in metadata:
            if metadata['feat_detect']:
                feature_ready = 'feat_least_depth' in metadata and 'feat_size' in metadata
            else:
                feature_ready = True
        else:
            feature_ready = False

        if not feature_ready:
            logger.warning('Quality metadata for features is not yet available.')
        else:
            logger.info('Quality metadata for features was found')

        # check the uncertainty metadata
        vert_uncert_ready = 'vert_uncert_fixed' in metadata and 'vert_uncert_vari' in metadata
        horiz_uncert_ready = 'horiz_uncert_fixed' in metadata and 'horiz_uncert_vari' in metadata

        if not vert_uncert_ready or not horiz_uncert_ready:
            logger.warning('Quality metadata for uncertainty is not yet available.')
        else:
            logger.info('Quality metadata for uncertainty was found')

        # check the coverage
        coverage_ready = 'complete_coverage' in metadata and 'bathymetry' in metadata

        if not coverage_ready:
            logger.warning('Coverage metadata is not yet available.')
        else:
            logger.info('Coverage metadata was found')

        return feature_ready and vert_uncert_ready and horiz_uncert_ready and coverage_ready

    def _date_metadata_ready(self, metadata):
        """
        Check the metadata to see if the required fields are populated.

        Parameters
        ----------
        metadata
            dictionary of metadata
        Returns
        ----------
            whether metadata has all date fields
        """

        return all(key in metadata for key in FuseProcessor._dates)

    def _score_metadata_ready(self, metadata):
        """
        Check the metadata to see if the required fields are populated.

        Parameters
        ----------
        metadata
            dictionary of metadata
        Returns
        ----------
            whether metadata has all score fields
        """

        return all(key in metadata for key in FuseProcessor._scores)


class ConfigParseError(Exception):
    pass
