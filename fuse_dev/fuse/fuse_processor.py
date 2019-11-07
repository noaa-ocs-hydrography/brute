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
        self._meta_obj = _mr.MetadataTable('ocs-vs-nbs01', 'metadata', self._config['metapath'], self._cols)
        self._set_data_reader()
        self._set_data_transform()
        self._set_data_interpolator()
        self._set_data_writer()
        self._db = None
        self._meta = {}  # initialize the metadata holder
        self.logger = _logging.getLogger('fuse')
        self.logger.setLevel(_logging.DEBUG)

    def _read_configfile(self, configuration_file: str):
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
            rawtmp = []
            for p in config['rawpaths']:
                rawtmp.append(_os.path.join(root, p))
            config['rawpaths'] = rawtmp
            # output paths
            config['outpath'] = _os.path.join(root, config['outpath'])
            # metapath
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
            'to_vert_datum': 'output vertical datum',
            'metapath': 'file path to CSV of metadata',
        }
        for required_config_key in required_config_keys:
            if required_config_key not in config_dict:
                raise ValueError(
                    f'no {required_config_keys[required_config_key]} ("{required_config_key}") found in configuration file')
        # check the paths
        for p in config_dict['rawpaths']:
            if not _os.path.isdir(p):
                raise ValueError(f'Invalid input path: {p}')
        if not _os.path.isdir(config_dict['outpath']):
            raise ValueError(f'Invalid output data path: {config_dict["outpath"]}')

    def _set_data_reader(self):
        """
        Use information from the config file to set the reader to use for
        converting the raw data to usable metadata and bathymetry.
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
            elif reader_type == 'bps':
                self._reader = _noaa.bps.BPSRawReader()
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
        self._raster_writer = ProcIO('gdal', self._raster_extension)
        self._point_writer = ProcIO('point', self._point_extension)

    def read(self, dataid: str) -> str:
        """
        Read survey bathymetry and metadata into useable forms.

        Parameters
        ----------
        dataid
            Filename of survey bathymetry.

        Returns
        ----------
        str
            input survey path
        """

        self._set_log(dataid)

        # get the metadata
        try:
            metadata = self._reader.read_metadata(dataid).copy()
        except (RuntimeError, TypeError) as e:
            self.logger.log(_logging.DEBUG, e)
            return None

        # get the config file information
        metadata = self._add_config_metadata(metadata)

        # check to see if the quality metadata is available.
        if not self._quality_metadata_ready(metadata):
            msg = f'Not all quality metadata was found during read.'
        else:
            msg = f'All quality metadata was found during read.'
        self.logger.log(_logging.DEBUG, msg)

        # write out the metadata and close the log
        self._meta_obj.extend(metadata)
        self._close_log()

        return metadata['from_path']

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

    def process(self, filename: str) -> str:
        """
        Do the datum transformtion and interpolation.

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
        if self._read_type == 'ehydro':
            meta_entry = self._reader.name_gen(_os.path.basename(filename), '', sfx=None)
        else:
            meta_entry = _os.path.splitext(_os.path.basename(filename))[0]
        metadata = self._get_stored_meta(meta_entry)
        self._set_log(meta_entry)
        metadata['read_type'] = self._read_type

        if self._datum_metadata_ready(metadata):
            # convert the bathy for the original data
            input_directory = _os.path.splitext(_os.path.basename(filename))[0]
            metadata['outpath'] = _os.path.join(self._config['outpath'], input_directory)
            metadata['new_ext'] = self._point_extension

            try:
                dataset, metadata, transformed = self._transform.reproject(filename, metadata)
                if self._read_type == 'ehydro' or self._read_type == 'bps':
                    outfilename = f"{metadata['outpath']}.{metadata['new_ext']}"
                    self._point_writer.write(dataset, outfilename)
                    metadata['to_filename'] = outfilename
                elif self._read_type == 'bag' and transformed:
                    metadata['to_filename'] = f"{metadata['outpath']}.{self._raster_extension}"
                    self._raster_writer.write(dataset, metadata['to_filename'])
            except (ValueError, RuntimeError, IndexError) as error:
                message = f' Transformation error: {error}'
                self.logger.warning(message)
                metadata['interpolate'] = 'False'

            self._meta_obj.extend(metadata)

            if 'interpolate' in metadata:
                interpolate = metadata['interpolate'].lower()

                if self._read_type == 'bag' and interpolate == 'true':
                    if ('support_files' not in metadata or len(metadata['support_files']) < 1):
                        interpolate = 'False'
                        self.logger.warning("No coverage files provided; no interpolation can occur")

                if interpolate == 'true':
                    meta_interp = metadata.copy()
                    meta_interp = self._transform.reproject_support_files(meta_interp, self._config['outpath'])

                    try:
                        root, filename = _os.path.split(meta_interp['outpath'])
                        base = _os.path.splitext(filename)[0]
                        meta_interp['from_filename'] = f'{base}.interpolated'

                        if self._config['interpolation_engine'] == 'raster':
                            output_filename = f"{_os.path.join(root, base)}_interp.{self._raster_extension}"
                        else:
                            resolution = float(self._config['to_resolution'])
                            resolution_string = f'{int(resolution)}m' if resolution >= 1 else f'{int(resolution * 100)}cm'
                            output_filename = f'{_os.path.join(root, base)}_{resolution_string}_interp.{self._raster_extension}'

                        meta_interp['to_filename'] = output_filename

                        dataset, meta_interp = self._interpolator.interpolate(dataset, meta_interp)
                        self._raster_writer.write(dataset, meta_interp['to_filename'])
                    except (ValueError, RuntimeError, IndexError) as error:
                        message = f' Interpolation error: {error}'
                        print(message)
                        self.logger.warning(message)

                    self._meta_obj.extend(meta_interp)
                    metadata.update(meta_interp)
                elif interpolate == 'false':
                    self.logger.log(_logging.DEBUG, f'{input_directory} - No interpolation required')
            else:
                del dataset
                raise ValueError('metadata has no "interpolate" value')
        else:
            self.logger.log(_logging.DEBUG, 'metadata is missing required datum transformation entries')

        self._close_log()
        return metadata['to_filename']

    def post(self, filename):
        """
        Make the data available for amalgamation.

        TODO: need to add checks to make sure the metadata is ready.
            Perhaps this should be added to the metadata object?
        """

        metadata = self._get_stored_meta(filename)
        self._set_log(filename)
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
            self._meta_obj.extend(metadata)
        self._close_log()

    def score(self, filename, date):
        """
        Provided a date, get the decayed quality metric and insert in the
        database, making the information available for amalgamation.
        """

        metadata = self._get_stored_meta(filename)
        self._set_log(filename)

        if metadata['posted'].upper() == 'TRUE':
            dscore = score.decay(metadata, date)
            if self._db == None:
                self._connect_to_db()
            procfile = metadata['to_filename']
            s57_meta = self._metadata_to_s57(metadata)
            s57_meta['dcyscr'] = dscore
            self._db.write(procfile, 'metadata', s57_meta)
            log = f'Posting new decay score of {dscore} to database.'
        else:
            log = 'Insertion of decay score failed.'

        self.logger.log(_logging.DEBUG, log)
        self._close_log()

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

    def _set_log(self, filename: str):
        """
        Set the global logger to the given filename.

        Parameters
        ----------
        filename
            filename to set logger to
        """

        root, extension = _os.path.splitext(filename)
        if extension == '.interpolated':
            filename = root
        log_filename = _os.path.join(_os.path.dirname(self._config['metapath']),
                                     f'{_os.path.splitext(_os.path.split(filename)[-1])[0]}.log')
        self._meta['logfilename'] = log_filename

        # remove handlers that might have existed from previous files
        for handler in self.logger.handlers:
            self.logger.removeHandler(handler)

        # create file handler for this filename
        file_handler = _logging.FileHandler(log_filename)
        file_handler.setLevel(_logging.DEBUG)
        file_handler.setFormatter(_logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(file_handler)

    def _close_log(self):
        """ Close the object logging file. """

        # remove handlers
        for handler in self.logger.handlers:
            self.logger.removeHandler(handler)

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
                self._meta = self._meta_obj.read_meta_record(f)
            else:
                # need to catch if this file is not in the metadata record yet here.
                raise KeyError(f'File not referenced in stored metadata: {f}')
            return self._meta
        except KeyError as e:
            return None

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

        return self._meta_obj.csv_to_s57(metadata)

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
        test_keys = FuseProcessor._datums.copy()
        # if geographic input remove the need for a zone key
        if metadata['from_horiz_type'] == 'geo' and 'from_horiz_key' in test_keys:
            idx = test_keys.index('from_horiz_key')
            test_keys.pop(idx)
        return all(key in metadata for key in test_keys if key != 'to_horiz_key')

    def _quality_metadata_ready(self, metadata):
        """
        Check the metadata to see if the required fields are populated.

        Parameters
        ----------
        metadata
            dictionary of metadata
        Returns
        ----------
            whether metadata has all quality fields
        """

        # check the feature metadata
        if 'feat_detect' in metadata:
            if metadata['feat_detect']:
                feature_ready = 'feat_least_depth' in metadata and 'feat_size' in metadata
            else:
                feature_ready = True
        else:
            feature_ready = False

        if not feature_ready:
            self.logger.log(_logging.DEBUG, 'Quality metadata for features is not yet available.')
        else:
            self.logger.log(_logging.DEBUG, 'Quality metadata for features was found')

        # check the uncertainty metadata
        vert_uncert_ready = 'vert_uncert_fixed' in metadata and 'vert_uncert_vari' in metadata
        horiz_uncert_ready = 'horiz_uncert_fixed' in metadata and 'horiz_uncert_vari' in metadata

        if not vert_uncert_ready or not horiz_uncert_ready:
            self.logger.log(_logging.DEBUG, 'Quality metadata for uncertainty is not yet available.')
        else:
            self.logger.log(_logging.DEBUG, 'Quality metadata for uncertainty was found')

        # check the coverage
        coverage_ready = 'complete_coverage' in metadata and 'bathymetry' in metadata

        if not coverage_ready:
            self.logger.log(_logging.DEBUG, 'Coverage metadata is not yet available.')
        else:
            self.logger.log(_logging.DEBUG, 'Coverage metadata was found')

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

        for key in FuseProcessor._dates:
            if key not in metadata:
                return False
        else:
            return True

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

        for key in FuseProcessor._scores:
            if key not in metadata:
                return False
        else:
            return True
