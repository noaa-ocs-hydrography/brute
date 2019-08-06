# -*- coding: utf-8 -*-
"""
fuse_ehydro.py

Created on Thu Jan 31 10:30:11 2019

@author: grice
"""

import logging as _logging
import fuse.fuse_processor as _fbc
import fuse.raw_read.usace as _usace

class FuseProcessor_eHydro(_fbc.FuseProcessor):
    """TODO write description"""
    _default_quality_metrics = {'complete_coverage': False,
                                'bathymetry': True,
                                'vert_uncert_fixed': 0.5,
                                'vert_uncert_vari': 0.1,
                                'horiz_uncert_fixed': 5.0,
                                'horiz_uncert_vari': 0.05,
                                'feat_detect': False,
                                }

    def __init__(self, config_filename):
        super().__init__(config_filename)

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
