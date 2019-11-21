# -*- coding: utf-8 -*-
"""
read_raw_cenan.py

Created on Fri Feb  1 16:35:18 2019

@author: grice

Read the various data sources available for a particular data stream such that
any available bathymetry or metadata can be accessed.
"""

import logging

from fuse.raw_read.usace.usace import USACERawReader


class CENANRawReader(USACERawReader):
    """An abstract raw data reader."""

    def __init__(self, logger: logging.Logger = None):
        if logger is None:
            logger = logging.getLogger('fuse')
        self.logger = logger

        super().__init__('CENAN', self.logger)

    def read_metadata(self, survey_folder: str) -> dict:
        """
        Function overwite of :func:`usace.USACERawReader.read_metadata` based
        on where the best metadata is for this district

        The CENAN metadata is retuned in order of precedence:
            1. The file name.
            2. The survey's ``.xyz`` header.
            3. The metadata pickle pulled from eHydro.

        Parameters
        ----------
        survey_folder
            folder path of the input ``.xyz`` data

        Returns
        -------
        dict
            complete metadata pulled from multiple sources
        """

        meta_supplement = {}
        meta_determine, filename = self._data_determination(meta_supplement, survey_folder)
        basexyzname, suffix = self.name_gen(filename, ext='.xyz')
        meta_xyz = self._parse_ehydro_xyz_header(basexyzname)
        meta_filename = self._parse_filename(filename)
        meta_pickle = self._parse_pickle(filename)
        meta_date = self._parse_start_date(filename, {**meta_pickle, **meta_xyz})
        meta_supplement = {**meta_determine, **meta_date, **meta_supplement}
        meta_defaults = self._cenan_defaults()
        meta_combined = {**meta_defaults, **meta_pickle, **meta_xyz, **meta_filename, **meta_supplement}
        meta_final = self._finalize_meta(meta_combined)
        if meta_final['interpolate']:
            meta_orig = meta_final.copy()
            meta_orig['interpolate'] = False
            meta_final['from_filename'] = f"{meta_orig['from_filename']}.interpolate"
            return [meta_orig, meta_final]
        else:
            return [meta_final]

    def _cenan_defaults(self):
        """
        Return default expectations for this disctrict.
        
        If a reader returns a value it should supersede these assuptions.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        dict
            The default metadata
        
        """
        meta = {}
        meta['from_vert_direction'] = 'height'
        return meta
