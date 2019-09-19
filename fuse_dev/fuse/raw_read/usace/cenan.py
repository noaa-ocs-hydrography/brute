# -*- coding: utf-8 -*-
"""
read_raw_cenan.py

Created on Fri Feb  1 16:35:18 2019

@author: grice

Read the various data sources available for a particular data stream such that
any available bathymetry or metadata can be accessed.
"""
from fuse.raw_read.usace.usace import USACERawReader


class CENANRawReader(USACERawReader):
    """An abstract raw data reader."""

    def __init__(self):
        """
        No init needed?
        """
        super().__init__('CENAN')

    def read_metadata(self, filename: str):
        """
        Function overwite of :func:`usace.USACERawReader.read_metadata` based
        on where the best metadata is for this district

        The CENAN metadata is retuned in order of precedence:
            1. The file name.
            2. The survey's ``.xyz`` header.
            3. The metadata pickle pulled from eHydro.

        Parameters
        ----------
        filename : str
            File path of the input ``.xyz`` data

        Returns
        -------
        dict :
            The complete metadata pulled from multiple sources

        """

        meta_supplement = {}
        basexyzname, suffix = self.name_gen(filename, ext='.xyz')
        meta_xyz = self._parse_ehydro_xyz_header(basexyzname)
        meta_filename = self._parse_filename(filename)
        meta_pickle = self._parse_pickle(filename)
        meta_date = self._parse_start_date(filename,
                                           {**meta_pickle, **meta_xyz})
        meta_determine = self._data_determination(meta_supplement, filename)
        meta_supplement = {**meta_determine, **meta_date, **meta_supplement}
        return {**meta_pickle, **meta_xyz, **meta_filename,
                **meta_supplement}
