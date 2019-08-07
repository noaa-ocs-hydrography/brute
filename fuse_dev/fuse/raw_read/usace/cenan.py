# -*- coding: utf-8 -*-
"""
read_raw_cenan.py

Created on Fri Feb  1 16:35:18 2019

@author: grice

Read the various data sources available for a particular data stream such that
any available bathymetry or metadata can be accessed.
"""

from . import usace


class CENANRawReader(usace.USACERawReader):
    """An abstract raw data reader."""

    def __init__(self):
        """
        No init needed?
        """
        usace.USACERawReader.__init__(self, version='CENAN')

    def read_metadata(self, infilename: str):
        """
        Read all available meta data.
        returns dictionary

        Parameters
        ----------
        infilename: str

        Returns
        -------

        """

        basexyzname, suffix = self.name_gen(infilename, ext='.xyz')
        meta_xyz = self._parse_ehydro_xyz_header(basexyzname)
        meta_filename = self._parse_filename(basexyzname)
        meta_pickle = self._parse_pickle(infilename)
        meta_date = self._parse_start_date(infilename,
                                           {**meta_pickle, **meta_xyz})
        return {**meta_pickle, **meta_filename, **meta_xyz,
                **meta_date}
