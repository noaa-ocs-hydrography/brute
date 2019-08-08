# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:36:08 2019

@author: Casiano.Koprowski
"""

from . import usace


class CENAERawReader(usace.USACERawReader):
    """An abstract raw data reader."""

    def __init__(self):
        """
        No init needed?
        """
        usace.USACERawReader.__init__(self, version='CENAE')

    def read_metadata(self, infilename: str):
        """
        Function overwite of :func:`usace.USACERawReader.read_metadata` based
        on where the best metadata is for this district

        The CENAE metadata is retuned in order of precedence:
            1. The survey's ``.xyz`` header.
            2. The survey's ``.xml``.
            3. The file name.
            4. The metadata pickle pulled from eHydro.

        Parameters
        ----------
        infilename : str
            File path of the input ``.xyz`` data

        Returns
        -------
        dict :
            The complete metadata pulled from multiple sources

        """
        basexyzname, suffix = self.name_gen(infilename, ext='.xyz')
        meta_xml = self._parse_usace_xml(infilename)
        meta_xyz = self._parse_ehydro_xyz_header(basexyzname)
        meta_filename = self._parse_filename(basexyzname)
        meta_pickle = self._parse_pickle(infilename)
        meta_date = self._parse_start_date(infilename,
                                           {**meta_pickle, **meta_xyz,
                                            **meta_xml})
        if suffix is not None and suffix.upper() in self.xyz_suffixes:
#            meta_xml['from_horiz_reolution'] = 3
#            self._check_grid(infilename)
        return {**meta_pickle, **meta_filename, **meta_xml, **meta_xyz,
                **meta_date}
