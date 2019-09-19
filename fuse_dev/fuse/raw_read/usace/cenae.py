# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:36:08 2019

@author: Casiano.Koprowski
"""

from fuse.raw_read.usace.usace import USACERawReader


class CENAERawReader(USACERawReader):
    def __init__(self):
        super().__init__('CENAE')

    def read_metadata(self, filename: str):
        """
        Function overwite of :func:`usace.USACERawReader.read_metadata` based
        on where the best metadata is for this district

        The CENAE metadata is retuned in order of precedence:
            1. The file name.
            2. The survey's ``.xyz`` header.
            3. The survey's ``.xml``.
            4. The metadata pickle pulled from eHydro.

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
        meta_xml = self._parse_usace_xml(filename)
        meta_xyz = self._parse_ehydro_xyz_header(basexyzname)
        meta_filename = self._parse_filename(filename)
        meta_pickle = self._parse_pickle(filename)
        meta_date = self._parse_start_date(filename,
                                           {**meta_pickle, **meta_xyz,
                                            **meta_xml})
        meta_determine = self._data_determination(meta_supplement, filename)
        meta_supplement = {**meta_determine, **meta_date, **meta_supplement}
        return {**meta_pickle, **meta_xml, **meta_xyz, **meta_filename,
                **meta_supplement}
