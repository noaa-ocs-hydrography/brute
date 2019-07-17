# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:49:10 2019

@author: Casiano.Koprowski
"""

import logging as _logging
import os as _os
import re as _re
import sys as _sys

import numpy as _np
from . import parse_usace_xml
from . import parse_usace_pickle

class Base:
    def __init__(self, version=None):
        self.version = version
        self.ussft2m = 0.30480060960121924  # US survey feet to meters
        self.xyz_suffixes = ('_A', '_FULL')
        self.xyz_files = {}

        self._logger = _logging.getLogger(f'fuse_{self.version}')

        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            ch.setLevel(_logging.DEBUG)
            self._logger.addHandler(ch)

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

        self.version = 'CENAE'

    def meta_xml(self, filename: str) -> dict:
        """
        Retrieves metadata for USACE E-Hydro files
        function returns metadata dictionary

        Parameters
        ----------
        filename: str


        Returns
        -------
        metadata : dict


        """
        xml_name = self.name_gen(filename, ext='xml', sfx=False)
        pickle_name = self.name_gen(filename, ext='pickle', sfx=False)
        pickle_dict = parse_usace_pickle.read_pickle(pickle_name, pickle_ext=True)
        return pickle_dict

    def name_gen(self, filename: str, ext: str = None, sfx: bool = True) -> str:
        """
        Returns the suffix of the a survey's xyz file and the file name
        for a different extension

        """
        filebase, fileext = _os.path.splitext(filename)
        suffix = None
        for item in self.xyz_suffixes:
            if _re.compile(f'{item}$', _re.IGNORECASE).search(filebase):
                suffix = item
        if suffix is not None and ext is not None:
            base = _re.sub(_re.compile(f'{suffix}$', _re.IGNORECASE), '', filebase) + f'.{ext}'
        elif suffix is None and ext is not None:
            base = filebase + f'.{ext}'
        else:
            base = filename

        if sfx:
            return base, suffix
        else:
            return base


