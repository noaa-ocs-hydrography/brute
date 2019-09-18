# -*- coding: utf-8 -*-
"""
read_raw_abstract.py

Created on Fri Feb  1 16:35:18 2019

@author: grice

Read the various data sources available for a particular data stream such that
any available bathymetry or metadata can be accessed.
"""

import numpy


class RawReader:
    """An abstract raw data reader."""

    def __init__(self, datatype: str):
        """
        Initialize the data reader with some knowledge of what types of data
        it is trying to read.
        """

        pass

    def read_metadata(self, infilename: str) -> dict:
        """
        Read all available meta data.

        Parameters
        ----------
        infilename
            filename of survey

        Returns
        -------
            dictionary of metadata
        """

        pass

    def read_bathymetry(self, infilename: str) -> numpy.array:
        """
        Read the bathymetry.

        Parameters
        ----------
        infilename
            filename of survey

        Returns
        -------
            array of XYZ points
        """

        pass
