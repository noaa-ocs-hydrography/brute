# -*- coding: utf-8 -*-
"""
read_raw_abstract.py

Created on Fri Feb  1 16:35:18 2019

@author: grice

Read the various data sources available for a particular data stream such that
any available bathymetry or metadata can be accessed.
"""

from abc import ABC, abstractmethod

import numpy


class RawReader(ABC):
    """An abstract raw data reader."""

    @abstractmethod
    def read_metadata(self, filename: str) -> dict:
        """
        Read all available meta data.

        Parameters
        ----------
        filename
            filename of survey

        Returns
        -------
        dict
            dictionary of metadata
        """

        pass

    @abstractmethod
    def read_bathymetry(self, filename: str) -> numpy.array:
        """
        Read the bathymetry.

        Parameters
        ----------
        filename
            filename of survey

        Returns
        -------
        numpy.array
            N x 3 array of XYZ points
        """

        pass
