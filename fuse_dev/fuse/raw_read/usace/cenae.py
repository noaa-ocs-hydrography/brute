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
