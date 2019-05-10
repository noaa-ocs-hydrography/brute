# -*- coding: utf-8 -*-
"""
BDB51.py

Created on Thu May  9 15:46:49 2019

@author: grice

Classes and methods for working with CARIS Bathy DataBASE 5.1.
"""

import subprocess

class bdb51:
    """
    A class for doing I/O with the CARIS Bathy DataBASE 5.1 server.
    
    This class spawns a conda environment containing the CARIS python API and
    communicates with that environment over TCP sockets.
    """
    def __init__(self, database_loc, database_name):
        """
        Instantiate the object and connect to the database referenced, waiting
        until the database responds.
        """
        self._start_bdb51_env(database_loc, database_name)
        
    def _start_bdb51_env(self, database_loc, database_name):
        """
        Instantiate an environment containing CARIS BDB51 object for talking
        to the database.
        
        Communications with the environment will happen over sockets.
        """
        pass
    
    def _bdb51_alive(self):
        """
        Check to see if the subprocess is still communicating and connected to
        the database.
        """
        
    def _send2bdb51(self, dataset, instruction):
        """
        Send the BDB environment instructions on where data is and what to do
        with it.
        """
        pass