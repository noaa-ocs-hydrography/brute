# -*- coding: utf-8 -*-
"""
fuse_abstract_obj

Created on Thu Jan 31 10:03:30 2019

@author: grice
"""

import configparser as _cp
import os as _os


class FuseProcessor:
    """The fuse object."""

    def __init__(self, configfilename: str = 'generic.config'):
        """
        Initialize with the metadata file to use and the horizontal and
        vertical datums of the workflow.

        Parameters
        ----------
        configfilename
            test
        """
        self._configfilename = configfilename
        self._config = self._read_configfile(configfilename)
        self.rawdata_path = self._config['rawpaths']
        self.procdata_path = self._config['outpath']
        self._meta = {}

    def _read_configfile(self, confile: str):
        """
        Read, parse, and return the configuration information in the provided
        file.  The actual format of this file is ....

        rawpaths
        outpath
        to_horiz_datum
        to_vert_datum
        metapath

        Parameters
        ----------
        confile :

        confile: str :


        Returns
        -------

        """

        config = {}
        config_file = _cp.ConfigParser()
        config_file.read(confile)
        sections = config_file.sections()
        for section in sections:
            for key in config_file[section]:
                if key == 'rawpaths':
                    rawpaths = []
                    raw = config_file[section][key].split(';')
                    for r in raw:
                        if _os.path.isdir(r):
                            rawpaths.append(r)
                        else:
                            raise ValueError(f'Invalid input path: {r}')
                    config[key] = rawpaths
                elif key == 'outpath':
                    if _os.path.isdir(config_file[section][key]):
                        config[key] = config_file[section][key]
                    else:
                        raise ValueError(f'Invalid input path: {config_file[section][key]}')
                else:
                    config[key] = config_file[section][key]

        if len(config) == 0:
            raise ValueError('Failed to read configuration file.')
        else:
            self._check_config(config)
        return config

    def _check_config(self, config_dict: dict):
        """
        Check to ensure at least the basic elements of the configuration are
        available.

        Parameters
        ----------
        config_dict :

        config_dict: dict :


        Returns
        -------

        """
        options = {'rawpaths': 'path to raw data',
                   'outpath': 'path to output data',
                   'to_horiz_datum': 'output horizontal datum',
                   'to_vert_datum': 'output vertical datum',
                   'metapath': 'metadata output',
                   }
        for key in options.keys():
            if key not in config_dict:
                raise ValueError(f'No {options[key]} found in configuration file.')

    def read(self, infilename: str):
        """
        The file name to use to read the bathymetry and metadata into useable
        forms.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        pass

    def process(self, infilename: str):
        """
        If the right metadata is available, perform any required datum
        transformations and interpolation.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        pass

    def upload(self, infilename: str):
        """
        If the metadata is checked and the data is interpolated and transformed
        as needed, upload for amalgamation.

        If no filename is provided try to upload all files from the metadata
        file.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        pass
