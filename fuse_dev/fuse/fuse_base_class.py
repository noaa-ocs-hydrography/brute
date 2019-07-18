# -*- coding: utf-8 -*-
"""
fuse_abstract_obj

Created on Thu Jan 31 10:03:30 2019

@author: grice
"""

import os as _os


class fuse_base_class:
    """The fuse object."""

    def __init__(self, configfilename: str = 'generic.config'):
        """
        Initialize with the metadata file to use and the horizontal and
        vertical datums of the workflow.

        Parameters
        ----------
        configfilename
            filename of config file
        """

        self._configfilename = configfilename
        self._config = self._read_configfile(configfilename)
        self.rawdata_path = self._config['rawpaths']
        self.procdata_path = self._config['outpath']
        self._meta = {}

    def _read_configfile(self, confile: str):
        """
        Read, parse, and return the configuration information in the provided file.
        The actual format of this file is ....
        
        rawpaths
        outpath
        to_horiz_datum
        to_vert_datum
        metapath

        Parameters
        ----------
        confile: str :
            

        Returns
        -------

        """

        config = {}

        with open(confile, 'r') as configfile:
            for line in configfile:
                if len(line) > 0:  # ignore lines with nothing
                    if line[0] != '#':  # ignore commented lines
                        stub, info = line.split('=')

                        # clean these up a bit
                        info = info.replace('\n', '')
                        stub = stub.rstrip()
                        info = info.rstrip().lstrip()

                        if stub == 'outpath':
                            if _os.path.isdir(info):
                                config[stub] = info
                            else:
                                raise ValueError(f'Invalid output folder: {info}')
                        elif stub == 'to_horiz_datum':
                            config[stub] = int(info)
                        elif stub == 'rawpaths':
                            rawpaths = []
                            raw = info.split(';')

                            for r in raw:
                                if _os.path.isdir(r):
                                    rawpaths.append(r)
                                else:
                                    raise ValueError(f'Invalid input path: {r}')

                            config[stub] = rawpaths
                        else:
                            config[stub] = info

        if len(config) > 0:
            self._check_config(config)
        else:
            raise ValueError('Failed to read configuration file.')

        return config

    def _check_config(self, config_dict: dict):
        """
        Check to ensure at least the basic elements of the configuration are
        available.

        Parameters
        ----------
        config_dict: dict :
            

        Returns
        -------

        """

        if 'rawpaths' not in config_dict:
            raise ValueError('No path to raw data found in configuration file.')

        if 'outpath' not in config_dict:
            raise ValueError('No path to output data found in configuration file.')

        if 'to_horiz_datum' not in config_dict:
            raise ValueError('No output horizontal datum found in configuration file.')

        if 'to_vert_datum' not in config_dict:
            raise ValueError('No output vertical datum found in configuration file.')

        if 'metapath' not in config_dict:
            raise ValueError('No metadata output location found in configuration file.')

    def read(self, infilename: str):
        """
        The file name to use to read the bathymetry and metadata into useable
        forms.

        Parameters
        ----------
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
        infilename: str :
            

        Returns
        -------

        """

        pass
