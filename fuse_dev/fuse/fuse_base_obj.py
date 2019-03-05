# -*- coding: utf-8 -*-
"""
fuse_abstract_obj

Created on Thu Jan 31 10:03:30 2019

@author: grice
"""

import os as _os

class fuse_obj:
    """
    The fuse object.
    """
    def __init__(self, configfilename = 'generic.config'):
        """
        Initialize with the metadata file to use and the horizontal and
        vertical datums of the workflow.
        """
        self._configfilename
        self._config = self._read_configfile(configfilename)
        self._meta = {}
        
    def _read_configfile(self, confile):
        """
        Read the config file and return the contents as a dictionary.
        """
        """
        Read, parse, and return the configuration information in the provided file.
        The actual format of this file is ....
        
        rawpaths
        outpath
        to_horiz_datum
        to_vert_datum
        metapath
        district
        """
        config = {}
        with open(confile, 'r') as configfile:
            for line in configfile:
                stub, info = line.split('=')
                # clean these up a bit
                info = info.replace('\n','')
                stub = stub.rstrip()
                info = info.rstrip().lstrip()
                if stub == 'outpath':
                    if _os.path.isdir(info):
                        config[stub] = info
                    else:
                        raise ValueError('Invalid output folder: ' + info)
                elif stub == 'to_horiz_datum':
                    config[stub] = int(info)
                elif stub == 'rawpaths':
                    rawpaths = []
                    raw = info.split(';')
                    for r in raw:
                        if _os.path.isdir(r):
                            rawpaths.append(r)
                        else:
                            raise ValueError('Invalid input path: ' + r)
                    config[stub] = rawpaths
                else:
                    config[stub] = info
        return config
        
    def read(self, infilename):
        """
        The file name to use to read the bathymetry and metadata into useable
        forms.
        """
        pass
    
    def process(self, infilename):
        """
        If the right metadata is available, perform any required datum
        transformations and interpolation.
        """
        pass
    
    def upload(self, infilename):
        """
        If the metadata is checked and the data is interpolated and transformed
        as needed, upload for amalgamation.
        
        If no filename is provided try to upload all files from the metadata
        file.
        """
        pass