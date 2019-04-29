# -*- coding: utf-8 -*-
"""
transform.py

@author: grice

Abstract datum transformation.
"""

from fuse.datum_transform import use_vdatum as uv

class transform:
    """
    An object for abstracting the datum transformation API.  This should allow
    for different transformation machines and versions.
    """
    def __init__(self, config):
        """
        Provided a key and a config file, get the method to use for datum
        conversion.
        """
        self._config = config
        self._setup()
        
    def _setup(self):
        """
        Set up and configure the transformation tools based on the information
        provided in the configruation file.
        """
        if 'vdatum_path' in self._config:
            self._engine = uv.vdatum(self._config)
        else:
            raise ValueError('No java path provided')
            
    def translate(self, infilename, metadata):
        """
        Run the specified transformation engine to translate the provided
        dataset.
        """
        self._meta = metadata
        in_fips = int(self._meta['from_fips'])
        in_verdat = self._meta['from_vert_datum']
        out_epsg = int(self._meta['to_horiz_datum'])
        out_verdat = self._meta['to_vert_datum']
        return self._engine.translate(infilename, 
                                      in_fips, 
                                      in_verdat, 
                                      out_epsg, 
                                      out_verdat)
        
        
