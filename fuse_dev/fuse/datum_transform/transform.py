# -*- coding: utf-8 -*-
"""
transform.py

@author: grice

Abstract datum transformation.
"""

from fuse.datum_transform import use_vdatum as uv


class DatumTransformer:
    """
    An object for abstracting the datum transformation API.  This should allow
    for different transformation machines and versions.

    Parameters
    ----------

    Returns
    -------

    """

    def __init__(self, config: dict, reader):
        """
        Provided a key and a config file, get the method to use for datum
        conversion.
        """

        self._config = config
        self._reader = reader
        self._setup()

    def _setup(self):
        """
        Set up and configure the transformation tools based on the information
        provided in the configruation file.

        Parameters
        ----------

        Returns
        -------

        """

        if 'vdatum_path' in self._config:
            self._engine = uv.VDatum(self._config, self._reader)
        else:
            raise ValueError('No vdatum path provided')

    def translate(self, infilename: str, metadata: dict):
        """
        Run the specified transformation engine to translate the provided
        dataset.

        Parameters
        ----------
        infilename :
            param metadata:
        infilename: str :
            
        metadata: dict :
            

        Returns
        -------

        """

        self._meta = metadata
        return self._engine.translate(infilename, metadata)
