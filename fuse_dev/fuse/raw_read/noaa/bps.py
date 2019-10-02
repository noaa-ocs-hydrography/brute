# -*- coding: utf-8 -*-
"""
bps.py

20191002 grice

V0.0.1 20191002

A reader for Bathy Point Store data as extracted from the NCEI database.

These files appear to be written with the following header:
    
    survey_id,
    lat,
    long,
    depth,
    sounding method,
    end_year,
    active_flag
    
A gross assumption is made about the quality of the data since no additional
metadata is available.
"""

import sys as _sys
import os as _os
import numpy as _np
import logging as _logging
from fuse.raw_read.raw_read import RawReader

class BPSRawReader(RawReader):
    """
    An object for handling Bathy Point Store data is the same fashion as other
    data types.
    """
    def __init__(self):
        """
        Initialize the logger for the read process.
        """
        self._logger = _logging.getLogger(f'fuse')
        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            ch.setLevel(_logging.DEBUG)
            self._logger.addHandler(ch)
    
    def read_metadata(self, filename: str) -> dict:
        """
        The bathy point store files contain limited to no meta data.  The
        values returned from this method are hard coded.

        Parameters
        ----------
        filename
            file to read

        Returns
        -------
            dictionary of metadata
        """
        
        metadata = {}
        # date info
        year = int(_np.loadtxt(filename, delimiter=',', skiprows=1,
                         usecols = (5), max_rows = 1))
        date = str(year) + '0101'
        metadata['end_date'] = date
        # datum info
        metadata['from_horiz_datum'] = 'Chart Datum'
        metadata['from_horiz_frame'] = 'NAD83'
        metadata['from_horiz_type'] = 'geo'
        metadata['from_horiz_units'] = 'deg'
        # metadata['from_horiz_key'] = ''
        metadata['from_vert_datum'] = 'Chart Datum'
        metadata['from_vert_key'] = 'mllw'
        metadata['from_vert_units'] = 'm'
        metadata['from_vert_direction'] = 'sounding'
        # path info
        base = _os.path.basename(filename)
        name, ext = _os.path.splitext(base)
        metadata['from_filename'] = base
        metadata['from_path'] = filename
        # source info
        metadata['agency'] = 'US'
        metadata['source_indicator'] = 'US,US,graph,' + name
        metadata['source_type'] = ext
        metadata['license'] = 'CC0'
        # processing info
        metadata['interpolate'] = False
        if year <= 1940:
            # CATZOC C quality
            metadata['from_horiz_unc'] = 500
            # metadata['from_horiz_resolution'] = 
            metadata['from_vert_unc'] = 'CATZOC C'
            metadata['complete_coverage'] = False
            metadata['bathymetry'] = False
            metadata['vert_uncert_fixed'] = 2.0
            metadata['vert_uncert_vari'] = 0.05
            metadata['horiz_uncert_fixed'] = 500
            metadata['horiz_uncert_vari'] = 0
            # metadata['feat_size'] = 
            metadata['feat_detect'] = False 
            # metadata['feat_least_depth'] = 
        else:
            # CATZOC B quality
            metadata['from_horiz_unc'] = 50
            # metadata['from_horiz_resolution'] = 
            metadata['from_vert_unc'] = 'CATZOC B'
            metadata['complete_coverage'] = False
            metadata['bathymetry'] = False
            metadata['vert_uncert_fixed'] = 1
            metadata['vert_uncert_vari'] = 0.02
            metadata['horiz_uncert_fixed'] = 50
            metadata['horiz_uncert_vari'] = 0
            # metadata['feat_size'] = 
            metadata['feat_detect'] = False
            # metadata['feat_least_depth'] = 
        return metadata

    def read_bathymetry(self, infilename: str) -> _np.array:
        """
        Returns a BagFile data object

        Parameters
        ----------
        infilename : str
            BAG filepath

        Returns
        -------
            BagFile object
        """
        xyz = _np.loadtxt(infilename, delimiter=',', skiprows=1,
                         usecols = (1,2,3))
        xyz[:,[0,1]] = xyz[:,[1,0]]
        return xyz