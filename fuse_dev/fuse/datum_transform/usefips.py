# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 15:26:20 2018

@author: grice
"""
import logging
import os as _os
import sys


def fips2wkt(fips: int, units: str = 'FEET', logger: logging.Logger = None) -> str:
    """
    Given an ESRI FIPS code, return the associated wkt string as found in the
    gdal module data file 'esri_StatePlane_extra.wkt'.
    
    units default to 'FEET', but 'METER' can also be provided.
    
    Only NAD83 codes without HARN are supported.

    Parameters
    ----------
    fips
        ESRI FIPS code
    units
        (Default value = 'FEET')
    logger
        logging object

    Returns
    -------
    str
        well-known text of FIPS
    """

    if logger is None:
        logger = logging.getLogger('fuse')

    # combine the fips code with units to get the ERSI code
    if units == 'FEET':
        esri_code = f'{fips}2,'
    elif units == 'METER':
        esri_code = f'{fips}1,'
    else:
        logger.warning("Unit type not recognized")
        esri_code = "-1"
    # search the file 'esri_StatePlane_extra.wkt' for the code
    # path to the gdal data
    gdal_data_path = _os.environ['GDAL_DATA'] if 'GDAL_DATA' in _os.environ \
        else _os.path.join(sys.exec_prefix, 'Library', 'share', 'gdal')

    wktfilename = _os.path.join(gdal_data_path, 'esri_StatePlane_extra.wkt')
    wktline = ''
    with open(wktfilename, 'r') as wktfile:
        for line in wktfile:
            if line.startswith(esri_code):
                wktline = line
                break
    split_idx = wktline.find(',') + 1
    wkt = wktline[split_idx:]
    return wkt
