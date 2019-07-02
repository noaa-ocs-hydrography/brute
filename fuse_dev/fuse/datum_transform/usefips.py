# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 15:26:20 2018

@author: grice
"""

import os as _os


def fips2wkt(fips: int, units: str = 'FEET'):
    """
    Given an ESRI FIPS code, return the associated wkt string as found in the
    gdal module data file 'esri_StatePlane_extra.wkt'.
    
    units default to 'FEET', but 'METER' can also be provided.
    
    Only NAD83 codes without HARN are supported.

    Parameters
    ----------
    fips :
        param units:  (Default value = 'FEET')
    fips: int :
        
    units: str :
         (Default value = 'FEET')

    Returns
    -------

    """

    # combine the fips code with units to get the ERSI code
    fipsstr = str(fips)
    if units == 'FEET':
        esri_code = f'{fipsstr}2,'
    elif units == 'METER':
        esri_code = f'{fipsstr}1,'
    else:
        print("Unit type not recognized")
        esri_code = "-1"
    # search the file 'esri_StatePlane_extra.wkt' for the code
    gdal_data_path = _os.environ['GDAL_DATA']  # path to the gdal data
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
