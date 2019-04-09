# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 15:26:20 2018

@author: grice
"""
import os as _os
import osgeo as _osgeo

def fips2wkt(fips, units = 'FEET'):
    """
    Given an ESRI FIPS code, return the associated wkt string as found in the
    gdal module data file 'esri_StatePlane_extra.wkt'.  
    
    units default to 'FEET', but 'METER' can also be provided.  
    
    Only NAD83 codes without HARN are supported.
    """
    # combine the fips code with units to get the ERSI code
    fipsstr = str(fips)
    if units == 'FEET':
        esri_code = fipsstr + '2,'
    elif units == 'METER':
        esri_code = fipsstr + '1,'
    else:
        print("Unit type not recognized")
        esri_code = "-1"
    # search the file 'esri_StatePlane_extra.wkt' for the code
    # hard link to file in Pydro directory until I can figure out how to find it
    osgeo_path = _osgeo.__path__[0]
    relpath = 'data/gdal/esri_StatePlane_extra.wkt'
    wktfilename = _os.path.join(osgeo_path,relpath)
    wktline = ''
    with open(wktfilename, 'r') as wktfile:
        for line in wktfile:
            if line.startswith(esri_code):
                wktline = line
                break
    split_idx = wktline.find(',') + 1
    wkt = wktline[split_idx:]
    return wkt