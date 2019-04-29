# -*- coding: utf-8 -*-
"""
use_vdatum.py
Created on Wed Aug 22 12:27:39 2018

@author: grice

Use VDatum for conversions. 
"""
__version__ = 'use_vdatum 0.0.1'

import subprocess as _subprocess
import os as _os
import numpy as _np
from osgeo import gdal,ogr,osr

class vdatum:
    """
    An object for working with VDatum.
    """
    def __init__(self, config):
        self._config = config
        self._setup()
        
    def _setup(self):
        """
        Set up the object based on the provided configuration.
        """
        if 'vdatum_path' in self._config:
            pth = self._config['vdatum_path']
            if _os.path.isdir(pth):
                self._vdatum_path = pth
            else:
                raise ValueError('Invalid vdatum folder: ' + pth)
        else:
            raise ValueError('No VDatum path provided')
        if 'java_path' in self._config:
            pth = self._config['java_path']
            if _os.path.isdir(pth):
                self._java_path = pth
            else:
                raise ValueError('Invalid java path: ' + pth)
        else:
            raise ValueError('No java path provided')
            
    def translate(self, infilename, in_fips, in_verdat, out_epsg, out_verdat):
        """
        NSRS2007 is assumed for the out EPSG code. 
        """
        return self._ehydro2gdal(infilename, 
                                 in_fips, 
                                 in_verdat, 
                                 out_epsg, 
                                 out_verdat)
    
    def _ehydro2gdal(self, infilename, in_fips, in_verdat, out_epsg, 
                    out_verdat = 'MLLW', out_dir = '', rm_files = True):
        """
        Combine ehydro2newdat and dat2gdal into one call.
        
        If out_dir is not specified the VDatum working files will be written to the
        currnet working directory.
        
        rm_files defaults to True, which removes the intermediary VDatum file.
        
        NSRS2007 is assumed.
        """
        if out_dir == '':
            out_dir = _os.getcwd()
        out_zone = self._epsg2zone(out_epsg)
        self._ehydro2newdat(infilename = infilename,
                      in_fips = in_fips, in_verdat = in_verdat,
                      out_zone = out_zone, out_verdat = out_verdat,
                      out_dir = out_dir)
        inpath, infile = _os.path.split(infilename)
        intermed_file = _os.path.join(out_dir, infile)
        
        gdal_ds = self._dat2gdal(intermed_file, out_epsg)
        if rm_files:
            _os.remove(intermed_file)
            _os.remove(intermed_file+'.log')
        return gdal_ds
    
    def _ehydro2newdat(self, infilename, in_fips, in_verdat, 
                      out_dir, out_zone, out_verdat = 'MLLW'):
        """
        Use VDatum to read an eHydro *.dat file, assuming US Survey Feet, and 
        convert it into UTM (meters) with the otherwise the specified vertical
        datums.  Vertical assumed to be positive down for both input and output.
        NAD83 is assumed for horizontal datums.  
        """
        ihorz = r'ihorz:NAD83:spc:us_ft:' + str(in_fips)
        ivert = ' ivert:' + in_verdat.lower() + ':us_ft:sounding'
        ohorz = ' ohorz:NAD83:utm:m:' + str(out_zone)
        overt = ' overt:' + out_verdat.lower()  + ':m:sounding'
        georef = ihorz + ivert + ohorz + overt
        java_str = _os.path.join(self._java_path,'java')
        shell = java_str + r' -jar vdatum.jar %s -file:txt:space,0,1,2:%s;%s' % (georef, infilename, out_dir)
        try:
            proc = _subprocess.Popen(shell, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, cwd=self._vdatum_path)
        except:
            print('Error executing: ' + shell + '\nat: ' + self._vdatum_path)
            raise
        try:
            (stdout, stderr) = proc.communicate()
        except:
            print(stdout)
            print(stderr)
            
    def _dat2gdal(self, infilename, dest_epsg):
        """
        Read a dat file and turn it into a GDAL point cloud.
        """
        points = _np.loadtxt(infilename, delimiter = ' ')
        # datum and unit conversion function declaration
        dest = osr.SpatialReference()
        dest.ImportFromEPSG(dest_epsg)
        # turn numpy points into ogr points in a gdal dataset
        dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = dataset.CreateLayer('pts', dest, geom_type=ogr.wkbPoint)
        for p in points:
            newp = ogr.Geometry(ogr.wkbPoint)
            newp.AddPoint(p[0], p[1], p[2])
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(newp)
            layer.CreateFeature(feature)
        return dataset
        
            
    def _zone2epsg(self, zone):
        """
        Assume the EPSG code for a UTM zone is 3707 + the zone number.  This should
        work for zones 1 through 19 for NSRS2007.
        http://spatialreference.org/ref/?search=nad83+utm+zone
        """
        return int(zone) + 3707
    
    def _epsg2zone(self, epsg):
        """
        Assume the EPSG code for a UTM zone is 3707 + the zone number.  This should
        work for zones 1 through 19 for NSRS2007.
        http://spatialreference.org/ref/?search=nad83+utm+zone
        """
        return epsg - 3707