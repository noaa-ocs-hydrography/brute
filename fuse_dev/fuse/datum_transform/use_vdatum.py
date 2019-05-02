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
from tempfile import TemporaryDirectory as tempdir

class vdatum:
    """
    An object for working with VDatum.
    """
    def __init__(self, config, reader):
        self._config = config
        self._reader = reader
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

    def translate(self, infilename, in_hordat, in_verdat, out_epsg, out_verdat):
        """
        Translate the provided filename from the provided in datums to the out
        datums and return a gdal object.
        
        NSRS2007 is assumed for the out EPSG code. 
        """
        outxyz =  self._translatexyz(infilename, in_hordat, in_verdat, 
                                out_epsg, out_verdat)
        out_gdal = self._xyz2gdal(outxyz, out_epsg, out_verdat)
        return out_gdal
    
    def _translatexyz(self, infilename, in_hordat, in_verdat, out_epsg, 
                     out_verdat):
        """
        were you planning on actually documenting anything?
        
        """
        # read the bathy and put it in a temp file for vdatum to read
        bathy = self._reader.read_bathymetry(infilename)
        d = tempdir()
        outfilename = _os.path.join(d.name,'outfile.txt')
        _np.savetxt(outfilename, bathy, delimiter = ',')
        # set up vdatum
        self._setup_vdatum(in_hordat, in_verdat, out_epsg, out_verdat)
        # run vdatum
        self._convert_file(outfilename, d.name)
        bathy = _np.loadtxt(outfilename, delimiter = ',')
        return bathy
        
    def _setup_vdatum(self, in_fips, in_verdat, out_epsg, out_verdat):
        """
        Setup the VDatum command line arguments to convert points.
        
        This method current assums US Survey Feet, and convert it into UTM 
        (meters) with the otherwise the specified vertical datums.  Vertical 
        assumed to be positive down for both input and output. NAD83 is assumed
        for horizontal datums.
        
        The output epsg code is converted to a NSRS2007 zone using a dumb
        conversion.
        """
        out_zone = self._epsg2zone(out_epsg)
        ihorz = r'ihorz:NAD83:spc:us_ft:' + str(in_fips)
        ivert = ' ivert:' + in_verdat.lower() + ':us_ft:sounding'
        ohorz = ' ohorz:NAD83:utm:m:' + str(out_zone)
        overt = ' overt:' + out_verdat.lower()  + ':m:sounding'
        georef = ihorz + ivert + ohorz + overt
        java_str = _os.path.join(self._java_path,'java')
        file_str = ' -file:txt:comma,0,1,2,skip0:{};{}'
        self._shell = java_str + ' -jar vdatum.jar ' + georef + file_str
        
    def _convert_file(self, vdinfilename, vdoutdir):
        """
        The provided file of xyz points will be converted and returned
        according to the from and to the datums provided through the 
        _setup_vdatum_point_converstion method.
        """
        command = self._shell.format(vdinfilename, vdoutdir)
        try:
            proc = _subprocess.Popen(command, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, cwd=self._vdatum_path)
        except:
            print('Error executing: ' + command + '\nat: ' + self._vdatum_path)
            raise
        try:
            (stdout, stderr) = proc.communicate()
            print(stdout)
        except:
            print(stdout)
            print(stderr)
            
    def _xyz2gdal(self, outxyz, out_epsg, out_verdat):
        """
        Convert from numpy xyz array to a gdal dataset.
        """
        # setup the gdal bucket
        dest = osr.SpatialReference()
        dest.ImportFromEPSG(out_epsg)
        dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = dataset.CreateLayer('pts', dest, geom_type=ogr.wkbPoint)
        for p in outxyz:
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