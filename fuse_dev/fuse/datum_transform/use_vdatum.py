# -*- coding: utf-8 -*-
"""
use_vdatum.py
Created on Wed Aug 22 12:27:39 2018

@author: grice

Use VDatum for conversions. 
"""

from typing import Tuple, List

__version__ = 'use_vdatum 0.0.1'

import logging as _logging
import os as _os
import subprocess as _subprocess
from tempfile import TemporaryDirectory as tempdir

import numpy as _np
from osgeo import gdal, ogr, osr


from_hdatum = ['from_horiz_frame',
                'from_horiz_type',
                'from_horiz_units',
                'from_horiz_key',
                ]

from_vdatum = ['from_vert_key',
               'from_vert_units',
               'from_vert_direction',
               ]

to_hdatum = ['to_horiz_frame',
             'to_horiz_type',
             'to_horiz_units',
             ]

to_vdatum = ['to_vert_key',
             'to_vert_units',
             'to_vert_direction',
             ]

class VDatum:
    """An object for working with VDatum."""

    def __init__(self, config: dict, reader):
        """

        Parameters
        ----------
        config
        reader
        """

        self._config = config
        self._reader = reader
        self._setup()
        self._logger = _logging.getLogger('fuse')

    def _setup(self):
        """Set up the object based on the provided configuration."""

        if 'vdatum_path' in self._config:
            pth = self._config['vdatum_path']
            if _os.path.isdir(pth):
                self._vdatum_path = pth
            else:
                raise ValueError(f'Invalid vdatum folder: {pth}')
        else:
            raise ValueError('No VDatum path provided')
        if 'java_path' in self._config:
            pth = self._config['java_path']
            if _os.path.isdir(pth):
                self._java_path = pth
            else:
                raise ValueError(f'Invalid java path: {pth}')
        else:
            raise ValueError('No java path provided')

    def translate(self, infilename: str, instructions: dict) -> gdal.Dataset:
        """
        Translate the provided filename from the provided in datums to the out
        datums and return a gdal object.
        
        NSRS2007 is assumed for the out EPSG code.

        Parameters
        ----------
        infilename :
            param in_hordat:

        instructions : dict :
            

        Returns
        -------

        """

        self._logger.log(_logging.DEBUG, 'Begin datum transformation')
        if not self._have_instructions(instructions):
            raise ValueError('The required fields for transforming datums are not available')
        outxyz, out_zone = self._translatexyz(infilename, instructions)
        out_verdat = instructions['to_vert_key'].upper()
        out_gdal = self._xyz2gdal(outxyz, out_zone, out_verdat)  # passing UTM zone instead of EPSG code
        self._logger.log(_logging.DEBUG, 'Datum transformation complete')
        return out_gdal

    def _translatexyz(self, infilename: str, instructions: dict) -> Tuple[
        _np.array, int]:
        """
        TODO write description

        Parameters
        ----------
        infilename :
            param in_hordat:
        in_verdat :
            param out_epsg:
        out_verdat :
            
        infilename: str :
            
        in_hordat: str :
            
        in_verdat: str :
            
        out_epsg: int :
            
        out_verdat: str :
            

        Returns
        -------

        """

        # read the bathy and put it in a temp file for vdatum to read
        bathy = self._reader.read_bathymetry(infilename)
        d = tempdir()
        outfilename = _os.path.join(d.name, 'outfile.txt')
        vd_dir = tempdir()
        vdfilename = _os.path.join(vd_dir.name, 'outfile.txt')
        vdlogfilename = _os.path.join(vd_dir.name, 'outfile.txt.log')
        _np.savetxt(outfilename, bathy, delimiter=',')
        # set up vdatum
        self._setup_vdatum(instructions)
        # run vdatum
        self._convert_file(outfilename, vd_dir.name)
        if 'to_horiz_key' in instructions:
            out_zone = instructions['to_horiz_key']
        else:
            # read out UTM Zone from VDatum log file
            with open(vdlogfilename, 'r') as vd_log:
                for line in vd_log.readlines():
                    if line.startswith('Zone:'):
                        Output = line[54:82]
                        Inputzone = line[27:53]
                        out_zone = Output.rstrip(' ')
        new_bathy = _np.loadtxt(vdfilename, delimiter=',')
        return new_bathy, out_zone

    def _setup_vdatum(self, instructions: dict):
        """
        Setup the VDatum command line arguments to convert points.
        
        This method current assums US Survey Feet, and convert it into UTM
        (meters) with the otherwise the specified vertical datums.  Vertical
        assumed to be positive down for both input and output. NAD83 is assumed
        for horizontal datums.
        
        The output epsg code is converted to a NSRS2007 zone using a dumb
        conversion.

        Parameters
        ----------
        instructions : dict :
            

        Returns
        -------

        """
        ihorz_vals = [instructions[v] for v in from_hdatum]
        ihorz = 'ihorz:' + ':'.join(ihorz_vals)
        ivert_vals = [instructions[v] for v in from_vdatum]
        ivert = 'ivert:'  + ':'.join(ivert_vals)
        ohorz_vals = [instructions[v] for v in to_hdatum]
        ohorz = 'ohorz:' + ':'.join(ohorz_vals)
        # having the output zone is optional
        if 'to_horiz_key' in instructions:
            ohorz = ohorz + ':' + instructions['to_horiz_key']
        else:
            ohorz = ohorz + ':'
        overt_vals = [instructions[v] for v in to_vdatum]
        overt = 'overt:' + ':'.join(overt_vals)
        georef = f'{ihorz} {ivert} {ohorz} {overt}'
        java_str = _os.path.join(self._java_path, 'java')
        file_str = ' -file:txt:comma,0,1,2,skip0:'
        self._shell = f'{java_str} -jar vdatum.jar {georef}{file_str}'

    def _convert_file(self, vdinfilename: str, vdoutdir: str):
        """
        The provided file of xyz points will be converted and returned
        according to the from and to the datums provided through the
        _setup_vdatum_point_converstion method.

        Parameters
        ----------
        vdinfilename :
            param vdoutdir:
        vdinfilename: str :
            
        vdoutdir: str :
            

        Returns
        -------

        """

        command = f'{self._shell}{vdinfilename};{vdoutdir}'
        self._logger.log(_logging.DEBUG, command)
        try:
            proc = _subprocess.Popen(command, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE, cwd=self._vdatum_path)
        except:
            print(f'Error executing: {command}\nat: {self._vdatum_path}')
            raise
        try:
            (output, outerr) = proc.communicate()
            output_str = output.decode('utf-8')
            outerr_str = outerr.decode('utf-8')
            self._logger.log(_logging.DEBUG, output_str)
            if len(outerr_str) > 0:
                self._logger.log(_logging.DEBUG, outerr_str)
            else:
                self._logger.log(_logging.DEBUG, 'No datum transformation errors reported')
        except:
            print(output)
            print(outerr)

    def _xyz2gdal(self, outxyz: List[Tuple[float, float, float]], out_zone: int, out_verdat: str) -> gdal.Dataset:
        """
        Convert from numpy xyz array to a gdal dataset.

        Parameters
        ----------
        outxyz :
            param out_zone:
        out_verdat :
            
        outxyz: List[Tuple[float :
            
        float :
            
        float]] :
            
        out_zone: int :
            
        out_verdat: str :
            

        Returns
        -------

        """

        # setup the gdal bucket
        dest = osr.SpatialReference()
        dest.SetWellKnownGeogCS('NAD83')
        if int(out_zone) < 0:  # if out_zone is positive it is in the northern hemisphere
            hemisphere = '0'  # 0 = South
        else:  # if out_zone is negative it is in the southern hemisphere
            hemisphere = '1'  # 1 =North#boolean test being passed to SetUTM
        dest.SetUTM(int(out_zone), int(hemisphere))
        dest.SetVertCS(out_verdat, out_verdat, 2000)
        dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = dataset.CreateLayer('pts', dest, geom_type=ogr.wkbPoint)
        for p in outxyz:
            newp = ogr.Geometry(ogr.wkbPoint)
            newp.AddPoint(p[0], p[1], p[2])
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(newp)
            layer.CreateFeature(feature)
        return dataset

    def _have_instructions(self, instructions):
        """Confirm the existance of required information for datum transformations."""
        ready = True
        req = from_hdatum + from_vdatum + to_hdatum + to_vdatum
        for key in req:
            if key not in instructions:
                ready = False
                break
        return ready