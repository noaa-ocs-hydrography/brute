# -*- coding: utf-8 -*-
"""
use_vdatum.py
Created on Wed Aug 22 12:27:39 2018

@author: grice

Use VDatum for conversions.
"""

__version__ = 'use_vdatum 0.0.1'

import logging as _logging
import os as _os
import subprocess as _subprocess
from tempfile import TemporaryDirectory as tempdir

import numpy as _np
from osgeo import gdal, ogr, osr

from_hdatum = [
    'from_horiz_frame',
    'from_horiz_type',
    'from_horiz_units',
    'from_horiz_key'
]

from_vdatum = [
    'from_vert_key',
    'from_vert_units',
    'from_vert_direction'
]

to_hdatum = [
    'to_horiz_frame',
    'to_horiz_type',
    'to_horiz_units'
]

to_vdatum = [
    'to_vert_key',
    'to_vert_units',
    'to_vert_direction'
]


class VDatum:
    """An object for working with VDatum."""

    def __init__(self, vdatum_path: str, java_path: str, reader):
        """
        Create a new object for using VDatum.

        Parameters
        ----------
        config
            dictionary of configuration
        reader
            reader object
        """

        self._reader = reader

        if _os.path.isfile(_os.path.join(vdatum_path, 'vdatum.jar')):
            self._vdatum_path = vdatum_path
        else:
            raise ValueError(f'Invalid vdatum folder: {vdatum_path}')

        if _os.path.isfile(_os.path.join(java_path, 'java.exe')):
            self._java_path = java_path
        else:
            raise ValueError(f'Invalid java path: {java_path}')

        self._logger = _logging.getLogger('fuse')

    def translate(self, filename: str, instructions: dict) -> gdal.Dataset:
        """
        Translate the provided filename from the provided in datums to the out
        datums and return a gdal object.

        NSRS2007 is assumed for the out EPSG code.

        Parameters
        ----------
        filename
            filename of data file
        instructions
            dictionary of metadata

        Returns
        -------
            GDAL point cloud dataset
        """

        self._logger.log(_logging.DEBUG, 'Begin datum transformation')
        if not _has_required_instructions(instructions):
            instructions['interpolate'] = 'False'
            raise ValueError('The required fields for transforming datums are not available')
        if instructions['read_type'].upper() == 'BAG':
            dataset, utm_zone = self.__translate_bag(filename, instructions)
        else:
            points, utm_zone = self.__translate_xyz(filename, instructions)
            vertical_datum = instructions['to_vert_key'].upper()
            # passing UTM zone instead of EPSG code
            dataset = self.__xyz2gdal(points, utm_zone, vertical_datum)
        self._logger.log(_logging.DEBUG, 'Datum transformation complete')
        return dataset

    def create(self, filename: str, instructions: dict) -> gdal.Dataset:
        """
        Get a GDAL point cloud dataset from the given file.

        Parameters
        ----------
        filename
            filename of data file
        instructions
            dictionary of metadata

        Returns
        -------
            GDAL point cloud dataset
        """

        if instructions['read_type'] == 'ehydro':
            return self.__read_points(filename, instructions)
        elif instructions['read_type'] == 'bag':
            return self.read_bathymetry(filename, instructions)
        else:
            raise ValueError('Reader type not implemented')

    def __read_points(self, filename: str, instructions: dict) -> gdal.Dataset:
        """
        Pass back a un-translated GDAL point cloud dataset

        Parameters
        ----------
        filename
            filename of data file
        instructions
            dictionary of metadata

        Returns
        -------
            GDAL point cloud dataset
        """

        points = self._reader.read_bathymetry(filename)

        if 'to_horiz_key' in instructions:
            utm_zone = int(instructions['from_horiz_key'])
        if 'from_vert_key' in instructions:
            vertical_datum = instructions['from_vert_key']

        return self.__xyz2gdal(points, utm_zone, vertical_datum)

    def __read_bag_bathy(self, filename: str, instructions: dict) -> gdal.Dataset:
        """
        Create a GDAL point cloud dataset from the given BAG file.

        Parameters
        ----------
        filename
            filename of BAG
        instructions
            dictionary of metadata

        Returns
        -------
            GDAL point cloud dataset
        """

        return self._reader.read_bathymetry(filename, instructions['to_vert_key'])

    def __translate_xyz(self, filename: str, instructions: dict) -> (_np.array, int):
        """
        Reproject XYZ from the given filename.

        Parameters
        ----------
        filename
            filename of XYZ points
        instructions
            dictionary of metadata

        Returns
        -------
            array of XYZ points and UTM zone
        """

        # read the points and put it in a temp file for VDatum to read
        points = self._reader.read_bathymetry(filename)
        original_directory = tempdir()
        output_filename = _os.path.join(original_directory.name, 'outfile.txt')
        reprojected_directory = tempdir()
        reprojected_filename = _os.path.join(reprojected_directory.name, 'outfile.txt')
        log_filename = _os.path.join(reprojected_directory.name, 'outfile.txt.log')
        # save point array to file
        _np.savetxt(output_filename, points, delimiter=',')
        # translate points with VDatum
        self.__setup_vdatum(instructions, 'points')
        self.__convert_file(output_filename, reprojected_directory.name)
        if 'to_horiz_key' in instructions:
            utm_zone = int(instructions['to_horiz_key'])
        else:
            # read out UTM Zone from VDatum log file
            with open(log_filename, 'r') as logfile:
                for line in logfile.readlines():
                    if line.startswith('Zone:'):
                        # input_zone = line[27:53].strip()
                        utm_zone = int(line[54:82].strip())
                        break
                else:
                    raise ValueError(f'no UTM zone found in file "{filename}"')
        return _np.loadtxt(reprojected_filename, delimiter=','), utm_zone

    def __translate_bag(self, filename: str, instructions: dict) -> (gdal.Dataset, int):
        """
        Reproject BAG bathy from the given filename.

        Parameters
        ----------
        filename
            filename of the BAG
        instructions
            dictionary of metadata

        Returns
        -------
            array of XYZ points and UTM zone
        """

        # create a gdal.Dataset and assign a temp file name for VDatum to read
        dataset = self._reader.read_bathymetry(filename, out_verdat=None)
        original_directory = tempdir()
        output_filename = _os.path.join(original_directory.name, 'outfile.tif')
        reprojected_directory = tempdir()
        reprojected_filename = _os.path.join(reprojected_directory.name, 'outfile.tif')
        log_filename = _os.path.join(reprojected_directory.name, 'outfile.tif.log')
        # save dataset to file
        driver = gdal.GetDriverByName('GTiff')
        driver.CreateCopy(output_filename, dataset)
        del driver
        # translate points with VDatum
        self.__setup_vdatum(instructions, 'geotiff')
        self.__convert_file(output_filename, reprojected_directory.name)
        if 'to_horiz_key' in instructions:
            utm_zone = int(instructions['to_horiz_key'])
        else:
            # read out UTM Zone from VDatum log file
            with open(log_filename, 'r') as logfile:
                for line in logfile.readlines():
                    if line.startswith('Zone:'):
                        # input_zone = line[27:53].strip()
                        utm_zone = int(line[54:82].strip())
                        break
                else:
                    raise ValueError(f'no UTM zone found in file "{filename}"')
        return gdal.Open(reprojected_filename), utm_zone

    def __setup_vdatum(self, instructions: dict, mode: str = 'points'):
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
        instructions
            dictionary of metadata
        """

        # having the output zone is optional
        local_to_hdatum = to_hdatum.copy()
        if 'to_horiz_key' in instructions:
            local_to_hdatum.append('to_horiz_key')
        self._shell = f'{_os.path.join(self._java_path, "java")} -jar vdatum.jar ' + \
                      f'ihorz:{":".join(instructions[key] for key in from_hdatum)} ' + \
                      f'ivert:{":".join(instructions[key] for key in from_vdatum)} ' + \
                      f'ohorz:{":".join(instructions[key] for key in local_to_hdatum)} ' + \
                      f'overt:{":".join(instructions[key] for key in to_vdatum)} '
        if mode == 'points':
            self._shell += f'-file:txt:comma,0,1,2,skip0:'
        elif mode == 'geotiff':
            self._shell += f'-file:geotiff:geotiff:'

    def __convert_file(self, filename: str, output_directory: str):
        """
        The provided file of xyz points will be converted and returned
        according to the from and to the datums provided through the
        _setup_vdatum_point_converstion method.

        Parameters
        ----------
        filename
            filename of XYZ points file
        output_directory
            output directory
        """

        command = f'{self._shell}{filename};{output_directory}'
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

    def __xyz2gdal(self, points: [(float, float, float)], utm_zone: int, vertical_datum: str) -> gdal.Dataset:
        """
        Get a GDAL point cloud dataset from XYZ points.

        Parameters
        ----------
        points
            XYZ points
        utm_zone
            desired UTM zone
        vertical_datum
            name of desired vertical datum

        Returns
        -------
            GDAL point cloud dataset
        """

        # setup the gdal bucket
        spatial_reference = osr.SpatialReference()
        spatial_reference.SetWellKnownGeogCS('NAD83')
        # positive UTM zone is in the northern hemisphere
        spatial_reference.SetUTM(abs(utm_zone), 1 if utm_zone > 0 else 0)
        spatial_reference.SetVertCS(vertical_datum, vertical_datum, 2000)
        dataset = gdal.GetDriverByName('Memory').Create('', 0, 0, 0, gdal.GDT_Unknown)
        layer = dataset.CreateLayer('pts', spatial_reference, geom_type=ogr.wkbPoint)
        for point in points:
            geometry = ogr.Geometry(ogr.wkbPoint)
            geometry.AddPoint(*point)
            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(geometry)
            layer.CreateFeature(feature)
        return dataset


def _has_required_instructions(instructions: dict) -> bool:
    """
    Confirm the existance of required information for datum transformations.

    Parameters
    ----------
    instructions
        dictionary of metadata
    """

    for key in from_hdatum + from_vdatum + to_hdatum + to_vdatum:
        if key not in instructions:
            return False
    else:
        return True
