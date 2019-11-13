# -*- coding: utf-8 -*-
"""
use_vdatum.py
Created on Wed Aug 22 12:27:39 2018

@author: grice

Use VDatum for conversions.
"""

from fuse.datum_transform.use_gdal import _xyz_to_gdal, spatial_reference_from_metadata

__version__ = 'use_vdatum 0.0.1'

import logging as _logging
import os as _os
import subprocess as _subprocess
from tempfile import TemporaryDirectory

import numpy as _np
from osgeo import gdal

FROM_HDATUM = [
    'from_horiz_frame',
    'from_horiz_type',
    'from_horiz_units'
]

FROM_VDATUM = [
    'from_vert_key',
    'from_vert_units',
    'from_vert_direction'
]

TO_HDATUM = [
    'to_horiz_frame',
    'to_horiz_type',
    'to_horiz_units'
]

TO_VDATUM = [
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

        TODO rename to reproject

        Parameters
        ----------
        filename
            filename of data file
        instructions
            dictionary of metadata

        Returns
        -------
        gdal.Dataset
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
            dataset = _xyz_to_gdal(points, utm_zone, vertical_datum)
        self._logger.log(_logging.DEBUG, 'Datum transformation complete')
        return dataset

    def __translate_xyz(self, filename: str, instructions: dict) -> (_np.array, int):
        """
        Reproject XYZ from the given filename.

        TODO rename to __reproject_xyz

        Parameters
        ----------
        filename
            filename of XYZ points
        instructions
            dictionary of metadata

        Returns
        -------
        numpy.array, int
            N x 3 array of XYZ points and UTM zone
        """

        # read the points and put it in a temp file for VDatum to read
        points = self._reader.read_bathymetry(filename)
        points = self.__filter_xyz(points, instructions)

        # save point array to file
        points_file_directory = TemporaryDirectory()
        points_filename = _os.path.join(points_file_directory.name, 'outfile.txt')
        _np.savetxt(points_filename, points, delimiter=',')

        # translate points with VDatum
        reprojected_directory = TemporaryDirectory()
        self.__setup_vdatum(instructions, 'points')
        self.__convert_file(points_filename, reprojected_directory.name)
        reprojected_filename = _os.path.join(reprojected_directory.name, 'outfile.txt')
        log_filename = _os.path.join(reprojected_directory.name, 'outfile.txt.log')

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

        if not _os.path.isfile(reprojected_filename):
            raise RuntimeError(f'VDatum was unable to reproject survey from {filename} to {reprojected_filename}')

        return _np.loadtxt(reprojected_filename, delimiter=','), utm_zone

    def __filter_xyz(self, xyz: _np.array, instructions: dict) -> _np.array:
        """
        Filter the given geographic XYZ points to exclude points outside the UTM zone of the given spatial reference frame.
        If the provided data is not geographic, or the given frame is not a UTM zone, no filter is applied.

        Parameters
        ----------
        xyz
            N x 3 array of XYZ points
        instructions
            dictionary of metadata defining geographic frame

        Returns
        ----------
        numpy.array
            N x 3 array of XYZ points, filtered to only include the given UTM zone
        """

        if instructions['from_horiz_type'] == 'geo' and instructions['to_horiz_type'] == 'utm':
            srs = spatial_reference_from_metadata(instructions)
            c_meridian = srs.GetProjParm(gdal.osr.SRS_PP_CENTRAL_MERIDIAN)
            west = c_meridian - 3
            east = c_meridian + 3
            x = xyz[:, 0]
            xyz = xyz[(x > west) & (x < east), :]

        return xyz

    def __translate_bag(self, filename: str, instructions: dict) -> (gdal.Dataset, int):
        """
        Reproject BAG bathy from the given filename.

        TODO rename to __reproject_bag

        Parameters
        ----------
        filename
            filename of the BAG
        instructions
            dictionary of metadata

        Returns
        -------
        numpy.array, int
            N x 3 array of XYZ points and UTM zone
        """

        # create a gdal.Dataset and assign a temp file name for VDatum to read
        dataset = self._reader.read_bathymetry(filename, out_verdat=None)
        original_directory = TemporaryDirectory()
        output_filename = _os.path.join(original_directory.name, 'outfile.tif')
        reprojected_directory = TemporaryDirectory()
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

        # having the input zone is optional if the input type is geographic
        local_from_hdatum = FROM_HDATUM.copy()
        if instructions['from_horiz_type'] != 'geo':
            local_from_hdatum.append('from_horiz_key')

        # having the output zone is optional
        local_to_hdatum = TO_HDATUM.copy()
        if 'to_horiz_key' in instructions:
            local_to_hdatum.append('to_horiz_key')

        self._shell = f'{_os.path.join(self._java_path, "java")} -jar {_os.path.join(self._vdatum_path, "vdatum.jar")} ' + \
                      f'ihorz:{":".join(instructions[key] for key in local_from_hdatum)} ' + \
                      f'ivert:{":".join(instructions[key] for key in FROM_VDATUM)} ' + \
                      f'ohorz:{":".join(instructions[key] for key in local_to_hdatum)} ' + \
                      f'overt:{":".join(instructions[key] for key in TO_VDATUM)} '

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
            output, outerr = proc.communicate()
            self._logger.log(_logging.DEBUG, output.decode('utf-8'))
            if len(outerr) > 0:
                self._logger.log(_logging.DEBUG, outerr.decode('utf-8'))
            else:
                self._logger.log(_logging.DEBUG, 'No datum transformation errors reported')
        except:
            print(output)
            print(outerr)


def _has_required_instructions(instructions: dict) -> bool:
    """
    Confirm the existance of required information for datum transformations.

    Parameters
    ----------
    instructions
        dictionary of metadata
    """

    for key in FROM_HDATUM + FROM_VDATUM + TO_HDATUM + TO_VDATUM:
        if key not in instructions:
            return False
    else:
        return True
