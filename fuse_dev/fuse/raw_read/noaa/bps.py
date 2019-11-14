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

import logging as _logging
import lxml.html as _html
import numpy as _np
import os as _os
import sys as _sys

from glob import glob as _glob

from fuse.raw_read.raw_read import RawReader
from fuse.raw_read import parse_file_pickle

h93_to_meta = {'Surv Id': 'survey',
               'Strt Yr': 'start_date',
               'End Yr': 'end_date',
               'Horizontal Datum of Records': 'from_horiz_datum',
               'Vertical Datum': 'from_vert_datum',
               }

horiz_datum = {
        'NORTH AMERICAN DATUM OF 1983': 'NAD83',
        'NORTH AMERICAN DATUM 1983': 'NAD83'
        }

vert_datum = {
    'MEAN LOW WATER SPRINGS': 'MLWS',
    'MEAN LOWER LOW WATER SPRINGS': 'MLLWS',
    'MEAN SEA LEVEL': 'MSL',
    'LOWER LOW WATER': 'LLW',
    'MEAN LOW WATER': 'MLW',
    'INDIAN SPRING LOW WATER': 'ISLW',
    'MEAN LOWER LOW WATER': 'MLLW',
    'MEAN HIGH WATER': 'MHW',
    'MEAN HIGH WATER SPRINGS': 'MHWS',
    'MEAN HIGHER HIGH WATER': 'MHHW',
    'LOWEST ASTRONOMICAL TIDE': 'LAT',
    'INTERNATION GREAT LAKES DATUM': 'IGLD',
    'LOWER LOW WATER LARGE TIDE': 'LLWLT',
    'HIGHER HIGH WATER LARGE TIDE': 'HHWLT',
    'HIGHEST ASTRONOMICAL TIDE': 'HAT',
    'HUDSON RIVER DATUM': 'HRD',  # Adding this for the Hudson River Datum
    }

bps_spec = r"https://www.ngdc.noaa.gov/mgg/dat/geodas/docs/hyd93.htm"

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

    def read_metadata(self, survey_folder: str) -> dict:
        """

        Parameters
        ----------
        filename
            file to read

        Returns
        -------
            dictionary of metadata
        """
        filename = self._data_determination(survey_folder)
        meta_h93 = self._parse_htm(filename)
        meta_pickle = parse_file_pickle.read_pickle(filename)
        meta_default = self._defualt_meta(filename)
        meta_combined = {**meta_default, **meta_h93, **meta_pickle}
        meta_final = self._finalize_meta(meta_combined)
        return meta_final

    def _xyz_precedence(self, file_list: [str]) -> (str, bool):
        xyz_scores = {}
        file_list = [xyz for xyz in file_list if _os.path.splitext(xyz)[1].lower() in ('.xyz', '.a93')]
        for xyz in file_list:
            ext = _os.path.splitext(xyz)[1].lower()
            if ext == '.a93':
                xyz_scores[2] = xyz
            elif ext == '.xyz':
                xyz_scores[1] = xyz

        max_score = max(list(xyz_scores.keys()))
        return xyz_scores[max_score]

    def _data_determination(self, survey_folder: str) -> str:
        """
        Evaluates files present to determine file passed back for metadata reading

        Parameters
        ----------
        survey_folder
            Complete filepath of the input data

        Returns
        -------
        str
            The filename of the data chosen data source

        """
        data_files = _glob(_os.path.join(survey_folder, '*xyz', '*a93'))

        return self._xyz_precedence(data_files)

    def _parse_htm(self, filename: str) -> dict:
        fileroot, filename = _os.path.split(filename)
        filebase, fileext = _os.path.splitext(filename)
        file_htm = _os.path.join(fileroot, f"{filebase}_h93.htm")

        metadata = {}

        if _os.path.isfile(file_htm):
            body = _html.parse(file_htm).getroot()[1]
            for table in body[1:-2]:
                columns = range(len([header for header in table[0]]))
                for index in columns:
                    try:
                        header = table[0][index].text.strip() if table[0][index].text is not None else None
                        value = table[1][index][0].text.strip()
                        if header is None or header.lower() in ('', 'seq', 'c', 'code'):
                            pass
                        elif header in ('Surv Id'):
                            metadata[h93_to_meta[header]] = value
                        elif header in ('Strt Yr', 'End Yr'):
                            metadata[h93_to_meta[header]] = f"{value}0101"
                            metadata['year'] = int(value)
                        elif header in ('Horizontal Datum of Records'):
                            metadata[h93_to_meta[header]] = value
                            if value.upper() in horiz_datum:
                                metadata['from_vert_key'] = horiz_datum[value.upper()]
                            else:
                                self._logger.warning(f"Unrecognized horizontal datum: {value}")
                        elif header in ('Vertical Datum'):
                            metadata[h93_to_meta[header]] = value
                            if value.upper() in vert_datum:
                                metadata['from_vert_key'] = vert_datum[value.upper()]
                            else:
                                self._logger.warning(f"Unrecognized vertical datum: {value}")
                        else:
                            metadata[header] = value
                    except IndexError as e:
                        self._logger.warning(f"{e}: {index}")
        else:
            self._logger.warning(f"Unable to find {file_htm} to read {filename} metadata")
        return metadata


    def _defualt_meta(self, filename) -> dict:
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
        metadata['from_filename'] = name
        metadata['from_path'] = filename
        # source info
        metadata['agency'] = 'US'
        metadata['source_indicator'] = 'US,US,graph,' + name
        metadata['source_type'] = ext
        metadata['license'] = 'CC0'
        # processing info
        metadata['interpolate'] = False

        return metadata

    def _finalize_meta(self, metadata):
        """
        Update the metadata to standard values.

        Parameters
        ----------
        meta
            The combined metadata dictionary to be updated with standard
            values.

        Returns
        -------
        dict
            The final metadata for return to the metadata requesting method
        """
        if 'year' in metadata:
            if metadata['year'] >= 1940:
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
                metadata['feat_size'] = 9999
                metadata['feat_detect'] = False
                metadata['feat_least_depth'] = False
            else:
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
                metadata['feat_size'] = 9999
                metadata['feat_detect'] = False
                metadata['feat_least_depth'] = False
        else:
            self._logger.warning(f"Could not finalize metadata due to missing 'year' assignment")
        return metadata

    def _format_a93_z(self, z_string: str) -> str:
        negative = False
        if '-' in z_string:
            z_string = z_string[1:]
            negative = True
        if len(z_string) == 6 and z_string.count('9') == 6:
            z_string = z_string
        elif len(z_string) > 2:
            z_string = f"{z_string[:2]}.{z_string[2:]}"
        if negative:
            z_string = f"-{z_string}"
        return z_string

    def _parse_a93(self, infilename: str) -> _np.array:
        """
        Opens and parse a93 data to return points marked as soundings

        Parameters
        ----------
        infilename : str
            Data filepath

        Returns
        -------
        _np.array
            xyz data
        """
        xyz = []
        a93 = open(infilename, mode='r')
        non_sounding_found = False
        line = 0
        for row in a93:
            snum = row[:8].strip()
            y = float(f"{row[8:11]}.{row[11:17]}".strip())
            x = float(f"{row[17:21]}.{row[21:26]}".strip())
            pre_format_z = f"{row[27:33]}".strip()
            z = float(self._format_a93_z(pre_format_z))
            depth_code, carto_code = float(row[33].strip()), float(row[34:].strip())
            if z == 999999 or carto_code != 711:
                non_sounding_found = True
                _logging.warning(f"Non-sounding data. Line {line} --> {x, y, z, carto_code}")
            else:
                xyz.append([x, y, z])
            line += 1
        if non_sounding_found:
            _logging.info(f"For more info on non-sounding points, see APPENDIX A: {bps_spec}")
        a93.close()
        return _np.array(xyz)

    def _parse_xyz(self, infilename: str) -> _np.array:
        """
        Opens and filters xyz data to only return active soundings

        Parameters
        ----------
        infilename : str
            Data filepath

        Returns
        -------
        _np.array
            xyz data

        """
        xyz = _np.loadtxt(infilename, delimiter=',', skiprows=1,
                          usecols=(1, 2, 3, 5))
        xyz[:, [0, 1]] = xyz[:, [1, 0]]
        active_bool = (xyz[:, 3] == 0)
        active = xyz[active_bool == 0, :]
        return active[:, [0, 1, 2]]

    def _data_lookup(self, infilename: str) -> [str]:
        data_files = [infilename]
        base, ext = _os.path.splitext(infilename)
        if ext == '.a93':
            xyz = f"{base}.xyz"
            if _os.path.isfile(xyz):
                data_files.append(xyz)
        elif ext == '.xyz':
            a93 = f"{base}.a93"
            if _os.path.isfile(a93):
                data_files.append(a93)
        return data_files

    def _data_check(self, infilename: str) -> _np.array:
        """
        Reads in all available data files and returns the data
        """

        data_files = self._data_lookup(infilename)

        data_hold = {}
        for data in data_files:
            ext = _os.path.splitext(data)[1]
            if ext == '.a93':
                data_hold[ext] = self._parse_a93(data)
            elif ext == '.xyz':
                data_hold[ext] = self._parse_xyz(data)

        return data_hold

    def read_bathymetry(self, infilename: str) -> _np.array:
        """
        Returns a numpy array

        Parameters
        ----------
        infilename : str
            data filepath

        Returns
        -------
        _np.array
            numpy array
        """


        return self._data_check(infilename)
