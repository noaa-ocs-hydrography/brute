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

from fuse.raw_read.raw_read import RawReader

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

        Parameters
        ----------
        filename
            file to read

        Returns
        -------
            dictionary of metadata
        """

        meta_h93 = self._parse_htm(filename)
        meta_default = self._defualt_meta(filename)
        meta_combined = {**meta_default, **meta_h93}
        meta_final = self._finalize_meta(meta_combined)
        return meta_final

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
                         usecols = (1, 2, 3, 5))
        xyz[:, [0, 1]] = xyz[:, [1, 0]]
        active_bool = (xyz[:, 3] == 0)
        active = xyz[active_bool == 0, :]
        return active[:, [0, 1, 2]]
