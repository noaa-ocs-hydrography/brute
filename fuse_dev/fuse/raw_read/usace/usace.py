# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:49:10 2019

@author: Casiano.Koprowski
"""

import logging as _logging
import os as _os
import pickle as _pickle
import re as _re
import sys as _sys
from datetime import datetime as _datetime
from glob import glob as _glob
from typing import Union
from xml.etree.ElementTree import parse as _parse

import numpy as _np
from fuse.datum_transform import usefips as _usefips
from fuse.raw_read.raw_read import RawReader

from fuse.raw_read import parse_file_pickle
from fuse.raw_read.usace import parse_usace_xml

_ehydro_quality_metrics = {
    'complete_coverage': False,
    'bathymetry': True,
    'vert_uncert_fixed': 0.5,
    'vert_uncert_vari': 0.1,
    'horiz_uncert_fixed': 5.0,
    'horiz_uncert_vari': 0.05,
    'feat_detect': False,
}

class USACERawReader(RawReader):
    def __init__(self, district: str = None):
        self.district = district
        self.survey_feet_per_meter = 0.30480060960121924  # US survey feet to meters
        self.xyz_suffixes = ('_A', '_FULL')
        self.xyz_files = {}

        self._logger = _logging.getLogger(f'fuse')

        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            ch.setLevel(_logging.DEBUG)
            self._logger.addHandler(ch)

    def read_metadata(self, survey_folder: str) -> dict:
        """
        Read all available meta data.

        The USACE metadata is retuned in order of precedence:
            1. The survey's ``.xml``.
            2. The file name.
            3. The survey's ``.xyz`` header.
            4. The metadata pickle pulled from eHydro.

        If the data are to be interpolated two dictionaries are returned, one
        that represents the interpolated dataset and one that is not.

        Parameters
        ----------
        survey_folder
            Folder path of the input ``.xyz`` data

        Returns
        -------
        list
            The complete metadata pulled from multiple sources within the
            survey as a dict within a list.
        """

        meta_supplement = {}
        meta_determine, filename = self._data_determination(meta_supplement, survey_folder)
        meta_xml = self._parse_usace_xml(filename)
        meta_xyz = self._parse_ehydro_xyz_header(filename)
        meta_filename = self._parse_filename(filename)
        meta_pickle = self._parse_pickle(filename)
        meta_date = self._parse_start_date(filename, {**meta_pickle, **meta_xyz, **meta_xml})
        meta_supplement = {**meta_determine, **meta_date, **meta_supplement}
        meta_combined = {**meta_pickle, **meta_xyz, **meta_filename, **meta_xml, **meta_supplement}
        meta_final = self._finalize_meta(meta_combined)
        if meta_final['interpolate']:
            meta_orig = meta_final.copy()
            meta_orig['interpolate'] = False
            meta_final['from_filename'] = f"{meta_orig['from_filename']}.interpolate" 
            return [meta_orig, meta_final]
        else:
            return [meta_final]

    def read_bathymetry(self, filename: str) -> _np.array:
        """
        Read the bathymetry and return an array of the xyz points.

        Parameters
        ----------
        filename
            filename of XYZ file

        Returns
        -------
        numpy.array
            A 2D array containing the complete list of points found in the
            input file
        """

        bathy = self._parse_ehydro_xyz_bathy(filename)
        return bathy

    def read_bathymetry_by_point(self, filename: str) -> _np.array:
        """
        Get the next point in the given bathymetry.

        Parameters
        ----------
        filename
            Complete filepath of the input data

        Yields
        ------
        numpy.array
            XYZ point
        """

        for point in self.read_bathymetry(filename):
            yield point

    def _parse_pickle(self, filename: str) -> dict:
        """
        Retrieves metadata for USACE E-Hydro files
        function returns metadata dictionary

        Parameters
        ----------
        filename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata returned via this method
        """

        root = _os.path.split(filename)[0]
        pickle_name = self.name_gen(filename, ext='.pickle', sfx=False)
        if _os.path.isfile(pickle_name):
            pickle_dict = parse_file_pickle.read_pickle(pickle_name, pickle_ext=True)
            if 'poly_name' in pickle_dict:
                pickle_dict['support_files'] = [_os.path.join(root, pickle_dict['poly_name'])]
        else:
            pickle_dict = {}
        return pickle_dict

    def _parse_usace_xml(self, filename: str):
        """
        Read all available meta data.
        returns dictionary

        Parameters
        ----------
        filename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata assigned via this method
        """

        xmlfilename = self.name_gen(filename, ext='.xml', sfx=False)
        if _os.path.isfile(xmlfilename):
            with open(xmlfilename, 'r') as xml_file:
                xml_txt = xml_file.read()
            xmlbasename = _os.path.basename(xmlfilename)
            xml_data = parse_usace_xml.XMLMetadata(xml_txt, filename=xmlbasename)
            if xml_data.version == 'USACE_FGDC':
                meta_xml = xml_data._extract_meta_USACE_FGDC()
            elif xml_data.version == 'ISO-8859-1':
                meta_xml = xml_data._extract_meta_USACE_ISO()
                if 'ISO_xml' not in meta_xml:
                    meta_xml = xml_data._extract_meta_USACE_FGDC(override='Y')
            else:
                meta_xml = xml_data.convert_xml_to_dict2()
            ext_dict = xml_data.extended_xml_fgdc()
            ext_dict = parse_usace_xml.ext_xml_map_enddate(ext_dict)
            meta_xml = parse_usace_xml.xml_SPCSconflict_flag(meta_xml)
        else:
            ext_dict = {}
            meta_xml = {}
        meta_xml['from_path'] = filename
        meta_xml['from_filename'] = _os.path.basename(filename)
        return {**meta_xml, **ext_dict}

    def _parse_start_date(self, infilename: str, metadata: dict) -> dict:
        """
        Reads the start data metadata avaiable via formatted filename.

        Parameters
        ----------
        infilename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata assigned via this method
        """

        start = {}
        file_date = None
        xml_date = None
        pickle_date = None
        filename_dict = self._parse_filename(infilename)

        if 'start_date' in filename_dict:
            file_date = filename_dict['start_date']

        if 'start_date' in metadata:
            xml_date = metadata['start_date']

        if 'SURVEYDATESTART' in metadata:
            pickle_date = metadata['SURVEYDATESTART']

        if xml_date is not None:
            start = {'start_date': xml_date}
        elif file_date is not None:
            start = {'start_date': file_date}
        elif pickle_date is not None:
            start = {'start_date': pickle_date}

        return start

    def name_gen(self, filename: str, ext: str = None, sfx: bool = True) -> str:
        """
        Returns the suffix of the a survey's xyz file and the file name
        for a different extension

        If the survey's xyz file name contains ``_A.xyz`` or ``_FULL.xyz`` this
        suffix is removed from the 'base' name of the file.

        Parameters
        ----------
        filename
            Input file
        ext
            New extention to be applied to the base
        sfx
            Whether or not the function passes back the suffix found in
            `filename`

        Returns
        -------
        type :
            ``tuple(base, suffix)`` if ``sfx`` is ``True``;
            ``base`` if ``sfx`` is ``False``
        """

        filebase, fileext = _os.path.splitext(filename)
        suffix = None

        for item in self.xyz_suffixes:
            if _re.compile(f'{item}$', _re.IGNORECASE).search(filebase):
                suffix = item

        if suffix is not None and ext is not None:
            if fileext.lower() == ext.lower():
                ext = fileext
            base = _re.sub(_re.compile(f'{suffix}$', _re.IGNORECASE), '', filebase) + f'{ext}'
        elif suffix is None and ext is not None:
            if fileext.lower() == ext.lower():
                ext = fileext
            base = filebase + f'{ext}'
        else:
            base = filename

        if sfx:
            return base, suffix
        else:
            return base

    def _xyz_precedence(self, file_list: [str]) -> (str, bool):
        xyz_scores = {}
        file_list = [xyz for xyz in file_list if _os.path.splitext(xyz)[1].lower() == '.xyz']
        for xyz in file_list:
            xyz_upper = xyz.upper()
            if '_FULL.' in xyz_upper:
                xyz_scores[3] = xyz
            elif '_A.' in xyz_upper:
                xyz_scores[2] = xyz
            else:
                xyz_scores[1] = xyz

        max_score = max(list(xyz_scores.keys()))
        return xyz_scores[max_score], True if max_score != 3 else False

    def _data_determination(self, meta_dict: dict, survey_folder: str) -> dict:
        """
        Determines cerain metadata values based on the known quality of the
        data.

        Currently, this function checks the ``.xyz`` files for ``_A`` or
        ``_FULL`` suffixes.  If one of these are found
        ``metadict['interpolate']`` is set to False; True otherwise.

        Parameters
        ----------
        meta_dict
            Dictionary to add values to
        survey_folder
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata assigned via this method
        """

        if _os.path.isfile(survey_folder):
            survey_folder = _os.path.dirname(survey_folder)

        xyz_files = _glob(_os.path.join(survey_folder, '*xyz'))
        if len(xyz_files) < 1:
            raise RuntimeError('No files to process in {survey_folder}')
        filename, meta_dict['interpolate'] = self._xyz_precedence(xyz_files)

        meta_dict['file_size'] = self._size_finder(filename)
        meta_dict['from_filename'] = self.name_gen(_os.path.split(filename)[1], '', sfx=None)
        meta_dict['from_path'] = filename

        return meta_dict, filename

    def _size_finder(self, filepath: Union[str, _os.PathLike]) -> int:
        """
        Returns the rounded size of a file in MB as an integer

        Parameters
        ----------
        filepath
            path of file

        Returns
        -------
        int
            size in MB
        """

        return int(_np.round(_os.path.getsize(filepath) / 1000))

    def _check_grid(self, infilename):
        #        data = self._parse_ehydro_xyz_bathy(infilename)
        ...

    def _parse_ehydro_xyz_header(self, filename: str, meta_source: str = 'xyz', default_meta: str = '') -> dict:
        """
        Parse an USACE eHydro file for the available meta data.

        Default metadata (values predetermined for the file but not in the file)
        can be stored at the location defined by 'default_meta' as a pickled
        dicitonary.  If no path is provided the dictionary in the same folder as
        the data but in the file 'default.pkl' will be used.  If this file does
        not exist no default metadata will be loaded.  If the same keyword for
        the metadata exists both in the file metadata and in the default location,
        the file metadata will take precidence.

        Parameters
        ----------
        infilename
            Complete filepath of the input data
        meta_source
            Choice of ``xyz`` or ``xml`` (Default == 'xyz')
        default_meta
            (Default == ''); Optional filepath of a default metadata ``dict``
            stored as a ``.pkl`` file

        Returns
        -------
        dict
            The metadata found via this method
        """
        meta = {}
        for file_ext in ('.xyz', '_A.xyz', '_FULL.xyz'):
            infilename, suffix = self.name_gen(filename, ext=file_ext)
            if not _os.path.isfile(infilename):
                continue
            name_sections = ('projid', 'uniqueid', 'subprojid', 'start_date',
                             'statuscode', 'optional')
            name_meta = self._parse_filename(infilename)
            sorind = '_'.join([name_meta[key] for key in name_sections if key in name_meta])
            name_meta['source_indicator'] = f'US,US,graph,{sorind}'
            if meta_source == 'xyz':
                file_meta = self._parse_xyz_header(infilename, mode=self.district)
            elif meta_source == 'xml':
                file_meta = self._parse_ehydro_xml(infilename)
            default_meta = self._load_default_metadata(infilename, default_meta)
            merged_meta = {**default_meta, **name_meta, **file_meta}
            if 'from_horiz_unc' in merged_meta:
                if merged_meta['from_horiz_units'] == 'US Survey Foot':
                    val = self.survey_feet_per_meter * float(merged_meta['from_horiz_unc'])
                    merged_meta['horiz_uncert'] = val
            if 'from_vert_unc' in merged_meta:
                if merged_meta['from_vert_units'] == 'US Survey Foot':
                    val = self.survey_feet_per_meter * float(merged_meta['from_vert_unc'])
                    merged_meta['vert_uncert_fixed'] = val
                    merged_meta['vert_uncert_vari'] = 0
            # merged_meta['script_version'] = __version__
            meta = {**meta, **merged_meta}
        return meta

    def _parse_filename(self, infilename: str) -> dict:
        """
        Parse the provided infilename for the channel project code, unique id,
        subproject code, survey acquistion start date, the survey code, and
        optional field and return a dictionary of these fields.  The dictionary
        contains the following keys:
            projid : the project id code
            uniqueid : the unique identifier, and can represent items such as
                        separate projects, map scale, or CWIS code.
            subprojid : the subproject within the channel project
            startdate : the start date of acquisition
            statuscode : survey description code, such as for a condition survey
            optional : this is the contents of the condition field
            from_path : this is named to match other scripts downstream
            from_filename : this is also named to match other file downstream

        Parameters
        ----------
        infilename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata found via this method
        """

        base = _os.path.basename(infilename)
        name, ext = _os.path.splitext(base)
        splitname = name.split('_')

        split_sections = ('projid', 'uniqueid', 'subprojid', 'start_date',
                          'statuscode', 'optional')

        if len(splitname) >= 5:
            meta = {
                'from_path': infilename,
                'from_filename': base,
                'projid': splitname[0],
                'uniqueid': splitname[1],
                'subprojid': splitname[2],
                'start_date': splitname[3],
                'statuscode': splitname[4],
            }
            if len(splitname) > 5:
                option = splitname[5]
                if len(splitname) > 6:
                    for n in range(6, len(splitname)):
                        option += f'_{splitname[n]}'
                meta['optional'] = option
        elif len(splitname) >= 2 and len(splitname) < 5:
            meta = {
                'from_path': infilename,
                'from_filename': base,
            }
            for index in range(len(splitname)):
                meta[split_sections[index]] = splitname[index]
        else:
            print(f'{name} appears to have a nonstandard naming convention.')
        return meta

    def _parse_xyz_header(self, filename: str, mode: str = None) -> dict:
        """
        Parse the xyz file header for meta data and return a dictionary.  The
        key words used to search are
            - NOTES.
            - PROJECT_NAME.
            - SURVEY_NAME.
            - DATES_OF_SURVEY.

        Parameters
        ----------
        filename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata found via this method
        """

        header = []
        metalist = []
        # get the header
        with open(filename, 'r') as xyz_file:
            for line in xyz_file.readlines():
                if len(line.strip()) == 0:
                    continue
                elif self._is_header(line):
                    header.append(line.strip())
                else:
                    break

        if mode == 'CENAE':
            for line in header:
                metalist.append(self._parse_assignment(line))
        else:
            # search the header for lines starting with the key words
            for line in header:
                if line.startswith('NOTES'):
                    metalist.append(self._parse_note(line))
                elif line.startswith('PROJECT SPECIFIC NOTES'):
                    metalist.append(self._parse_note(line))
                elif line.startswith('PROJECT_NAME'):
                    metalist.append(self._parse_projectname(line))
                elif line.startswith('SURVEY_NAME'):
                    metalist.append(self._parse_surveyname(line))
                elif line.startswith('DATES_OF_SURVEY'):
                    metalist.append(self._parse_surveydates(line))
            # bring all the dictionaries together
        meta = {}
        for m in metalist:
            meta = {**meta, **m}
        return meta

    def _is_header(self, line: str) -> bool:
        """
        Test if a line contains anything other than numbers it is a meta data
        line.

        ``line`` used is determined by :func:`_parse_ehydro_xyz_bathy`

        Parameters
        ----------
        line
            string of text row

        Returns
        -------
        bool
            True, if no characters match; False otherwise
        """

        return _re.search('[a-zA-Z]', line) is not None

    def _parse_note(self, line: str) -> dict:
        """
        Parse the notes line.

        ``line`` used is determined by :func:`_parse_xyz_header`

        Parameters
        ----------
        line
            A string identified to contain metadata

        Returns
        -------
        dict
            The metadata found via this method
        """

        metadata = {}
        # find the horizontal datum information.
        zone_idx = line.find('ZONE')
        zone_len = line[zone_idx:].find('.')
        horiz_datum = line[zone_idx:zone_idx + zone_len]
        if len(horiz_datum) > 0:
            fips = horiz_datum.split()[1]
            fips = fips.rstrip(',')
            metadata['from_horiz_key'] = fips
            metadata['from_wkt'] = _usefips.fips2wkt(fips)
            horiz_units = horiz_datum.split(',')[1]
            if horiz_units.strip().upper() in ('US SURVEY FEET', 'U.S. SURVEY FEET', 'FEET'):
                metadata['from_horiz_units'] = 'US Survey Foot'
            elif horiz_units.strip().upper() in ('INTL FOOT'):
                metadata['from_horiz_units'] = 'Intl Foot'
            else:
                metadata['from_horiz_units'] = horiz_units.strip()
            metadata['from_horiz_datum'] = horiz_datum
        # find the vertical datum information
        if line.find('MEAN LOWER LOW WATER') > 0:
            metadata['from_vert_datum'] = 'MEAN LOWER LOW WATER'
            metadata['from_vert_key'] = 'MLLW'
        elif line.find('MLLW') > 0:
            metadata['from_vert_datum'] = 'MEAN LOWER LOW WATER'
            metadata['from_vert_key'] = 'MLLW'
        elif line.find('MEAN LOW WATER') > 0:
            metadata['from_vert_datum'] = 'MEAN LOWER WATER'
            metadata['from_vert_key'] = 'MLW'
        vert_units_tags = ['NAVD88', 'NAVD1988', 'NAVD 1988']
        for tag in vert_units_tags:
            vert_units_end = line.find(tag)
            if vert_units_end >= 0:
                vert_units_end += len(tag)
                break
            else:
                vert_units_end = 0
        vert_units_start = vert_units_end - line[vert_units_end::-1].find('>krb<')
        vert_units = line[vert_units_start + 1:vert_units_end]
        if vert_units.find('FEET') > 0:
            metadata['from_vert_units'] = 'US Survey Foot'
        return metadata

    def _parse_projectname(self, line: str) -> dict:
        """
        Parse the project name line.

        ``line`` used is determined by :func:`_parse_xyz_header`

        Parameters
        ----------
        line
            A string identified to contain metadata

        Returns
        -------
        dict
            The metadata found via this method
        """

        name = line.split('=')[-1]
        name = name.strip('\n')
        metadata = {'projectname': name}
        return metadata

    def _parse_surveyname(self, line: str) -> dict:
        """
        Parse the survey name line.

        ``line`` used is determined by :func:`_parse_xyz_header`

        Parameters
        ----------
        line
            A string identified to contain metadata

        Returns
        -------
        dict
            The metadata found via this method
        """

        name = line.split('=')[-1]
        name = name.strip('\n')
        metadata = {'surveyname': name}
        return metadata

    def _parse_surveydates(self, line: str) -> dict:
        """
        Parse the project dates line.

        ``line`` used is determined by :func:`_parse_xyz_header`

        Parameters
        ----------
        line
            A string identified to contain metadata

        Returns
        -------
        dict
            The metadata found via this method
        """

        metadata = {}
        datestr = line.split('=')[-1]
        datestr = datestr.strip('\n')
        if datestr.find('-') > 0:
            delim = '-'
        else:
            delim = ' to '
        dateout = datestr.split(delim)
        start_date = self._xyztext2date(dateout[0])
        if start_date is not None:
            metadata['start_date'] = start_date
        if len(dateout) == 1:
            pass
        elif len(dateout) == 2:
            metadata['end_date'] = self._xyztext2date(dateout[1])
        else:
            print('ambiguous date found!')
        return metadata

    def _xyztext2date(self, textdate: str) -> str:
        """
        Take the date as provided in a text string as "day month year" as in
        "20 March 2017" and return the format "YearMonthDay" as in "20170320".

        Parameters
        ----------
        textdate
            Date string as "DD Month YYYY"

        Returns
        -------
        str
            Date String as "YYYYMMDD"
        """

        try:
            date = _datetime.strptime(textdate, '%d %B %Y')
            numdate = date.strftime('%Y%m%d')
            return numdate
        except:
            return None

    def _parse_assignment(self, line: str) -> dict:
        meta = {}
        try:
            key, value = line.split('==')
        except ValueError as e:
            _logging.debug(_logging.DEBUG, f'ValueError: {e}')
            return meta
        if key == 'Horizontal_Datum':
            meta['from_horiz_datum'] = value.strip()
            if ',' in value and '-' in value:
                fips = value.split(',')[1]
                fips = fips.strip().split()[0].split('-')[1]
            elif ',' not in value and '-' in value:
                fips = value.split(' ')
                fips = [segment.strip() for segment in fips if '-' in segment][0]
                fips = fips.strip().split()[0].split('-')[1]
            elif _re.compile(r'[0-9]{4}').search(value.strip()):
                if ',' in value:
                    value = _re.sub(',', '', value)
                fips = value.split(' ')
                fips = [segment.strip() for segment in fips if _re.compile(r'[0-9]{4}').search(segment.strip())][0]
            try:
                fips = _re.sub(r'\D', '', fips)
                meta['from_horiz_key'] = fips
            except NameError:
                _logging.debug(_logging.DEBUG, f"Unable to parse 'from_horiz_key' from: {value}")
            try:
                meta['from_wkt'] = _usefips.fips2wkt(int(fips))
            except ValueError as e:
                _logging.debug(_logging.DEBUG, f'ValueError: {e}')
                return meta
        elif key == 'Distance_Units':
            if value.strip().upper() in ('US SURVEY FEET', 'U.S. SURVEY FEET', 'FEET'):
                meta['from_horiz_units'] = 'US Survey Foot'
            elif value.strip().upper() in ('INTL FOOT'):
                meta['from_horiz_units'] = 'Intl Foot'
            else:
                meta['from_horiz_units'] = value.strip()
        elif key == 'Vertical_Datum':
            meta['from_vert_datum'] = value.strip()
            upper_key = value.strip().upper()
            if upper_key in ('MEAN LOWER LOW WATER', 'MLLW'):
                meta['from_vert_key'] = 'MLLW'

        elif key == 'Depth_Units':
            if value.strip().upper() in ('US SURVEY FEET', 'U.S. SURVEY FEET', 'FEET'):
                meta['from_vert_units'] = 'US Survey Foot'
            else:
                meta['from_vert_units'] = value.strip()
        return meta

    def _load_default_metadata(self, infilename: str, default_meta: str):
        """
        Given the file name for data and a default metadata file (containing a
        picked dictionary), look for the default file.  If that files does not
        exist, look for a file named 'default.pkl' in the same directory as the
        provided file name.

        Parameters
        ----------
        infilename
            Complete filepath of the input data
        default_meta
            Complete filepath of the input default values

        Returns
        -------
        dict
            The metadata found via this method
        """

        if len(default_meta) == 0:
            path, infile = _os.path.split(infilename)
            default_meta = _os.path.join(path, 'default.pkl')
        if _os.path.exists(default_meta):
            with open(default_meta, 'rb') as metafile:
                meta = _pickle.load(metafile)
        else:
            meta = {}
        return meta

    def _parse_ehydro_xml(self, infilename: str) -> dict:
        """
        Parse the eHydro XML file as provided by Wilmington, Charleston, and
        Norfolk Districts.

        Parameters
        ----------
        infilename : str
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata found via this method

        """
        xml_meta = self._parse_xml(infilename)
        text_meta = self._parse_xml_text(infilename)
        meta_out = {**xml_meta, **text_meta}
        return meta_out

    def _parse_xml(self, infilename: str) -> dict:
        """
        Parse the xml portion of the xml file

        Parameters
        ----------
        infilename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata found via this method
        """

        xml_meta = {}
        tree = _parse(infilename)
        root = tree.getroot()
        val = root.findtext('.//altdatum')
        if val is not None:
            xml_meta['from_vert_key'] = val
        val = root.findtext('.//altunits')
        if val is not None:
            xml_meta['from_vert_units'] = val
        val = root.findtext('.//abstract')
        if val is not None:
            vals = val.split(' ')
            for n, v in enumerate(vals):
                if v == 'dates':
                    break
            date_str = vals[n + 1]
            start, end = date_str.split(',')
            xml_meta['start_date'] = start.replace('-', '')
            xml_meta['end_date'] = end.replace('-', '')
        return xml_meta

    def _parse_xml_text(self, infilename: str) -> dict:
        """
        Pase the text portion of the xml file

        Parameters
        ----------
        infilename
            Complete filepath of the input data

        Returns
        -------
        dict
            The metadata found via this method
        """

        txt_meta = {}
        txt_keys = {
            'Implied_Vertical_Accuracy': 'from_vert_unc',
            'Implied_Horizontal_Accuracy': 'from_horiz_unc',
            'Horizontal_Zone': 'from_horiz_datum',
            'Horizontal_Datum': 'from_horiz_datum',
            'Units': 'from_horiz_units'
        }
        keys = txt_keys.keys()
        with open(infilename, 'r') as metafile:
            for line in metafile:
                for key in keys:
                    if line.startswith(key):
                        meta_key = txt_keys[key]
                        txt_meta[meta_key] = line
        for key in txt_meta:
            if key in ('from_vert_acc', 'from_horiz_acc'):
                line = txt_meta[key]
                val = line.split()[-2]
                txt_meta[key] = val
            elif key == 'from_horiz_datum':
                line = txt_meta[key]
                val = line.split(':')[-1]
                val = val.lstrip(' ')
                val = val.rstrip('\n')
                txt_meta[key] = val
                fips = _re.search('\d{4}', val)
            elif key == 'from_horiz_units':
                line = txt_meta[key]
                val = line.split(':')[-1]
                val = val.strip(' ')
                txt_meta[key] = val
        if fips is not None:
            txt_meta['fips'] = int(fips.group())
        return txt_meta

    def _finalize_meta(self, meta):
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
        # if the data is coming from this reader these should be true.
        meta['read_type'] = 'ehydro'
        meta['posted'] = False
        meta['interpolated'] = 'False'
        meta['from_horiz_frame'] = 'NAD83'
        meta['from_horiz_type'] = 'spc'
        # make sure the dates are useful
        if 'end_date' not in meta and 'start_date' in meta:
            meta['end_date'] = meta['start_date']
        elif ('end_date' in meta and meta['end_date'] == '') and 'start_date' in meta:
            meta['end_date'] = meta['start_date']
        # convert from text to standard values
        if 'from_horiz_units' in meta:
            if meta['from_horiz_units'].upper() in ('US SURVEY FOOT'):
                meta['from_horiz_units'] = 'us_ft'
            elif meta['from_horiz_units'].upper() in ('INTL FOOT'):
                meta['from_horiz_units'] = 'ft'
            else:
                raise ValueError(f'Input datum units are unknown: {meta["from_horiz_units"]}')
        if 'from_vert_key' in meta:
            meta['from_vert_key'] = meta['from_vert_key'].lower()
        if 'from_vert_units' in meta:
            if meta['from_vert_units'].upper() == 'US SURVEY FOOT':
                meta['from_vert_units'] = 'us_ft'
            else:
                raise ValueError(f'Input datum units are unknown: {meta["from_vert_units"]}')
        # add the default quality metrics
        final_meta = {**_ehydro_quality_metrics, **meta}
        return final_meta

    def _parse_ehydro_xyz_bathy(self, filename: str) -> _np.array:
        """
        Read the best available point bathymetry for the district.

        This method assumes the provided file name is the XYZ file.

        Parameters
        ----------
        filename
            Complete filepath of the input data

        Returns
        -------
        numpy.array
            The xyz values as a 2d array
        """

        points = []

        # get the header
        with open(filename, 'r') as input_file:
            message = f'{filename}: parser '
            comma = False
            for line in input_file.readlines():
                if line not in ('\n', '\x1a', '') and not self._is_header(line):
                    if ',' in line:
                        comma = True
                        row = [float(entry.strip()) for entry in line.split(',')]
                        if len(row) == 3:
                            points.append(row)
                    else:
                        row = [float(entry) for entry in line.split()]
                        if len(row) == 3:
                            points.append(row)
            if comma:
                message += f'found comma-delimited file "{filename}"'
            else:
                message += f'found whitespace-delimited file (tab or space) "{filename}"'
            self._logger.log(_logging.DEBUG, message)

        points = _np.asarray(points)
        return points
