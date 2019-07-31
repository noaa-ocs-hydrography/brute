# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:49:10 2019

@author: Casiano.Koprowski
"""

import os as _os
import pickle as _pickle
import re as _re
import sys as _sys
import logging as _logging
from datetime import datetime as _datetime
from xml.etree.ElementTree import parse as _parse

import numpy as _np
from . import parse_usace_xml
from . import parse_usace_pickle
from fuse.datum_transform import usefips as _usefips


class Base:
    def __init__(self, version=None):
        self.version = version
        self.ussft2m = 0.30480060960121924  # US survey feet to meters
        self.xyz_suffixes = ('_A', '_FULL')
        self.xyz_files = {}

        self._logger = _logging.getLogger(f'fuse')

        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            ch.setLevel(_logging.DEBUG)
            self._logger.addHandler(ch)

    def read_metadata(self, infilename: str):
        """
        Read all available meta data.
        returns dictionary

        Parameters
        ----------
        infilename: str

        Returns
        -------

        """
        meta_pickle = self._parse_pickle(infilename)
        meta_xml = self._parse_usace_xml(infilename)
        meta_xml['poly_name'] = meta_pickle['poly_name']
        return meta_xml

    def read_bathymetry(self, infilename: str):
        """
        Read the bathymetry and return an array of the xyz points.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        xyz = self._parse_ehydro_xyz_bathy(infilename)
        return xyz

    def read_bathymetry_by_point(self, infilename: str):
        """
        Read the bathymetry and return point by point.

        Parameters
        ----------
        infilename :

        infilename: str :


        Returns
        -------

        """

        xyz = self.read_bathymetry(infilename)
        for n in xyz:
            yield n

    def _parse_pickle(self, filename: str) -> dict:
        """
        Retrieves metadata for USACE E-Hydro files
        function returns metadata dictionary

        Parameters
        ----------
        filename: str


        Returns
        -------
        metadata : dict


        """
        pickle_name = self.name_gen(filename, ext='pickle', sfx=False)
        pickle_dict = parse_usace_pickle.read_pickle(pickle_name, pickle_ext=True)
        pickle_keys = parse_usace_pickle.dict_keys(pickle_dict)
        return pickle_dict

    def name_gen(self, filename: str, ext: str = None, sfx: bool = True) -> str:
        """
        Returns the suffix of the a survey's xyz file and the file name
        for a different extension

        """
        filebase, fileext = _os.path.splitext(filename)
        suffix = None
        for item in self.xyz_suffixes:
            if _re.compile(f'{item}$', _re.IGNORECASE).search(filebase):
                suffix = item
        if suffix is not None and ext is not None:
            base = _re.sub(_re.compile(f'{suffix}$', _re.IGNORECASE), '', filebase) + f'.{ext}'
        elif suffix is None and ext is not None:
            base = filebase + f'.{ext}'
        else:
            base = filename

        if sfx:
            return base, suffix
        else:
            return base

    def _parse_ehydro_xyz_header(self, infilename: str, meta_source: str = 'xyz', default_meta: str = '') -> dict:
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
        infilename :
            param meta_source:  (Default value = 'xyz')
        default_meta :
            Default value = '')
        infilename: str :

        meta_source: str :
             (Default value = 'xyz')
        default_meta: str :
             (Default value = '')

        Returns
        -------

        """

        name_meta = self._parse_filename(infilename)
        if meta_source == 'xyz':
            file_meta = self._parse_xyz_header(infilename)
        elif meta_source == 'xml':
            file_meta = self._parse_ehydro_xml(infilename)
        default_meta = self._load_default_metadata(infilename, default_meta)
        merged_meta = {**default_meta, **name_meta, **file_meta}
        if 'from_horiz_unc' in merged_meta:
            if merged_meta['from_horiz_units'] == 'US Survey Foot':
                val = self.ussft2m * float(merged_meta['from_horiz_unc'])
                merged_meta['horiz_uncert'] = val
        if 'from_vert_unc' in merged_meta:
            if merged_meta['from_vert_units'] == 'US Survey Foot':
                val = self.ussft2m * float(merged_meta['from_vert_unc'])
                merged_meta['vert_uncert_fixed'] = val
                merged_meta['vert_uncert_vari'] = 0
        sorind = f"{name_meta['projid']}_{name_meta['uniqueid']}_{name_meta['subprojid']}_{name_meta['start_date']}_" + \
                 f"{name_meta['statuscode']}"
        merged_meta['source_indicator'] = f'US,US,graph,{sorind}'
        # merged_meta['script_version'] = __version__
        return merged_meta

    def _parse_filename(self, infilename):
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
        infilename :


        Returns
        -------

        """
        base = _os.path.basename(infilename)
        name, ext = _os.path.splitext(base)
        splitname = name.split('_')

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
        else:
            print(f'{name} appears to have a nonstandard naming convention.')
        return meta

    def _parse_xyz_header(self, infilename):
        """
        Parse the xyz file header for meta data and return a dictionary.  The
        key words used to search are
            NOTES
            PROJECT_NAME
            SURVEY_NAME
            DATES_OF_SURVEY

        Parameters
        ----------
        infilename :


        Returns
        -------

        """
        header = []
        metalist = []
        # get the header
        with open(infilename, 'r') as infile:
            for line in infile.readlines():
                if line == '\n':
                    continue
                elif self._is_header(line):
                    header.append(line)
                else:
                    break
        # search the header for lines starting with the key words
        for line in header:
            if line.startswith('NOTES'):
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

    def _is_header(self, line):
        """
        Test if a line contains anything other than numbers it is a meta data line.

        Parameters
        ----------
        line :


        Returns
        -------

        """
        pattern = '[a-zA-Z]'
        if _re.search(pattern, line) is None:
            return False
        else:
            return True

    def _parse_note(self, line):
        """
        Parse the notes line.

        Parameters
        ----------
        line :


        Returns
        -------

        """
        metadata = {}
        # find the horizontal datum information.
        zone_idx = line.find('ZONE')
        zone_len = line[zone_idx:].find('.')
        horiz_datum = line[zone_idx:zone_idx + zone_len]
        if len(horiz_datum) > 0:
            fips = horiz_datum.split()[1]
            fips = fips.rstrip(',')
            metadata['from_fips'] = fips
            metadata['from_wkt'] = _usefips.fips2wkt(fips)
            horiz_units = horiz_datum.split(',')[1]
            if horiz_units.lstrip(' ') == 'US SURVEY FEET':
                metadata['from_horiz_units'] = 'US Survey Foot'
            else:
                metadata['from_horiz_units'] = horiz_units.lstrip(' ')
            metadata['from_horiz_datum'] = horiz_datum
        else:
            metadata['from_wkt'] = 'unknown'
            metadata['from_horiz_units'] = 'unknown'
            metadata['from_horiz_datum'] = 'unknown'
        # find the vertical datum information
        if line.find('MEAN LOWER LOW WATER') > 0:
            metadata['from_vert_key'] = 'MLLW'
        elif line.find('MLLW') > 0:
            metadata['from_vert_key'] = 'MLLW'
        elif line.find('MEAN LOW WATER') > 0:
            metadata['from_vert_key'] = 'MLW'
        else:
            metadata['vert_key'] = 'Unknown'
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
        metadata['from_vert_datum'] = vert_units
        if vert_units.find('FEET') > 0:
            metadata['from_vert_units'] = 'US Survey Foot'
        else:
            metadata['from_vert_units'] = 'unknown'
        return metadata

    def _parse_projectname(self, line):
        """
        Parse the project name line.

        Parameters
        ----------
        line :


        Returns
        -------

        """
        name = line.split('=')[-1]
        name = name.strip('\n')
        metadata = {'projectname': name}
        return metadata

    def _parse_surveyname(self, line):
        """
        Parse the survey name line.

        Parameters
        ----------
        line :


        Returns
        -------

        """
        name = line.split('=')[-1]
        name = name.strip('\n')
        metadata = {'surveyname': name}
        return metadata

    def _parse_surveydates(self, line):
        """
        Parse the project dates line.

        Parameters
        ----------
        line :


        Returns
        -------

        """
        metadata = {}
        datestr = line.split('=')[-1]
        datestr = datestr.strip('\n')
        if datestr.find('-') > 0:
            delim = '-'
        else:
            delim = ' to '
        dateout = datestr.split(delim)
        metadata['start_date'] = self._xyztext2date(dateout[0])
        if len(dateout) == 1:
            metadata['end_date'] = 'unknown'
        elif len(dateout) == 2:
            metadata['end_date'] = self._xyztext2date(dateout[1])
        else:
            print('ambiguous date found!')
        return metadata

    def _xyztext2date(self, textdate):
        """
        Take the date as provided in a text string as "day month year" as in
        "20 March 2017" and return the format "YearMonthDay" as in "20170320".

        Parameters
        ----------
        textdate :


        Returns
        -------

        """
        try:
            date = _datetime.strptime(textdate, '%d %B %Y')
            numdate = date.strftime('%Y%m%d')
            return numdate
        except:
            return 'unknown'

    def _load_default_metadata(self, infilename, default_meta):
        """
        Given the file name for data and a default metadata file (containing a
        picked dictionary), look for the default file.  If that files does not
        exist, look for a file named 'default.pkl' in the same directory as the
        provided file name.

        Parameters
        ----------
        infilename :
            param default_meta:
        default_meta :


        Returns
        -------

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

    def _parse_usace_xml(self, infilename):
        """
        Read all available meta data.
        returns dictionary

        Parameters
        ----------
        infilename: str

        """
        xmlfilename = self.name_gen(infilename, ext='xml', sfx=False)
        if _os.path.isfile(xmlfilename):
            with open(xmlfilename, 'r') as xml_file:
                xml_txt = xml_file.read()
            xmlbasename = _os.path.basename(xmlfilename)
            xml_data = parse_usace_xml.XML_Meta(xml_txt, filename=xmlbasename)
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
        meta_xml['from_path'] = infilename
        meta_xml['from_filename'] = _os.path.basename(infilename)
        return {**meta_xml, **ext_dict}

    def _parse_ehydro_xml(self, infilename):
        """
        Parse the eHydro XML file as provided by Wilmington, Charleston, and
        Norfolk Districts.

        Parameters
        ----------
        infilename :


        Returns
        -------

        """
        xml_meta = self._parse_xml(infilename)
        text_meta = self._parse_xml_text(infilename)
        meta_out = {**xml_meta, **text_meta}
        return meta_out

    def _parse_xml(self, infilename):
        """
        Parse the xml portion of the xml file

        Parameters
        ----------
        infilename :


        Returns
        -------

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

    def _parse_xml_text(self, infilename):
        """
        Pase the text portion of the xml file

        Parameters
        ----------
        infilename :


        Returns
        -------

        """
        txt_meta = {}
        txt_keys = {'Implied_Vertical_Accuracy': 'from_vert_unc',
                    'Implied_Horizontal_Accuracy': 'from_horiz_unc',
                    'Horizontal_Zone': 'from_horiz_datum',
                    'Units': 'from_horiz_units'}
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
                val = val.lstrip(' ')
                val = val.rstrip('\n')
                txt_meta[key] = val
        if fips is not None:
            txt_meta['fips'] = int(fips.group())
        return txt_meta

    def _parse_ehydro_xyz_bathy(self, infilename):
        """
        Read the best available point bathymetry for the district.

        This method assumes the provided file name is the XYZ file.

        Parameters
        ----------
        infilename :


        Returns
        -------

        """
        bathy = []
        # get the header
        with open(infilename, 'r') as infile:
            for line in infile.readlines():
                if line == '\n':
                    continue
                elif line == '\x1a':
                    continue
                elif self._is_header(line):
                    pass
                else:
                    if ',' in line:
                        output = f'Comma delimited file found: {infilename}'
                        self._logger.log(_logging.DEBUG, output)
                        raise ValueError(output)
                    else:
                        bathy.append([float(x) for x in line.split(' ')])
        bathy = _np.asarray(bathy)
        return bathy
