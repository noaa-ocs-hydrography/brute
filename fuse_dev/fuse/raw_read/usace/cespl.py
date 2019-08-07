# -*- coding: utf-8 -*-
"""
Edited by Juliet Kinney
extract_ehydro_meta_class_CESAJ.py

This script takes as an input the filename and path of an USACE e-Hydro .xyz
text file and pull the metadata from the xyz header and file name.
Additionally, it passes the matching .xml file and utilizes
the parse_usace_xml script to pull out the relevant metadata if it
is either FGDC or ISO FGDC USACE metadata format.

update 4/5/19
major update April 2, 2019
update July 12,2019 adding in call to pickle reader need to commit 8/7/2019
"""
__version__ = 'FUSE'
import os as os
import pickle as _pickle
import re as _re

_ussft2m = 0.30480060960121924  # US survey feet to meters
import dateutil.parser as parser
from datetime import datetime
import numpy as _np

try:
    import fuse.raw_read.usace.parse_usace_xml as p_usace_xml
except:
    try:
        from . import parse_usace_xml as p_usace_xml
    except:
        print('importing fuse.raw_read.usace.parse_usace_xml as p_usace_xml did not work')
try:
    import fuse.raw_read.usace.parse_usace_pickle as parse_usace_pickle
except:
    try:
        from . import parse_usace_pickle as parse_usace_pickle
    except:
        print('importing fuse.raw_read.usace.parse_usace_pickle as parse_usace_pickle  did not work')


##-----------------------------------------------------------------------------


class CESPLRawReader:
    """
    This class passes back bathymetry
    & a metadata dictionary from the e-Hydro files
    
    Parameters
    ----------
    
    Returns
    -------
    
    """


    def read_metadata(self, infilename: str) -> dict:
        """
        Read all available meta data.
        returns dictionary
        
        Parameters
        ----------
        infilename: str :
        
        
        Returns
        -------
        dict
        
        """
        version = 'CESPL'
        self.version = version
        return retrieve_meta_for_Ehydro_out_onefile(infilename)

    def read_bathymetry_dat(self, infilename: str) -> _np.array:
        """
        Read the bathymetry from the .dat file. The dat file is less precise,
        but had no header and is in a standardized format
        
        Parameters
        ----------
        infilename :
        
        
        Returns
        -------
        xyz
        
        """
        # get the dat file for CESAJ# Jacksonville
        stub, ext = os.path.splitext(infilename)
        bathyfilename = f'{stub}.dat'
        """
        F-strings provide a way to embed expressions inside string literals, using a minimal syntax.
        It should be noted that an f-string is really an expression evaluated at run time, not a constant
         value. In Python source code, an f-string is a literal string, prefixed with f, which contains 
        expressions inside braces. The expressions are replaced with their values.
        https://realpython.com/python-f-strings/
        one can also include expressions within the quoted strings, The expressions in an f-string are evaluated in left-to-right order. This is detectable only if the expressions have side effects:
        https://www.python.org/dev/peps/pep-0498/
        """

        xyz = _np.loadtxt(bathyfilename, delimiter=' ')
        return xyz

    def read_bathymetry(self, infilename: str) -> _np.array:
        """
        Read the bathymetry from the xyz files, this tells it to not include
        the header when reading the file
        
        Note: The high resolution multibeam files are available as .xyz on E-Hydro
        
        Parameters
        ----------
        infilename: str :
        
        Returns
        -------
        xyz
        
        """
        version = 'CESAJ'
        self.version = version
        first_instance, commas_present = _start_xyz(infilename)
        print(infilename)# remove later
        if first_instance != '':
            xyz = _np.loadtxt(infilename, delimiter=' ', skiprows=first_instance, usecols=(0, 1, 2))
        else:
            xyz = _np.loadtxt(infilename, delimiter=' ',
                              usecols=(0, 1, 2))  # ignoring anything after the first 3 columns on import
            # other option from np.loadtxt(infilename, converters={4:datestr2num})
        return xyz


# ------------------------------------------------------------------------------

def return_surveyid(filenamepath: str, ex_string: str) -> str:
    """
    strip end of filename off
    surveybasename =return_surveyid(filenamepath, ex_string)
    
    Parameters
    ----------
    filenamepath: str :
        param ex_string:
    ex_string: str :
    
    Returns
    -------
    surveybasename: str:
    
    """
    basename = os.path.basename(filenamepath)
    surveybasename = basename.rstrip(ex_string)
    return surveybasename


# ------------------------------------------------------------------------------

def retrieve_meta_for_Ehydro_out_onefile(filename: str) -> dict:
    """
    retrieve metadata for USACE E-Hydro files
    function returns metadata dictionary
    
    input is filename of .xyz file with path
    
    Parameters
    ----------
    filename: str :
    
    Returns
    -------
    dict:
    
    """
    # next if pull the subset of the table in the dataframe related to the list of files passed to it.
    merged_meta = {}
    merge2 = {}
    f = filename
    basename = os.path.basename(f)
    e_t = XYZMetaReader(f)
    # xml pull here.
    xmlfilename = get_xml_match(f)
    if os.path.isfile(xmlfilename):
        with open(xmlfilename, 'r') as xml_file:
            xml_txt = xml_file.read()
        xmlbasename = os.path.basename(xmlfilename)
        xml_data = p_usace_xml.XMLMetadata(xml_txt, filename=xmlbasename)
        if xml_data.version == 'USACE_FGDC':
            meta_xml = xml_data._extract_meta_USACE_FGDC()
        elif xml_data.version == 'ISO-8859-1':
            meta_xml = xml_data._extract_meta_USACE_ISO()
            if 'ISO_xml' not in meta_xml:
                meta_xml = xml_data._extract_meta_USACE_FGDC(
                    override='Y')  # xml_data._extract_meta_ISOlabel_USACE_FGDC()
        else:
            meta_xml = xml_data.convert_xml_to_dict2()
        ext_dict = xml_data.extended_xml_fgdc()
        ext_dict = p_usace_xml.ext_xml_map_enddate(ext_dict)
        meta_xml = p_usace_xml.xml_SPCSconflict_flag(meta_xml)
    else:
        ext_dict = {}
        meta_xml = {}
    meta = e_t.parse_ehydro_xyz(f, meta_source='xyz', version='CESAJ', default_meta='')
    meta['special_handling'] = _check_special_handling(basename)#special handling is saved with text meta as it has to do with the text file
    # bringing ehydro table attributs(from ehydro REST API)saved in pickle during ehydro_move #empty dictionary place holder for future ehydro table ingest (make come from imbetween source TBD)
    meta_from_ehydro = {}
    
    e_pick = eHydroPickleReader(xmlfilename)
    meta_from_ehydro = e_pick._read_pickle()  # to handle files
    meta_from_ehydro = e_pick._when_use_pickle(meta_xml)
    meta_from_ehydro = e_pick._when_use_pickle_startdate(meta_xml)
    
    list_keys_empty = []
    combined_row = {}
    subset_row = {}
    subset_no_overlap = {}
    subset_dict = {}
    for key in meta:
        if meta[key] in ('unknown', ''):
            list_keys_empty.append(key)
        else:
            subset_row[key] = meta[key]
            # non blank columns only
            if key in meta_xml:
                if meta[key] == meta_xml[key]:
                    combined_row[key] = meta[key]
                    # only make list within cell if values from different sources are different
                else:
                    combined_row[key] = f'{meta[key]} , {meta_xml[key]}'
            else:
                subset_no_overlap[key] = meta[key]
    for key in ext_dict:
        if ext_dict[key] in ('unknown', '') or ext_dict[key] is None:
            list_keys_empty.append(key)
        else:
            if key in meta_xml:
                if meta_xml[key] == ext_dict[key]:
                    combined_row[key] = ext_dict[key]
                    # only make list within cell if values from different sources are different
                else:
                    combined_row[key] = f'{ext_dict[key]} , {meta_xml[key]}'
            else:
                subset_dict[key] = ext_dict[key]
    merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row}  # this one excluded 'unknown' keys, and
    # in merging sources from the text file and xml it will show any values that do not match as a list.
    merged_meta = {**meta, **meta_from_ehydro, **meta_xml}  # this method overwrites
    try:
        merged_meta = check_date_order(merged_meta, merged_meta)
    except:
        err_file = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\CESAJ\metadata\Error_file_if_date_checkfail.txt"
        with open(err_file, 'a') as error:
            error.write(f'{f} : extra dict END DATE SEARCH call fail \n')

    return merged_meta


###----------------------------------------------------------------------------
class EhydroPickleReader(object):
    """
    Reading in ehydro pickle file, looking for conflicts with other metadata inputs
    and determining when /how to pass on metadata attributes 
    """
    
    def __init__(self, infilename: str):
        """
        Pass filename that matches the pickle file you want to match 
        but with any extension
        (Here we tend to pass the xmlfilename as it already has been matched
        in the cases of _A.xyz etc., but one could use a .xyz file)
          
          
        Parameters
        ----------
        infilename: str:
        
        
        Returns
        -------
        self.filename = infilename: str:
        
        """
        self.filename = infilename
    
    def _read_pickle(self) -> dict:
        """
        Read in picklefile that ehydro_move creates from the E-Hydro REST API
        table attributes.
        
        
        Parameters
        ----------
        self.infilename: str:
        
        
        Returns
        -------
        pickle_meta: dict:
        
        """
        print(f'reading in pickle based on: {self.filename}')  # making sure pickle passing is working
        pickle_meta = parse_usace_pickle.read_pickle(self.filename)
        self.meta_from_ehydro = pickle_meta
        return pickle_meta

    def _Check_for_SPCSconflicts(self, meta_xml):
        """
        Checking to see if the SPCS codes conflict between sources
        
        Parameters
        ----------
        meta_xml: dict:
        self.meta_from_ehydro: dict
        
        Returns
        -------
        no_SPCS_conflict: str:
        no_SPCS_conflict_withpickle: str:
        meta_from_ehydro: dict:
        
        """
        
        meta_from_ehydro = self.meta_from_ehydro
        
        no_SPCS_conflict = ''
        no_SPCS_conflict_withpickle = ''
        if 'SPCS_conflict_XML' in meta_from_ehydro:
            if meta_from_ehydro['SPCS_conflict_XML'] != '':
                no_SPCS_conflict = 'False'
            else:
                no_SPCS_conflict = 'True'
        
        if 'SOURCEPROJECTION' in meta_from_ehydro:
            if 'from_fips' in meta_xml:
                meta_xml = p_usace_xml.xml_SPCSconflict_otherspcs(meta_xml,
                                                                  f"{p_usace_xml.SOURCEPROJECTION_dict, meta_from_ehydro['SOURCEPROJECTION']}")
                if p_usace_xml.convert_tofips(p_usace_xml.SOURCEPROJECTION_dict,
                                              meta_from_ehydro['SOURCEPROJECTION']) == \
                        meta_xml['from_fips']:
                    no_SPCS_conflict_withpickle = 'True'
                else:
                    no_SPCS_conflict_withpickle = 'False'
            if meta_xml['SPCS_conflict_XML_other'] != '':
                no_SPCS_conflict = 'False'
                # We know for CEMVN thath this will conflict with some of the SPCS values but have a method that works.
                # this way we pass on that there are conflicts but do not raise a flag unless the final from_fips disagrees
        
        meta_from_ehydro['no_SPCS_conflict_withpickle'] = no_SPCS_conflict_withpickle
        self.meta_from_ehydro
        return no_SPCS_conflict, no_SPCS_conflict_withpickle, meta_from_ehydro

    def _when_use_pickle(self, meta_xml):
        """
        If there is no SPCS code in the xml, use the pickle/ REST API SPCS code
        
        Additional check to see if their is a conflict. District specific rules on conflict resolution may need to apply.
        1st assumption is that the REST API has the correct SPCS code according to E-Hydro team. (John McKenzie) and reinterated by
        District contacts thus far (as of June 2019) base on E-hydro upload procedures.
        
        Parameters
        ----------
        meta_xml: dict:
        self.meta_from_ehydro : dict:
        
        
        Returns
        -------
        meta_from_ehydro: dict:
        
        
        """
        meta_from_ehydro = self.meta_from_ehydro
        if 'SOURCEPROJECTION' in meta_from_ehydro:
            if 'from_FIPS' in meta_xml:
                # run check for conflict
                no_SPCS_conflict, no_SPCS_conflict_withpickle = self._Check_for_SPCSconflicts(meta_xml,
                                                                                              meta_from_ehydro)
                if no_SPCS_conflict_withpickle == 'False':
                    meta_from_ehydro['from_fips'] = p_usace_xml.convert_tofips(p_usace_xml.SOURCEPROJECTION_dict,
                                                                               meta_from_ehydro['SOURCEPROJECTION'])
            else:
                meta_from_ehydro['from_fips'] = p_usace_xml.convert_tofips(p_usace_xml.SOURCEPROJECTION_dict,
                                                                           meta_from_ehydro['SOURCEPROJECTION'])
        self.meta_from_ehydro = meta_from_ehydro
        return meta_from_ehydro

    def _when_use_pickle_startdate(self, meta_xml):
        """
        if xml_meta is blank and if meta does not have information use pickle data for date
        next: Check survey start & end date against filename and other locations
        
        Parameters
        ----------
        meta_xml :
        self: # uses meta_from_ehydro :
        
        
        Returns
        """
        meta_from_ehydro = self.meta_from_ehydro
        if meta_from_ehydro:  # check if dictionary empty
            if meta_xml:  # check if dictionary empty
                print(meta_from_ehydro['SURVEYDATEEND'])
                # Check survey start & end date against filename and other locations
            else:  # if xml_meta is blank and if meta does not have information use pickle data:
                meta_from_ehydro['start_date'] = meta_from_ehydro['SURVEYDATESTART']
                # "SURVEYDATESTART"
                # "SURVEYDATEEND"
        return meta_from_ehydro


###----------------------------------------------------------------------------
class XYZMetaReader(object):
    """
    Extract both information from the filename as well as from the text file's header
    """
    def __init__(self, preloadeddata, version='', filename=''):
        """
        xyz file (the ascii text file) handler for metadata parsing  gets initiated here
        Parameters
        ----------
        preloadeddata
        version
        filename
        """
        self.filename = preloadeddata
        if filename != "" or None:
            self.filename_1 = filename
            self.errorfile = os.path.join(os.path.dirname(filename),
                                          'TEST_extract_ehdyro_meta_class_CESAJ_ErrorFile1.txt')
        else:
            self.errorfile = os.path.join(os.path.dirname(filename),
                                          'Default_extract_ehdyro_meta_class_CESAJ_error.txt')
    def parse_ehydro_xyz(self, infilename: str, meta_source: str = 'xyz', version: str = 'CESAJ',
                         default_meta: str = '') -> dict:  # need to change version to None
        """
        'CESAJ' Jacksonville USACE
        
        Parameters
        ----------
        infilename :
            param meta_source:  (Default value = 'xyz')
        version :
            Default value = 'CESAM')
        default_meta :
            Default value = '')
        meta_source :
             (Default value = 'xyz')
        
        Returns
        -------
        
        """
        """
        Parse an USACE eHydro file for the available meta data.
        
        Default metadata (values predetermined for the file but not in the file)
        can be stored at the location defined by 'default_meta' as a pickled
        dicitonary.  If no path is provided the dictionary in the same folder as
        the data but in the file 'default.pkl' will be used.  If this file does
        not exist no default metadata will be loaded.  If the same keyword for
        the metadata exists both in the file metadata and in the default location,
        the file metadata will take precidence.
        """
        name_meta = self.parse_ehydro_filename(infilename)
        if 'start_date' in name_meta:
            name_meta['filename_date'] = name_meta['start_date']
        if meta_source == 'xyz':
            file_meta = self.parse_xyz_header(infilename, version)
        # elif meta_source == 'xml':
        #    file_meta = parse_ehydro_xml(infilename)
        default_meta = self.load_default_metadata(infilename, default_meta)
        merged_meta = {**default_meta, **name_meta, **file_meta}
        if 'from_horiz_unc' in merged_meta:
            if merged_meta['from_horiz_units'] == 'US Survey Foot':
                val = _ussft2m * float(merged_meta['from_horiz_unc'])
                merged_meta['horiz_uncert'] = val
        if 'from_vert_unc' in merged_meta:
            if merged_meta['from_vert_units'] == 'US Survey Foot':
                val = _ussft2m * float(merged_meta['from_vert_unc'])
                merged_meta['vert_uncert_fixed'] = val
                merged_meta['vert_uncert_vari'] = 0
        sorind = f"{name_meta['projid']}_{name_meta['uniqueid']}_{name_meta['subprojid']}_{name_meta['start_date']}_" + \
                 f"{name_meta['statuscode']}"
        merged_meta['source_indicator'] = f'US,US,graph,{sorind}'
        merged_meta['script_version'] = __version__
        return merged_meta

    def parse_ehydro_filename(self, infilename):
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
        base = os.path.basename(infilename)
        name, ext = os.path.splitext(base)
        splitname = name.split('_')
        if len(splitname) >= 4:
            meta = {
                'from_path': infilename,
                'from_filename': base,
                'projid': splitname[0],
                'uniqueid': splitname[1],
                'subprojid': splitname[2],
                'start_date': splitname[3],
            }
            if len(splitname) > 4:
                meta['statuscode'] = splitname[4]
                if len(splitname) > 5:
                    option = splitname[5]
                    if len(splitname) > 6:
                        for n in range(6, len(splitname)):
                            option += f'_{splitname[n]}'
                    meta['optional'] = option
            else:
                meta['statuscode'] = ''
        else:
            print(f'{name} appears to have a nonstandard naming convention.')
        return meta

    def parse_xyz_header(self, infilename, version=None):
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
            param version:  (Default value = None)
        version :
             (Default value = None)
        
        Returns
        -------
        
        """
        header = []
        metalist = []
        more_metalist = []
        # get the header
        if version == 'CESAJ':
            with open(infilename, 'r') as infile:
                for line in infile.readlines():
                    if line == '\n':
                        continue
                    elif _is_header2(line):
                        header.append(line)
                    else:
                        break
                for line in infile.readlines():
                    if line.startswith('Survey_Number=='):
                        metalist.append(_parse_Survey_Number(line))
                    if line.startswith('Survey_Type=='):
                        metalist.append(_parse_Survey_Type(line))
        try:
            for m in metalist:
                meta = {**meta, **m}
            return meta
        except:
            meta = {}
            errorfile = self.errorfile
            with open(errorfile, 'a') as metafail:
                metafail.write(f'{infilename}\n')
            return meta

    def load_default_metadata(self, infilename, default_meta):
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
            path, infile = os.path.split(infilename)
            default_meta = os.path.join(path, 'default.pkl')
        if os.path.exists(default_meta):
            with open(default_meta, 'rb') as metafile:
                meta = _pickle.load(metafile)
        else:
            meta = {}
        return meta


##-----------------------------------------------------------------------------

def get_xml(filename):
    """
    input USACE .xyz/.XYZ filename or any last extension and return .xml
    xmlname = get_xml(filename) this makes this friendlier to .ppxyz files
    for instance
    
    Parameters
    ----------
    filename :
        
    
    Returns
    -------
    
    """
    basef = filename.rpartition('.')[0]
    xml_name = f'{basef}.xml'
    return xml_name


def get_xml_xt(filename, extension):
    """
    input USACE text filename and ending to chop to get to basename
    output will be the .xml file name
    (_A.xyz for instance or _FULL.XYZ are examples of extensions)
    xmlname = get_xml_xt(filename, extension)
    
    Parameters
    ----------
    filename :
        param extension:
    extension :
    
    
    Returns
    -------
    
    """
    end_len = len(extension)
    if filename[-end_len:].upper() == extension:
        basef = filename[:-end_len]
    else:
        basef = filename
    xml_name = f'{basef}.xml'
    return xml_name


def get_xml_match(f):
    """
    input USACE .xyz/.XYZ filename or any last extension and return .xml
    it will try to match the non-full survey to the full density survey
    inorder to use the matching xml
    
    Parameters
    ----------
    f :
        
    
    Returns
    -------
    
    """
    ext_list = ['_FULL.XYZ', '_A.XYZ', '.PPXYZ']
    for extension in ext_list:
        if extension in f.upper():
            xmlfilename = get_xml_xt(f, extension)
        else:
            xmlfilename = get_xml(f)
    return xmlfilename


##-----------------------------------------------------------------------------
def _check_special_handling(basename):
    """
    Doing a check if the xyz file type is full resolution or may have
    other special handling flags that should be passed
    
    Parameters
    ----------
    basename :
        
    
    Returns
    -------
    """
    special_handling = ''
    if basename.find('.ppxyz') > 0:
        special_handling = 'ppxyz'
    full_res = ['_A.xyz', '_A.XYZ', '_FULL.xyz', '_FULL.XYZ']
    for ext_full in full_res:
        if basename.find(ext_full) > 0:
            special_handling = 'FullRES'
    return special_handling


##-----------------------------------------------------------------------------

def _start_xyz(infilename):
    """
    looks for the first line of the xyz data after the header
    returns the row number of first line of data
    
    Parameters
    ----------
    infilename :
    
    Returns
    -------
    
    """
    first_instance = ''
    numberofrows = []
    commas_present = ''
    pattern_coordinates = '[\d][\d][\d][\d][\d][\d]'  # at least six digits# should be seven then . plus two digits
    with open(infilename, 'r') as infile:
        for (index1, line) in enumerate(infile):
            if _re.match(pattern_coordinates, line) is not None:
                numberofrows.append(index1)
                if line.find(',') > 0:
                    commas_present = ','
        first_instance = numberofrows[0]
    return first_instance, commas_present


def _is_header2(line, version=None):
    """
    looks at header
    
    Parameters
    ----------
    line :
        param version:  (Default value = None)
    version :
         (Default value = None)
    
    Returns
    -------
    
    """
    if version is None:
        version = ''
    if version == 'CESAJ':
        pattern_coordinates = '[\d][\d][\d][\d][\d][\d]'  # at least six digits# should be seven then . plus two digits
        if _re.match(pattern_coordinates, line) is not None:
            return False
        else:
            return True
    elif version == '':
        pattern = '[a-zA-Z]'
        if _re.search(pattern, line) is None:
            return False
        else:
            return True
    else:
        pattern = '[^0-9]'  # anything except 0-9
        if _re.match(pattern, line) is None:  # does the string start with this pattern?
            return False
        else:
            return True


def _parse_projectname(line):
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


def _parse_note(line):
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
        horiz_units = horiz_datum.split(',')[1]
        if horiz_units.lstrip(' ') == 'US SURVEY FEET':
            metadata['from_horiz_units'] = 'US Survey Foot'
        else:
            metadata['from_horiz_units'] = horiz_units.lstrip(' ')
        metadata['from_horiz_datum'] = horiz_datum
    else:
        metadata['from_horiz_units'] = 'unknown'
        metadata['from_horiz_datum'] = 'unknown'
    # find the vertical datum information
    if line.find('MEAN LOWER LOW WATER') >= 0:
        metadata['from_vert_key'] = 'MLLW'
    elif line.find('MLLW') >= 0:
        metadata['from_vert_key'] = 'MLLW'
    elif line.find('MEAN LOW WATER') >= 0:
        metadata['from_vert_key'] = 'MLW'
    else:
        metadata['vert_key'] = 'unknown'
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
    if vert_units.find('FEET') >= 0:
        metadata['from_vert_units'] = 'US Survey Foot'
    else:
        metadata['from_vert_units'] = 'unknown'
    return metadata


def _parse_surveyname(line):
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


def _parse_surveydates(line):
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
    if datestr.find('-') >= 0:
        delim = '-'
    else:
        delim = ' to '
    dateout = datestr.split(delim)
    metadata['start_date'] = _xyztext2date(dateout[0])
    if len(dateout) == 1:
        metadata['end_date'] = 'unknown'
    elif len(dateout) == 2:
        metadata['end_date'] = _xyztext2date(dateout[1])
    else:
        print('ambiguous date found!')
    return metadata


def _xyztext2date(textdate):
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
        date = datetime.strptime(textdate, '%d %B %Y')
        numdate = date.strftime('%Y%m%d')
        return numdate
    except:
        try:
            date = datetime.strptime(textdate, '%d\%B\%Y')
            numdate = date.strftime('%Y%m%d')
            return numdate
        except:
            return 'unknown'


def _parse_sounding_frequency(line):
    """
    parse sounding frequency.
    Note: LOW & HIGH are usually settings for the
    single beam in New Orleans
    400kHz seems to be their multibeam.
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    name = line.split('SOUNDING_FREQUENCY==')[-1].strip('\n')
    metadata = {'sounding_frequency': name}
    return metadata


def _parse_survey_type(line):
    """
    returns survey type
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    name = line.split('SURVEY_TYPE==')[-1]
    name = name.strip('\n')
    metadata = {'text: survey_type': name}
    return metadata


def _parse_survey_crew(line):
    """
    returns survey crew
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    name = line.split('SURVEY_CREW==')[-1]
    name = name.strip('\n')
    metadata = {'survey_crew': name}
    return metadata


def _parse_sea_condition(line):
    """
    sea conditions
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    name = line.split('SEA_CONDITION==')[-1]
    name = name.strip('\n')
    metadata = {'sea_condition': name}
    return metadata


def _parse_vessel_name(line):
    """
    vessel name
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    name = line.split('VESSEL_NAME==')[-1]
    name = name.strip('\n')
    metadata = {'vessel_name': name}
    return metadata


def _parse_LWRP_(line):
    """
    Checks to see if its in Low Water Reference Plane
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    name = line.split('LWRP==')[-1]
    name = name.split('LWRP=')[-1]
    name = name.strip('\n')
    name = name.strip('\t')
    name = name.rstrip(' ')
    if name == 'N/A':
        metadata = {'LWRP': ''}
    elif name == 'NA':
        metadata = {'LWRP': ''}
    elif name == 'N/A0':
        metadata = {'LWRP': ''}
    else:
        metadata = {'LWRP': name}
        metadata = {'from_vert_datum': name}
        metadata = {'from_vert_key': 'LWRP'}
        print('Data in LWRP')
    return metadata


def _parse_Gage_Reading(line, allcap1):
    """
    Looks for the water level Gage
    
    Parameters
    ----------
    line :
        param allcap1:
    allcap1 :
    
    Returns
    -------
    
    """
    if allcap1 == 1:
        name = line.split('GAGE_READING==')[-1]
        name = name.strip('\n')
        metadata = {'GAGE_READING': name}
    if allcap1 == 2:
        name = line.split('Gage_Reading==')[-1]
        name = name.strip('\n')
        metadata = {'GAGE_READING': name}
    return metadata


def _parse_sound_velocity(line):
    """
    Looks for Sound Velocity
    
    Parameters
    ----------
    line :
    
    Returns
    -------
    
    """
    name = line.split('SOUND VELOCITY')[-1]
    name = name.strip('\n')
    metadata = {'sound_velocity': name}
    return metadata


def _parse_Ranges(line):
    """
    looks at ranges
    
    Parameters
    ----------
    line :
    
    Returns
    -------
    
    """
    name = line.split('Range:')[-1]
    name = name.strip('\n')
    metadata = {'Range': name}
    return metadata


def _is_RTK(line):
    """
    pulls RTK line
    
    Parameters
    ----------
    line :
    
    
    Returns
    -------
    
    """
    pattern_coordinates = '[RTK]'  # at least six digits# should be seven then . plus two digits
    if _re.findall(pattern_coordinates, line) is not None:
        return False
    else:
        return True


def _is_RTK_Tide(line):
    """
    looks for RTK Tide
    
    Parameters
    ----------
    line :
    
    
    Returns
    -------
    
    """
    if _re.findall('[VRS RTK TIDES]', line) is not None:
        return False
    else:
        return True


def _parse_Survey_Type(line):
    """
    parse Survey Type
    
    Parameters
    ----------
    line :
        
    
    Returns
    -------
    
    """
    metadata = {'Survey_Type==': line.split('Survey_Type==')[1]}
    return metadata


def _parse_Survey_Number(line):
    """
    parse Survey Number
    
    Parameters
    ----------
    line :
    
    
    Returns
    -------
    metadata
    
    """
    metadata = {'Survey_Number==': line.split('Survey_Number==')[1]}
    return metadata


##-----------------------------------------------------------------------------
def check_date_order(m, mm):
    """
    ingest dates from e-hydro file name, and xml if available
    do a date check.
    
    Parameters
    ----------
    m :
        param mm:
    mm :
    
    
    Returns
    -------
    
    """
    date_list = []  # date_list = [begdate, enddate,filename_date]
    if 'begdate' in m:
        # parser.parse(text_date, dayfirst=False)
        if m['begdate'] != '' and m['begdate'] is not None:
            est_begdate, ans1 = check_date_format_hasday(m['begdate'])
            if ans1 == 'yes':
                begdate = datetime.date(datetime.strptime(m['begdate'], '%Y%m%d'))
                date_list.append(begdate)
                # m['start_date'] = m['begdate']
    if 'enddate' in m:
        if m['enddate'] != '' and m['enddate'] is not None:
            est_enddate, ans1 = check_date_format_hasday(m['begdate'], int(
                '30'))  # may want better logic here other date modules can handle this better once situation flagged
            if ans1 == 'yes':
                enddate = datetime.date(datetime.strptime(m['enddate'], '%Y%m%d'))
                date_list.append(enddate)
    filename_date = datetime.date(datetime.strptime(mm['filename_date'], '%Y%m%d'))
    date_list.append(filename_date)
    if 'daterange' in m:
        try:  # catches ValueError returns if daterange ends up being unexpected
            next_date = check_abst_date(mm['filename_date'], m['daterange'])
            for day in next_date:
                date_list.append(day)
        except:
            print('unusual date format')
    date_list.sort()
    date_list2 = []
    for d in date_list:
        date_list2.append(datetime.strftime(d, '%Y%m%d'))
    m['start_date'] = date_list2[0]
    m['end_date'] = date_list2[-1]
    return m


##-----------------------------------------------------------------------------
def check_abst_date(filename_date, daterange):
    """
    check_abst_date(filename_date, daterange)
    Expecting values from:
    #filename_date = m['filename_date']
    #dateramge = xml_meta['daterange']
    
    Parameters
    ----------
    filename_date :
        param daterange:
    daterange :
    
    
    Returns
    -------
    
    """
    next_date = []
    mnum = ''
    XX = datetime.strptime(filename_date, '%Y%m%d')  # Create default value based on filename date
    if daterange != '' and daterange is not None:
        dates = []
        if '&' in daterange:
            dates = daterange.split('&')
        elif 'thru' in daterange:
            dates = daterange.split('thru')
        elif 'through' in daterange:
            dates = daterange.split('through')
        elif '-' in daterange:
            dates = daterange.split('-')
        if len(dates) > 0:
            date1 = parser.parse(dates[-1], parser.parserinfo(dayfirst=True), default=XX)
        if len(dates) > 1:
            # next_date.append(date1)
            splitters = ['&', '-', 'thru', 'through']
            dates2 = []
            for split1 in splitters:
                if split1 in dates[-1]:
                    date1 = parser.parse(dates[-1].split(split1)[-1], parser.parserinfo(dayfirst=True), default=XX)
                for days in dates:
                    if split1 in days:
                        # recalculate date1!
                        dates2.append(days.split(split1))
            for numday, day_ in enumerate(dates):
                if len(dates2) > 1:
                    if numday < len(dates) - 1:
                        # test for number list#
                        for m in months:
                            if m in day_.lower():
                                mnum = months_d[m]
                                temp_date = date1.replace(month=int(mnum))  # changed month
                                # day_.replace(m,'')
                                # month is keyword for funciton, as is day, year
                    if numday < len(dates) - 1:
                        for day_ in dates2[numday]:
                            m1 = months_bynum_d[temp_date.strftime('%m')]
                            day_ = day_.lower().replace(m1, '')
                            if ',' in day_:
                                days = day_.split(',')
                                for day_ in days:  # try to reduce to calendar day integers for input
                                    day_ = day_.strip('on').replace('from', '').replace('of', '').strip()
                                    if day_ != '':
                                        if len(mnum) > 0:
                                            next_date.append(temp_date.replace(day=int(day_)))
                                        else:
                                            next_date.append(date1.replace(day=int(day_)))
                            else:
                                day_ = day_.strip('on').replace('from', '').replace('of', '').strip()
                                if len(mnum) > 0:
                                    t = temp_date.replace(day=int(day_))
                                    next_date.append(t)
                                else:
                                    next_date.append(date1.replace(day=int(day_)))
                    else:  # last section
                        for i, day_ in enumerate(dates2[numday]):
                            if i < len(dates2[-1]) - 1:
                                if ',' in day_:
                                    days = day_.split(',')
                                    for day_ in days:  # try to reduce to calendar day integers for input
                                        day_ = day_.replace('on', '').replace('from', '').replace('of', '').strip()
                                        if day_ != '':
                                            next_date.append(date1.replace(day=int(day_)))
                                else:
                                    next_date.append(date1.replace(
                                        day=int(day_.replace('on', '').replace('from', '').replace('of', '').strip())))
                            else:
                                next_date.append(date1)
    next_date = check_datelist(next_date)  # convert from datetime to date format
    return next_date


def check_datelist(next_date):
    """
    checks date list
    
    Parameters
    ----------
    next_date :
    
    
    Returns
    -------
    dateonly_list
    
    """
    dateonly_list = []
    for day in next_date:
        day = datetime.date(day)
        dateonly_list.append(day)
    return dateonly_list


def check_date_format_hasday(date_string, b_or_e=None):
    """
    check_date_format_hasday
    
    Parameters
    ----------
    date_string :
        param b_or_e:  (Default value = None)
    b_or_e :
         (Default value = None)
    
    Returns
    -------
    
    """
    pattern_missing_valid_day = '[\d][\d][\d][\d][\d][\d][0][0]'
    if b_or_e is None:
        day = '1'
    elif type(b_or_e) == int:
        day = str(b_or_e)
    else:
        day = '1'
    if _re.match(pattern_missing_valid_day, date_string) is not None:
        args = _re.search(pattern_missing_valid_day, date_string)
        end_position = args.endpos
        d_l = list(date_string)
        d_l[end_position - 1] = day  # '1' default value
        date_st1 = "".join(d_l)
        ans1 = 'no'
    else:
        date_st1 = date_string
        ans1 = 'yes'
    return date_st1, ans1


months = ['jan', 'january', 'feb', 'february', 'mar', 'march', 'apr', 'april', 'may', 'jun', 'june', 'jul', 'july',
          'aug', 'august', 'sep', 'september', 'oct', 'october', 'nov', 'november', 'dec',
          'december']  # Check for other months
months_d = {'jan': '01', 'january': '01', 'feb': '02', 'february': '02', 'mar': '03', 'march': '03', 'apr': '04',
            'april': '04', 'may': '05', 'jun': '06', 'june': '06', 'jul': '07', 'july': '07', 'aug': '08',
            'august': '08', 'sep': '09', 'september': '09', 'oct': '10', 'october': '10', 'nov': '11', 'november': '11',
            'dec': '12', 'december': '12'}
months_bynum_d = {val1: key1 for (key1, val1) in months_d.items()}  # swap keys and values to new dictionary
##-----------------------------------------------------------------------------
