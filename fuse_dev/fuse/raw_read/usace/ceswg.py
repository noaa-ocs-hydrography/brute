# -*- coding: utf-8 -*-
"""
Edited by Juliet Kinney
extract_ehydro_meta_class_CESWG.py

Galveston uses MLLW in the name of some files, but not all?
Many are in MLLW even if its not in the filename.

This script takes as an input the filename and path of an USACE e-Hydro .xyz 
text file and pull the metadata from the xyz header and file name.
Additionally, it passes the matching .xml file and utilizes 
the parse_usace_xml script to pull out the relevant metadata if it
is either FGDC or ISO FGDC USACE metadata format.

update 4/17,4/5/19
major update April 2, 2019

"""
__version__ = 'FUSE'
import os as os
import pickle as _pickle
import re as _re

_ussft2m = 0.30480060960121924 # US survey feet to meters
from datetime import datetime
import numpy as _np 
try:
    import fuse.raw_read.usace.parse_usace_xml as p_usace_xml
except:
    try:
        from . import parse_usace_xml as p_usace_xml
    except:
        print('importing fuse.raw_read.usace.parse_usace_xml as p_usace_xml did not work') 
##-----------------------------------------------------------------------------

class read_raw:
    """
    This class passes back bathymetry
    & a metadata dictionary from the e-Hydro files

    Parameters
    ----------

    Returns
    -------

    """
    
    def read_metadata(self, infilename):
        """
        Read all available meta data.
        returns dictionary

        Parameters
        ----------
        infilename :
            

        Returns
        -------

        """
        version='CESWG'
        self.version = version
        return retrieve_meta_for_Ehydro_out_onefile(infilename)#
    
    def read_bathymetry_dat(self, infilename):
        """
        Read the bathymetry from the .dat file. The dat file is less precise,
        but had no header and is in a standardized format

        Parameters
        ----------
        infilename :
            

        Returns
        -------

        """
        # get the dat file for CESWG# Galveston
        stub, ext = os.path.splitext(infilename)
        bathyfilename = stub + '.dat'
        xyz = _np.loadtxt(bathyfilename, delimiter = ' ')
        self.xyz
        return xyz

    def read_bathymetry(self, infilename):
        """
        Read the bathymetry from the xyz files, this tells it to not include
        the header when reading the file
        
        Note: The high resolution multibeam files are available as .xyz on E-Hydro

        Parameters
        ----------
        infilename :
            

        Returns
        -------

        """
        version='CESWG'
        self.version = version
        first_instance, commas_present = _start_xyz(infilename)
        if first_instance != '':
            if commas_present == ',':
                xyz =  _np.loadtxt(infilename, delimiter = ',', skiprows = first_instance, usecols=(0,1,2))
            else:
                xyz = _np.loadtxt(infilename, delimiter = ' ', skiprows = first_instance, usecols=(0,1,2))
            #xyz = _np.loadtxt(infilename, delimiter = ',', skiprows = first_instance)
        #else:
            #xyz = _np.loadtxt(infilename, delimiter = ',')
        return xyz        

#------------------------------------------------------------------------------
def return_surveyid(filenamepath, ex_string):
    """
    strip end of filename off
    surveybasename =return_surveyid(filenamepath, ex_string)

    Parameters
    ----------
    filenamepath :
        param ex_string:
    ex_string :
        

    Returns
    -------

    """
    basename = os.path.basename(filenamepath)
    surveybasename = basename.rstrip(ex_string)
    return surveybasename    
#------------------------------------------------------------------------------

def retrieve_meta_for_Ehydro_out_onefile(filename):
    """
    retrieve metadata for USACE E-Hydro files
    function returns metadata dictionary
    
    input is filename of .xyz file with path

    Parameters
    ----------
    filename :
        

    Returns
    -------

    """
    #next if pull the subset of the table in the dataframe related to the list of files passed to it.
    merged_meta = {}
    merge2 = {}      
    f = filename
    basename = os.path.basename(f)
    ex_string1 = '*_A.xyz'
    ex_string2 = '*_FULL.xyz'
    ex_string3 = '*_FULL.XYZ'
    ex_string4 = '*_A.XYZ'
    basename = os.path.basename(basename)
    basename = return_surveyid(basename, ex_string1)
    basename = return_surveyid(basename, ex_string2)
    basename = return_surveyid(basename, ex_string3)
    basename = return_surveyid(basename, ex_string4)
    basename = basename.rstrip('.XYZ')  
    basename = basename.rstrip('.xyz')    
    #empty dictionary place holder for future ehydro table ingest (make come from imbetween source TBD)
    meta_from_ehydro={}
    e_t = Extract_Txt(f)
    # xml pull here.
    xmlfilename = get_xml_match(f)
    if os.path.isfile(xmlfilename):
        with open(xmlfilename, 'r') as xml_file:
            xml_txt = xml_file.read()
        xmlbasename = os.path.basename(xmlfilename)
        xml_data = p_usace_xml.XML_Meta(xml_txt, filename = xmlbasename)
        if xml_data.version == 'USACE_FGDC':
            meta_xml = xml_data._extract_meta_USACE_FGDC()#CEMVN()
        elif xml_data.version == 'ISO-8859-1':
            meta_xml = xml_data._extract_meta_USACE_ISO()
            if 'ISO_xml' not in meta_xml:
                meta_xml = xml_data._extract_meta_USACE_FGDC(override = 'Y')#xml_data._extract_meta_ISOlabel_USACE_FGDC()
        else:
            meta_xml = xml_data.convert_xml_to_dict2()
        ext_dict = xml_data.extended_xml_fgdc()
        ext_dict =  p_usace_xml.ext_xml_map_enddate(ext_dict)
    else:
        ext_dict = {}
        meta_xml = {}
    meta = e_t.parse_ehydro_xyz(f, meta_source = 'xyz', version='CESWG', default_meta = '')#
    list_keys_empty =[]
    combined_row = {}
    subset_row = {}
    subset_no_overlap = {}
    subset_dict ={}
    for key in meta:
        if meta[key] == 'unknown' or  meta[key] == '':
            list_keys_empty.append(key)
        else:
            subset_row[key] = meta[key]
            #non blank columns only
            if key in meta_xml:
                if meta[key] == meta_xml[key]:
                    combined_row[key] = meta[key]
                    """
                    only make list within cell if values from different sources are different
                    """
                else:
                    combined_row[key] = meta[key] + ' , ' + meta_xml[key]
            else:
                subset_no_overlap[key] = meta[key]
    for key in ext_dict:
        if ext_dict[key] == 'unknown' or  ext_dict[key] == '' or ext_dict[key] == None:
            list_keys_empty.append(key)
        else:
            if key in meta_xml:
                if meta_xml[key] == ext_dict[key]:
                    combined_row[key] = ext_dict[key]
                    """
                    only make list within cell if values from different sources are different
                    """
                else:
                    combined_row[key] = ext_dict[key] + ' , ' + meta_xml[key]
            else:
                subset_dict[key] = ext_dict[key]
    merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row } #this one excluded 'unknown' keys, and 
    #in merging sources from the text file and xml it will show any values that do not match as a list.
    merged_meta = { **meta, **meta_from_ehydro,**meta_xml }#this method overwrites
    merged_meta = check_date_order(merged_meta, merged_meta)
    return merged_meta

###---------------------------------------------------------------------------- 
class Extract_Txt(object):
    """Extract both information from the filename as well as from the text file's header"""
    def __init__(self, preloadeddata, version = '', filename = ''):
        self.filename = preloadeddata
        if filename != "" or None:
            self.filename_1 = filename
            self.errorfile = os.path.dirname(filename) + 'TEST_extract_ehdyro_meta_class_CESWG_ErrorFile1.txt'
        else:
            self.errorfile = os.path.dirname(filename) + 'Default_extract_ehdyro_meta_class_CESWG_error.txt'       

    def parse_ehydro_xyz(self, infilename, meta_source = 'xyz', version= 'CESWG', default_meta = ''):#need to change version to None
        """
        'CESWG'

        Parameters
        ----------
        infilename :
            param meta_source:  (Default value = 'xyz')
        version :
            Default value = 'CESWG')
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
        #elif meta_source == 'xml':
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
        sorind = (name_meta['projid'] + '_' + 
                  name_meta['uniqueid'] + '_' + 
                  name_meta['subprojid'] + '_' + 
                  name_meta['start_date'] + '_' + 
                  name_meta['statuscode'])
        merged_meta['source_indicator'] = 'US,US,graph,' + sorind
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
                    'from_path' : infilename,
                    'from_filename' : base,
                    'projid' : splitname[0],
                    'uniqueid' : splitname[1],
                    'subprojid' : splitname[2],
                    'start_date' : splitname[3],
                    }
            if len(splitname) >4:
                meta['statuscode']=splitname[4]
                if len(splitname) > 5:
                    option = splitname[5]
                    if len(splitname) > 6:
                        for n in range(6, len(splitname)):
                            option = option + '_' + splitname[n]
                    meta['optional'] = option
            else:
                meta['statuscode']=''               
        else:
            print(name + ' appears to have a nonstandard naming convention.')
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
        # get the header
        if version=='CESWG':
            with open(infilename, 'r') as infile:
                for line in infile.readlines():
                    if line == '\n':
                        continue
                    elif _is_header2(line):
                        header.append(line)
                    else:
                        break
                    #do header check:
                    #within header:
            for line in header:
                if line.startswith('notes_chart== 1.'):
                    metalist.append(_parse_note(line))
                    metalist.append(_parse_notes_chart(line))
                    #do something with this line metaline.append()
                    #tokens = line.split('\n')#tokens= ['value', 'value2', etc]
                elif line.startswith('Notes_chart== 1.'):
                    metalist.append(_parse_note(line))
                    metalist.append(_parse_notes_chart(line))
                elif line.startswith('notes_chart=='):
                    metalist.append(_parse_note(line))
                    metalist.append(_parse_notes_chart(line))
                if line.startswith('ProcessedBy=='):
                    metalist.append(_parse_processedBy(line))      
                if line.startswith('CheckedBy=='):                
                    metalist.append(_parse_CheckedBy(line))
                if line.startswith('ReviewedBy=='):
                    metalist.append(_parse_ReviewedBy(line))
            meta = {}
        try:
            for m in metalist:
                meta = {**meta, **m}
            return meta
        except:
            meta = {}
            errorfile = self.errorfile 
            with open(errorfile,'a') as metafail:
                metafail.write(infilename + '\n')
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
            default_meta = os.path.join(path,'default.pkl')
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
    xml_name = basef + '.xml'
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
    if filename[-end_len:] == extension:
        basef =  filename[:-end_len]
    else:
        basef = filename
    xml_name = basef + '.xml'
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
    if '_A.xyz' in f:
        xmlfilename = get_xml_xt(f,'_A.xyz')
    elif '_FULL.xyz' in f:
        xmlfilename = get_xml_xt(f,'_FULL.xyz')
    elif '_FULL.XYZ' in f:
        xmlfilename = get_xml_xt(f,'_FULL.XYZ')
    elif '_A.XYZ' in f:
        xmlfilename = get_xml_xt(f,'_A.XYZ')
    else:
        xmlfilename = get_xml(f)
    return xmlfilename
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
    pattern_coordinates = '[\d][\d][\d][\d][\d][\d]'#at least six digits# should be seven then . plus two digits
    with open(infilename, 'r') as infile:
        for (index1, line) in enumerate (infile):
            if _re.match(pattern_coordinates, line) is not None:
                numberofrows.append(index1)
                if line.find(',')>0:
                    commas_present=','
        first_instance = numberofrows[0]
        return first_instance, commas_present
    return first_instance, commas_present

def _is_header2(line, version = None):
    """
    

    Parameters
    ----------
    line :
        param version:  (Default value = None)
    version :
         (Default value = None)

    Returns
    -------

    """
    if version == None:
        version = ''
    if version == 'CESWG':
        pattern_coordinates = '[\d][\d][\d][\d][\d][\d]'#at least six digits# should be seven then . plus two digits
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
        pattern ='[^0-9]'#anything except 0-9
        if _re.match(pattern, line) is None:#does the string start with this pattern?
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
    metadata = {'projectname' : name}
    return metadata

def _parse_notes_chart(line):
    """
    

    Parameters
    ----------
    line :
        

    Returns
    -------

    """
    lines = line.split('\\n')
    metadata = {}
    metadata ['notes_chart']= line
    for aline in lines:
        if aline != '':
            if aline.find('ALL ELEVATIONS SHOWN ARE REFERENCED') >= 0:                
                if aline.find('FEET') >= 0:
                    metadata['from_vert_units'] = 'US Survey Foot'
                else:
                    metadata['from_vert_units']= aline.split('ALL ELEVATIONS SHOWN ARE REFERENCED')[-1]
            if aline.find('COORDINATES ARE REFERENCED TO') >= 0:
                metadata['horiz_sys'] = aline.split('COORDINATES ARE REFERENCED TO')[-1]
            if aline.find('COORDINATE SYSTEM') >= 0:
                metadata['COORDINATE SYSTEM'] = aline.split('COORDINATE SYSTEM')[-1]              
            if aline.find('SURVEY VESSEL') >= 0:
                metadata['SURVEY VESSEL'] = aline.split('SURVEY VESSEL:')[-1]
            if aline.find('SURVEY DATE:') >= 0:
                metadata['SURVEY DATE'] = aline.split(':')[-1]
                #metadata =_parse_surveydates(aline.split(':')[-1])
            if aline.find('SURVEYED BY:') >= 0:
                metadata['SURVEYED_BY'] = aline.split(':')[-1]#metadata = _split_at_colon(key, line)
            if aline.find('FREQUENCY SOUNDINGS') >= 0:
                metadata['FREQUENCY SOUNDINGS'] = aline.split(':')[-1]
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
#    zone_idx = line.find('ZONE')
#    zone_len = line[zone_idx:].find('.')
#    horiz_datum = line[zone_idx:zone_idx + zone_len]
#    if len(horiz_datum) > 0:
#        fips = horiz_datum.split()[1]
#        fips = fips.rstrip(',')
#        metadata['from_fips'] = fips
#        horiz_units = horiz_datum.split(',')[1]
#        if horiz_units.lstrip(' ') == 'US SURVEY FEET':
#            metadata['from_horiz_units'] = 'US Survey Foot'
#        else:
#            metadata['from_horiz_units'] = horiz_units.lstrip(' ')
#        metadata['from_horiz_datum'] = horiz_datum
#    else:
#        metadata['from_horiz_units'] = 'unknown'
#        metadata['from_horiz_datum'] = 'unknown'
    # find the vertical datum information
    if line.find('MEAN LOWER LOW WATER') >= 0:
        metadata['from_vert_key'] = 'MLLW'
        metadata['script: from_vert_key'] = 'MLLW'
    elif line.find('MLLW') >= 0:
        metadata['from_vert_key'] = 'MLLW'
        metadata['script: from_vert_key'] = 'MLLW'
    elif line.find('MEAN LOW WATER') >= 0:
        metadata['from_vert_key'] = 'MLW'
        metadata['script: from_vert_key'] = 'MLW'
    else:
        metadata['vert_key'] = 'unknown'
#    vert_units_tags = ['NAVD88','NAVD1988','NAVD 1988']
#    for tag in vert_units_tags:
#        vert_units_end = line.find(tag) 
#        if vert_units_end >= 0:
#            vert_units_end += len(tag)
#            break
#        else:
#            vert_units_end = 0
#    vert_units_start = vert_units_end - line[vert_units_end::-1].find('>krb<')
#    vert_units = line[vert_units_start+1:vert_units_end]
#    metadata['from_vert_datum'] = vert_units
#    if vert_units.find('FEET') >= 0:
#        metadata['from_vert_units'] = 'US Survey Foot'
#    else:
#        metadata['from_vert_units'] = 'unknown'
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
    metadata = {'surveyname' : name}
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
        numdate=date.strftime('%Y%m%d')
        return numdate
    except:
        try:
            date = datetime.strptime(textdate, '%d\%B\%Y')
            numdate=date.strftime('%Y%m%d')
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
    metadata = {'sounding_frequency' : name}
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
    metadata = {'text: survey_type' : name}
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
    metadata = {'survey_crew' : name}
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
    metadata = {'sea_condition' : name}
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
    metadata = {'vessel_name' : name}
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
        metadata = {'LWRP' : ''}
    elif name == 'NA':
        metadata = {'LWRP' : ''}
    elif name == 'N/A0':
        metadata = {'LWRP' : ''}
    else:
        metadata = {'LWRP' : name}
        metadata = {'from_vert_datum':name}
        metadata = {'from_vert_key' : 'LWRP'}
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
        metadata = {'GAGE_READING' : name}
    if allcap1 == 2:
        name = line.split('Gage_Reading==')[-1]
        name = name.strip('\n')
        metadata = {'GAGE_READING' : name}
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
    metadata = {'sound_velocity' : name}
    return metadata

def _parse_Ranges(line):
    """
    

    Parameters
    ----------
    line :
        

    Returns
    -------

    """
    name = line.split('Range:')[-1]
    name = name.strip('\n')
    metadata = {'Range' : name}
    return metadata

def _is_RTK(line):
    """
    

    Parameters
    ----------
    line :
        

    Returns
    -------

    """
    pattern_coordinates = '[RTK]'#at least six digits# should be seven then . plus two digits
    if _re.findall(pattern_coordinates, line) is not None:
        return False
    else:
        return True
        
def _is_RTK_Tide(line):
    """
    

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

def _parse_processedBy(line):
    """
    

    Parameters
    ----------
    line :
        

    Returns
    -------

    """
    metadata = {'ProcessedBy' : line}
    return metadata

def _parse_CheckedBy(line):
    """
    

    Parameters
    ----------
    line :
        

    Returns
    -------

    """
    metadata = {'CheckedBy':line.split('CheckedBy==')[1]}
    return metadata

def _parse_ReviewedBy(line):
    """
    

    Parameters
    ----------
    line :
        

    Returns
    -------

    """
    metadata = {'ReviewedBy':line.split('ReviewedBy==')[1]}
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
    date_list = []#date_list = [begdate, enddate,filename_date]
    if 'begdate' in m:
        #parser.parse(text_date, dayfirst=False)
        if  m['begdate'] != '' and  m['begdate'] != None:
            est_begdate, ans1 = check_date_format_hasday(m['begdate'])
            if ans1 == 'yes':
                begdate = datetime.date(datetime.strptime(m['begdate'],'%Y%m%d'))
                date_list.append(begdate)
            #begdate = datetime.date(datetime.strptime(m['begdate'],'%Y%m%d'))
            #date_list.append(begdate)
            #m['start_date'] = m['begdate']
    if 'enddate' in m:
        if  m['enddate'] != '' and  m['enddate'] != None:
            est_enddate, ans1 = check_date_format_hasday(m['begdate'],int('30'))#may want better logic here other date modules can handle this better once situation flagged
            if ans1 == 'yes':
                enddate = datetime.date(datetime.strptime(m['enddate'],'%Y%m%d'))
                date_list.append(enddate)
    filename_date = datetime.date(datetime.strptime(mm['filename_date'],'%Y%m%d'))
    date_list.append(filename_date)
    if 'daterange' in m:
        next_date = check_abst_date(mm['filename_date'], m['daterange'])
        for day in next_date:
            date_list.append(day)
    #if 'daterange' in m:
    #    next_date = check_abst_date(mm['filename_date'], m['daterange'])
    #    for day in next_date:
    #        date_list.append(day)
    date_list.sort()
    date_list2 = []
    for d in date_list:
        date_list2.append(datetime.strftime(d,'%Y%m%d'))
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
    XX=datetime.strptime(filename_date,'%Y%m%d')#Create default value based on filename date
    if daterange != '' and daterange != None:
        dates=[]
        if '&' in daterange:
            dates = daterange.split('&')
        elif 'thru' in daterange:
            dates = daterange.split('thru')
        elif 'through' in daterange:
            dates = daterange.split('through')
        elif '-' in daterange:
            dates = daterange.split('-')
        if len(dates) > 0:
            date1 = parser.parse(dates[-1],parser.parserinfo(dayfirst=True), default=XX)
        if len(dates) >1:
            #next_date.append(date1)
            splitters =['&', '-', 'thru', 'through']
            dates2=[]
            for split1 in splitters:
                if split1 in dates[-1]:
                    date1 = parser.parse(dates[-1].split(split1)[-1],parser.parserinfo(dayfirst=True), default=XX)
                for days in dates:
                    if split1 in days:
                        #recalculate date1!
                            dates2.append(days.split(split1))
            for numday, day_ in enumerate(dates):
                if len(dates2)>1:
                    if numday< len(dates)-1:
                            #test for number list#
                            for m in months:
                                if m in day_.lower():
                                    mnum=months_d[m]
                                    temp_date = date1.replace(month=int(mnum))#changed month
                                    #day_.replace(m,'')
                                    #month is keyword for funciton, as is day, year
                    if numday< len(dates)-1:
                        for day_ in dates2[numday]:
    
                                m1 = months_bynum_d[temp_date.strftime('%m')]
                                day_ = day_.lower().replace(m1, '')
                                if ',' in day_:
                                    days = day_.split(',')
                                    for day_ in days:#try to reduce to calendar day integers for input
                                        day_ = day_.strip('on').replace('from', '').replace('of', '').strip()
                                        if day_ != '':
                                            if len(mnum)>0:
                                                next_date.append(temp_date.replace(day=int(day_)))
                                            else:
                                                next_date.append(date1.replace(day=int(day_)))
                                else:
                                    day_ = day_.strip('on').replace('from', '').replace('of', '').strip()
                                    if len(mnum)>0:
                                        t=temp_date.replace(day=int(day_))
                                        next_date.append(t)
                                    else:
                                        next_date.append(date1.replace(day=int(day_)))                                        
                    else:#last section
                        for i, day_ in enumerate(dates2[numday]):
                            if i< len(dates2[-1])-1:
                                if ',' in day_:
                                    days = day_.split(',')
                                    for day_ in days:#try to reduce to calendar day integers for input
                                        day_ = day_.replace('on', '').replace('from', '').replace('of', '').strip()
                                        if day_ != '':
                                            next_date.append(date1.replace(day=int(day_)))
                                else:
                                     next_date.append(date1.replace(day=int(day_.replace('on', '').replace('from', '').replace('of', '').strip())))
                            else:
                                next_date.append(date1)

    next_date = check_datelist(next_date)#convert from datetime to date format
    return next_date

def check_datelist(next_date):
    """
    

    Parameters
    ----------
    next_date :
        

    Returns
    -------

    """
    dateonly_list =[]
    for day in next_date:
        day = datetime.date(day)
        dateonly_list.append(day)
    return dateonly_list

def check_date_format_hasday(date_string, b_or_e =None):
    """
    

    Parameters
    ----------
    date_string :
        param b_or_e:  (Default value = None)
    b_or_e :
         (Default value = None)

    Returns
    -------

    """
    pattern_missing_valid_day='[\d][\d][\d][\d][\d][\d][0][0]'

    if b_or_e == None:
        day = '1'
    elif type(b_or_e) == int:
        day = str(b_or_e)
    else:
        day = '1'
    if _re.match(pattern_missing_valid_day, date_string) is not None:
        args= _re.search(pattern_missing_valid_day, date_string)
        end_position=args.endpos
        d_l = list(date_string)
        d_l[end_position-1]=day#'1' default value
        date_st1 ="".join(d_l)
        ans1= 'no'
    else:
        date_st1 = date_string
        ans1 ='yes'
    return date_st1, ans1

months=['jan', 'january', 'feb', 'february', 'mar', 'march', 'apr', 'april', 'may', 'jun', 'june', 'jul', 'july', 'aug', 'august', 'sep', 'september', 'oct', 'october', 'nov', 'november', 'dec', 'december']#Check for other months
months_d={'jan':'01',  'january':'01',  'feb':'02',  'february':'02',  'mar':'03',  'march':'03',  'apr':'04',  'april':'04',  'may':'05',  'jun':'06',  'june':'06',  'jul':'07',  'july':'07',  'aug':'08',  'august':'08',  'sep':'09',  'september':'09',  'oct':'10',  'october':'10',  'nov':'11',  'november':'11',  'dec':'12',  'december':'12'}
months_bynum_d={val1:key1 for (key1, val1) in months_d.items()}#swap keys and values to new dictionary
##-----------------------------------------------------------------------------