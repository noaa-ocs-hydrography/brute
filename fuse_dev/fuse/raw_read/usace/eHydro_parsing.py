# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 17:26:38 2019

jkinney 
eHydro_Parsing

Taken from:
    extract_ehydro_meta.py

    Created on Wed Aug  8 10:11:36 2018
    
    @author: grice
    
    Parse the eHydro xyz header using the NY region's files as an example.  The
    metadata is returned as a dictionary with keys that match the meta2csv naming
    convention where required.
"""

import re as _re

def _is_header(line):
    """
    Test if a line contains anything other than numbers it is a meta data line.
    ###THIS METHOD DOES NOT ALWAYS WORK!!!!!##### THIS MAY WORK IN SOME DISTRICTS MOST OF THE TIME, BUT NOT ALL SCENARIOS###
    EX:
    3932544.19,265747.33,-29.83,R01MS04_0035.pst,-1.90
    from: SW_04_SWP_20181005_CS.xyz
    """
    pattern = '[a-zA-Z]'
    if _re.search(pattern, line) is None:
        return False
    else:
        return True


def _is_header2(line, version = None):
    if version == None:
        v=0
    if version == 'CEMVN':
        pattern_coordinates = '[\d\][\d\][\d\][\d\][\d\][\d\]'#at least six digits# should be seven then . plus two digits
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
           # m = _re.match(pattern, line)
            #X = m.span()#returns the start and end position of the match#for match the starts is always 0
            #print(X)
            return True

##pattern ='[\d][][][][][][\d\]

def _parse_note(line):
    """
    Parse the notes line.
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
        metadata['from_wkt'] = _fips2wkt(fips)
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
        metadata['vert_key'] = 'unknown'
    vert_units_tags = ['NAVD88','NAVD1988','NAVD 1988']
    for tag in vert_units_tags:
        vert_units_end = line.find(tag) 
        if vert_units_end >= 0:
            vert_units_end += len(tag)
            break
        else:
            vert_units_end = 0
    vert_units_start = vert_units_end - line[vert_units_end::-1].find('>krb<')
    vert_units = line[vert_units_start+1:vert_units_end]
    metadata['from_vert_datum'] = vert_units
    if vert_units.find('FEET') > 0:
        metadata['from_vert_units'] = 'US Survey Foot'
    else:
        metadata['from_vert_units'] = 'unknown'
    return metadata

def _parse_projectname(line):
    """
    Parse the project name line.
    """
    name = line.split('=')[-1]
    name = name.strip('\n')
    metadata = {'projectname' : name}
    return metadata

def _parse_surveyname(line):
    """
    Parse the survey name line.
    """
    name = line.split('=')[-1]
    name = name.strip('\n')
    metadata = {'surveyname' : name}
    return metadata

def _parse_survey_condition(line):
    name = line.split('VESSEL_NAME==')[-1]
    name = name.strip('\n')
    metadata = {'vessel_name' : name}
    return metadata

    
def _parse_LWRP_(line):
    name = line.split('LWRP==')[-1]
    name = name.strip('\n')
    if name == 'N/A':
        metadata = {'LWRP' : ''}
        #pass
    elif name =="":
        metadata = {'LWRP' : ''}
    else:
        metadata = {'LWRP' : name}
        metadata = {'script: from_vert_datum':name}
        metadata = {'script: from_vert_key' : 'LWRP'}
        print('Data in LWRP')
    return metadata 

def _parse_sea_condition(line):
    name = line.split('SEA_CONDITION==')[-1]
    name = name.strip('\n')
    metadata = {'sea_condition' : name}
    return metadata     
        
def _parse_survey_type(line):
    name = line.split('SURVEY_TYPE==')[-1]
    name = name.strip('\n')
    metadata = {'text: survey_type' : name}
    return metadata
            
def _parse_survey_crew(line):
    name = line.split('SURVEY_CREW==')[-1]
    name = name.strip('\n')
    metadata = {'survey_crew' : name}
    return metadata     
        
def _parse_Gage_Reading(line, allcap1):
    if allcap1 == 1:
        name = line.split('GAGE_READING==')[-1]
        name = name.strip('\n')
        metadata = {'GAGE_READING' : name}
    if allcap1 == 2:
        name = line.split('Gage_Reading==')[-1]
        name = name.strip('\n')
        metadata = {'GAGE_READING' : name}
    return metadata  
            
def _parse_vessel_name(line):
    name = line.split('SEA_CONDITION==')[-1]
    name = name.strip('\n')
    metadata = {'survey_condition' : name}
    return metadata

def _parse_sound_velocity(line):
    name = line.split('SOUND VELOCITY')[-1]
    name = name.strip('\n')
    metadata = {'sound_velocity' : name}

def _is_RTK(line, version = None):
    if version == None:
        v=0
    if version == 'CEMVN':
        pattern_coordinates = '[RTK]'#at least six digits# should be seven then . plus two digits
        if _re.findall(pattern_coordinates, line) is not None:
            return False
        else:
            return True
        
def _is_RTK_Tide(line, version = None):
    if version == None:
        v=0
    if version == 'CEMVN':
        if _re.findall('[VRS RTK TIDES]', line) is not None:
            return False
        else:
            return True

def _parse_sounding_frequency(line, version = None):
    if version == None:
        v = 0
    if version == 'CEMVN':
        v = 1 
    name = line.split('SOUNDING_FREQUENCY==')[-1].strip('\n')
    metadata = {'sounding_frequency' : name}
    return metadata


def _parse_surveydates(line):
    """
    Parse the project dates line.
    """
    metadata = {}
    datestr = line.split('=')[-1]
    datestr = datestr.strip('\n')
    if datestr.find('-') > 0:
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
    """
    try:
        date = _datetime.strptime(textdate, '%d %B %Y')
        numdate=date.strftime('%Y%m%d')
        return numdate
    except:
        return 'unknown'

def convert_tofips(SOURCEPROJECTION_dict, SOURCEPROJECTION):
    """
    FIPS = SOURCEPROJECTION_dict[SOURCEPROJECTION]
    """
    FIPS = SOURCEPROJECTION_dict[SOURCEPROJECTION]
    return FIPS

SOURCEPROJECTION_dict = {}
SOURCEPROJECTION_dict = {
     'Alabama East' : '0101',
     'Alabama West' : '0102',
     'Alabama_West' : '0102',
     'Alaska 1' : '5001',
     'Alaska 2' : '5002',
     'Alaska 3' : '5003',
     'Alaska 4' : '5004',
     'Alaska 5' : '5005',
     'Alaska 6' : '5006',
     'Alaska 7' : '5007',
     'Alaska 8' : '5008',
     'Alaska 9' : '5009',
     'Alaska 10' : '5010',
     'California I' : '0401',
     'California II' : '0402',
     'California III' : '0403',
     'California IV' : '0404',
     'California V' : '0405',
     'California VI' : '0406',
     'Connecticut' : '0600',
     'Delaware' : '0700',
     'Florida East' : '0901',
     'Florida North' : '0903',
     'Florida West' : '0902',
     'Florida_North' : '0903',
     'Georgia East' : '1001',
     'Georgia West' : '1002',
     'Hawaii 1' : '5101',
     'Hawaii 2' : '5102',
     'Hawaii 3' : '5103',
     'Hawaii 4' : '5104',
     'Hawaii 5' : '5105',
     'Illinois East' : '1201',
     'Illinois West' : '1202',
     'Illinois_East' : '1201',
     'Illinois_West' : '1202',
     'Indiana East' : '1301',
     'Indiana West' : '1302',
     'Iowa_North' : '1401',
     'Iowa_South' : '1402',
     'Kentucky North' : '1601',
     'Kentucky South' : '1602',
     'Louisiana North' : '1701',
     'Louisiana South' : '1702',
     'Maine East' : '1801',
     'Maine West' : '1802',
     'Maryland' : '1900',
     'Massachusetts Island' : '2002',
     'Massachusetts Mainland' : '2001',
     'Michigan North' : '2111',
     'Michigan Central' : '2112',
     'Michigan South' : '2113',
     'Minnesota_Central' : '2202',
     'Minnesota Central' : '2202',
     'Minnesota_North' : '2201',
     'Minnesota North' : '2201',
     'Minnesota_South' : '2203',
     'Minnesota South' : '2203',
     'Mississippi East' : '2301',
     'Mississippi_East' : '2301',
     'Mississippi_West' : '2302',
     'Mississippi West' : '2302',
     'Missouri West' : '2403',
     'Missouri Central' : '2402',
     'Missouri East' : '2401',
     'Missouri_East' : '2401',
     'New Hampshire' : '2800',
     'New Jersey' : '2900',
     'New York Long Island' : '3104',
     'New_Jersey' : '2900',
     'New York Central' : '3102',
     'New York West' : '3103',
     'New_York_East' : '3101',
     'New York East' : '3101',
     'New_York_Long_Island' : '3104',
     'North Carolina' : '3200',
     'Ohio North' : '3401',
     'Ohio_South' : '3402',
     'Oregon North' : '3601',
     'Oregon South' : '3602',
     'Puerto Rico Virgin Islands' : '5200',
     'Rhode Island' : '3800',
     'South Carolina' : '3900',
     'Texas North' : '4201',
     'Texas North Central' : '4202',
     'Texas Central' : '4203',
     'Texas South' : '4205',
     'Texas South Central' : '4204',
     'Virginia North' : '4501',
     'Virginia South' : '4502',
     'Washington North' : '4601',
     'Washington South' : '4602',
     'West_Virginia_North' : '4701',
     'West_Virginia_South' : '4702',
     'Wisconsin Central' : '4802',
     'Wisconsin North' : '4801',
     'Wisconsin South' : '4803'}
#Please note all of the above are in State Plane, and the corresponding FIPS code
#
    
#def convert_tofips(self, SOURCEPROJECTION):
#    """
#    FIPS = convert_tofips(SOURCEPROJECTION)
#    """
#    FIPS = self.SOURCEPROJECTION_dict(SOURCEPROJECTION)
#    return FIPS