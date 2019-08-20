# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:57:56 2019

@author: jkinney
"""

#Review of Expected formats
#exoected format YYYYMMDD
#SURDAT
#SORDAT
#SUREND
#AGENCY US # final product will come from NOAA

#SORIND# 'US,US,graph,' + sorind

#1 or 0 

#TECSOU


def row2s57(row):
    """
    Convert the expanded column names used in the csv to an S57 name and value.
    """
    s57row = {}
    # remap the keys
    for key in row:
        if key in _field_map:
            s57row[_field_map[key]] = row[key]
            if row[key] == 'TRUE' or row[key] == 'True':
                s57row[_field_map[key]] = 0
            if row[key] == 'FALSE' or row[key] == 'False':
                s57row[_field_map[key]] = 1
    # enforce additional required formating
    if 'VERDAT' in s57row:
        s57row['VERDAT'] = _vert_datum[s57row['VERDAT']]
    if 'HORDAT' in s57row:
        h = s57row['HORDAT']
        for name in _horz_datum:
            if name in h:
                s57row['HORDAT'] = _horz_datum[name]
    if 'SUREND' in s57row:
        s57row['SORDAT'] = s57row['SUREND']
    return s57row

_vert_datum = {
        'MLWS'  : '1',
        'MLLWS' : '2',
        'MSL'   : '3',
        'LLW'   : '4',
        'MLW'   : '5',
        'ISLW'  : '8',
        'MLLW'  : '12',
        'MHW'   : '16',
        'MHWS'  : '17',
        'MHHW'  : '21',
        'LAT'   : '23',
        'LOC'   : '24',
        'IGLD'  : '25',
        'LLWLT' : '27',
        'HHWLT' : '28',
        'HAT'   : '30',
        'Unknown': '701',
        'Other' : '703',
        'HRD'   : '24',  # Adding this for the Hudson River Datum
        }
        
_horz_datum = {
        'WGS72' : '1',
        'WGS84' : '2',
        'WGS_1984' : '2',
        'NAD27' : '74',
        'NAD83' : '75',
        'North_American_Datum_1983' : '75',
        'Local' : '131',
        }

Units = [ 'U.S. Survey Feet.', 'US Survey Foot'
]
"""
Notes: CEMVN conventions
#from_horiz_units
 U.S. Survey Feet.
script: from_vert_units
US Survey Foot
Horizontal_Units
 U.S. Survey Feet.


"""

 
    """    
    #S-57 Definitions long 
    Code:Meaning
    """
    vert_datumS57_longname = {        
            '1':'Mean low water springs',
            '2':'Mean lower low water springs',
            '3':'Mean sea level',
            '4':'Lowest low water',
            '5':'Mean low water',
            '6':'Lowest low water springs',
            '7':'Approximate mean low water springs',
            '8':'Indian spring low water',
            '9':'Low water springs',
            '10':'Approximate lowest astronomical tide',
            '11':'Nearly lowest low water',
            '12':'Mean lower low water',
            '13':'Low water',
            '14':'Approximate mean low water',
            '15':'Approximate mean lower low water',
            '16':'Mean high water',
            '17':'Mean high water springs',
            '18':'High water',
            '19':'Approximate mean sea level',
            '20':'High water springs',
            '21':'Mean higher high water',
            '22':'Equinoctial spring low water',
            '23':'Lowest astronomical tide',
            '24':'Local datum',
            '25':'International Great Lakes Datum 1985',
            '26':'Mean water level',
            '27':'Lower low water large tide',
            '28':'Higher high water large tide',
            '29':'Nearly highest high water',
            '30':'Highest astronomical tide (HAT)',
            }
    """		
    BPS vertical datum codes
    code     Description
    """
    BPS_vert_datum = {     
            '0':'Undetermined vertical datum',
            '1':'Mean Sea Level',
            '2':'Mean Low Water',
            '3':'Mean Low Water Springs',
            '4':'Mean Lower Low Water',
            '5':'Mean Lower Water Springs',
            '6':'Lowest Normal Low Water',
            '7':'Lowest Low Water',
            '8':'Indian Spring Low Water',
            '9':'Great Lakes Low Water',
            '10':'Low Lake Level - Lake Champlain',
            '11':'Normal Pool Level - Cayuga and Seneca Lakes',
            '12':'Hudson River Datum',
            '13':'Normal Lake Level - Franklin D. Roosevelt Lake',
            '14':'Sacramento River - Sacramento to Old Ferry',
            '15':'Columbia River Datum',
            '16':'Local Low Water',
            '17':'Gulf Coast Low Water',
            '18':'Low Water Datum 600.0 ft IGLD-1955 Lake Superior',
            '19':'Low Water Datum 576.8 ft IGLD-1955 L Michigan,Huron',
            '20':'Low Water Datum 571.7 ft IGLD-1955 Lake St. Clair',
            '21':'Low Water Datum 568.6 ft IGLD-1955 Lake Erie',
            '22':'Low Water Datum 242.8 ft IGLD-1955 Lake Ontario',
            '23':'Mystic River Datum',
            '24':'Mean High Water',
            '25':'Low Water Datum 601.1 ft IGLD-1985 Lake Superior',
            '26':'Low Water Datum 577.5 ft IGLD-1985 L Michigan,Huron',
            '27':'Low Water Datum 572.3 ft IGLD-1985 Lake St. Clair',
            '28':'Low Water Datum 569.2 ft IGLD-1985 Lake Erie',
            '29':'Low Water Datum 243.3 ft IGLD-1985 Lake Ontario',
            '30':'LW Reference Plane 1993',
            '31':'Gulf Coast LWD',
            '32':'Lake Level above MSL',
            '33':'Lake Washington and Lake Union MLD',
            '40':'Other',
            '34':'Columbia River above MSL 1929',
            '35':'Low Water Datum',
            '36':'Low Water Datum IGLD-1955',
            '37':'Low Water Datum IGLD-1985',
            '38':'Low Water Datum at Ordinary Springtides',
            '39':'Lake Washington Low Water Datum',
            }
    
    BPStoS57_VertDat={
            '1':'3',#Mean Sea Level
            '2':'5',#Mean Low Water
            '3':'1',#	Mean Low Water Springs	 
            '4':'12',#Mean Lower Low Water
            '24':'16',#Mean High Water
            '12':'24',#Hudson River Datum	-> Local Datum,#USACE
            '15':'24',#Columbia River Datum-> Local Datum,#USACE
            '30':'24',#Low Water Reference Plane 1993-> Local Datum (Mississippi River),#USACE
            '14':'24',#'Sacramento River - Sacramento to Old Ferry',#USACE
            '33':'24',#'Lake Washington and Lake Union MLD',#USACE,#(Seattle Puget Sound above locks)
            '39':'24',#'Lake Washington Low Water Datum',#USACE
            '10':'24',#'Low Lake Level - Lake Champlain',#NY/VT Border
            '32':'24',#'Lake Level above MSL',#Places like Lake Tahoe are reported in lowest lake level above MSL,#California
            '11':'24',#'Normal Pool Level - Cayuga and Seneca Lakes',#NY (See chart 14786, Eerie Canal, and connections to Lake Champlain on Hudson River, Lake Cayuga & Seneca)
            '13':'24',#'Normal Lake Level - Franklin D. Roosevelt Lake',#Franklin D. Roosevelt Lake is part of the Columbia River in Washington State, that is above a large dam(Coulee Dam). Normal Lake level (from chart) is 1288.6ft above mean sea level
            '23':'24',#'Mystic River Datum',#Mystic & Malden River Datums (Chart says 6.2ft above MLLW), off of Boston Harbor, Massachusetts
            '37':'25',#Low Water Datum IGLD-1985-> IGLD 1985
            }
    BPS_GreatLakeDatums={
            '18':'25',#Low Water Datum 600.0 ft IGLD-1955 Lake Superior		 --> IGLD 1985
            '19':'25',
            '20':'25',
            '21':'25',
            '22':'25',
            '25':'25',
            '26':'25',
            '27':'25',
            '28':'25',
            '29':'25',
            '36':'25',
            '37':'25',#Low Water Datum IGLD-1985-> IGLD 1985
            }
    
            #    CALCULATED_VERTICAL_DATUM_CODE			
            #Code     Description     'VERDAT Code'     VERDAT     Description
            #0	Undetermined vertical datum		
            #1	Mean Sea Level     3	     Mean Sea Level
            #2	Mean Low Water     5     Mean Low Water
            #3	Mean Low Water Springs	     1     mean Low water Springs
            #4	Mean Lower Low Water     12     	Mean lower low water
            #5	Mean Lower Water Springs		
            #6	Lowest Normal Low Water		
            #7	Lowest Low Water		
            #8	Indian Spring Low Water		
            #9	Great Lakes Low Water		
            #10	Low Lake Level - Lake Champlain		
            #11	Normal Pool Level - Cayuga and Seneca Lakes		
            #12	Hudson River Datum		
            #13	Normal Lake Level - Franklin D. Roosevelt Lake		
            #14	Sacramento River - Sacramento to Old Ferry		
            #15	Columbia River Datum		
            #16	Local Low Water		
            #17	Gulf Coast Low Water		
            #18	Low Water Datum 600.0 ft IGLD-1955 Lake Superior		
            #19	Low Water Datum 576.8 ft IGLD-1955 L Michigan,Huron		
            #20	Low Water Datum 571.7 ft IGLD-1955 Lake St. Clair		
            #21	Low Water Datum 568.6 ft IGLD-1955 Lake Erie		
            #22	Low Water Datum 242.8 ft IGLD-1955 Lake Ontario		
            #23	Mystic River Datum		
            #24	Mean High Water     16     Mean High water
            #25	Low Water Datum 601.1 ft IGLD-1985 Lake Superior		
            #26	Low Water Datum 577.5 ft IGLD-1985 L Michigan,Huron		
            #27	Low Water Datum 572.3 ft IGLD-1985 Lake St. Clair		
            #28	Low Water Datum 569.2 ft IGLD-1985 Lake Erie		
            #29	Low Water Datum 243.3 ft IGLD-1985 Lake Ontario		
            #30	LW Reference Plane 1993		
            #31	Gulf Coast LWD		
            #32	Lake Level above MSL		
            #33	Lake Washington and Lake Union MLD		
            #40	Other		
            #34	Columbia River above MSL 1929		
            #35	Low Water Datum		
            #36	Low Water Datum IGLD-1955		
            #37	Low Water Datum IGLD-1985		
            #38	Low Water Datum at Ordinary Springtides		
            #39	Lake Washington Low Water Datum		