"""

01112019 BDB_meta_match_plus_Table.py 
BDB_Loadmeta_matchfromcsv_BPS_plusTable_A1.py
breaking up into subscripts (BPS_Loadmeta_BPS.py) to make it easier to read
and differentiate if different file formats needed to be added, and BPS vs BAG metadata and calculations

BDB_Loadmeta_matchfromcsv_BPS_plusTable_A1.py
03072018
20180220

v0.0.1 20180129
update 01102019 for caris 5.1

This script is designed to find a list of csv files from the Bathy Point Store and pull in their associated metadata from
a CSV list and load it into Caris BathyDatabase.

Will need to use csv module in standard python install to do a better job. For now just trying to make it work.
"text,text" as value exists in the export and we will need to deal with that for some column string values.
Putting in plus 1 to move over a column this non-column break comma

Juliet Kinney first round
20180220

I you pass your BDB metatadata BACK OUT OF THE DATABASE it translates it to the description. 
It stays in the 'value' form within the script in the version that you need to load it before you set and retrieve the values.

For example:
metadata['TECSOU'] = '1' (in the main script)
then when passed to the other script its read as
metadata['TECSOU'] = "found by echo-sounder" 

Adding pull of SUREND metadata from spreadsheet
"""

# called in script that calls these.
#import caris.bathy.db as bdb
#import os
#from glob import glob
#import xml.etree.ElementTree as et  
#from datetime import datetime
#from caris.bathy.db import *
#import caris

###import caris.coverage#NEW!! with 5.1 probably not needed till later for coverage actions
##import parse_ncei_meta_HSMDB_c2 as pncei #parse_ncei_meta_gr.py# called in script that calls these.# importing python file in the same directory of that name, trying to give function an alias that is shorter to call and bring in caris batch function


#moving to other script
#import CalculatePOSACC as C_P
#import DataQual1 as dqual#DataQual1.py
#import BringInSuppression as Sup
import S101_fvals_ as S_f
#import BPS_qual1 as B_Q1

#import load_BPS_metadataFeb2_Ed as lbm # importing python file in the same directory of that name, trying to give function an alias that is shorter to call and bring in caris batch function
#------------------------------------------------------------------------------
"""
REQUIRES that the BathyPointStore Table format recieved from Jason/Coung is defined in the calling script.
Examples:
BPStableversion=2#DEFINE BPS TABLE VERSION RECIEVED FROM NCEI. 1:NY format, 2: Format as recieved for Mississippi 
load_meta_BPS(database, csv_file, metadata, Man_Table, directory_path, File_NameBPS)

"""
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def load_meta_BPS(db, csv_file, metadata, Man_Table, directory_path, File_NameBPS):
    """
    Upload the BPS at file_path into the database with the metadata captured 
    in metadata_path.  If metadata_path is empty then the metadata is assumed
    to be in the same location as the BPS file and with the same name.
    
    The required metadata fields include
        resolution
        end survey date
        catzoc = From NCEI XML?
        source authority = ???, only applies to army corps?
        dton = False
    """    
    try:

        metadata=load_csv_metadata(db, csv_file, metadata, directory_path, File_NameBPS)#
    except:
        print('error csv metadata')
    return  metadata #CATZOC #, OBJNAM, Information

#
#------------------------------------------------------------------------------
def load_ncei_metadata_xml(sid, bpath, metadata, nceipath):
    """
    Load the metadata from the NCEI xml into BDB using the load_ncei_metadata
    module, which was originally based on "parse_bag_meta".  Only the SORDAT
    information used.
    """
    # read the xml from the diretory with a similar name
    try:
        nceimeta = sid + '_hsmdb.xml'
        #nceimeta = sid + '.xml'#
        #check format
        #Try this if you are loading a different way
        ###head, tail = os.path.split(bpath)
        # assuming the bag name and the NCEI xml start with the survey number.
        ###xml_name = tail[:6] + '.xml'
        ###xml_path = os.path.join(head,xml_name)
        xml_path = os.path.join(nceipath, nceimeta)#xml_path = os.path.join(bpath, nceimeta)
        xml_dict = pncei.extract_s57_dict(xml_path)#need to call the ncei script not the metascript#pbm.extract_s57_dict(xml_path)

    except:
        print ("error getting NCEI metadata for testing other naming convention " + bpath)
        try:
            nceimeta = sid + '.xml'#
            xml_path = os.path.join(nceipath, nceimeta)#xml_path = os.path.join(bpath, nceimeta)
            xml_dict = pncei.extract_s57_dict(xml_path)#need to call the ncei script not the metascript#pbm.extract_s57_dict(xml_path)
        except:
            print('neither H#####.xml or H#####_hsmdb.xml worked')

    # write them into the provided metadata object.
    try:
        if 'SORDAT' in xml_dict:
            metadata['SORDAT'] = xml_dict['SORDAT']
        if 'SURSTA' in xml_dict:
            metadata['SURSTA'] = xml_dict['SURSTA']
        if 'SUREND' in xml_dict:
            metadata['SUREND'] = xml_dict['SUREND']
    except:
        print ("Error writing to NCEI xml metadata object")
    return metadata
#
#------------------------------------------------------------------------------

def load_csv_metadata(db, csv_file, metadata, file_path, File_NameBPS):
    """
    This method was developed during training to demonstrate loading metadata in
    a specific format.  It is not intended to be used for any specific future
    purpose, but is maintained here for archival purposes.
    """
    #print (csv_file)
    #print (file_path)
    if BPStableversion:
        if BPStableversion == 2:
            try:
                SurveyList = []
                bname= File_NameBPS.split('_')#use '.' if its just the .csv from the BPS store
                Sname = bname[0]# first column that should equal survey name
                SurveyList.append(Sname)
                print('Survey name? ' + Sname)
                found = False # adding in some debugging flags with this #False
                with open(csv_file, 'r') as f:#'r' read only mode
                    for line in f: #lets us loop through every line in the file
                        tokens = line.split('\t')#tokens= ['value', 'value2', etc]
                        tokens = [tt.replace('"','') for tt in tokens]#remove double quotes]#
                        try:#try catch statement to help us skip over problem files
                            BPS_internalid = tokens[0]
                            if tokens[1] == Sname: 
                                found = True
                                if found ==True:
                                    print('found a match')
                                Information = tokens[1]#what is survey ID test load in here
                                metadata['INFORM'] = Information
                                SURVEY_ID = tokens[1] 
                                START_YEAR = tokens [14]#11+1]#moves it over a column
                                END_YEAR = tokens[15]#[12+1]
                                SCALE = tokens[16]#[13+1]
                                SOURCE_INSTITUTION = tokens[8]#5]
                                SURVEY_TYPE_CODE = tokens[3]#2]
                                CALCULATED_VERTICAL_DATUM_CODE = tokens[24]#18+1]
                                CALCULATED_HORIZ_DATUM_CODE = tokens[21]#15+1]
                                SOUNDING_METHOD = tokens[28]#22+1]
                                metadata['CSCALE'] = SCALE #match bdb metadata field with coresponding defined token column in BPS header
                                metadata['SORIND'] = 'US,US,Graph,' + SURVEY_ID
                                metadata['autrty'] = '1' #'False'#False =1, undefined, or True =0
                                metadata['dgrton'] = '1' #'False'#False =1, undefined, or True =0
                                metadata['DUNITS'] = '1' #1 is for meters
                                if END_YEAR == "":
                                    print('no end date found')
                                else:
                                    metadata['SORDAT'] = END_YEAR + '1231'#Dateformat needs day and month not just year, so this part adds the last day of the year to end date
                                    metadata['SUREND'] = END_YEAR + '1231'
                                if START_YEAR == "":
                                    print('no start date found')
                                else: 
                                    metadata['SURSTA'] = START_YEAR + '0101'#Dateformat needs day and month not just year, so this part adds the first day of the year to start date
                                metadata['OBJNAM'] = SURVEY_ID + '_BathyPointStore'
                                metadata['AGENCY'] = 'US'
                                metadata['SURATH'] = SOURCE_INSTITUTION
                                try:
                                    metadata['VERDAT'],  vertdat_other = outputverdat(CALCULATED_VERTICAL_DATUM_CODE)
                                    #vertdat, vertdat_other = outputverdat(CALCULATED_VERTICAL_DATUM_CODE)
                                    if vertdat_other is not None:
                                        print('local datum is listed in db as vd_oth')
                                        metadata['vd_oth'] = vertdat_other
                                except:
                                    print('may be in a local datum')
                                #metadata['VERDAT'] = '12' #12 is MLLW #could incorporate look up table between BPS CALCULATED_VERTICAL_DATUM_CODE and S-57, however all values were MLLW
                                metadata['s_ftyp'] = 'xyz'
                                metadata['s_orig'] = 'BPS'
                                if SOUNDING_METHOD == 'Lead Line assumed':
                                    metadata['TECSOU'] = '5'
                                else:
                                    if SOUNDING_METHOD == 'Digital Echo Sounder w/ Graphical Record assumed':
                                        metadata['TECSOU'] = '1'
                                    else:
                                        print('TECSOU not found in BPS')                            
                                if CALCULATED_HORIZ_DATUM_CODE == 'NOS31':
                                    metadata['HORDAT'] = '75' #need to incorporate look up table between BPS CALCULATED_HORIZ_DATUM_CODE and S-57
                                else:
                                    metadata['HORDAT'] = '75'
                                    with open(File_NameBPS,'a') as error:
                                        error.write(File_NameBPS + ' HORDAT issue \n')
                                    try:
                                        a= metadata['r_note'] 
                                        metadata['r_note'] = a + ',' + 'hdatum code' #CALCULATED_HORIZ_DATUM_CODE + ' hdatum in BPS'
                                        #CALCULATED_HORIZ_DATUM_CODE + ' hdatum in BPS'
                                    except:
                                        metadata['r_note'] = 'hdatum code' 
                                        print('hdatum')
                                print ('metadata for ' + Sname + ' found')                
                        except:
                            None
                            print ('exception')
                if  not found:#if !found:
                    print ('no metadata found from text  for {}', format(os.path.basename(File_NameBPS)))
            except:
                print ("error getting metadata for " + file_path)

        elif BPStableversion == 1:         
            try:
                #csv_file = r'C:\Training\NOAA\WGOM_BIS_Metadata.csv'
                #csv_file = r'N:\Data\arcgis\Export_HCELL_Final_ed.csv'
                ####BPSlist1=glob.glob("N:/Data/ncei/BPSs/*.csv")
                SurveyList = []
                    #BPSlist1=bb
                ####for nn,filefromBPSlist in enumerate(BPSlist1):
                    ####bb2=os.path.basename(BPSlist1[nn])               
                    ####bname = bb2.split('_')# split the BPS name by underscore i.e. parse into a comma delimited list
                #bb2=os.path.basename(file_path)
                bname= File_NameBPS.split('_')#use '.' if its just the .csv from the BPS store
                Sname = bname[0]# first column that should equal survey name
                SurveyList.append(Sname)
                print('Survey name? ' + Sname)
                    ####Sname = bname[0]# first column that should equal survey name
                    ####SurveyList.append(Sname)
                    ####print('Survey name? ' + Sname)
                    
                found = False # adding in some debugging flags with this #False
                with open(csv_file, 'r') as f:#'r' read only mode
                    for line in f: #lets us loop through every line in the file
                        tokens = line.split(',')#tokens= ['value', 'value2', etc]
                        #will let us split the file based on commas ',' , using tokens
                        #token = tokens.strip('"')
                        tokens = [tt.replace('"','') for tt in tokens]#
                        try:#try catch statement to help us skip over problem files
                            #python stackoverflow basename
                            ####BPSfiles = glob(os.path.join(file_path,'*.BPS'))
                            ####if tokens[3] == SurveyList:
                            #print(tokens[0])
                            if tokens[0] == Sname: 
                                found = True
                                #print(tokens[2])
                                #index 19 in file #ExtVessel field in Paul's data
                                Information = tokens[0]#what is survey ID test load in here
                                #CATZOC = tokens[3]
                                #OBJNAM = os.path.basename(file_path)
                                ##metadata = surface.get_metadata()
                                #metadata['CATZOC'] = CATZOC #tokens[3]# should be CATZOC if available
                                ##metadata['OBJNAM'] = os.path.basename(file_path)
                                metadata['INFORM'] = Information
                                #metadata['SORIND'] = 'US,US,Graph,' + SourceID
                                #surface.set_metadata(metadata)
                                SURVEY_ID = tokens[0] 
                                START_YEAR = tokens [11+1]#moves it over a column
                                END_YEAR = tokens[12+1]
                                SCALE = tokens[13+1]
                                SOURCE_INSTITUTION = tokens[5]
                                SURVEY_TYPE_CODE = tokens[2]
                                CALCULATED_VERTICAL_DATUM_CODE = tokens[18+1]
                                CALCULATED_HORIZ_DATUM_CODE = tokens[15+1]
                                SOUNDING_METHOD = tokens[22+1]
                                #metadata['OBJNAM'] = os.path.basename(file_path)
              
                                metadata['CSCALE'] = SCALE #match bdb metadata field with coresponding defined token column in BPS header
                                metadata['SORIND'] = 'US,US,Graph,' + SURVEY_ID
                                metadata['autrty'] = '1' #'False'#False =1, undefined, or True =0
                                metadata['dgrton'] = '1' #'False'#False =1, undefined, or True =0
                                metadata['DUNITS'] = '1' #1 is for meters
                                if END_YEAR == "":
                                    print('no end date found')
                                else:
                                    metadata['SORDAT'] = END_YEAR + '1231'#Dateformat needs day and month not just year, so this part adds the last day of the year to end date
                                    metadata['SUREND'] = END_YEAR + '1231'
                                if START_YEAR == "":
                                    print('no start date found')
                                else: 
                                    metadata['SURSTA'] = START_YEAR + '0101'#Dateformat needs day and month not just year, so this part adds the first day of the year to start date
                                metadata['OBJNAM'] = SURVEY_ID + '_BathyPointStore'
                                metadata['AGENCY'] = 'US'
                                metadata['SURATH'] = SOURCE_INSTITUTION
                                #metadata['VERDAT'] = '12' #12 is MLLW #could incorporate look up table between BPS CALCULATED_VERTICAL_DATUM_CODE and S-57, however all values were MLLW
                                try:
                                    metadata['VERDAT']=outputverdat(CALCULATED_HORIZ_DATUM_CODE)
                                except:
                                    print('may be in a local datum')
                                metadata['s_ftyp'] = 'xyz'
                                metadata['s_orig'] = 'BPS'
                                if SOUNDING_METHOD == 'Lead Line assumed':
                                    metadata['TECSOU'] = '5'
                                else:
                                    if SOUNDING_METHOD == 'Digital Echo Sounder w/ Graphical Record assumed':
                                        metadata['TECSOU'] = '1'
                                    else:
                                        print('TECSOU not found in BPS')                            
        ##                        metadata['TECSOU'] = 
        ##
        ##                        TECSOU: SOUNDING_METHOD
        ##
        ##                        5: Lead Line assumed
        ##                        1: Digital Echo Sounder w/ Graphical Record assumed
        
        ##                        PARAMETERS_SURVEYED_CODE
        ##                        SF
        ##                        S
        ##
        ##                        Quality_Code
        ##                        0 SUSPECT                                                      DEPTH
        ##                        1 ORIGINAL                                                     DEPTH
        ##                        2 DUPLICATED - CO-LOCATED SOUNDINGS                            DEPTH
        ##                        3 DUPLICATED OF SHALLOW  -  CO-LOCATED SOUNDINGS               DEPTH
        ##                        4 SELECTED DEPTH - CO-LOCATED SOUNDINGS                        DEPTH
        ##                        5 NON-SELECTED DEPTH - CO-LOCATED SOUNDINGS                    DEPTH
        ##                        6 UN-DETERMINED DEPTH - CO-LOCATED SOUNDINGS                   DEPTH
        ##                        7 CORRECTION OF ORIGINAL DEPTH - OUTLIERS                      POSITION
        ##                        8 UN-DETERMINED DEPTH - OUTLIERS                               POSITION
        ##                        9 CORRECTION OF ORIGINAL DEPTH - NOS ERRORS                    DEPTH
        ##                        10 UN-DETERMINED DEPTH - NO SMOOTHSHEET - NOS ERRORS            DEPTH
        ##                        11 SUSPECTED SURFACE OBJECT                                     DEPTH
        ##                        12 POSITIONAL ERROR - OUTLIERS FROM CONVEX HULL                 POSITION
        ##                        14 SURVEYS WITH LOW SOUNDING COUNT                              POSITION
        ##                        15 SURVEYS WITH SOUNDINGS IN WRONG LOCATION                     POSITION
        ##                        16 POSITIONAL ERROR - SOUNDINGS CONFLICT W/ GENERATED POLYGON   POSITION
        ##                        17 SUPERSEDED SOUNDINGS                                         POSITION
        ##                        18 NON-CHARTED SMOOTHSHEET SOUNDING                             POSITION
        ##                        19 SUPERSEDED SOUNDINGS - HCELL PROCESS                         POSITION
        ##
        
        
        ##                        SURVEY_TYPE_CODE
        ##                                       0 Unknown                                                     
        ##                        1 Basic Hydrographic Survey                                   
        ##                        2 D survey                                                    
        ##                        3 EEZ Bathymetric Survey                                      
        ##                        4 HYDROGRAPHIC SURVEY                                         
        ##                        5 Hydro Field Examination                                     
        ##                        6 Navigable Area Hydro Survey                                 
        ##                        7 Reconnaissance survey                                       
        ##                        8 Side Scan Sonar Survey                                      
        ##                        9 Wire-drag field examination                                 
        ##                        10 Wire-drag survey                                            
        ##                        11 CANADIAN HYDROGRAPHIC SERVICE  
        ##
        ##                        S-57
        ##                        UNDEFINED: 'Value is undefined'
        ##                        UNKNOWN: 'Value is unknown'
        ##                        1: 'reconnasiance sketch survey'
        ##                        2: 'controlled survey'
        ##                        3: '{unsurveyed}'
        ##                        4: 'examination survey'
        ##                        5: 'passage survey'
        ##                        6: 'remotely sensed'
                                
                                ##      metadata['SURTYP'] = #need to incorporate look up table between BPS SURVEY_TYPE_CODE and S-57
                                if CALCULATED_HORIZ_DATUM_CODE == 'NOS31':
                                    metadata['HORDAT'] = '75' #need to incorporate look up table between BPS CALCULATED_HORIZ_DATUM_CODE and S-57
                                else:
                                    metadata['HORDAT'] = '75'
                                    with open(File_NameBPS,'a') as error:
                                        error.write(File_NameBPS + ' HORDAT issue \n')
                                    try:
                                        a= metadata['r_note'] 
                                        metadata['r_note'] = a + ',' + 'hdatum code' #CALCULATED_HORIZ_DATUM_CODE + ' hdatum in BPS'
                                        #CALCULATED_HORIZ_DATUM_CODE + ' hdatum in BPS'
                                    except:
                                        metadata['r_note'] = 'hdatum code' 
                                        print('hdatum')
        
              
                                ##      metadata['DRVAL1'] = #will need to be sorted and extracted from all sounding values - Sounding Table Column 4 (CURRENT_DEPTH) - or taken from gridded surface
                                ##      metadata['DRVAL2'] = #will need to be sorted and extracted from all sounding values - Sounding Table Column 4 (CURRENT_DEPTH) - or taken from gridded surface
                                ##      metadata['srfres'] = #possible grid resolution deermined on scale value, or this is autopopulated from surface once grided
                                ##      
                                ##      metadata['TECSOU'] = #Survey type code include side scan and wire drag. sounding method includes lead line or digital echosounder.  string to number list conversion from SOUNDING_METHOD
                                ##      
                                ##      metadata['CATZOC'] = #possible calculated value based on other fields, starting with END_YEAR and SOUNDING_METHOD
                                ##      metadata['SOUACC'] = #there is a quality code value in the sounding table.  possibly contribute to this assesment?
                                ##      metadata['POSACC'] = #field POSITION_DETERMINATION lists horizontal positioning technique and could be used to evaluate 
          
                                print ('metadata for ' + Sname + ' found')
        
        
                        except:
                            None
                            print ('exception')
                if  not found:#if !found:
                    print ('no metadata found from text  for {}', format(os.path.basename(File_NameBPS)))
            except:
                print ("error getting metadata for " + file_path)
    return metadata#CATZOC #, OBJNAM, Information
#------------------------------------------------------------------------------
def outputverdat(CALCULATED_VERTICAL_DATUM_CODE):
    vert_datum = {
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
            }
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
    print('Vertical Datum as listed in BPS as: ' + BPS_vert_datum[CALCULATED_VERTICAL_DATUM_CODE])
    try:
        vertdat = BPStoS57_VertDat[CALCULATED_VERTICAL_DATUM_CODE]
        print('Is converted to S57 as: ' + vert_datumS57_longname[vertdat] + ', code: ' + vertdat)
        if vertdat == '24':
            vertdat_other =BPS_vert_datum[CALCULATED_VERTICAL_DATUM_CODE]
        else:
            vertdat_other = None
    except:
        print('try one of the other sets of conversion dictionaries')
        #verdat = BPS_GreatLakeDatums[CALCULATED_VERTICAL_DATUM_CODE]#some great lake datum
    return vertdat, vertdat_other
#------------------------------------------------------------------------------

def plusTable(metadata, Man_Table, File_NameBPS):
    """
    plusTable is loading in data from NCEI metadata that was exported from XML and hand manipulated into a readable format.
    at this point the example table contains TECSOU information
    """
    print('trying to add in TABLE metadata')
    try:
        SurveyList = []
        T = []
        bname= File_NameBPS.split('_')#use '.' if its just the .csv from the BPS store
        Sname = bname[0]# first column that should equal survey name
        SurveyList.append(Sname)
        print('Survey name? ' + Sname)            
        found = False # adding in some debugging flags with this #False
        with open(Man_Table, 'r') as f:#'r' read only mode
            for line in f: #lets us loop through every line in the file
                tokens = line.split(',')#tokens= ['value', 'value2', etc]
                #will let us split the file based on commas ',' , using tokens
                #token = tokens.strip('"')
                tokens = [tt.replace('"','') for tt in tokens]# if necessary
                try:#try catch statement to help us skip over problem files
                    if tokens[0] == Sname: 
                        found = True
                        Information = tokens[0]#what is survey ID test load in here
                        #metadata['INFORM'] = Information
                        SURVEY_ID = tokens[0] 
                        END_YEAR = tokens[1]
                        FILE_TYPE = tokens[2]
                        Leadline = tokens[3]
                        MBES = tokens[4]
                        SSS = tokens[5]
                        SB = tokens[6]
                        DR_SUREND = tokens[15]
                        print(DR_SUREND)
                        #Lidar = tokens[8]
                        try:
                            d=datetime.strptime(END_YEAR, '%m/%d/%Y')
                            #print(d)
                            numdat1=d.strftime('%Y%m%d')
                            metadata['SUREND'] = numdat1
                            metadata['SORDAT'] = numdat1


                        except:
                            print('date end in spreadsheet problem')
                            with open(errorfile1,'a') as error:
                                error.write(File_NameBPS + ' date end in spreadsheet problem \n')
                        try:
                            if DR_SUREND !="":
                                try:
                                    d=datetime.strptime(DR_SUREND, '%m/%d/%Y')
                                    numdat1=d.strftime('%Y%m%d')
                                    metadata['SUREND'] = numdat1
                                    metadata['SORDAT'] = numdat1
                                except:
                                    print('DR Survey End issue')
                                    with open(errorfile1,'a') as error:
                                        error.write(File_NameBPS + ' DR Survey End issue \n')
                        except:
                            print('DR time issue')
                        try:
                            if MBES == 'MBES':
                                T.append('3')
                            #else:
                            #    print ('not multibeam')
                            if SSS == 'SSS':
                                T.append('2')
                            #else:
                            #    print('not sidescan')
                            if SB == 'SB':
                                T.append('1')
                            if Leadline == 'Leadline':
                                T.append('5')
                        except:
                            print('not found')
                        l = len(T)
                        if l == 2:
                            TECSOU = T[0] + ',' + T[1]
                        elif l == 3:
                            TECSOU = T[0] + ',' + T[1] + ',' + T[2]
                        elif l == 1:
                            TECSOU = T[0]
                        else:
                            if l>-1:
                                print('complicated TECSOU')
                            else:
                                print('no TECSOU found')
                        
                        #else:
                        #    print('not single beam')
                        #if Lidar == 'LI':
                        #    TECSOU.append = 7
                        #else:
                        #    print('not Lidar')
                        print(TECSOU)
                        metadata['TECSOU'] = TECSOU
                        #metadata['STATUS'] = '3'
                        #print(metadata['TECSOU'])
                    
#                            T = 1
#                            T = 2
#                            T = 3
#                            T = 4
#                            T = 5
#                            T = 6
#                            T = 7
#                            T = 8
#                            T = 9
#                            T = 10
#                            T = 11
#                            T = 12
#                            T = 13
#                            T = 14
#                    1	found by echo-sounder
#                    2	found by side scan sonar
#                    3	found by multi-beam
#                    4	found by diver
#                    5	found by lead-line
#                    6	swept by wire-drag
#                    7	found by laser
#                    8	swept by vertical acoustic system
#                    9	found by electromagnetic sensor
#                    10	photogrammetry
#                    11	satellite imagery
#                    12	found by levelling
#                    13	swept by side-scan sonar
#                    14	computer generated
#                    TECSOU S-57 definitions
#                    1 found by echo-sounder: the depth was determined by using an instrument that determines depth of water by measuring the time interval between emission of a sonic or ultrasonic signal and return of its echo from the bottom. (adapted from IHO Dictionary, S-32, 1547)
#                    2 found by side-scan-sonar: the depth was computed from a record produced by active sonar in which fixed acoustic beams are directed into the water perpendicularly to the direction of travel to scan the bottom and generate a record of the bottom configuration. (adapted from IHO Dictionary, S-32, 4710)
#                    3 found by multi-beam: the depth was determined by using a wide swath echo sounder that uses multiple beams to measure depths directly below and transverse to the ship's track. (adapted from IHO Dictionary, S-32, 3339)
#                    4 found by diver: the depth was determined by a person skilled in the practice of diving. (adapted from IHO Dictionary, S-32, 1422)
#                    5 found by lead-line: the depth was determined by using a line, graduated with attached marks and fastened to a sounding lead. (adapted from IHO Dictionary, S-32, 2698)
#                    6 swept by wire-drag: the given area was determined to be free from navigational dangers to a certain depth by towing a buoyed wire at the desired depth by two launches, or a least depth was identified using the same technique. (adapted from IHO Dictionary, S-32, 5248, 6013)
#                    7 found by laser: the depth was determined by using an instrument that measures distance by emitting timed pulses of laser light and measuring the time between emission and reception of the reflected pulses. (adapted from IHO Dictionary, S-32, 2763)
#                    8 swept by vertical acoustic system: the given area has been swept using a system comprised of multiple echo sounder transducers attached to booms deployed from the survey vessel.
#                    9 found by electromagnetic sensor: the depth was determined by using an instrument that compares electromagnetic signals. (adapted from IHO Dictionary, S-32, 1571)
#                    10 photogrammetry: the depth was determined by applying mathematical techniques to photographs. (adapted from IHO Dictionary, S-32, 3791)
#                    11 satellite imagery: the depth was determined by using instruments placed aboard an artificial satellite. (adapted from IHO Dictionary, S-32, 4509)
#                    12 found by levelling: the depth was determined by using levelling techniques to find the elevation of the point relative to a datum. (adapted from IHO Dictionary, S-32, 2741)
#                    13 swept by side-scan-sonar: the given area was determined to be free from navigational dangers to a certain depth by towing a side-scan-sonar. (adapted from IHO Dictionary, S-32, 5248, 4710) [415.2]
#                    14 computer generated: the sounding was determined from a bottom model constructed using a computer.

                except:
                    None
                    print ('exception with table')
        if  not found:#if !found:
            print ('no metadata found in table for {}', format(os.path.basename(File_NameBPS)))      
    except:
        print ("error getting metadata for " + File_NameBPS)
    #print (metadata)
    try:
        metadata, vertACC, vertACCPerct = S_f.S101_fvals(metadata)#needs to be after plusTable's TECSOU, but before other columsn as it will overwrite with manual entries.
#        try:
#            print(metadata['flcvrg'])
#        except:
#            None
    except:
        #None
        with open(errorfile1,'a') as error:
            error.write(File_NameBPS + ' S101_fvals did not work \n')
        print('S101_fvals did not work')
        
    try:
        with open(Man_Table, 'r') as f:#'r' read only mode
            for line in f: #lets us loop through every line in the file
                tokens = line.split(',')#tokens= ['value', 'value2', etc]
                #will let us split the file based on commas ',' , using tokens
                #token = tokens.strip('"')
                tokens = [tt.replace('"','') for tt in tokens]# if necessary
                try:#try catch statement to help us skip over problem files
                    if tokens[0] == Sname: 
                        found = True
                        #Information = tokens[0]#what is survey ID test load in here
                        #metadata['INFORM'] = Information
                        #SURVEY_ID = tokens[0] 
                        MCDPriority = tokens[8]#MCD Priority
#                        if tokens[9]!= "": #if value is not a blank entry than save value to metadata
#                            metadata['flcvrg'] = tokens[9]#flcvrg
#                        if tokens[10]!= "":
#                            metadata['f_dtct'] = tokens[10]#f_dtct
#                        if tokens[11]!= "":
#                            metadata['f_lstd'] = tokens[11]#f_lstd
#                        if tokens[12]!= "":
#                            metadata['f_size'] = tokens[12]#f_size
#                        if tokens[13]!= "":
#                            metadata['flbath'] = tokens[13]#flbath
#                        try:
#                            DR_SUREND = tokens[15]
#                            d=datetime.strptime(END_YEAR, '%m/%d/%Y')
#                            #print(d)
#                            numdat1=d.strftime('%Y%m%d')
#                            metadata['SUREND'] = numdat1
#                            metadata['SORDAT'] = numdat1
#                            if DR_SUREND !="":
#                                try:
#                                    d=datetime.strptime(DR_SUREND, '%m/%d/%Y')
#                                    numdat1=d.strftime('%Y%m%d')
#                                    metadata['SUREND'] = numdat1
#                                    metadata['SORDAT'] = numdat1
#                                except:
#                                    print('DR Survey End issue')
#
#                        except:
#                            print('date end in spreadsheet problem')
                        
                        if tokens[9]!= "": #if value is not a blank entry than save value to metadata
                            flcvrg = tokens[9]#flcvrg
                            if flcvrg == "TRUE":  
                                metadata['flcvrg'] = '0'#'True'                                                                       
                            elif flcvrg == "FALSE":
                                metadata['flcvrg'] = '1'#'False' 
                        if tokens[10]!= "":
                            f_dtct = tokens[10]#f_dtct
                            if f_dtct == "TRUE":
                                metadata['f_dtct'] = '0'#'True'
                            elif f_dtct == "FALSE":
                                metadata['f_dtct'] = '1'#'False' 
                        if tokens[11]!= "":
                            f_lstd = tokens[11]#f_lstd
                            if f_lstd == "TRUE":
                                metadata['f_lstd']  = '0'#'True'
                            elif f_lstd == "FALSE":
                                metadata['f_lstd'] = "1"#'False'
                        if tokens[12]!= "":
                            metadata['f_size'] = tokens[12]#f_size
                            print(metadata['f_size'])
                            #This doesn't need any conversion
                        try:
                            if tokens[13]!= "":
                                flbath = tokens[13]#flbath
                                print(flbath)
                                if flbath == "TRUE":
                                    metadata['flbath'] = '0'#'True'
                                elif flbath == "FALSE":
                                    metadata['flbath'] = '1'#'False'
                                else:
                                    print('flbath unsure')
                        except:
                            print('flbath err')
                        Notes1 = tokens[14]#notes
                        """
                        if tokens[9]!= "":
                             metadata['flcvrg']
                        #if value is not a blank entry than save value to metadata
                        and repeating for  #f_dtct, f_lstd, f_size and flbath
                        """
                except:
                    print('manual load of flcvrg, f_dtct, f_lstd, f_size and flbath did not work')
                    with open(errorfile1,'a') as error:
                        error.write(File_NameBPS + ' manual load of flcvrg, f_dtct, f_lstd, f_size and flbath did not work \n')

    except:
        with open(errorfile1,'a') as error:
             error.write(File_NameBPS + ' table err \n')
        #None
    return metadata, vertACC, vertACCPerct# TECSOU #, vertACC, vertACCPerct

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def convert_date1(inputdate):
    dateformat1 = datetime.strptime(inputdate, '%Y%m%d')
    return dateformat1                    
#------------------------------------------------------------------------------


