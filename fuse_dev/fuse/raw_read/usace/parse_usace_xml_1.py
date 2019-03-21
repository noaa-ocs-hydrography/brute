"""
parse_usace_xml_1.py
Edit 1/3/2018 J. Kinney modifying for HSMDB format xml from Jason Baillio (FGDC not ISO).


20180302 G.Rice

V0.0.1 Last Updated 20180312

This is a collection of tools for working with metadata from NCEI. Created
specifically to extract the metadata for building a bathymetric database,
this collection of method attempt to both serve extraction of the meta data in
a general sense, and also for specific S57 needs.

"""
import re
import xml as xml
from xml.etree import ElementTree as et 
import logging as log
from os import path

try:
    import dateutil.parser as parser
except:
    pass


#log = logging.getLogger(__name__)

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
        
horz_datum = {
        'WGS72' : '1',
        'WGS84' : '2',
        'NAD27' : '74',
        'NAD83' : '75',
        'Local' : '131'}
       
def extract_s57_dict(xmlfilename):
    """
    Open the filename provided and extract the s57 dictionary from an object of
    this module's Meta class type.
    """
    with open(xmlfilename, 'r') as xml_file:
        xml_txt = xml_file.read()
    xmlbasename = path.basename(xmlfilename)
    #we want to be able to figure out the data source as E-Hydro at this point or not
    # and pass that on to the class function
    xml_data = XML_Meta(xml_txt, filename = xmlbasename)
    s57_dict = xml_data.get_s57_dict()
    return s57_dict


        
class XML_Meta(object):
    """ 
    Helper class to manage BAG xml metadata. This class takes an xml string and
    parses the string based on a dictionary (named 'source') with keys that
    name the item to extract and the tree path string to find the data.
    Different versions of the source dictionary can be set based on the version
    of the BAG being parsed.  All extracted data is placed in a dictionary
    (named 'data').
    """

    def __init__(self, meta_xml, version = '', filename = ''):
        """
        Provided a BAG xml string for parsing, a tree will be created and the
        name speace parsed from the second line.  Values are then extracted
        based on the source dictionary.  If this dictionary
        """
        self.filename = filename
        self.xml_tree = et.fromstring(meta_xml)
        self.ns = parse_namespace(meta_xml)
        self.xml_txt = meta_xml
        

        if len(version) > 0:
            self.version = float(version)
        else:
            self.version = self._guess_version()
            print(version)
            
        self._set_format()
        self.get_fields()
    
    def _guess_version(self):
        """
        Try to guess the version of the bag for the purpose of parsing the
        fields.  The guess is based off comparing the name spaces.
        """
        version_1 = {
                    'gco': 'http://www.isotc211.org/2005/gco',
                    'gmd': 'http://www.isotc211.org/2005/gmd',
                    'gmi': 'http://www.isotc211.org/2005/gmi',
                    'gml': 'http://www.opengis.net/gml/3.2',
                    'srv': 'http://www.isotc211.org/2005/srv',
                    'xlink': 'http://www.w3.org/1999/xlink',
                    'gmx': 'http://www.isotc211.org/2005/gmx',
                    'gss': 'http://www.isotc211.org/2005/gss',
                    'gts': 'http://www.isotc211.org/2005/gts',
                    'gsr': 'http://www.isotc211.org/2005/gsr',
                    }

        version_HSMDB = {'HSMDB'}
        version_USACE_FGDC = {'USACE_FGDC'}
        if self.ns == version_1:
            return 1.0
        elif self.xml_txt.startswith('<?xml version="1.0" encoding="ISO-8859-1"?>\n'):
            print('ISO-8859-1 xml version')
            #xml_version = 'ISO-8859-1'
            return 'ISO-8859-1'
        elif self.xml_tree.tag == 'metadata':
            print(version_USACE_FGDC)#:
            return 'USACE_FGDC'
            print('FGDC format not ISO, USACE example')
        elif self.xml_tree.tag == 'HSMDB':
            """
            need to add if loop to check that its source is E-Hydro
            self.source needs to be defined above first
            """
            print(version_HSMDB)#:
            return 'HSMDB'
            print('FGDC format not ISO, HSMDB example')
        else:
            return -1.0
            print('We do not have a template for this verion yet!')

    def _set_format(self):
        """
        Set the locations of the desired data types based on the version of the
        bag.
        """
        version = self.version
        
        if version =='USACE_FGDC':
            self.source = {}
            try:
               my_etree_dict = self.convert_xml_to_dict2() 
               if my_etree_dict['metstdv']:
                   Metadataformat = my_etree_dict['metstdv']
                   print(Metadataformat)
                   self.metadataformat = Metadataformat
            except:
                print('unexpected issue with assumed USACE FGDC format parsing')
                #try:
                    #my_etree_dict = self.convert_xml_to_dict()
                    #this version just pulls children and grandchildren generically,
                    #it is not looking for specific structure
                #except:
                    #print('unexpected issue with assumed USACE FGDC format parsing')
        elif version == 'ISO-8859-1':
            self.source = {}
            try:
               my_etree_dict = self.convert_xml_to_dict() 
               if my_etree_dict['metstdv']:
                   Metadataformat = my_etree_dict['metstdv']
                   print(Metadataformat)
                   self.metadataformat = Metadataformat
            except:
                print('unexpected issue with assumed USACE ISO 88591 FGDC format parsing')

                        
        elif version =='HSMDB':
            self.source = {}
            try:
                my_etree_dict= self.convert_xml_to_dict()
                self.source['SORDAT'] = my_etree_dict['SURVEY__DATE_SURVEY_END']
                #self.source['SURATH'] = my_etree_dict['']
                self.source['SUREND'] = my_etree_dict['SURVEY__DATE_SURVEY_END']#'.//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:temporalElement/gmd:EX_TemporalExtent/gmd:extent/gml:TimePeriod/gml:endPosition' #updated
                self.source['SURSTA'] = my_etree_dict['SURVEY__DATE_SURVEY_BEGIN']#'.//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:temporalElement/gmd:EX_TemporalExtent/gmd:extent/gml:TimePeriod/gml:beginPosition' # updated
                #TECSOU information is in the Survey Instruments list
                #self.source['TECSOU'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:instrument/gmi:MI_Instrument/gmi:type/gco:CharacterString' # updated
                #self.source['DATUM2'] =  my_etree_dict['SURVEY__VERTICAL_DATUM_SURVEYED']
    
                self.source['DATUM'] =  my_etree_dict['SURVEY__VERTICAL_DATUM']
                self.source['survey'] = my_etree_dict['SURVEY__NUMBER']#'.//gmd:fileIdentifier/gco:CharacterString' #updated
                self.source['planam'] = my_etree_dict['PROJECT__PROJECT_NUMBER']#'.//gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:platform/gmi:MI_Platform/gmi:identifier/gmd:RS_Identifier/gmd:code/gco:CharacterString' # updated#
                ##self.source['sensor'] = #'.//gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:instrument/gmi:MI_Instrument/gmi:description/gco:CharacterString' #updated
                #my_etree_dict['SURVEY__VERTICAL_UNITS']
                #my_etree_dict['SURVEY__SUBLOCALITY']
                #my_etree_dict['SURVEY__SCALE']                
                my_etree_dict['SURVEY_INSTRUMENT__TYPE']
                #my_etree_dict['SURVEY__LOCALITY']
                #my_etree_dict['SURVEY__HORIZONTAL_DATUM_SURVEYED']
                #my_etree_dict['SURVEY__HORIZONTAL_RESOLUTION']
                #my_etree_dict['SURVEY__HORIZONTAL_UNITS']
                #my_etree_dict['SURVEY__HORIZONTAL_DATUM']
                #my_etree_dict['SURVEY__FIELD_UNIT']
                #my_etree_dict['SURVEY__FIELD_NUMBER']
                
                #my_etree_dict['SURVEY__DATA_MODIFIED_BY']
                
                #my_etree_dict['SURVEY_STATE__STATE_ID']
                #my_etree_dict['SURVEY_STATE__STATE']
                
##                my_etree_dict['SURVEY_PRODUCTS__PUBLISH_DATE']
##                my_etree_dict['SURVEY_PRODUCTS__FORMAT_ID']
##                my_etree_dict['SURVEY_PRODUCTS__DATA_FILE']
##                my_etree_dict['SURVEY_DATA_CONTENT__DATA_CONTENT_ID']
##                my_etree_dict['SURVEY_DATA_CONTENT__DATA_CONTENT']
##                my_etree_dict['SURVEY_CHIEF__CHIEF']
##                my_etree_dict['SURVEY_BOUNDS__TOP_BOUND']
##                my_etree_dict['SURVEY_BOUNDS__RIGHT_BOUND']
##                my_etree_dict['SURVEY_BOUNDS__MIN_DEPTH']
##                my_etree_dict['SURVEY_BOUNDS__MAX_DEPTH']
##                my_etree_dict['SURVEY_BOUNDS__LEFT_BOUND']
##                my_etree_dict['SURVEY_BOUNDS__BOTTOM_BOUND']
                #my_etree_dict['PROJECT__PURPOSE']
                
                #my_etree_dict['SURVEY__PROJECTION_NAME']
                #my_etree_dict['PROJECT__PROJECT_ID']
                #my_etree_dict['SURVEY__UUID']
#                my_etree_dict['SURVEY_STATS__DAYS_OF_ACQUISITION']
#                my_etree_dict['SURVEY_STATS__LINEAR_NM_BOTTOM_DRAG']
#                my_etree_dict['SURVEY_STATS__LINEAR_NM_CHAINDRAG']
#                my_etree_dict['SURVEY_STATS__LINEAR_NM_LIDAR']
#                my_etree_dict['SURVEY_STATS__LINEAR_NM_MULTIBEAM']
#                my_etree_dict['SURVEY_STATS__LINEAR_NM_SIDESCAN']
#                my_etree_dict['SURVEY_STATS__LINEAR_NM_VERTICALBEAM']
#                my_etree_dict['SURVEY_STATS__NUMBER_OF_POSITIONS']
#                my_etree_dict['SURVEY_STATS__PERCENT_LIDAR_COVERAGE']
#                my_etree_dict['SURVEY_STATS__PERCENT_MULTIBEAM_COVERAGE']
#                my_etree_dict['SURVEY_STATS__PERCENT_SIDESCAN_COVERAGE']
#                my_etree_dict['SURVEY_STATS__PERCENT_VERTICALBEAM_COVERAGE']
#                my_etree_dict['SURVEY_STATS__SELECTED_SOUNDINGS']
#                my_etree_dict['SURVEY_STATS__SQUARE_NM_BOTTOM_DRAG']
#                my_etree_dict['SURVEY_STATS__SQUARE_NM_LIDAR']
#                my_etree_dict['SURVEY_STATS__SQUARE_NM_MULTIBEAM']
#                my_etree_dict['SURVEY_STATS__SQUARE_NM_SIDESCAN']
#                my_etree_dict['SURVEY_STATS__SQUARE_NM_VERTICALBEAM']
#                my_etree_dict['SURVEY_STATS__TOTAL_LINEAR_NM']
#                my_etree_dict['SURVEY_STATS__TOTAL_LINEAR_NM_CROSSLINES']
#                my_etree_dict['SURVEY_STATS__TOTAL_SQUARE_NM']
#                my_etree_dict['SURVEY_TIDE_STATION__REFERENCE']
#                my_etree_dict['SURVEY_TIDE_STATION__STATUS']
#                my_etree_dict['SURVEY_TIDE_STATION__TIDE_STATION_ID']
#                my_etree_dict['SURVEY_TIDE_STATION__TIDE_STATION_NAME']
#                my_etree_dict['SURVEY_TIDE_STATION__TYPE']
#                my_etree_dict['SURVEY__CONTRACTOR_ID']
#                my_etree_dict['SURVEY__CONTRACT_NUMBER']                
#                my_etree_dict['SURVEY__TIDE_EPOCH']
#                my_etree_dict['SURVEY__TIDE_NOTE']
#                my_etree_dict['SURVEY__TIDE_ZONES']
            except:
                print('unexpected issue with HSMDB format parsing')     
        else:           
            version = float(self.version)
            
            if version == 1.0:
                self.source = {}
                self.source['SORDAT'] ='.//gmd:identificationInfo/gmd:MD_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date' #updated 
                self.source['SURATH'] = './/gmd:contact/gmd:CI_ResponsibleParty/gmd:organisationName/gco:CharacterString' # Updated
                self.source['SUREND'] = './/gmd:identificationInfo/gmd:MD_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:temporalElement/gmd:EX_TemporalExtent/gmd:extent/gml:TimePeriod/gml:endPosition' #updated
                self.source['SURSTA'] = './/gmd:identificationInfo/gmd:MD_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:temporalElement/gmd:EX_TemporalExtent/gmd:extent/gml:TimePeriod/gml:beginPosition' # updated
                self.source['TECSOU'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:instrument/gmi:MI_Instrument/gmi:type/gco:CharacterString' # updated
                self.source['DATUM'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:authority/gmd:CI_Citation/gmd:title/gco:CharacterString'
                #self.source['DATUM2'] =  './/' #updated
                self.source['survey'] = './/gmd:fileIdentifier/gco:CharacterString' #updated
                self.source['planam'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:platform/gmi:MI_Platform/gmi:identifier/gmd:RS_Identifier/gmd:code/gco:CharacterString' # updated#
                self.source['sensor'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:instrument/gmi:MI_Instrument/gmi:description/gco:CharacterString' #updated
               
            else:
                log.warning("verison not compatible")

    ##def convert_xml_to_dict(self):
    #def convert_xml_to_dict(my_etree):            
    #    my_etree_dict={}            
    #   for ch in self.xml_tree:
    #    #for ch in my_etree:
    #        #print(ch)#list of children
    #        grandchildren = ch.getchildren()
    #        for ss in grandchildren:
    #            #print(ss)
    #            print(ss.tag, ss.text)
    #            tag1=ss.tag
    #            key=tag1
    #            value=ss.text
    #            my_etree_dict[key]=value
    #    return my_etree_dict
    #previous version
    
    def convert_xml_to_dict(self):
        """
        This version exports out multiple entries into a comma delimited list within the dictionary
        it is just appended within the dictionary
        The method may be modified if needed.
        This approach is just going through the first few layers of the etree
        and pulling them out, it is not for more deeply nested files
        it is meant as a generic pull when more specific paths have yet
        to be developed
        """           
        my_etree_dict1={}        
        for ch in self.xml_tree:
            #print(ch)#list of children
            grandchildren = ch.getchildren()
            for ss in grandchildren:
                #print(ss)
                print(ss.tag, ss.text)
                tag1=ss.tag
                key=tag1
                value1=ss.text
                if value1 is None:
                    value1 = ''
                if my_etree_dict1:
                    if bool(my_etree_dict1.get(key)) == False:
                        #if this key is not populated for the dictionary yet, then populate it
                        #print('False, there is no key for this yet')
                        my_etree_dict1[key] = value1
                        #print('populated')
                    else:#if this key already exists, then append to it
                        SI = my_etree_dict1[key]
                        my_etree_dict1[key] = SI + ',' + value1                   
                else:#if no dictionary exists yet populate
                    my_etree_dict1[key] = value1 
            #my_etree_dict[key]=value
            self.my_etree_dict1 = my_etree_dict1
        return my_etree_dict1
    
    def convert_xml_to_dict2(self):
        """
        This version exports out  entries into a dictionary using the dictionary
        xml_path_to_baseattribute
        defined below for USACE FGDC data
        The method may be modified if needed.
        
        my_etree_dict1={}

        for key in xml_path_to_baseattribute:
            my_etree_dict1[xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[8:])[0].text
        self.my_etree_dict1 = my_etree_dict1
        """           
        my_etree_dict1={}
        len_root_name_to_remove = len(self.xml_tree.tag)        
        for key in xml_path_to_baseattribute:
            my_etree_dict1[xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[len_root_name_to_remove:])[0].text
            #my_etree_dict1[xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[8:])[0].text
            #editing path to add ./ and then remove root name ('metadata'), the first 8 characters.
            #my_etree_dict1[xml_path_to_baseattribute[key]] = xml_data.xml_tree.findall('./' + key[8:])[0].text
            #pathlist = key.split('/')
        self.my_etree_dict1 = my_etree_dict1
        return my_etree_dict1
    
    def convert_xml_to_dict_ISO_FGDC(self):
        """
        This version exports out  entries into a dictionary using the dictionary
        xml_path_to_baseattribute
        defined below for USACE FGDC data
        The method may be modified if needed.
        
        my_etree_dict1={}

        for key in xml_path_to_baseattribute:
            my_etree_dict1[xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[8:])[0].text
        self.my_etree_dict1 = my_etree_dict1
        """           
        my_etree_dict1={}
        len_root_name_to_remove = len(self.xml_tree.tag)   
        vertdatum = {'metadata/spref/vertdef/altsys/altdatum':'altdatum'}
        for key in vertdatum:#iso_xml_path_to_baseattribute:
            #if isinstance(self.xml_tree.findall('./'+ key[len_root_name_to_remove:]), list) == True: #check if list
            #    my_etree_dict1[iso_xml_path_to_baseattribute[key]]  = self.xml_tree.findall('./'+ key[len_root_name_to_remove:][0])
            #else:
            my_etree_dict1[vertdatum[key]]  = self.xml_tree.find('./'+ key[len_root_name_to_remove:]).text
            my_etree_dict1['script: from_vert_key'] = my_etree_dict1[vertdatum[key]]
            #zz = self.xml_tree.find('./'+ key[len_root_name_to_remove:]
            #if type(zz) == xml.etree.ElementTree.Element
            #    print('true')                
            # my_etree_dict1[iso_xml_path_to_baseattribute[key]]  = self.xml_tree.find('./'+ key[len_root_name_to_remove:]).text
            #my_etree_dict1[iso_xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[len_root_name_to_remove:])[0].text
            #my_etree_dict1[xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[8:])[0].text
            #editing path to add ./ and then remove root name ('metadata'), the first 8 characters.
            #my_etree_dict1[xml_path_to_baseattribute[key]] = xml_data.xml_tree.findall('./' + key[8:])[0].text
            #pathlist = key.split('/')
        self.my_etree_dict1 = my_etree_dict1
        return my_etree_dict1
    
    def find_Instruments(self):
        """
        This method just takes out the Survey Instruments into a dictionary, 
        rather than loading the full attribute dictionary
        """            
        Survey_Instruments = {}
        my_etree=self.xml_tree            
        for S_INST in my_etree.iter('SURVEY_INSTRUMENT'):
            S_INST_ch  = S_INST.getchildren()
            print(S_INST_ch)
            for si in S_INST_ch:
                si_key = si.tag
                si_value = si.text
                if si_value is None:
                    si_value = ''
                print(si_key, si_value)
                try:
                    if Survey_Instruments:
                        if bool(Survey_Instruments.get(si_key)) == False:
                            #if this key is not populated for the dictionary yet, then populate it
                            #print('False, there is no key for this yet')
                            Survey_Instruments[si_key] = si_value
                            #print('populated')
                        else:#if this key already exists, then append to it
                            SI = Survey_Instruments[si_key]
                            Survey_Instruments[si_key] = SI + ',' + si_value                    
                    else:#if no dictionary exists yet populate
                        Survey_Instruments[si_key] = si_value 
                except:
                    print('problem with loop')
        return Survey_Instruments
            
    def get_fields(self):
        """
        Using the field available for the version type, get the data for those
        fields.
        """
        self.data = {}
        if self.version == 1.0:
            if 'filename' in self.source:
                self._read_file_name()
            if 'SORDAT' in self.source:
                self._read_SORDAT()
            if 'SURATH' in self.source:
                self._read_survey_authority()
            if 'SURSTA' in self.source:
                self._read_survey_start_date()
            if 'SUREND' in self.source:
                self._read_survey_end_date()
            if 'TECSOU' in self.source:
                self._read_tecsou()
            if 'DATUM' in self.source:
                self._read_datum()
            if 'survey' in self.source:
                self._read_survey_name()
            if 'planam' in self.source:
                self._read_planam()
            if 'sensor' in self.source:
                self._read_sensor_desc()
        elif self.version == 'HSMDB':
            for key in self.source.keys():
                self.data[key]=self.source[key]
            self.convert_dateformat_HSMDB()    
            """
            #converts format of SORDAT, SUREND & SURSTA
            YYYYMMDD
            populates self.data[]
            """
            
    def get_s57_dict(self):
        """
        Convert the object dictionary 'data' keys to the desired S57 keys and
        return the dicitonary.
        """
        s57 = {}
        for key in self.data.keys():
            if key == 'filename':
                label = self.data[key]
                p = re.compile('[A-Z][0-9][0-9][0-9][0-9][0-9]')
                s = p.findall(label)
                if len(s) == 1:
                    s57['OBJNAM'] = s[0]
                else:
                    s57['OBJNAM'] = label
            elif key == 'SORDAT':
                datestring = self.data[key]
                newds = datestring.replace('-','')
                no_time_ds = newds.split('T')[0]
                s57[key] = no_time_ds
            elif key == 'SURATH':
                s57[key] = self.data[key]
            elif key == 'SURSTA':
                datestring = self.data[key]
                newds = datestring.replace('-','')
                no_time_ds = newds.split('T')[0]
                s57[key] = no_time_ds
            elif key == 'SUREND':
                datestring = self.data[key]
                newds = datestring.replace('-','')
                no_time_ds = newds.split('T')[0]
                s57[key] = no_time_ds
            elif key == 'VERDAT':
                if self.data[key] in vert_datum:
                    s57[key] = vert_datum[self.data[key]]
                else:
                    s57[key] = vert_datum['Other']
            elif key == 'HORDAT':
                if self.data[key] in horz_datum:
                    s57[key] = horz_datum[self.data[key]]
                else:
                    print (self.data[key])
                    s57[key] = horz_datum['Local']
        if self.version=='HSMDB':
            s57['OBJNAM'] = self.data['survey']
        return s57
    
    #--------------------------------------------------------------------------
    def _extract_meta_CEMVN(self):
        if self.version =='USACE_FGDC':
            meta_xml = self.convert_xml_to_dict2()
            meta = parsing_xml_FGDC_attributes_s57(meta_xml)
            try:
                m = convert_meta_to_input(meta)
            except:
                print('still debugging')
                m={}
            meta_all_fields = {**meta_xml, **meta, **m}
        return meta_all_fields
    
    def _extract_meta_USACE_ISO(self):
        if self.version == 'ISO-8859-1':
            meta_xml = self.convert_xml_to_dict_ISO_FGDC()#convert_xml_to_dict2()
            #meta = parsing_xml_ISO_FGDC_attributes_s57(meta_xml)
            meta={}
            try:
                m = convert_meta_to_input(meta)
            except:
                print('still debugging')
                m={}
            meta_all_fields = {**meta_xml, **meta, **m}
        return meta_all_fields


            
    def _read_file_name(self):
        """ 
        Read the source file name and store it in the object 'data' dictionary
        with the key 'filename'.
        """

        try:
            ret = self.xml_tree.find(self.source['filename'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the survey name string")
            return

        try:
            self.data['filename'] = ret.text
        except Exception as e:
            log.warning("unable to read the survey name attribute: %s" % e)
            return
    #--------------------------------------------------------------------------
    def _read_SORDAT(self):
        """
        Reads a date, but what it means exactly needs to be researched...
        """
        if self.version==1.0:
            try:
                ret = self.xml_tree.find(self.source['SORDAT'],
                                          namespaces=self.ns)
            except:
                log.warning("unable to read the SORDAT date string")
                return
    
            try:
                text_date = ret.text
            except Exception as e:
                log.warning("unable to read the SORDAT date string: %s" % e)
                return
    
            tm_date = None
            try:
                parsed_date = parser.parse(text_date)
                tm_date = parsed_date.strftime('%Y%m%d')#('%Y-%m-%dT%H:%M:%SZ')
            except Exception:
                log.warning("unable to handle the date string: %s" % text_date)
        elif self.version=='HSMDB':
            date1 =  parser.parse(self.source['SORDAT'])
            tm_date = date1.strftime('%Y%m%d')#('%Y-%m-%dT%H:%M:%SZ')
        if tm_date is None:
            self.data['SORDAT'] = text_date
        else:
            self.data['SORDAT'] = tm_date
            
    def convert_dateformat_HSMDB(self):
        if self.version=='HSMDB':
            self.datekeys=['SURSTA', 'SORDAT', 'SUREND']
            for key in self.datekeys: 
                if key in self.source:
                    date1 =  parser.parse(self.source[key])
                    self.data[key]=date1.strftime('%Y%m%d')

        #for key in xml_data.datekeys: 
        #    if key in xml_data.source:
        #        print(key)
        #        date1 =  parser.parse(xml_data.source[key])
        #        xml_data.data[key]=date1.strftime('%Y%m%d')
                
    def _read_survey_authority(self):
        """ 
        Read the survey authority name and store it in the object 'data'
        dictionary with the key 'SURATH'.
        """
        if self.version==1.0:
            try:
                ret = self.xml_tree.find(self.source['SURATH'],
                                          namespaces=self.ns)
            except:
                log.warning("unable to read the survey authority name string")
                return
    
            try:
                if ret is not None:
                    self.data['SURATH'] = ret.text
            except Exception as e:
                log.warning("unable to read the survey authority name attribute: %s" % e)
                return
            
        elif self.version=='HSMDB':
            try:
                self.data['SURATH']=self.source['SURATH']
            except Exception as e:
                log.warning("unable to read the survey authority name attribute: %s" % e)
                return
        
    def _read_survey_start_date(self):
        """
        Read the survey start date store it in the object 'data'
        dictionary with the key 'SURSTA'.
        """

        try:
            rets = self.xml_tree.find(self.source['SURSTA'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the survey start date string")
            return

        if rets is not None:
            try:
                text_start_date = rets.text
            except Exception as e:
                log.warning("unable to read the survey start date string: %s" % e)
                return
            
            tms_date = None
            try:
                parsed_date = parser.parse(text_start_date)
                tms_date = parsed_date.strftime('%Y%m%d')#('%Y-%m-%dT%H:%M:%SZ')
            except Exception:
                log.warning("unable to handle the survey start string: %s" % text_start_date)
    
            if tms_date is None:
                self.data['SURSTA'] = text_start_date
            else:
                self.data['SURSTA'] = tms_date
            
    def _read_survey_end_date(self):
        """
        Read the survey end date and store it in the object 'data'
        dictionary with the key 'SUREND'.
        """
        try:
            rete = self.xml_tree.find(self.source['SUREND'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the survey end date string")
            return

        if rete is not None:
            try:
                text_end_date = rete.text
            except Exception as e:
                log.warning("unable to read the survey end date string: %s" % e)
                return
            
            tme_date = None
            try:
                parsed_date = parser.parse(text_end_date)
                tme_date = parsed_date.strftime('%Y%m%d')#('%Y-%m-%dT%H:%M:%SZ')
            except Exception:
                log.warning("unable to handle the survey end string: %s" % text_end_date)
    
            if tme_date is None:
                self.data['SUREND'] = text_end_date
            else:
                self.data['SUREND'] = tme_date
            
    def _read_tecsou(self):
        """
        Read tehcnology used for sounding the seafloor during the described
        survey.
        """

        try:
            ret = self.xml_tree.findall(self.source['TECSOU'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the TECSOU name string")
            return

        try:
            if len(ret) > 0:
                self.data['TECSOU'] = []
                for r in ret:
                    self.data['TECSOU'].append(r.text)
        except Exception as e:
            log.warning("unable to read the TECSOU attribute: %s" % e)
            return
            
    def _read_datum(self):
        """ 
        Read the survey vertical datum and store it in the object 'data'
        dictionary with the key 'VERDAT'.
        """

        try:
            ret = self.xml_tree.findall(self.source['DATUM'],
                                      namespaces=self.ns)                     
        except:
            log.warning("unable to read the survey datum name string(s)")
            return

        try:
            if len(ret) > 0:
                self.data['DATUM'] = []
                for r in ret:
                    self.data['DATUM'].append(r.text)
                
        except Exception as e:
            log.warning("unable to read the survey datum name attribute: %s" % e)
            return
        
    def _read_survey_name(self):
        """ 
        Read the survey name.
        """

        try:
            ret = self.xml_tree.find(self.source['survey'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the survey name string")
            return

        try:
            if ret is not None:
                self.data['survey'] = ret.text
        except Exception as e:
            log.warning("unable to read the survey name attribute: %s" % e)
            return
    
    def _read_planam(self):
        """ 
        Read the name of the survey platform.
        """

        try:
            ret = self.xml_tree.find(self.source['planam'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the survey platform name string")
            return

        try:
            if ret is not None:
                self.data['planam'] = ret.text
        except Exception as e:
            log.warning("unable to read the survey platform name attribute: %s" % e)
            return
        
    def _read_sensor_desc(self):
        """ 
        Read a description of the survey sensor.
        """

        try:
            ret = self.xml_tree.findall(self.source['sensor'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the sensor description name string")
            return

        try:
            if len(ret) > 0:
                self.data['sensor'] = []
                for r in ret:
                    self.data['sensor'].append(r.text)
        except Exception as e:
            log.warning("unable to read the sensor descriptioin name attribute: %s" % e)
            return
#------------------------------------------------------------------------------

def parse_namespace(meta_str):
    """
    Catch the xml and read the second line assuming it is the namespace
    information.  Return a namespace dictionary for use when parsing.
    """
    xmlns_loc = meta_str.find("xmlns:")
    xmlns_start = meta_str.rfind("<", 0, xmlns_loc)
    xmlns_end = meta_str.find(">", xmlns_start)
    ns_info = meta_str[xmlns_start: xmlns_end + 1]

    namespace = {}
    vals = ns_info.split('xmlns')
    for v in vals:
        # deal with case of missing '\n' by looking for the ':' following 'xmlns'.
        if v.find('xsi')<0:#here we are looking for the xsi section that has multiple "= signs"
            if v[0] == ':':
                tmp = v[1:]
                tmp = tmp.split('>')[0]
                name, info = tmp.split('=')
                site = info.split('"')[1]
                namespace[name] = site
        else:#this handles most cases
            if v[0] == ':':
                tmp = v[1:]
                tmp = tmp.split('>')[0]
                xsi_info = tmp.split('=')
                #site = xsi_info.split('"')[1]
                #namespace[name] = site
    return namespace

#------------------------------------------------------------------------------
def check_firstline(meta_xml):
    xml_version = ''
    if meta_xml.startswith('<?xml version="1.0" encoding="ISO-8859-1"?>\n'):
        print('ISO-8859-1 xml version')
        xml_version = 'ISO-8859-1'
        #return 'ISO-8859-1'
    return xml_version
       

#------------------------------------------------------------------------------
xml_path_to_baseattribute = {
        'metadata/idinfo/citation/citeinfo/origin':'origin',
        'metadata/idinfo/citation/citeinfo/pubdate':'pubdate',
        'metadata/idinfo/citation/citeinfo/title':'title',
        'metadata/idinfo/citation/citeinfo/geoform':'geoform',
        'metadata/idinfo/citation/citeinfo/pubinfo/pubplace':'pubplace',
        'metadata/idinfo/citation/citeinfo/pubinfo/publish':'publish',
        'metadata/idinfo/citation/citeinfo/onlink':'onlink',
        'metadata/idinfo/descript/abstract':'abstract',
        'metadata/idinfo/descript/purpose':'purpose',
        'metadata/idinfo/timeperd/timeinfo/sngdate/caldate':'caldate',
        'metadata/idinfo/timeperd/current':'current',
        'metadata/idinfo/status/progress':'progress',
        'metadata/idinfo/status/update':'update',
        'metadata/idinfo/spdom/bounding/westbc':'westbc',
        'metadata/idinfo/spdom/bounding/eastbc':'eastbc',
        'metadata/idinfo/spdom/bounding/northbc':'northbc',
        'metadata/idinfo/spdom/bounding/southbc':'southbc',
        'metadata/idinfo/keywords/theme/themekt':'themekt',
        'metadata/idinfo/keywords/theme/themekey':'themekey',
        'metadata/idinfo/keywords/place/placekt':'placekt',
        'metadata/idinfo/keywords/place/placekey':'placekey',
        'metadata/idinfo/accconst':'accconst',
        'metadata/idinfo/useconst':'useconst',
        'metadata/idinfo/ptcontac/cntinfo/cntorgp/cntorg':'cntorg',
        'metadata/idinfo/ptcontac/cntinfo/cntorgp/cntper':'cntper',
        'metadata/idinfo/ptcontac/cntinfo/cntpos':'cntpos',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/addrtype':'addrtype',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/address':'address',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/city':'city',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/state':'state',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/postal':'postal',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/country':'country',
        'metadata/idinfo/ptcontac/cntinfo/cntvoice':'cntvoice',
        'metadata/idinfo/ptcontac/cntinfo/cntfax':'cntfax',
        'metadata/idinfo/ptcontac/cntinfo/cntemail':'cntemail',
        'metadata/idinfo/datacred':'datacred',
        'metadata/idinfo/native':'native',
        'metadata/dataqual/attracc/attraccr':'attraccr',
        'metadata/dataqual/logic':'logic',
        'metadata/dataqual/complete':'complete',
        'metadata/dataqual/posacc/horizpa/horizpar':'horizpar',
        'metadata/dataqual/posacc/vertacc/vertaccr':'vertaccr',
        'metadata/dataqual/lineage/procstep/procdesc':'procdesc',
        'metadata/dataqual/lineage/procstep/procdate':'procdate',
        'metadata/spdoinfo/direct':'direct',
        'metadata/spdoinfo/ptvctinf/sdtsterm/sdtstype':'sdtstype',
        'metadata/spdoinfo/ptvctinf/sdtsterm/ptvctcnt':'ptvctcnt',
        'metadata/spref/horizsys/planar/gridsys/gridsysn':'gridsysn',
        'metadata/spref/horizsys/planar/gridsys/spcs/spcszone':'spcszone',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/stdparll':'stdparll',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/longcm':'longcm',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/latprjo':'latprjo',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/feast':'feast',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/fnorth':'fnorth',
        'metadata/spref/horizsys/planar/planci/plance':'plance',
        'metadata/spref/horizsys/planar/planci/coordrep/absres':'absres',
        'metadata/spref/horizsys/planar/planci/coordrep/ordres':'ordres',
        'metadata/spref/horizsys/planar/planci/plandu':'plandu',
        'metadata/spref/horizsys/geodetic/horizdn':'horizdn',
        'metadata/spref/horizsys/geodetic/ellips':'ellips',
        'metadata/spref/horizsys/geodetic/semiaxis':'semiaxis',
        'metadata/spref/horizsys/geodetic/denflat':'denflat',
        'metadata/eainfo/detailed/enttyp/enttypl':'enttypl',
        'metadata/eainfo/detailed/enttyp/enttypd':'enttypd',
        'metadata/eainfo/detailed/enttyp/enttypds':'enttypds',
        'metadata/eainfo/detailed/attr':'attr',
        'metadata/eainfo/overview/eaover':'eaover',
        'metadata/eainfo/overview/eadetcit':'eadetcit',
        'metadata/distinfo/distrib/cntinfo/cntorgp/cntorg':'cntorg',
        'metadata/distinfo/distrib/cntinfo/cntorgp/cntper':'cntper',
        'metadata/distinfo/distrib/cntinfo/cntpos':'cntpos',
        'metadata/distinfo/distrib/cntinfo/cntaddr/addrtype':'addrtype',
        'metadata/distinfo/distrib/cntinfo/cntaddr/address':'address',
        'metadata/distinfo/distrib/cntinfo/cntaddr/city':'city',
        'metadata/distinfo/distrib/cntinfo/cntaddr/state':'state',
        'metadata/distinfo/distrib/cntinfo/cntaddr/postal':'postal',
        'metadata/distinfo/distrib/cntinfo/cntaddr/country':'country',
        'metadata/distinfo/distrib/cntinfo/cntvoice':'cntvoice',
        'metadata/distinfo/distrib/cntinfo/cntfax':'cntfax',
        'metadata/distinfo/distrib/cntinfo/cntemail':'cntemail',
        'metadata/distinfo/distliab':'distliab',
        'metadata/distinfo/stdorder/digform/digtinfo/formname':'formname',
        'metadata/distinfo/stdorder/digform/digtopt/onlinopt/computer/networka/networkr':'networkr',
        'metadata/distinfo/stdorder/fees':'fees',
        'metadata/metainfo/metd':'metd',
        'metadata/metainfo/metc/cntinfo/cntperp/cntper':'cntper',
        'metadata/metainfo/metc/cntinfo/cntperp/cntorg':'cntorg',
        'metadata/metainfo/metc/cntinfo/cntpos':'cntpos',
        'metadata/metainfo/metc/cntinfo/cntaddr/addrtype':'addrtype',
        'metadata/metainfo/metc/cntinfo/cntaddr/address':'address',
        'metadata/metainfo/metc/cntinfo/cntaddr/city':'city',
        'metadata/metainfo/metc/cntinfo/cntaddr/state':'state',
        'metadata/metainfo/metc/cntinfo/cntaddr/postal':'postal',
        'metadata/metainfo/metc/cntinfo/cntaddr/country':'country',
        'metadata/metainfo/metc/cntinfo/cntvoice':'cntvoice',
        'metadata/metainfo/metc/cntinfo/cntfax':'cntfax',
        'metadata/metainfo/metc/cntinfo/cntemail':'cntemail',
        'metadata/metainfo/metstdn':'metstdn',
        'metadata/metainfo/metstdv':'metstdv'}

#------------------------------------------------------------------------------
iso_xml_path_to_baseattribute = {
		'metadata/idinfo/citation/citeinfo/origin':'origin',
		'metadata/idinfo/citation/citeinfo/pubdate':'pubdate',
		'metadata/idinfo/citation/citeinfo/title':'title',
		'metadata/idinfo/descript/abstract':'abstract',
		'metadata/idinfo/descript/purpose':'purpose',
		'metadata/idinfo/timeperd/timeinfo/sngdate/caldate':'caldate',
		'metadata/idinfo/timeperd/timeinfo/current':'current',
		'metadata/idinfo/status/progress':'progress',
		'metadata/idinfo/status/update':'update',
		'metadata/idinfo/spdom/bounding/westbc':'westbc',
		'metadata/idinfo/spdom/bounding/eastbc':'eastbc',
		'metadata/idinfo/spdom/bounding/northbc':'northbc',
		'metadata/idinfo/spdom/bounding/southbc':'southbc',
		'metadata/idinfo/keywords/themekt':'themekt',
		'metadata/idinfo/accconst':'accconst',
		'metadata/idinfo/useconst':'useconst',
		'metadata/idinfo/ptcontac/cntinfo/cntperp/cntper':'cntper',
		'metadata/idinfo/ptcontac/cntinfo/cntperp/cntorg':'cntorg',
		'metadata/idinfo/ptcontac/cntinfo/cntaddr/addrtype':'addrtype',
		'metadata/idinfo/ptcontac/cntinfo/cntaddr/address':'address',
		'metadata/idinfo/ptcontac/cntinfo/cntaddr/city':'city',
		'metadata/idinfo/ptcontac/cntinfo/cntaddr/state':'state',
		'metadata/idinfo/ptcontac/cntinfo/cntaddr/postal':'postal',
		'metadata/idinfo/ptcontac/cntinfo/cntvoice':'cntvoice',
		'metadata/spref/horizsys/planar/gridsys/gridsysn':'gridsysn',
		'metadata/spref/horizsys/planar/planci/plance':'plance',
		'metadata/spref/horizsys/planar/planci/absres':'absres',
		'metadata/spref/horizsys/planar/planci/ordres':'ordres',
		'metadata/spref/vertdef/altsys/altdatum':'altdatum',
		'metadata/spref/vertdef/altsys/altres':'altres',
		'metadata/spref/vertdef/altsys/altunits':'altunits',
		'metadata/spref/vertdef/altsys/altenc':'altenc',
		'metadata/spref/vertdef/depthsys/depthdn':'depthdn',
		'metadata/spref/vertdef/depthsys/depthres':'depthres',
		'metadata/spref/vertdef/depthsys/depthdu':'depthdu',
		'metadata/spref/vertdef/depthsys/depthem':'depthem',
		'metadata/metainfo/metd':'metd',
		'metadata/metainfo/metc/cntinfo/cntperp/cntper':'cntper',
		'metadata/metainfo/metc/cntinfo/cntperp/cntorg':'cntorg',
		'metadata/metainfo/metc/cntinfo/cntaddr/addrtype':'addrtype',
		'metadata/metainfo/metc/cntinfo/cntaddr/address':'address',
		'metadata/metainfo/metc/cntinfo/cntaddr/city':'city',
		'metadata/metainfo/metc/cntinfo/cntaddr/state':'state',
		'metadata/metainfo/metc/cntinfo/cntaddr/postal':'postal',
		'metadata/metainfo/metc/cntinfo/cntvoice':'cntvoice',
		'metadata/metainfo/metstdn':'metstdn',
		'metadata/metainfo/metstdv':'metstdv',
		'metadata/ellips':'ellips',
		'metadata/#text':'#text',
        }
#------------------------------------------------------------------------------

xmlbasename_to_index1 = {
        'origin':'0',
        'pubdate':'1',
        'title':'2',
        'geoform':'3',
        'pubplace':'4',
        'publish':'5',
        'onlink':'6',
        'abstract':'7',
        'purpose':'8',
        'caldate':'9',
        'current':'10',
        'progress':'11',
        'update':'12',
        'westbc':'13',
        'eastbc':'14',
        'northbc':'15',
        'southbc':'16',
        'themekt':'17',
        'themekey':'18',
        'placekt':'19',
        'placekey':'20',
        'accconst':'21',
        'useconst':'22',
        'cntorg':'23',
        'cntper':'24',
        'cntpos':'25',
        'addrtype':'26',
        'address':'27',
        'city':'28',
        'state':'29',
        'postal':'30',
        'country':'31',
        'cntvoice':'32',
        'cntfax':'33',
        'cntemail':'34',
        'datacred':'35',
        'native':'36',
        'attraccr':'37',
        'logic':'38',
        'complete':'39',
        'horizpar':'40',
        'vertaccr':'41',
        'procdesc':'42',
        'procdate':'43',
        'direct':'44',
        'sdtstype':'45',
        'ptvctcnt':'46',
        'gridsysn':'47',
        'spcszone':'48',
        'stdparll':'49',
        'longcm':'50',
        'latprjo':'51',
        'feast':'52',
        'fnorth':'53',
        'plance':'54',
        'absres':'55',
        'ordres':'56',
        'plandu':'57',
        'horizdn':'58',
        'ellips':'59',
        'semiaxis':'60',
        'denflat':'61',
        'enttypl':'62',
        'enttypd':'63',
        'enttypds':'64',
        'attr':'65',
        'eaover':'66',
        'eadetcit':'67',
        'cntorg':'68',
        'cntper':'69',
        'cntpos':'70',
        'addrtype':'71',
        'address':'72',
        'city':'73',
        'state':'74',
        'postal':'75',
        'country':'76',
        'cntvoice':'77',
        'cntfax':'78',
        'cntemail':'79',
        'distliab':'80',
        'formname':'81',
        'networkr':'82',
        'fees':'83',
        'metd':'84',
        'cntper':'85',
        'cntorg':'86',
        'cntpos':'87',
        'addrtype':'88',
        'address':'89',
        'city':'90',
        'state':'91',
        'postal':'92',
        'country':'93',
        'cntvoice':'94',
        'cntfax':'95',
        'cntemail':'96',
        'metstdn':'97',
        'metstdv':'98'}

#------------------------------------------------------------------------------
xml_path_to_index1 ={
        'metadata/idinfo/citation/citeinfo/origin':'0',
        'metadata/idinfo/citation/citeinfo/pubdate':'1',
        'metadata/idinfo/citation/citeinfo/title':'2',
        'metadata/idinfo/citation/citeinfo/geoform':'3',
        'metadata/idinfo/citation/citeinfo/pubinfo/pubplace':'4',
        'metadata/idinfo/citation/citeinfo/pubinfo/publish':'5',
        'metadata/idinfo/citation/citeinfo/onlink':'6',
        'metadata/idinfo/descript/abstract':'7',
        'metadata/idinfo/descript/purpose':'8',
        'metadata/idinfo/timeperd/timeinfo/sngdate/caldate':'9',
        'metadata/idinfo/timeperd/current':'10',
        'metadata/idinfo/status/progress':'11',
        'metadata/idinfo/status/update':'12',
        'metadata/idinfo/spdom/bounding/westbc':'13',
        'metadata/idinfo/spdom/bounding/eastbc':'14',
        'metadata/idinfo/spdom/bounding/northbc':'15',
        'metadata/idinfo/spdom/bounding/southbc':'16',
        'metadata/idinfo/keywords/theme/themekt':'17',
        'metadata/idinfo/keywords/theme/themekey':'18',
        'metadata/idinfo/keywords/place/placekt':'19',
        'metadata/idinfo/keywords/place/placekey':'20',
        'metadata/idinfo/accconst':'21',
        'metadata/idinfo/useconst':'22',
        'metadata/idinfo/ptcontac/cntinfo/cntorgp/cntorg':'23',
        'metadata/idinfo/ptcontac/cntinfo/cntorgp/cntper':'24',
        'metadata/idinfo/ptcontac/cntinfo/cntpos':'25',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/addrtype':'26',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/address':'27',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/city':'28',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/state':'29',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/postal':'30',
        'metadata/idinfo/ptcontac/cntinfo/cntaddr/country':'31',
        'metadata/idinfo/ptcontac/cntinfo/cntvoice':'32',
        'metadata/idinfo/ptcontac/cntinfo/cntfax':'33',
        'metadata/idinfo/ptcontac/cntinfo/cntemail':'34',
        'metadata/idinfo/datacred':'35',
        'metadata/idinfo/native':'36',
        'metadata/dataqual/attracc/attraccr':'37',
        'metadata/dataqual/logic':'38',
        'metadata/dataqual/complete':'39',
        'metadata/dataqual/posacc/horizpa/horizpar':'40',
        'metadata/dataqual/posacc/vertacc/vertaccr':'41',
        'metadata/dataqual/lineage/procstep/procdesc':'42',
        'metadata/dataqual/lineage/procstep/procdate':'43',
        'metadata/spdoinfo/direct':'44',
        'metadata/spdoinfo/ptvctinf/sdtsterm/sdtstype':'45',
        'metadata/spdoinfo/ptvctinf/sdtsterm/ptvctcnt':'46',
        'metadata/spref/horizsys/planar/gridsys/gridsysn':'47',
        'metadata/spref/horizsys/planar/gridsys/spcs/spcszone':'48',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/stdparll':'49',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/longcm':'50',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/latprjo':'51',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/feast':'52',
        'metadata/spref/horizsys/planar/gridsys/spcs/lambertc/fnorth':'53',
        'metadata/spref/horizsys/planar/planci/plance':'54',
        'metadata/spref/horizsys/planar/planci/coordrep/absres':'55',
        'metadata/spref/horizsys/planar/planci/coordrep/ordres':'56',
        'metadata/spref/horizsys/planar/planci/plandu':'57',
        'metadata/spref/horizsys/geodetic/horizdn':'58',
        'metadata/spref/horizsys/geodetic/ellips':'59',
        'metadata/spref/horizsys/geodetic/semiaxis':'60',
        'metadata/spref/horizsys/geodetic/denflat':'61',
        'metadata/eainfo/detailed/enttyp/enttypl':'62',
        'metadata/eainfo/detailed/enttyp/enttypd':'63',
        'metadata/eainfo/detailed/enttyp/enttypds':'64',
        'metadata/eainfo/detailed/attr':'65',
        'metadata/eainfo/overview/eaover':'66',
        'metadata/eainfo/overview/eadetcit':'67',
        'metadata/distinfo/distrib/cntinfo/cntorgp/cntorg':'68',
        'metadata/distinfo/distrib/cntinfo/cntorgp/cntper':'69',
        'metadata/distinfo/distrib/cntinfo/cntpos':'70',
        'metadata/distinfo/distrib/cntinfo/cntaddr/addrtype':'71',
        'metadata/distinfo/distrib/cntinfo/cntaddr/address':'72',
        'metadata/distinfo/distrib/cntinfo/cntaddr/city':'73',
        'metadata/distinfo/distrib/cntinfo/cntaddr/state':'74',
        'metadata/distinfo/distrib/cntinfo/cntaddr/postal':'75',
        'metadata/distinfo/distrib/cntinfo/cntaddr/country':'76',
        'metadata/distinfo/distrib/cntinfo/cntvoice':'77',
        'metadata/distinfo/distrib/cntinfo/cntfax':'78',
        'metadata/distinfo/distrib/cntinfo/cntemail':'79',
        'metadata/distinfo/distliab':'80',
        'metadata/distinfo/stdorder/digform/digtinfo/formname':'81',
        'metadata/distinfo/stdorder/digform/digtopt/onlinopt/computer/networka/networkr':'82',
        'metadata/distinfo/stdorder/fees':'83',
        'metadata/metainfo/metd':'84',
        'metadata/metainfo/metc/cntinfo/cntperp/cntper':'85',
        'metadata/metainfo/metc/cntinfo/cntperp/cntorg':'86',
        'metadata/metainfo/metc/cntinfo/cntpos':'87',
        'metadata/metainfo/metc/cntinfo/cntaddr/addrtype':'88',
        'metadata/metainfo/metc/cntinfo/cntaddr/address':'89',
        'metadata/metainfo/metc/cntinfo/cntaddr/city':'90',
        'metadata/metainfo/metc/cntinfo/cntaddr/state':'91',
        'metadata/metainfo/metc/cntinfo/cntaddr/postal':'92',
        'metadata/metainfo/metc/cntinfo/cntaddr/country':'93',
        'metadata/metainfo/metc/cntinfo/cntvoice':'94',
        'metadata/metainfo/metc/cntinfo/cntfax':'95',
        'metadata/metainfo/metc/cntinfo/cntemail':'96',
        'metadata/metainfo/metstdn':'97',
        'metadata/metainfo/metstdv':'98'}
#------------------------------------------------------------------------------

def parsing_xml_FGDC_attributes_s57(meta_xml):
    """
    #PARSING XML attributes
    'Survey Type' = themekey.find('Condition Survey')
    
    #QC_checks
    #if expected results found ok, if not trigger more QC:
        logic.find("Horizontal_Positional_Accuracy_Explanation: Static Test")
        logic.find("Vertical_Accuracy_Explanation: Bar/Ball check")
    #QC horizontal coordinate system    
    if gridsysn =='State Plane Coordinate System 1983':
        print('expected')
    else:
        print('double check horizontal reference system')
    fipstr = spcszone
    #QC check against FIPS from table
    
    if plandu == 'Foot_US':
        horiz_units = plandu      
    #plandu = #horizontal units    
    if horizdn == 'D_North_American_1983':
        horizontal_datum = horizdn
    if ellips =='GRS_1980':
        ellipsoid_v = ellips
    # horizontal uncertainty check
    if  horizpar.find('RTK'):
        ?
        horiz_uncert = ?
    elif horizpar.find('DGPS, 1 Meter') == True:
        horiz_uncert='1'# (POSACC)    
    if vertaccr.find('Bar Checked to 0.1 foot') == True:
        vert_unc
        
        
    if procdesc.find():
    procdesc
        -> to coverage single beam (Cat B type coverage)
        'SonarSystem', 'SonarManufacturer',
         if =='Odom MKIII echosounder':
             'single beam'
             
   abstract = meta_xml['abstract']        
   if abstract.find('Survey Type: Single Beam Soundings') >= 0:
       TECSOU= 'single beam'
       
    -> 1= TECSOU, (NO)=1, 1= f_dict. 1= f_lstd, 9999=f_size, 1 =flbath, 1=flcvrg | abstract.find('Vertical Datum:) 
    -> find(Mean Lower Low Water) 
    ->(pass to function) 
    -> VERTDAT | abstract.find('State Plane Coordinate System (SPCS),
    """
    #------------------------------------------------------------------------------
    m={}    
    abstract = meta_xml['abstract']
    lines = abstract.split('\n')
    for line in lines:
        try:
            if line.find('Survey Type: ') >= 0:
                name = line.split('Survey Type:')[-1]
                m['survey_description'] = name
                if line.find('Survey Type: Single Beam Soundings') >= 0:
                   m['TECSOU']= '1'#'single beam'
                   m['f_dict'] = '1'
                   m['f_lstd'] = '1'
                   m['f_size'] = '9999'
                   m['flbath'] = '1'
                   m['flcvrg'] = '1' #where '1' = 'NO'
                elif line.find('Single Beam Soundings') >= 0:
                   m['TECSOU']= '1'#'single beam'
                   m['f_dict'] = '1'
                   m['f_lstd'] = '1'
                   m['f_size'] = '9999'
                   m['flbath'] = '1'
                   m['flcvrg'] = '1' #where '1' = 'NO'
                elif line.find('Multi Beam Soundings') >= 0:
                   m['TECSOU']= '3'#'multi beam'
                   m['f_dict'] = '1'
                   m['f_lstd'] = '1'
                   m['f_size'] = '9999'
                   m['flbath'] = '1'
                   m['flcvrg'] = '1' #where '1' = 'NO'
                else:
                   m['TECSOU']= ''#'multi beam'
                   m['f_dict'] = ''
                   m['f_lstd'] = ''
                   m['f_size'] = ''
                   m['flbath'] = ''
                   m['flcvrg'] = '' #where '1' = 'NO'
                    
            if  line.find('Horizontal Coordinate System:') >= 0:
                name = line.split('Horizontal Coordinate System:')[-1]
                m['horizontal_datum_i'] = name
    
            if  line.find('Coordinate System (SPCS)') >= 0:
                name = line.split('Coordinate System (SPCS),')[-1]
                name = name.split('. Distance units in')
                m['SPCS'] = name[0]
                m['Horizontal_Units'] = name[1]
                                  
            if  line.find('Vertical Datum:') >= 0:
                name = line.split('Vertical Datum:')[-1]
                m['Vertical Datum Description']= name
                if name.find('Soundings are shown in feet and indicate depths below Mean Lower Low Water') >= 0:
                    m['VERTDAT'] = 'MLLW'
                elif abstract.find('LWRP') >= 0:
                    m['VERTDAT'] = 'LWRP'
                elif abstract.find('Low Water Reference Plane 2007 (LWRP07)') >= 0:
                    m['VERTDAT'] = 'LWRP'
                elif abstract.find('MLG') >= 0:
                    m['VERTDAT'] = 'MLG'
                    print(name)
                elif abstract.find('depths below National Geodetic Vertical Datum or 1929 (NGVD29)') >= 0:
                    m['VERTDAT'] = 'NGVD29'
                if name.find('Soundings are shown in feet') >= 0:
                    m['script: from_vert_units'] = 'US Survey Foot'
        except:
            #Other way to split abstract, in case format changed over time
            print('issue parsing')
            if abstract.find('Survey Type: Single Beam Soundings') >= 0:
               m['TECSOU']= '1'#'single beam'
               m['f_dict'] = '1'
               m['f_lstd'] = '1'
               m['f_size'] = '9999'
               m['flbath'] = '1'
               m['flcvrg'] = '1' #where '1' = 'NO'
            else:
                m['TECSOU']= ''#'multi beam'
                m['f_dict'] = ''
                m['f_lstd'] = ''
                m['f_size'] = ''
                m['flbath'] = ''
                m['flcvrg'] = '' #where '1' = 'NO'
            if abstract.find('Vertical Datum:') >= 0:
                if abstract.find('Soundings are shown in feet and indicate depths below Mean Lower Low Water') >= 0:
                    m['VERTDAT'] = 'MLLW'                    
                elif abstract.find('LWRP') >= 0:
                    m['VERTDAT'] = 'LWRP'
                elif abstract.find('Low Water Reference Plane 2007 (LWRP07)') >= 0:
                    m['VERTDAT'] = 'LWRP'
                elif abstract.find('MLG') >= 0:
                    m['VERTDAT'] = 'MLG'
                    print(abstract)
                elif abstract.find('depths below National Geodetic Vertical Datum or 1929 (NGVD29)') >= 0:
                    m['VERTDAT'] = 'NGVD29'
                    print(abstract)
            else:
                m['VERTDAT'] = ''
    return m
                    
def convert_meta_to_input(m):                   
    m['from_vert_datum'] = m['Vertical Datum Description']
    #m['script: from_vert_units']
    m['from_horiz_datum'] = m['horizontal_datum_i'] + ',' + m['SPCS']
    m['from_horiz_units'] = m['Horizontal_Units']#may need to enforce some kind of uniform spelling etc. here
    m['script: from_vert_key'] = m['VERTDAT']
    m['script: feat_size'] = m['f_size']
    m['script: feat_detect'] = m['f_dict'] 
    m['script: feat_least_depth'] =m['f_lstd']
    m['script: complete_coverage'] = m['flcvrg']
    m['script: complete_bathymetry'] = m['flbath']
    return m
            

#------------------------------------------------------------------------------