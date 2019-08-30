"""
parse_ncei_meta.py
Edit 1/3/2018 J. Kinney modifying for HSMDB format xml from Jason Baillio (FGDC not ISO).


20180302 G.Rice

V0.0.1 Last Updated 20180312

This is a collection of tools for working with metadata from NCEI. Created
specifically to extract the metadata for building a bathymetric database,
this collection of method attempt to both serve extraction of the meta data in
a general sense, and also for specific S57 needs.

"""
import re
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
    xml_data = Meta(xml_txt, filename = xmlbasename)
    s57_dict = xml_data.get_s57_dict()
    return s57_dict

class Meta(object):
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

        if len(version) > 0:
            self.version = float(version)
        else:
            self.version = self._guess_version()
            
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

        version_HSMDB ={'HSMDB'}

        if self.ns == version_1:
            return 1.0
        elif self.xml_tree.tag == 'HSMDB':
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

        elif version =='HSMDB':
            self.source = {}
            try:
                my_etree_dict= convert_xml_to_dict(self.my_etree)
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
                
        return s57
        
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

    def _read_SORDAT(self):
        """
        Reads a date, but what it means exactly needs to be researched...
        """

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
            tm_date = parsed_date.strftime('%Y-%m-%dT%H:%M:%SZ')
        except Exception:
            log.warning("unable to handle the date string: %s" % text_date)

        if tm_date is None:
            self.data['SORDAT'] = text_date
        else:
            self.data['SORDAT'] = tm_date
            
    def _read_survey_authority(self):
        """ 
        Read the survey authority name and store it in the object 'data'
        dictionary with the key 'SURATH'.
        """

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
                tms_date = parsed_date.strftime('%Y-%m-%dT%H:%M:%SZ')
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
                tme_date = parsed_date.strftime('%Y-%m-%dT%H:%M:%SZ')
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
