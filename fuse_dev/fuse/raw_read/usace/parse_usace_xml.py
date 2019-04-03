"""
parse_usace_xml.py


This is a collection of tools for working with metadata from USACE e-hydro metadata. Created
specifically to extract the metadata for building a bathymetric database,
this collection of method attempt to both serve extraction of the meta data in
a general sense, and also for specific S57 needs.

J Kinney 
update April 3, 2019

"""
import re
from xml.etree import ElementTree as et 
import logging as log
from os import path

try:
    import dateutil.parser as parser
except:
    pass

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
    xml_data = XML_Meta(xml_txt, filename = xmlbasename)
    s57_dict = xml_data.get_s57_dict()
    return s57_dict
        
class XML_Meta(object):
    """ 
    Helper class to manage xml metadata. This class takes an xml string and
    parses the string based on a dictionary (named 'source') with keys that
    name the item to extract and the tree path string to find the data.
    Different versions of the source dictionary can be set based on the version
    of the metadata being parsed.  All extracted data is placed in a dictionary
    (named 'data').
    """

    def __init__(self, meta_xml, version = '', filename = ''):
        """
        Provided an xml string for parsing, a tree will be created and the
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
        Try to guess the version of the metadata for the purpose of parsing the
        fields.  The guess is based off comparing the name spaces and or first lines of xml
        files.
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
        else:
            return -1.0
            print('We do not have a template for this verion yet!')

    def _set_format(self):
        """
        Set the locations of the desired data types based on the version of the
        metadata.
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
                self.source['survey'] = './/gmd:fileIdentifier/gco:CharacterString' #updated
                self.source['planam'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:platform/gmi:MI_Platform/gmi:identifier/gmd:RS_Identifier/gmd:code/gco:CharacterString' # updated#
                self.source['sensor'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:instrument/gmi:MI_Instrument/gmi:description/gco:CharacterString' #updated               
            else:
                log.warning("verison not compatible")
   
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
            grandchildren = ch.getchildren()
            for ss in grandchildren:
                print(ss.tag, ss.text)
                tag1=ss.tag
                key=tag1
                value1=ss.text
                if value1 is None:
                    value1 = ''
                if my_etree_dict1:
                    if bool(my_etree_dict1.get(key)) == False:
                        #if this key is not populated for the dictionary yet, then populate it
                        my_etree_dict1[key] = value1
                    else:#if this key already exists, then append to it
                        SI = my_etree_dict1[key]
                        my_etree_dict1[key] = SI + ',' + value1                   
                else:#if no dictionary exists yet populate
                    my_etree_dict1[key] = value1 
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
            if self.xml_tree.findall('./' + key[len_root_name_to_remove:]):
                if self.xml_tree.findall('./'+ key[len_root_name_to_remove:]) is None:
                     my_etree_dict1[xml_path_to_baseattribute[key]] = ''
                elif len(self.xml_tree.findall('./' + key[len_root_name_to_remove:]))>0:
                    my_etree_dict1[xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[len_root_name_to_remove:])[0].text
            #editing path to add ./ and then remove root name ('metadata'), the first 8 characters in this case.
        self.my_etree_dict1 = my_etree_dict1
        return my_etree_dict1
    
    def convert_xml_to_dict_ISO_FGDC(self):
        """
        This version exports out entries into a dictionary using the dictionary
        iso_xml_path_to_baseattribute
        defined below for USACE ISO FGDC data
        The method may be modified if needed.
        
        my_etree_dict1={}

        for key in xml_path_to_baseattribute:
            my_etree_dict1[iso_xml_path_to_baseattribute[key]] = self.xml_tree.findall('./' + key[8:])[0].text
        self.my_etree_dict1 = my_etree_dict1
        
        vertical datum is returned in my_etree_dict1['from_vert_key']
        """           
        my_etree_dict1={}
        len_root_name_to_remove = len(self.xml_tree.tag)   
        vertdatum = {'metadata/spref/vertdef/altsys/altdatum':'altdatum'}
        for key in iso_xml_path_to_baseattribute:
            if self.xml_tree.findall('./'+ key[len_root_name_to_remove:]):
                if self.xml_tree.findall('./'+ key[len_root_name_to_remove:]) is None:
                    my_etree_dict1[iso_xml_path_to_baseattribute[key]]  = ''
                elif isinstance(self.xml_tree.findall('./'+ key[len_root_name_to_remove:]), list) == True: #check if list
                    if len(self.xml_tree.find('./'+ key[len_root_name_to_remove:]))>0:            
                        my_etree_dict1[iso_xml_path_to_baseattribute[key]]  = self.xml_tree.find('./'+ key[len_root_name_to_remove:])[0].text
                    else:
                        my_etree_dict1[iso_xml_path_to_baseattribute[key]]  = self.xml_tree.find('./'+ key[len_root_name_to_remove:]).text                
            else:
                my_etree_dict1[iso_xml_path_to_baseattribute[key]]  = ''
        for key in vertdatum:#
            if self.xml_tree.findall('./'+ key[len_root_name_to_remove:]):
                if isinstance(self.xml_tree.findall('./'+ key[len_root_name_to_remove:]), list) == True: #check if list                                    
                    if len(self.xml_tree.find('./'+ key[len_root_name_to_remove:]))>0:
                        if self.xml_tree.find('./'+ key[len_root_name_to_remove:]) is None:
                             my_etree_dict1['script: from_vert_key'] = ''
                             #Checks for NoneType object ('None')
                        else:
                            my_etree_dict1[vertdatum[key]]  = self.xml_tree.find('./'+ key[len_root_name_to_remove:]).text
                            my_etree_dict1['from_vert_key'] = my_etree_dict1[vertdatum[key]]
            else:
                my_etree_dict1['from_vert_key'] = ''
            my_etree_dict1['script: from_vert_key'] = my_etree_dict1['from_vert_key']
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
                            Survey_Instruments[si_key] = si_value
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
        """
        retrieves USACE e-hydro metadata that follows the FGDC format 
        and returns a dictionary
        """
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
        """
        retrieves USACE e-hydro metadata that follows the ISO-8859-1 
        FGDC format and returns a dictionary
        """
        if self.version == 'ISO-8859-1':
            meta_xml = self.convert_xml_to_dict_ISO_FGDC()#
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
                tm_date = parsed_date.strftime('%Y%m%d')
            except Exception:
                log.warning("unable to handle the date string: %s" % text_date)
        elif self.version=='HSMDB':
            date1 =  parser.parse(self.source['SORDAT'])
            tm_date = date1.strftime('%Y%m%d')
        if tm_date is None:
            self.data['SORDAT'] = text_date
        else:
            self.data['SORDAT'] = tm_date
            
                
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
                tme_date = parsed_date.strftime('%Y%m%d')
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
    return namespace

#------------------------------------------------------------------------------
def check_firstline(meta_xml):
    xml_version = ''
    if meta_xml.startswith('<?xml version="1.0" encoding="ISO-8859-1"?>\n'):
        print('ISO-8859-1 xml version')
        xml_version = 'ISO-8859-1'        
    return xml_version#returns 'ISO-8859-1'
       
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
        'metadata/metainfo/metstdv':'metstdv',
        'metadata/idinfo/citation/citeinfo/pubtime':'pubtime',
        'metadata/idinfo/keywords/theme':'theme',
        'metadata/idinfo/ptcontac':'ptcontac',
        'metadata/spref/horizsys/planar/mapproj/mapprojn':'mapprojn',
        'metadata/spref/horizsys/planar/mapproj/transmer/sfctrmer':'sfctrmer',
        'metadata/spref/horizsys/planar/mapproj/transmer/longcm':'longcm',
        'metadata/spref/horizsys/planar/mapproj/transmer/latprjo':'latprjo',
        'metadata/spref/horizsys/planar/mapproj/transmer/feast':'feast',
        'metadata/spref/horizsys/planar/mapproj/transmer/fnorth':'fnorth',
        'metadata/metainfo/metc/cntinfo/cntorgp/cntorg':'cntorg',
        'metadata/metainfo/metc/cntinfo/cntorgp/cntper':'cntper',
        'metadata/metainfo/mettc':'mettc',
        }
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
        }
#		'metadata/#text':'#text',
#------------------------------------------------------------------------------
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
                elif line.find('Single Beam Soundings') >= 0:
                   m['TECSOU']= '1'#'single beam'
                elif line.find('Multi Beam Soundings') >= 0:
                   m['TECSOU']= '3'#'multi beam'
                else:
                   m['TECSOU']= ''#'multi beam'/'single beam' etc.                    
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
                elif name.find('Soundings are shown in feet and are referenced to Mean Lower Low Water') >= 0:#CESAM
                    m['VERTDAT'] = 'MLLW'
                elif name.find('Values are based on the National Geodetic Vertical Datum (NGVD) of 1929') >=0:#CESAM
                    m['VERTDAT'] = 'NGVD29'
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

            else:
                m['TECSOU']= ''#'multi beam' 3
            if abstract.find('Vertical Datum:') >= 0:
                if abstract.find('Soundings are shown in feet and indicate depths below Mean Lower Low Water') >= 0:
                    m['VERTDAT'] = 'MLLW'
                elif name.find('Soundings are shown in feet and are referenced to Mean Lower Low Water') >= 0:#CESAM
                    m['VERTDAT'] = 'MLLW'                    
                elif abstract.find('LWRP') >= 0:
                    m['VERTDAT'] = 'LWRP'
                elif abstract.find('Low Water Reference Plane 2007 (LWRP07)') >= 0:
                    m['VERTDAT'] = 'LWRP'
                elif abstract.find('MLG') >= 0:
                    m['VERTDAT'] = 'MLG'
                    print(abstract)
                elif name.find('Values are based on the National Geodetic Vertical Datum (NGVD) of 1929') >=0:#CESAM
                    m['VERTDAT'] = 'NGVD29'
                elif abstract.find('depths below National Geodetic Vertical Datum or 1929 (NGVD29)') >= 0:
                    m['VERTDAT'] = 'NGVD29'
                    print(abstract)
            else:
                m['VERTDAT'] = ''
    return m
                    
def convert_meta_to_input(m):
    """
    m = convert_meta_to_input(m)
    maps dictionary keys to new keys
    """                   
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
    