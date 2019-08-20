"""
parse_bag_meta.py

20180222 G.Rice

V0.0.5 Last Updated 20180403

This is a collection of tools for working with BAG metadata.  Created
specifically to extract the bag metadata for building a bathymetric database,
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

__version__ = 'parse_bag_meta 0.0.4'
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
        'HRD'   : '24',  # Adding this for the Hudson River Datum
        }
        
horz_datum = {
        'WGS72' : '1',
        'WGS84' : '2',
        'WGS_1984' : '2',
        'NAD27' : '74',
        'NAD83' : '75',
        'North_American_Datum_1983' : '75',
        'Local' : '131',
        }

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
        post15 = {
                'bag': 'http://www.opennavsurf.org/schema/bag',
                'gco': 'http://www.isotc211.org/2005/gco',
                'gmd': 'http://www.isotc211.org/2005/gmd',
                'gmi': 'http://www.isotc211.org/2005/gmi',
                'gml': 'http://www.opengis.net/gml/3.2',
                'xlink': 'http://www.w3.org/1999/xlink',
                'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                }
        pre15 = {
                'smXML':"http://metadata.dgiwg.org/smXML",
                "xlink":"http://www.w3.org/1999/xlink",
                'xsi':"http://www.w3.org/2001/XMLSchema-instance", 
                'gml':"http://www.opengis.net/gml"
                }
        if self.ns == post15:
            return 1.5
        elif self.ns == pre15:
            return 1.0
        else:
            return -1.0

    def _set_format(self):
        """
        Set the locations of the desired data types based on the version of the
        bag.
        """
        version = float(self.version)
        
        if version >= 1.5:
            self.source = {}
            # these are the bag types
            self.source['rows_cols'] = './/gmd:spatialRepresentationInfo/gmd:MD_Georectified/gmd:axisDimensionProperties/gmd:MD_Dimension/gmd:dimensionSize/gco:Integer'
            self.source['resx_resy'] = './/gmd:spatialRepresentationInfo/gmd:MD_Georectified/gmd:axisDimensionProperties/gmd:MD_Dimension/gmd:resolution/gco:Measure'                          
            self.source['bbox'] = './/gmd:spatialRepresentationInfo/gmd:MD_Georectified/gmd:cornerPoints/gml:Point/gml:coordinates'
            self.source['wkt_srs'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString'
            self.source['lon_min'] = './/gmd:EX_GeographicBoundingBox/gmd:westBoundLongitude/gco:Decimal'
            self.source['lon_max'] = './/gmd:EX_GeographicBoundingBox/gmd:eastBoundLongitude/gco:Decimal'
            self.source['lat_min'] = './/gmd:EX_GeographicBoundingBox/gmd:southBoundLatitude/gco:Decimal'
            self.source['lat_max'] = './/gmd:EX_GeographicBoundingBox/gmd:northBoundLatitude/gco:Decimal'
            self.source['abstract'] = './/gmd:abstract/gco:CharacterString'
            self.source['date'] = './/gmd:CI_Date/gmd:date/gco:Date'
            self.source['unc_type'] = './/bag:verticalUncertaintyType/bag:BAG_VertUncertCode'
            self.source['z_min'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:verticalElement/gmd:EX_VerticalExtent/gmd:minimumValue/gco:Real'
            self.source['z_max'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:verticalElement/gmd:EX_VerticalExtent/gmd:maximumValue/gco:Real'
            self.source['filename'] = self.filename
            self.source['sourcename'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString'
            self.source['SORDAT'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date'
            self.source['SURATH'] = './/gmd:contact/gmd:CI_ResponsibleParty/gmd:organisationName/gco:CharacterString'
            self.source['SUREND'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date' 
            self.source['SURSTA'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:temporalElement/gmd:EX_TemporalExtent/gmd:extent/gml:TimePeriod/gml:beginPosition'
            self.source['VERDAT'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString'
            self.source['HORDAT'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString'
            self.source['planam'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:platform/gmi:MI_Platform/gmi:identifier/gmd:MD_Identifier/gmd:code/gco:CharacterString'
            self.source['sensor'] = './/gmi:instrument/gmi:MI_Instrument/gmi:identifier/gmd:MD_Identifier/gmd:code/gco:CharacterString'
       
        elif version < 1.5 and version > 0:
            self.source = {}
            self.source['filename'] = self.filename
            self.source['abstract'] = './/identificationInfo/smXML:BAG_DataIdentification/abstract'
            self.source['SORDAT'] = './/identificationInfo/smXML:BAG_DataIdentification/citation/smXML:CI_Citation/date/smXML:CI_Date/date'
            self.source['SURATH'] = './/identificationInfo/smXML:BAG_DataIdentification/citation/smXML:CI_Citation/citedResponsibleParty/smXML:CI_ResponsibleParty/organisationName'
            self.source['SUREND'] = './/identificationInfo/smXML:BAG_DataIdentification/citation/smXML:CI_Citation/date/smXML:CI_Date/date' 
            self.source['VERDAT'] = './/referenceSystemInfo/smXML:MD_CRS'
            self.source['HORDAT'] = './/referenceSystemInfo/smXML:MD_CRS'#/projection/smXML:RS_Identifier/code/ellipsoid/smXML:RS_Identifier/code/datum/smXML:RS_Identifier/code/projectionParameters/smXML:MD_ProjectionParameters/zone'
            self.source['resx_resy'] = './/spatialRepresentationInfo/smXML:MD_Georectified/axisDimensionProperties/smXML:MD_Dimension/resolution/smXML:Measure/smXML:value'
            
        else:
            log.warning("verison not compatible")
            
    def get_fields(self):
        """
        Using the field available for the version type, get the data for those
        fields.
        """
        self.data = {}
        if 'rows_cols' in self.source:
            self._read_rows_and_cols()
        if 'resx_resy' in self.source:
            self._read_res_x_and_y()
        if 'bbox' in self.source:
            self._read_corners_sw_and_ne()
        if 'wkt_srs' in self.source:
            self._read_wkt_prj()
        if 'lon_min' in self.source:
            self._read_bbox()
        if 'abstract' in self.source:
            self._read_abstract()
        if 'date' in self.source:
            self._read_date()
        if 'unc_type' in self.source:
            self._read_uncertainty_type()
        if 'z_min' in self.source:
            self._read_depth_min_max()
        if 'sourcename' in self.source:
            self._read_source_name()
        if 'filename' in self.source:
            self.data['filename'] = self.source['filename']
        if 'SORDAT' in self.source:
            self._read_SORDAT()
        if 'SURATH' in self.source:
            self._read_survey_authority()
        if 'SURSTA' in self.source:
            self._read_survey_start_date()
        if 'SUREND' in self.source:
            self._read_survey_end_date()
        if 'VERDAT' in self.source:
            self._read_vertical_datum()
        if 'HORDAT' in self.source:
            self._read_horizontal_datum()
        if 'planam' in self.source:
            self._read_platform_name()
        if 'sensor' in self.source:
            self._read_sensor_types()
            
    def get_s57_dict(self):
        """
        Convert the object dictionary 'data' keys to the desired S57 keys and
        return the dicitonary.
        """
        s57 = {}
        for key in self.data.keys():
            if key == 'z_min':
                s57['DRVAL1'] = str(round(self.data[key],1))
            elif key == 'z_max':
                s57['DRVAL2'] = str(round(self.data[key],1))
            elif key == 'filename':
                s57['OBJNAM'] = self.data['filename'][:-8] # this gets rid of _bag.xml
                label = self.data[key]
                p = re.compile('_MB_|_VB_')
                s = p.findall(label)
                if len(s) > 0:
                    if s[0] == '_MB_':
                        s57['TECSOU'] = '3'
                    elif s[0] == '_VB__':
                        s57['TECSOU'] = '1'
            elif key == 'abstract':
                s57['uniqid'] = self.data['abstract']
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
                    # print (self.data[key])
                    s57[key] = horz_datum['Local']
        # try to find a source for SORIND
        p = re.compile('[A-Z][0-9][0-9][0-9][0-9][0-9]')
        if 'filename' in self.data:
            s = p.findall(self.data['filename'])
            if len(s) == 1:
                s57['SORIND'] = 'US,US,graph,'+ s[0]
        elif 'abstract' in self.data and len(s) == 0:
            s = p.findall(self.data['abstract'])
            if len(s) == 1:
                s57['SORIND'] = 'US,US,graph,'+ s[0]
        elif 'sourcename' in self.data and len(s) == 0:
            s = p.findall(self.data['sourcename'])
            if len(s) == 1:
                s57['SORIND'] = 'US,US,graph,'+ s[0]
        return s57
    
    def _read_rows_and_cols(self):
        """ attempts to read rows and cols info """

        try:
            ret = self.xml_tree.findall(self.source['rows_cols'],
                                        namespaces=self.ns)
        except:
            log.warning("unable to read rows and cols")
            return

        try:
            self.data['rows'] = int(ret[0].text)
            self.data['cols'] = int(ret[1].text)

        except Exception as e:
            log.warning("unable to read rows and cols: %s" % e)
            return

    def _read_res_x_and_y(self):
        """ attempts to read resolution along x- and y- axes """

        try:
            ret = self.xml_tree.findall(self.source['resx_resy'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read res x and y")
            return

        try:
            self.data['res_x'] = float(ret[0].text)
            self.data['res_y'] = float(ret[1].text)

        except Exception as e:
            log.warning("unable to read res x and y: %s" % e)
            return

    def _read_corners_sw_and_ne(self):
        """ attempts to read corners SW and NE """

        try:
            ret = self.xml_tree.find(self.source['bbox'],
                                      namespaces=self.ns).text.split()
        except:
            log.warning("unable to read corners SW and NE")
            return

        try:
            self.data['southwest_corner'] = [float(c) for c in ret[0].split(',')]
            self.data['northeast_corner'] = [float(c) for c in ret[1].split(',')]

        except Exception as e:
            log.warning("unable to read corners SW and NE: %s" % e)
            return

    def _read_wkt_prj(self):
        """ attempts to read the WKT projection string """

        try:
            ret = self.xml_tree.find(self.source['wkt_srs'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the WKT projection string")
            return

        try:
            self.data['wkt_srs'] = ret.text

        except Exception as e:
            log.warning("unable to read the WKT projection string: %s" % e)
            return

    def _read_bbox(self):
        """ attempts to read the bounding box values """

        try:
            ret_x_min = self.xml_tree.find(self.source['lon_min'],
                                            namespaces=self.ns)
            ret_x_max = self.xml_tree.find(self.source['lon_max'],
                                            namespaces=self.ns)
        except:
            log.warning("unable to read the bbox's longitude values")
            return

        try:
            self.data['lon_min'] = float(ret_x_min.text)
            self.data['lon_max'] = float(ret_x_max.text)
        except Exception as e:
            log.warning("unable to read the bbox's longitude values: %s" % e)
            return

        try:
            ret_y_min = self.xml_tree.find(self.source['lat_min'],
                                            namespaces=self.ns)
            ret_y_max = self.xml_tree.find(self.source['lat_max'],
                                            namespaces=self.ns)
        except:
            log.warning("unable to read the bbox's latitude values")
            return

        try:
            self.data['lat_min'] = float(ret_y_min.text)
            self.data['lat_max'] = float(ret_y_max.text)
        except Exception as e:
            log.warning("unable to read the bbox's latitude values: %s" % e)
            return

    def _read_abstract(self):
        """ attempts to read the abstract string """

        try:
            ret = self.xml_tree.find(self.source['abstract'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the abstract string")
            return

        try:
            self.data['abstract'] = ret.text
        except Exception as e:
            log.warning("unable to read the abstract string: %s" % e)
            return

    def _read_date(self):
        """ attempts to read the date string """

        try:
            ret = self.xml_tree.find(self.source['date'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the date string")
            return

        try:
            text_date = ret.text
        except Exception as e:
            log.warning("unable to read the date string: %s" % e)
            return

        tm_date = None
        try:
            parsed_date = parser.parse(text_date)
            tm_date = parsed_date.strftime('%Y-%m-%d')
        except Exception:
            pass
            # log.warning("unable to handle the date string: %s" % text_date)

        if tm_date is None:
            self.data['date'] = text_date
        else:
            self.data['date'] = tm_date

    def _read_uncertainty_type(self):
        """ attempts to read the uncertainty type """

        try:
            ret = self.xml_tree.find(self.source['unc_type'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the uncertainty type string")
            return

        try:
            self.data['unc_type'] = ret.text
        except Exception as e:
            log.warning("unable to read the uncertainty type attribute: %s" % e)
            return
        
    def _read_depth_min_max(self):
        """
        Read the depth min and max values and store them in the object 'data'
        dictionary.  It should be noted that these are the min and max values 
        for the elevation cells.  If this is a varible resolution surface than 
        these may not be accurate.
        """

        try:
            ret_z_min = self.xml_tree.find(self.source['z_min'],
                                            namespaces=self.ns)
            ret_z_max = self.xml_tree.find(self.source['z_max'],
                                            namespaces=self.ns)
        except:
            log.warning("unable to read the depth min and max values")
            return

        try:
            if ret_z_min is not None:
                self.data['z_min'] = float(ret_z_min.text)
            if ret_z_max is not None:
                self.data['z_max'] = float(ret_z_max.text)
        except Exception as e:
            log.warning("unable to read the depth min and max values: %s" % e)
            return
        
    def _read_source_name(self):
        """ 
        Read the source file name and store it in the object 'data' dictionary
        with the key 'filename'.
        """

        try:
            ret = self.xml_tree.find(self.source['sourcename'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the survey name string")
            return

        try:
            self.data['sourcename'] = ret.text
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
            tm_date = parsed_date.strftime('%Y-%m-%d')
        except Exception:
            pass
            #log.warning("unable to handle the date string: %s" % text_date)

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
                pass
                #log.warning("unable to handle the survey start string: %s" % text_start_date)
    
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
                pass
                #log.warning("unable to handle the survey end string: %s" % text_end_date)
    
            if tme_date is None:
                self.data['SUREND'] = text_end_date
            else:
                self.data['SUREND'] = tme_date
            
    def _read_vertical_datum(self):
        """ 
        Read the survey vertical datum and store it in the object 'data'
        dictionary with the key 'VERDAT'.
        """

        try:
            ret = self.xml_tree.findall(self.source['VERDAT'],
                                      namespaces=self.ns)
            for r in ret:
                # first look to see if this is the wrong branch of < v1.5
                val = r.find('projection')
                if val is not None:
                    continue
                else:
                    # check if this is the other (right) branch of < v1.5
                    val = r.find('datum')
                    if val is not None:
                        val = r.find('datum/smXML:RS_Identifier/code',self.ns).text
                    else:
                        # by default, we know we are in >v1.5
                        datum_str = r.text
                        if datum_str[:4] == 'VERT':
                            vals = datum_str.split('"')
                            val = vals[1]
                        else:
                            continue
                        
        except:
            log.warning("unable to read the survey vertical datum name string")
            return

        try:
            self.data['VERDAT'] = val
        except Exception as e:
            log.warning("unable to read the survey vertical datum name attribute: %s" % e)
            return
        
    def _read_horizontal_datum(self):
        """ 
        Read the survey horizontal datum and store it in the object 'data'
        dictionary with the key 'HORDAT'.
        """

        try:
            ret = self.xml_tree.findall(self.source['HORDAT'],
                                      namespaces=self.ns)
            datum = None
            for r in ret:
                # first look to see if this is the right branch of < v1.5
                val = r.find('projection')
                if val is not None:
                    datum = r.find('datum/smXML:RS_Identifier/code',self.ns).text
                else:
                    # check if this is the wrong branch of < v1.5
                    val = r.find('datum')
                    if val is not None:
                        continue
                    else:
                        # by default, we know we are in >v1.5
                        datum_str = r.text
                        if datum_str[:4] == 'PROJ':
                            datum_sections = datum_str.split(',')
                            for d in datum_sections:
                                loc = d.find('DATUM')
                                if loc > -1:
                                    datum = d.split('"')[1]
                                    #print (datum)
                        else:
                            continue
        except:
            log.warning("unable to read the survey horizontal datum name string")
            return

        try:
            if datum is not None:
                self.data['HORDAT'] = datum
        except Exception as e:
            log.warning("unable to read the survey horizontal datum name attribute: %s" % e)
            return
        
    def _read_platform_name(self):
        """ 
        Read the platform name and store it in the object 'data' dictionary
        with the key 'planam'.
        """

        try:
            ret = self.xml_tree.findall(self.source['planam'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the platform name string")
            return

        try:
            if len(ret) > 0:
                self.data['planam'] = []
                for r in ret:
                    self.data['planam'].append(r.text)
        except Exception as e:
            log.warning("unable to read the platform name attribute: %s" % e)
            return
        
    def _read_sensor_types(self):
        """ 
        Read the sensor used and store it in the object 'data' dictionary
        with the key 'sensor'.
        """

        try:
            ret = self.xml_tree.findall(self.source['sensor'],
                                      namespaces=self.ns)
        except:
            log.warning("unable to read the sensor name string")
            return

        try:
            if len(ret) > 0:
                self.data['sensor'] = []
                for r in ret:
                    self.data['sensor'].append(r.text)
        except Exception as e:
            log.warning("unable to read the sensor name attribute: %s" % e)
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
        if v[0] == ':':
            tmp = v[1:]
            tmp = tmp.split('>')[0]
            name, info = tmp.split('=')
            site = info.split('"')[1]
            namespace[name] = site
    return namespace
    