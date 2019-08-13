"""
parse_bag_meta.py

20180222 G.Rice

V0.0.5 Last Updated 20180403

This is a collection of tools for working with BAG metadata.  Created
specifically to extract the bag metadata for building a bathymetric database,
this collection of method attempt to both serve extraction of the meta data in
a general sense, and also for specific S57 needs.

"""

import logging as _logging
import os as _os
import re as _re
import sys as _sys
import tables as _tb
import datetime as _datetime
from glob import glob as _glob

from osgeo import gdal as _gdal
from osgeo import osr as _osr
from xml.etree import ElementTree as _et

_gdal.UseExceptions()

__version__ = 'BAG'

vert_datum = {
    'MLWS': '1',
    'MLLWS': '2',
    'MSL': '3',
    'LLW': '4',
    'MLW': '5',
    'ISLW': '8',
    'MLLW': '12',
    'MHW': '16',
    'MHWS': '17',
    'MHHW': '21',
    'LAT': '23',
    'LOC': '24',
    'IGLD': '25',
    'LLWLT': '27',
    'HHWLT': '28',
    'HAT': '30',
    'Unknown': '701',
    'Other': '703',
    'HRD': '24',  # Adding this for the Hudson River Datum
}

horz_datum = {
    'WGS72': '1',
    'WGS84': '2',
    'WGS_1984': '2',
    'NAD27': '74',
    'NAD83': '75',
    'North_American_Datum_1983': '75',
    'Local': '131',
}


class BAGRawReader:
    """
    Helper class to manage BAG xml metadata. This class takes an xml string and
    parses the string based on a dictionary (named 'source') with keys that
    name the item to extract and the tree path string to find the data.
    Different versions of the source dictionary can be set based on the version
    of the BAG being parsed.  All extracted data is placed in a dictionary
    (named 'data').
    """

    def __init__(self):
        """
        Provided a BAG xml string for parsing, a tree will be created and the
        name speace parsed from the second line.  Values are then extracted
        based on the source dictionary.  If this dictionary
        """
        self.data = {}

        self._logger = _logging.getLogger(f'fuse')

        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            ch.setLevel(_logging.DEBUG)
            self._logger.addHandler(ch)

    def read_metadata(self, infilename: str, version=None) -> dict:
        """

        Parameters
        ----------
        infilename : str

        Returns
        -------
        dict

        """
        try:
            meta_gdal, version = self._parse_bag_gdal(infilename)
            meta_xml = self._parse_bag_xml(infilename, version=version)
            meta_support = self._known_meta(infilename)
            return {**meta_support, **meta_xml, **meta_gdal}
        except ValueError as e:
            print(f'{e}')


    def _parse_bag_gdal(self, infilename: str) -> dict:
        """
        Parses metadata available via gdal and returns them as a dict

        Parameters
        ----------
        infilename : str
            Input file path

        Returns
        -------
        dict
            A dictionary object containing found metadata

        """

        metadata = {}
        try:
            bag_file = _gdal.Open(infilename)
            metadata = {**bag_file.GetMetadata()}
            version = float(metadata['BagVersion'][:-2])
            metadata['rows'], metadata['cols'] = bag_file.RasterYSize, bag_file.RasterXSize
            geotransform = bag_file.GetGeoTransform()
            metadata['wkt_srs'] = bag_file.GetProjectionRef()
            spacial_ref = _osr.SpatialReference(wkt=metadata['wkt_srs'])
            if spacial_ref.IsProjected:
                metadata['horiz_datum'] = spacial_ref.GetAttrValue('projcs')
            metadata['horiz_frame'] = spacial_ref.GetAttrValue('geogcs')
#            metadata['res_x'], metadata['res_y'] = geotransform[1], geotransform[5]
        except RuntimeError as e:
            raise ValueError(f'{e}')
        return metadata, version

    def _parse_bag_xml(self, infilename: str, version=None) -> dict:
        """
        Parses metadata available via the file's xml and returns them as a dict

        Parameters
        ----------
        infilename : str
            Input file path

        Returns
        -------
        dict
            A dictionary object containing found metadata

        """
        try:
            with _tb.open_file(infilename, mode='r') as bag_file:
                meta_read = [str(x, 'utf-8', 'ignore') for x in bag_file.root.BAG_root.metadata.read()]
                meta_xml = ''.join(meta_read)

                encode_val = 0

                for x in meta_xml:
                    if meta_xml[encode_val] == '>':
                        meta_xml = meta_xml[encode_val:]
                        break
                    else:
                        encode_val += 1
                start_val = 0

                for x in meta_xml:
                    if meta_xml[start_val] == '<':
                        meta_xml = meta_xml[start_val:]
                        break
                    else:
                        start_val += 1

            self.xml_tree = _et.XML(meta_xml)
            self.namespace = self._assign_namspace(version=version)
            self.bag_format = self._set_format(infilename, version)
            self.get_fields()

            return self.data
        except _tb.HDF5ExtError as e:
            raise ValueError(f'{e}')

    def _known_meta(self, infilename: str) -> dict:
        """
        Identifies known metadata and returns them as a dict

        Parameters
        ----------
        infilename : str
            Input file path

        Returns
        -------
        dict
            A dictionary object containing found metadata

        """
        meta = {}
        coverage = self._find_coverage(infilename)
        root, name = _os.path.split(infilename)
        meta['from_filename'] = name
        meta['from_path'] = infilename

        return {**coverage, **meta}

    def _find_coverage(self, infilename: str) -> dict:
        """
        Identifies supporting coverage files and returns them as a dict

        Parameters
        ----------
        infilename : str
            Input file path

        Returns
        -------
        dict
            A dictionary object containing found metadata

        """
        meta = {}
        root, filename = _os.path.split(infilename)
        dir_files = [name for name in _os.listdir(root) if _os.path.isfile(_os.path.join(root, name))]
        meta['support_files'] = [support_file for support_file in dir_files if _os.path.splitext(support_file)[1].lower() in ('.tiff', '.tif', '.tfw', '.gpkg')]
        if len(meta['support_files']) > 0:
            meta['interpolate'] = True
        return meta

    def _assign_namspace(self, version=None, xml=None):
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
            'smXML': "http://metadata.dgiwg.org/smXML",
            "xlink": "http://www.w3.org/1999/xlink",
            'xsi': "http://www.w3.org/2001/XMLSchema-instance",
            'gml': "http://www.opengis.net/gml"
        }

        if version is not None:
            if version >= 1.5:
                namespace = post15
            elif version < 1.5:
                namespace = pre15
            return namespace
        elif version is None and xml != None:
            return parse_namespace(xml)



    def _set_format(self, infilename: str, version: float) -> object:
        """
        Set the locations of the desired data types based on the version of the
        bag.

        Parameters
        ----------
        infilename : str
        version : float
        """

        source = {'filename': infilename}

        if version >= 1.5:
            # these are the bag types
            source[
                'rows_cols'] = './/gmd:spatialRepresentationInfo/gmd:MD_Georectified/gmd:axisDimensionProperties/gmd:MD_Dimension/gmd:dimensionSize/gco:Integer'
            source[
                'resx_resy'] = './/gmd:spatialRepresentationInfo/gmd:MD_Georectified/gmd:axisDimensionProperties/gmd:MD_Dimension/gmd:resolution/gco:Measure'
            source[
                'bbox'] = './/gmd:spatialRepresentationInfo/gmd:MD_Georectified/gmd:cornerPoints/gml:Point/gml:coordinates'
            source[
                'wkt_srs'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString'
            source['lon_min'] = './/gmd:EX_GeographicBoundingBox/gmd:westBoundLongitude/gco:Decimal'
            source['lon_max'] = './/gmd:EX_GeographicBoundingBox/gmd:eastBoundLongitude/gco:Decimal'
            source['lat_min'] = './/gmd:EX_GeographicBoundingBox/gmd:southBoundLatitude/gco:Decimal'
            source['lat_max'] = './/gmd:EX_GeographicBoundingBox/gmd:northBoundLatitude/gco:Decimal'
            source['abstract'] = './/gmd:abstract/gco:CharacterString'
            source['date'] = './/gmd:CI_Date/gmd:date/gco:Date'
            source['unc_type'] = './/bag:verticalUncertaintyType/bag:BAG_VertUncertCode'
            source[
                'z_min'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:verticalElement/gmd:EX_VerticalExtent/gmd:minimumValue/gco:Real'
            source[
                'z_max'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:verticalElement/gmd:EX_VerticalExtent/gmd:maximumValue/gco:Real'
            source[
                'sourcename'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:title/gco:CharacterString'
            source[
                'SORDAT'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date'
            source['SURATH'] = './/gmd:contact/gmd:CI_ResponsibleParty/gmd:organisationName/gco:CharacterString'
            source[
                'SUREND'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:citation/gmd:CI_Citation/gmd:date/gmd:CI_Date/gmd:date/gco:Date'
            source[
                'SURSTA'] = './/gmd:identificationInfo/bag:BAG_DataIdentification/gmd:extent/gmd:EX_Extent/gmd:temporalElement/gmd:EX_TemporalExtent/gmd:extent/gml:TimePeriod/gml:beginPosition'
            source[
                'VERDAT'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString'
            source[
                'HORDAT'] = './/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString'
            source[
                'planam'] = './/gmi:acquisitionInformation/gmi:MI_AcquisitionInformation/gmi:platform/gmi:MI_Platform/gmi:identifier/gmd:MD_Identifier/gmd:code/gco:CharacterString'
            source[
                'sensor'] = './/gmi:instrument/gmi:MI_Instrument/gmi:identifier/gmd:MD_Identifier/gmd:code/gco:CharacterString'

        elif version < 1.5 and version > 0:
            source['abstract'] = './/identificationInfo/smXML:BAG_DataIdentification/abstract'
            source[
                'SORDAT'] = './/identificationInfo/smXML:BAG_DataIdentification/citation/smXML:CI_Citation/date/smXML:CI_Date/date'
            source[
                'SURATH'] = './/identificationInfo/smXML:BAG_DataIdentification/citation/smXML:CI_Citation/citedResponsibleParty/smXML:CI_ResponsibleParty/organisationName'
            source[
                'SUREND'] = './/identificationInfo/smXML:BAG_DataIdentification/citation/smXML:CI_Citation/date/smXML:CI_Date/date'
            source['VERDAT'] = './/referenceSystemInfo/smXML:MD_CRS'
            source[
                'HORDAT'] = './/referenceSystemInfo/smXML:MD_CRS'  # /projection/smXML:RS_Identifier/code/ellipsoid/smXML:RS_Identifier/code/datum/smXML:RS_Identifier/code/projectionParameters/smXML:MD_ProjectionParameters/zone'
            source[
                'resx_resy'] = './/spatialRepresentationInfo/smXML:MD_Georectified/axisDimensionProperties/smXML:MD_Dimension/resolution/smXML:Measure/smXML:value'

        else:
            _logging.warning("verison not compatible")

        return source

    def get_fields(self):
        """
        Using the field available for the version type, get the data for those
        fields.
        """
        self.data = {}
        if 'rows_cols' in self.bag_format:
            self._read_rows_and_cols()
        if 'resx_resy' in self.bag_format:
            self._read_res_x_and_y()
        if 'bbox' in self.bag_format:
            self._read_corners_sw_and_ne()
        if 'wkt_srs' in self.bag_format:
            self._read_wkt_prj()
        if 'lon_min' in self.bag_format:
            self._read_bbox()
        if 'abstract' in self.bag_format:
            self._read_abstract()
        if 'date' in self.bag_format:
            self._read_date()
        if 'unc_type' in self.bag_format:
            self._read_uncertainty_type()
        if 'z_min' in self.bag_format:
            self._read_depth_min_max()
        if 'sourcename' in self.bag_format:
            self._read_source_name()
        if 'filename' in self.bag_format:
            self.data['filename'] = self.bag_format['filename']
        if 'SORDAT' in self.bag_format:
            self._read_SORDAT()
        if 'SURATH' in self.bag_format:
            self._read_survey_authority()
        if 'SURSTA' in self.bag_format:
            self._read_survey_start_date()
        if 'SUREND' in self.bag_format:
            self._read_survey_end_date()
        if 'VERDAT' in self.bag_format:
            self._read_vertical_datum()
        if 'HORDAT' in self.bag_format:
            self._read_horizontal_datum()
        if 'planam' in self.bag_format:
            self._read_platform_name()
        if 'sensor' in self.bag_format:
            self._read_sensor_types()

    def get_s57_dict(self):
        """
        Convert the object dictionary 'data' keys to the desired S57 keys and
        return the dicitonary.
        """
        s57 = {}
        for key in self.data.keys():
            if key == 'z_min':
                s57['DRVAL1'] = str(round(self.data[key], 1))
            elif key == 'z_max':
                s57['DRVAL2'] = str(round(self.data[key], 1))
            elif key == 'filename':
                s57['OBJNAM'] = self.data['filename'][:-8]  # this gets rid of _bag.xml
                label = self.data[key]
                p = _re.compile('_MB_|_VB_')
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
                newds = datestring.replace('-', '')
                no_time_ds = newds.split('T')[0]
                s57[key] = no_time_ds
            elif key == 'SURATH':
                s57[key] = self.data[key]
            elif key == 'SURSTA':
                datestring = self.data[key]
                newds = datestring.replace('-', '')
                no_time_ds = newds.split('T')[0]
                s57[key] = no_time_ds
            elif key == 'SUREND':
                datestring = self.data[key]
                newds = datestring.replace('-', '')
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
        p = _re.compile(r'[A-Z]([0-9]{5})')
        if 'filename' in self.data:
            s = p.findall(self.data['filename'])
            if len(s) == 1:
                s57['SORIND'] = 'US,US,graph,' + s[0]
        elif 'abstract' in self.data and len(s) == 0:
            s = p.findall(self.data['abstract'])
            if len(s) == 1:
                s57['SORIND'] = 'US,US,graph,' + s[0]
        elif 'sourcename' in self.data and len(s) == 0:
            s = p.findall(self.data['sourcename'])
            if len(s) == 1:
                s57['SORIND'] = 'US,US,graph,' + s[0]
        return s57

    def _read_rows_and_cols(self):
        """ attempts to read rows and cols info """

        try:
            ret = self.xml_tree.findall(self.bag_format['rows_cols'],
                                        namespaces=self.namespace)
        except:
            _logging.warning("unable to read rows and cols")
            return

        try:
            self.data['rows'] = int(ret[0].text)
            self.data['cols'] = int(ret[1].text)

        except Exception as e:
            _logging.warning("unable to read rows and cols: %s" % e)
            return

    def _read_res_x_and_y(self):
        """ attempts to read resolution along x- and y- axes """

        try:
            ret = self.xml_tree.findall(self.bag_format['resx_resy'],
                                        namespaces=self.namespace)
        except:
            _logging.warning("unable to read res x and y")
            return

        try:
            self.data['res_x'] = float(ret[0].text)
            self.data['res_y'] = float(ret[1].text)

        except Exception as e:
            _logging.warning("unable to read res x and y: %s" % e)
            return

    def _read_corners_sw_and_ne(self):
        """ attempts to read corners SW and NE """

        try:
            ret = self.xml_tree.find(self.bag_format['bbox'],
                                     namespaces=self.namespace).text.split()
        except:
            _logging.warning("unable to read corners SW and NE")
            return

        try:
            self.data['soutwest_corner'] = [float(c) for c in ret[0].split(',')]
            self.data['northeast_corner'] = [float(c) for c in ret[1].split(',')]

        except Exception as e:
            _logging.warning("unable to read corners SW and NE: %s" % e)
            return

    def _read_wkt_prj(self):
        """ attempts to read the WKT projection string """

        try:
            ret = self.xml_tree.find(self.bag_format['wkt_srs'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the WKT projection string")
            return

        try:
            self.data['wkt_srs'] = ret.text

        except Exception as e:
            _logging.warning("unable to read the WKT projection string: %s" % e)
            return

    def _read_bbox(self):
        """ attempts to read the bounding box values """

        try:
            ret_x_min = self.xml_tree.find(self.bag_format['lon_min'],
                                           namespaces=self.namespace)
            ret_x_max = self.xml_tree.find(self.bag_format['lon_max'],
                                           namespaces=self.namespace)
        except:
            _logging.warning("unable to read the bbox's longitude values")
            return

        try:
            self.data['lon_min'] = float(ret_x_min.text)
            self.data['lon_max'] = float(ret_x_max.text)
        except Exception as e:
            _logging.warning("unable to read the bbox's longitude values: %s" % e)
            return

        try:
            ret_y_min = self.xml_tree.find(self.bag_format['lat_min'],
                                           namespaces=self.namespace)
            ret_y_max = self.xml_tree.find(self.bag_format['lat_max'],
                                           namespaces=self.namespace)
        except:
            _logging.warning("unable to read the bbox's latitude values")
            return

        try:
            self.data['lat_min'] = float(ret_y_min.text)
            self.data['lat_max'] = float(ret_y_max.text)
        except Exception as e:
            _logging.warning("unable to read the bbox's latitude values: %s" % e)
            return

    def _read_abstract(self):
        """ attempts to read the abstract string """

        try:
            ret = self.xml_tree.find(self.bag_format['abstract'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the abstract string")
            return

        try:
            self.data['abstract'] = ret.text
        except Exception as e:
            _logging.warning("unable to read the abstract string: %s" % e)
            return

    def _read_date(self):
        """ attempts to read the date string """

        try:
            ret = self.xml_tree.find(self.bag_format['date'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the date string")
            return

        try:
            text_date = ret.text
        except Exception as e:
            _logging.warning("unable to read the date string: %s" % e)
            return

        tm_date = None
        try:
            parsed_date = parser.parse(text_date)
            tm_date = parsed_date.strftime('%Y-%m-%d')
        except Exception:
            pass
            # _logging.warning("unable to handle the date string: %s" % text_date)

        if tm_date is None:
            self.data['date'] = text_date
        else:
            self.data['date'] = tm_date

    def _read_uncertainty_type(self):
        """ attempts to read the uncertainty type """

        try:
            ret = self.xml_tree.find(self.bag_format['unc_type'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the uncertainty type string")
            return

        try:
            self.data['unc_type'] = ret.text
        except Exception as e:
            _logging.warning("unable to read the uncertainty type attribute: %s" % e)
            return

    def _read_depth_min_max(self):
        """
        Read the depth min and max values and store them in the object 'data'
        dictionary.  It should be noted that these are the min and max values
        for the elevation cells.  If this is a varible resolution surface than
        these may not be accurate.
        """

        try:
            ret_z_min = self.xml_tree.find(self.bag_format['z_min'],
                                           namespaces=self.namespace)
            ret_z_max = self.xml_tree.find(self.bag_format['z_max'],
                                           namespaces=self.namespace)
        except:
            _logging.warning("unable to read the depth min and max values")
            return

        try:
            if ret_z_min is not None:
                self.data['z_min'] = float(ret_z_min.text)
            if ret_z_max is not None:
                self.data['z_max'] = float(ret_z_max.text)
        except Exception as e:
            _logging.warning("unable to read the depth min and max values: %s" % e)
            return

    def _read_source_name(self):
        """
        Read the source file name and store it in the object 'data' dictionary
        with the key 'filename'.
        """

        try:
            ret = self.xml_tree.find(self.bag_format['sourcename'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the survey name string")
            return

        try:
            self.data['sourcename'] = ret.text
        except Exception as e:
            _logging.warning("unable to read the survey name attribute: %s" % e)
            return

    def _read_SORDAT(self):
        """
        Reads a date, but what it means exactly needs to be researched...
        """

        try:
            ret = self.xml_tree.find(self.bag_format['SORDAT'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the SORDAT date string")
            return

        try:
            text_date = ret.text
        except Exception as e:
            _logging.warning("unable to read the SORDAT date string: %s" % e)
            return

        tm_date = None
        try:
            parsed_date = parser.parse(text_date)
            tm_date = parsed_date.strftime('%Y-%m-%d')
        except Exception:
            pass
            # _logging.warning("unable to handle the date string: %s" % text_date)

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
            ret = self.xml_tree.find(self.bag_format['SURATH'],
                                     namespaces=self.namespace)
        except:
            _logging.warning("unable to read the survey authority name string")
            return

        try:
            self.data['SURATH'] = ret.text
        except Exception as e:
            _logging.warning("unable to read the survey authority name attribute: %s" % e)
            return

    def _read_survey_start_date(self):
        """
        Read the survey start date store it in the object 'data'
        dictionary with the key 'SURSTA'.
        """

        try:
            rets = self.xml_tree.find(self.bag_format['SURSTA'],
                                      namespaces=self.namespace)
        except:
            _logging.warning("unable to read the survey start date string")
            return

        if rets is not None:
            try:
                text_start_date = rets.text
            except Exception as e:
                _logging.warning("unable to read the survey start date string: %s" % e)
                return

            tms_date = None
            try:
                parsed_date = parser.parse(text_start_date)
                tms_date = parsed_date.strftime('%Y-%m-%dT%H:%M:%SZ')
            except Exception:
                pass
                # _logging.warning("unable to handle the survey start string: %s" % text_start_date)

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
            rete = self.xml_tree.find(self.bag_format['SUREND'],
                                      namespaces=self.namespace)
        except:
            _logging.warning("unable to read the survey end date string")
            return

        if rete is not None:
            try:
                text_end_date = rete.text
            except Exception as e:
                _logging.warning("unable to read the survey end date string: %s" % e)
                return

            tme_date = None
            try:
                parsed_date = parser.parse(text_end_date)
                tme_date = parsed_date.strftime('%Y-%m-%dT%H:%M:%SZ')
            except Exception:
                pass
                # _logging.warning("unable to handle the survey end string: %s" % text_end_date)

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
            ret = self.xml_tree.findall(self.bag_format['VERDAT'],
                                        namespaces=self.namespace)
            for r in ret:
                # first look to see if this is the wrong branch of < v1.5
                val = r.find('projection')
                if val is not None:
                    continue
                else:
                    # check if this is the other (right) branch of < v1.5
                    val = r.find('datum')
                    if val is not None:
                        val = r.find('datum/smXML:RS_Identifier/code', self.namespace).text
                    else:
                        # by default, we know we are in >v1.5
                        datum_str = r.text
                        if datum_str[:4] == 'VERT':
                            vals = datum_str.split('"')
                            val = vals[1]
                        else:
                            continue

        except:
            _logging.warning("unable to read the survey vertical datum name string")
            return

        try:
            self.data['VERDAT'] = val
        except Exception as e:
            _logging.warning("unable to read the survey vertical datum name attribute: %s" % e)
            return

    def _read_horizontal_datum(self):
        """
        Read the survey horizontal datum and store it in the object 'data'
        dictionary with the key 'HORDAT'.
        """

        try:
            ret = self.xml_tree.findall(self.bag_format['HORDAT'],
                                        namespaces=self.namespace)
            datum = None
            for r in ret:
                # first look to see if this is the right branch of < v1.5
                val = r.find('projection')
                if val is not None:
                    datum = r.find('datum/smXML:RS_Identifier/code', self.namespace).text
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
                                    # print (datum)
                        else:
                            continue
        except:
            _logging.warning("unable to read the survey horizontal datum name string")
            return

        try:
            if datum is not None:
                self.data['HORDAT'] = datum
        except Exception as e:
            _logging.warning("unable to read the survey horizontal datum name attribute: %s" % e)
            return

    def _read_platform_name(self):
        """
        Read the platform name and store it in the object 'data' dictionary
        with the key 'planam'.
        """

        try:
            ret = self.xml_tree.findall(self.bag_format['planam'],
                                        namespaces=self.namespace)
        except:
            _logging.warning("unable to read the platform name string")
            return

        try:
            if len(ret) > 0:
                self.data['planam'] = []
                for r in ret:
                    self.data['planam'].append(r.text)
        except Exception as e:
            _logging.warning("unable to read the platform name attribute: %s" % e)
            return

    def _read_sensor_types(self):
        """
        Read the sensor used and store it in the object 'data' dictionary
        with the key 'sensor'.
        """

        try:
            ret = self.xml_tree.findall(self.bag_format['sensor'],
                                        namespaces=self.namespace)
        except:
            _logging.warning("unable to read the sensor name string")
            return

        try:
            if len(ret) > 0:
                self.data['sensor'] = []
                for r in ret:
                    self.data['sensor'].append(r.text)
        except Exception as e:
            _logging.warning("unable to read the sensor name attribute: %s" % e)
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
