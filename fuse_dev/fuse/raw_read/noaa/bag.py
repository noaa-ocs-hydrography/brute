"""
parse_bag_meta.py

20180222 G.Rice

V0.0.5 Last Updated 20180403

This is a collection of tools for working with BAG metadata.  Created
specifically to extract the bag metadata for building a bathymetric database,
this collection of method attempt to both serve extraction of the meta data in
a general sense, and also for specific S57 needs.

"""

import datetime as _datetime
import logging as _logging
import os as _os
import re as _re
import sys as _sys
from glob import glob as _glob

import lxml.etree as _et
import numpy as _np
import rasterio as _rio
import tables as _tb
from fuse.raw_read import parse_file_pickle
from fuse.raw_read.raw_read import RawReader
from numba import jit
from osgeo import gdal as _gdal
from osgeo import osr as _osr
from scipy.ndimage import binary_closing

try:
    import dateutil.parser as _parser
except ModuleNotFoundError as e:
    _logging.error(f'{e.__class__.__name__} - {e}')

_gdal.UseExceptions()

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

_ns = {
    'bag': 'http://www.opennavsurf.org/schema/bag',
    'gco': 'http://www.isotc211.org/2005/gco',
    'gmd': 'http://www.isotc211.org/2005/gmd',
    'gmi': 'http://www.isotc211.org/2005/gmi',
    'gml': 'http://www.opengis.net/gml/3.2',
    'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
}

_ns2 = {
    'gml': 'http://www.opengis.net/gml',
    'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
    'smXML': 'http://metadata.dgiwg.org/smXML',
}


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


class BAGRawReader(RawReader):
    """
    Helper class to manage BAG xml metadata. This class takes an xml string and
    parses the string based on a dictionary (named 'source') with keys that
    name the item to extract and the tree path string to find the data.
    Different versions of the source dictionary can be set based on the version
    of the BAG being parsed.  All extracted data is placed in a dictionary
    (named 'data').
    """

    def __init__(self, out_file_location):
        """
        Provided a BAG xml string for parsing, a tree will be created and the
        name speace parsed from the second line.  Values are then extracted
        based on the source dictionary.  If this dictionary
        """
        self.data = {}
        self.out_file_location = out_file_location
        self._logger = _logging.getLogger('fuse.read')

        if len(self._logger.handlers) == 0:
            ch = _logging.StreamHandler(_sys.stdout)
            ch.setLevel(_logging.DEBUG)
            self._logger.addHandler(ch)

    def read_metadata(self, filename: str) -> dict:
        """
        read metadata from file

        Parameters
        ----------
        filename
            file to read

        Returns
        -------
            dictionary of metadata
        """

        try:
            meta_gdal, bag_version = self._parse_bag_gdal(filename)
            meta_xml = self._parse_bag_xml(filename, bag_version=bag_version)
            meta_support = self._known_meta(filename)
            meta_pickle = parse_file_pickle.read_pickle(filename)
            meta_combined = {**meta_pickle, **meta_xml, **meta_gdal, **meta_support}
            meta_final = self._finalize_meta(meta_combined)
            return meta_final
        except ValueError as error:
            _logging.warning(f'Error in reading metadata for {filename}: {error}')
            return {}

    def read_bathymetry(self, infilename: str, out_verdat: str = None) -> _gdal.Dataset:
        """
        Returns a BagFile data object

        Parameters
        ----------
        infilename : str
            BAG filepath

        Returns
        -------
            BagFile object
        """

        bag_file = Open(infilename)
        dataset = BagToGDALConverter(out_verdat)
        dataset.bag2gdal(bag_file)
        del bag_file

        return dataset.dataset

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
            bag_version = float(metadata['BagVersion'][:-2])
            metadata['shape'] = (bag_file.RasterYSize, bag_file.RasterXSize)
            #            geotransform = bag_file.GetGeoTransform()
            #            metadata['res'] = (geotransform[1], geotransform[5])
            metadata['from_horiz_datum'] = bag_file.GetProjectionRef()
            spacial_ref = _osr.SpatialReference(wkt=metadata['from_horiz_datum'])

            if spacial_ref.GetAttrValue('projection') == 'Transverse_Mercator':
                metadata['from_horiz_type'] = 'UTM'

            metadata['from_horiz_key'] = spacial_ref.GetUTMZone()
            datum = spacial_ref.GetAttrValue('datum')

            if datum.lower() in ('north_american_datum_1983', 'north american datum 1983', 'nad83'):
                metadata['from_horiz_frame'] = 'NAD83'
            elif metadata['from_horiz_datum'].lower() in ('wgs_1984', 'wgs84'):
                metadata['from_horiz_frame'] = 'WGS84'

            if spacial_ref.GetAttrValue('unit').lower() in ('meter', 'meters', 'metre', 'm'):
                metadata['from_horiz_units'] = 'm'
            elif spacial_ref.GetAttrValue('unit').lower() in ('feet', 'ft'):
                metadata['from_horiz_units'] = 'f'

            # metadata['elevation'] = self._read_gdal_array(bag_file, 1)
            # metadata['uncertainty'] = self._read_gdal_array(bag_file, 2)
            # metadata['nodata'] = 1000000.0

            del bag_file
        except RuntimeError as e:
            raise ValueError(f'{e}')
        return metadata, bag_version

    def _parse_bag_xml(self, infilename: str, bag_version=None) -> dict:
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
            self.namespace = self._assign_namspace(bag_version=bag_version)
            self.bag_format = self._set_format(infilename, bag_version)
            self.get_fields()
            if 'from_vert_datum' not in self.data or 'from_vert_key' not in self.data:
                for datum in vert_datum.keys():
                    if datum in infilename:
                        self.data['from_vert_datum'], self.data['from_vert_key'] = datum, datum
                        break

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
        basename = _os.path.splitext(_os.path.basename(infilename))[0]
        meta['from_filename'] = basename
        meta['from_path'] = infilename
        meta['file_size'] = self._size_finder(infilename)

        return {**coverage, **meta}

    def _size_finder(self, filepath: str) -> int:
        """
        Returns the rounded size of a file in MB as an integer

        Parameters
        ----------
        filepath : str, os.Pathlike

        Returns
        -------
        int

        """

        return int(_np.round(_os.path.getsize(filepath) / 1000))

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
        snum_reg = _re.compile(r'[A-Z]([0-9]{5})')
        snum = snum_reg.search(filename).group()
        dir_files = [_os.path.join(root, name) for name in _os.listdir(root) if
                     (_os.path.isfile(_os.path.join(root, name)) and snum in name)]
        meta['support_files'] = [support_file for support_file in dir_files if
                                 _os.path.splitext(support_file)[1].lower() in ('.tiff', '.tif', '.gpkg')]
        exts = [_os.path.splitext(support_file)[1].lower() for support_file in dir_files if
                _os.path.splitext(support_file)[1].lower() in ('.tiff', '.tif', '.gpkg')]
        return meta

    def _finalize_meta(self, meta):
        if ('from_vert_datum' not in meta or meta['from_vert_datum'] == '') and 'from_vert_key' not in meta:
            for datum in vert_datum.keys():
                if datum in meta['from_filename'] and datum == meta['from_filename'].split('_')[3]:
                    meta['from_vert_datum'], meta['from_vert_key'] = datum, datum
                    break
        elif 'from_vert_datum' in meta and 'from_vert_key' not in meta:
            for datum in vert_datum.keys():
                if datum in meta['from_filename'] and datum == meta['from_filename'].split('_')[3]:
                    meta['from_vert_key'] = datum
                    break
        # this should be moved to the reader.
        meta['read_type'] = 'bag'
        # translate from the reader to common metadata keys for datum transformations
        if 'from_vert_direction' not in meta:
            meta['from_vert_direction'] = 'height'
        if 'from_vert_units' not in meta:
            meta['from_vert_units'] = 'm'
        meta['interpolated'] = 'False'
        meta['posted'] = False
        return meta

    def _assign_namspace(self, bag_version=None, xml=None):
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

        if bag_version is not None:
            if bag_version >= 1.5:
                namespace = post15
            elif bag_version < 1.5:
                namespace = pre15
            return namespace
        elif bag_version is None and xml != None:
            return parse_namespace(xml)

    def _set_format(self, infilename: str, bag_version: float) -> object:
        """
        Set the locations of the desired data types based on the version of the
        bag.

        Parameters
        ----------
        infilename : str
        bag_version : float
        """

        source = {'filename': infilename}

        if bag_version >= 1.5:
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
            source['date_stamp'] = './/gmd:dateStamp/gco:Date'
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

        elif bag_version < 1.5 and bag_version > 0:
            source['abstract'] = './/identificationInfo/smXML:BAG_DataIdentification/abstract'
            source['date_stamp'] = './/dateStamp'
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
        if 'date_stamp' in self.bag_format:
            self._read_date_stamp()
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
            elif key == 'source_date':
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

    def _read_gdal_array(self, bag_obj, band: int) -> _np.array:
        """
        Returns a numpy.array of a GDAL object

        This function uses a :obj:`gdal.Dataset` object and band number to pull
        and read the appropriate array from the object

        Parameters
        ----------
        bag_obj : gdal.Dataset
            TODO write description
        band : int
            raster band number

        Returns
        -------
        type
            numpy.array

        """

        return bag_obj.GetRasterBand(band).ReadAsArray()

    def _read_rows_and_cols(self):
        """ attempts to read rows and cols info """

        try:
            ret = self.xml_tree.findall(self.bag_format['rows_cols'],
                                        namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read res x and y: {e}")
            return

        try:
            self.data['shape'] = (int(ret[1].text), int(ret[0].text))

        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read rows and cols: {e}")
            return

    def _read_res_x_and_y(self):
        """ attempts to read resolution along x- and y- axes """

        try:
            ret = self.xml_tree.findall(self.bag_format['resx_resy'],
                                        namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read res x and y: {e}")
            return

        try:
            self.data['res'] = (float(ret[0].text), -float(ret[1].text))

        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read res x and y: {e}")
            return

    def _read_corners_sw_and_ne(self):
        """ attempts to read corners SW and NE """

        try:
            ret = self.xml_tree.find(self.bag_format['bbox'],
                                     namespaces=self.namespace).text.split()
        except _et.Error as e:
            _logging.warning(f"unable to read corners SW and NE: {e}")
            return

        try:
            self.data['bounds'] = ([float(c) for c in ret[0].split(',')], [float(c) for c in ret[1].split(',')])

        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read corners SW and NE: {e}")
            return

    def _read_wkt_prj(self):
        """ attempts to read the WKT projection string """

        try:
            ret = self.xml_tree.find(self.bag_format['wkt_srs'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the WKT projection string: {e}")
            return

        try:
            self.data['from_horiz_datum'] = ret.text

        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the WKT projection string: {e}")
            return

    def _read_bbox(self):
        """ attempts to read the bounding box values """

        try:
            ret_x_min = self.xml_tree.find(self.bag_format['lon_min'],
                                           namespaces=self.namespace)
            ret_x_max = self.xml_tree.find(self.bag_format['lon_max'],
                                           namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the bbox's longitude values: {e}")
            return

        try:
            self.data['lon_min'] = float(ret_x_min.text)
            self.data['lon_max'] = float(ret_x_max.text)
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the bbox's longitude values: {e}")
            return

        try:
            ret_y_min = self.xml_tree.find(self.bag_format['lat_min'],
                                           namespaces=self.namespace)
            ret_y_max = self.xml_tree.find(self.bag_format['lat_max'],
                                           namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the bbox's latitude values: {e}")
            return

        try:
            self.data['lat_min'] = float(ret_y_min.text)
            self.data['lat_max'] = float(ret_y_max.text)
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the bbox's latitude values: {e}")
            return

    def _read_abstract(self):
        """ attempts to read the abstract string """

        try:
            ret = self.xml_tree.find(self.bag_format['abstract'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the abstract string: {e}")
            return

        try:
            self.data['abstract'] = ret.text
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the abstract string: {e}")
            return

    def _read_date(self):
        """ attempts to read the date string """

        try:
            ret = self.xml_tree.find(self.bag_format['date'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the date string: {e}")
            return

        if ret is not None:
            try:
                text_date = ret.text
            except (ValueError, IndexError, AttributeError) as e:
                _logging.warning(f"unable to read the date string: {e}")
                return

            tm_date = None
            try:
                parsed_date = _parser.parse(text_date)
                tm_date = parsed_date.strftime('%Y%m%d')
            except Exception:
                pass
                # _logging.warning("unable to handle the date string: %s" % text_date)

            if tm_date is None and text_date != '':
                self.data['date'] = text_date
            elif tm_date != '':
                self.data['date'] = tm_date
        else:
            _logging.warning(f"unable to read the `date` date string: {ret}")

    def _read_date_stamp(self):
        """ attempts to read the date_stamp string """

        try:
            ret = self.xml_tree.find(self.bag_format['date_stamp'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the date_stamp string: {e}")
            return

        if ret is not None:
            try:
                text_date = ret.text
            except (ValueError, IndexError, AttributeError) as e:
                _logging.warning(f"unable to read the date_stamp string: {e}")
                return

            tm_date = None
            try:
                parsed_date = _parser.parse(text_date)
                tm_date = parsed_date.strftime('%Y%m%d')
            except Exception:
                pass
                # _logging.warning("unable to handle the date string: %s" % text_date)

            if tm_date is None and text_date != '':
                self.data['date_stamp'] = text_date
            elif tm_date != '':
                self.data['date_stamp'] = tm_date
        else:
            _logging.warning(f"unable to read the `date_stamp` date string: {ret}")

    def _read_uncertainty_type(self):
        """ attempts to read the uncertainty type """

        try:
            ret = self.xml_tree.find(self.bag_format['unc_type'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the uncertainty type string: {e}")
            return

        try:
            self.data['unc_type'] = ret.text
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the uncertainty type attribute: {e}")
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
        except _et.Error as e:
            _logging.warning(f"unable to read the depth min and max values: {e}")
            return

        try:
            if ret_z_min is not None:
                self.data['z_min'] = float(ret_z_min.text)
            if ret_z_max is not None:
                self.data['z_max'] = float(ret_z_max.text)
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the depth min and max values: %s" % e)
            return

    def _read_source_name(self):
        """
        Read the source file name and store it in the object 'data' dictionary
        with the key 'filename'.
        """

        try:
            ret = self.xml_tree.find(self.bag_format['sourcename'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the survey name string: {e}")
            return

        try:
            self.data['sourcename'] = ret.text
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the survey name attribute: {e}")
            return

    def _read_SORDAT(self):
        """
        Reads a date, but what it means exactly needs to be researched...
        """

        try:
            ret = self.xml_tree.find(self.bag_format['SORDAT'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the SORDAT date string: {e}")
            return

        if ret is not None:
            try:
                text_date = ret.text
            except (ValueError, IndexError, AttributeError) as e:
                _logging.warning(f"unable to read the SORDAT date string: {e}")
                return

            tm_date = None
            try:
                parsed_date = _parser.parse(text_date)
                tm_date = parsed_date.strftime('%Y%m%d')
            except Exception:
                _logging.warning("unable to handle the date string: %s" % text_date)
                pass

            if tm_date is None and text_date != '':
                self.data['source_date'] = text_date
            elif tm_date != '':
                self.data['source_date'] = tm_date
        else:
            _logging.warning(f"unable to read the SORDAT date string: {ret}")

    def _read_survey_authority(self):
        """
        Read the survey authority name and store it in the object 'data'
        dictionary with the key 'SURATH'.
        """

        try:
            ret = self.xml_tree.find(self.bag_format['SURATH'],
                                     namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the survey authority name string: {e}")
            return

        try:
            self.data['agency'] = ret.text
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the survey authority name attribute: {e}")
            return

    def _read_survey_start_date(self):
        """
        Read the survey start date store it in the object 'data'
        dictionary with the key 'SURSTA'.
        """

        try:
            rets = self.xml_tree.find(self.bag_format['SURSTA'],
                                      namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the survey start date string: {e}")
            return

        if rets is not None:
            try:
                text_start_date = rets.text
            except (ValueError, IndexError, AttributeError) as e:
                _logging.warning(f"unable to read the survey start date string: {e}")
                return

            tms_date = None
            try:
                parsed_date = _parser.parse(text_start_date)
                tms_date = parsed_date.strftime('%Y%m%d')
            except Exception:
                # pass
                _logging.warning("unable to handle the survey start string: %s" % text_start_date)

            if tms_date is None and text_start_date != '':
                self.data['start_date'] = text_start_date
            elif tms_date != '':
                self.data['start_date'] = tms_date
        else:
            _logging.warning(f"unable to read the survey start date string: {rets}")

    def _read_survey_end_date(self):
        """
        Read the survey end date and store it in the object 'data'
        dictionary with the key 'SUREND'.
        """
        try:
            rete = self.xml_tree.find(self.bag_format['SUREND'],
                                      namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the survey end date string: {e}")

        if rete is not None:
            try:
                text_end_date = rete.text
            except (ValueError, IndexError, AttributeError) as e:
                _logging.warning(f"unable to read the survey end date string: {e}")

            tme_date = None
            try:
                parsed_date = _parser.parse(text_end_date)
                tme_date = parsed_date.strftime('%Y%m%d')
            except Exception:
                # pass
                _logging.warning("unable to handle the survey end string: %s" % text_end_date)

            if tme_date is None and text_end_date != '':
                self.data['end_date'] = text_end_date
            elif tme_date != '':
                self.data['end_date'] = tme_date
        else:
            _logging.warning(f"unable to read the survey end date string: {rete}")

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

        except _et.Error as e:
            _logging.warning(f"unable to read the survey vertical datum name string: {e}")
            return

        try:
            if val.lower() == 'unknown':
                raise ValueError(f'Invalid vertical datum assignment --> {val}')
            elif val.lower() in ('mean_lower_low_water', 'mean lower low water', 'mllw', 'mllw depth'):
                self.data['from_vert_key'] = 'MLLW'
            elif val.lower() in ('hudson river datum', 'hrd'):
                self.data['from_vert_key'] = 'HRD'
            self.data['from_vert_datum'] = val
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"Unable to read the survey vertical datum name attribute: {e}")
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
                        else:
                            continue
        except _et.Error as e:
            _logging.warning(f"unable to read the survey horizontal datum name string: {e}")
            return

        try:
            if datum is not None:
                if datum.lower() in ('north_american_datum_1983', 'north american datum 1983', 'nad83'):
                    self.data['from_horiz_frame'] = 'NAD83'
                elif datum.lower() in ('wgs_1984', 'wgs84'):
                    self.data['from_horiz_frame'] = 'WGS84'
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the survey horizontal datum name attribute: {e}")
            return

    def _read_platform_name(self):
        """
        Read the platform name and store it in the object 'data' dictionary
        with the key 'planam'.
        """

        try:
            ret = self.xml_tree.findall(self.bag_format['planam'],
                                        namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the platform name string: {e}")
            return

        try:
            if len(ret) > 0:
                self.data['planam'] = []
                for r in ret:
                    self.data['planam'].append(r.text)
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the platform name attribute: {e}")
            return

    def _read_sensor_types(self):
        """
        Read the sensor used and store it in the object 'data' dictionary
        with the key 'sensor'.
        """

        try:
            ret = self.xml_tree.findall(self.bag_format['sensor'],
                                        namespaces=self.namespace)
        except _et.Error as e:
            _logging.warning(f"unable to read the sensor name string: {e}")
            return

        try:
            if len(ret) > 0:
                self.data['sensor'] = []
                for r in ret:
                    self.data['sensor'].append(r.text)
        except (ValueError, IndexError, AttributeError) as e:
            _logging.warning(f"unable to read the sensor name attribute: {e}")
            return


class BAGSurvey(BAGRawReader):
    """


    """
    survey_meta = []

    def __init__(self, out_file_location):
        self.logger = _logging.getLogger('data.read')
        BAGRawReader.__init__(self, out_file_location)

    def read_metadata(self, survey_folder: str) -> [dict]:
        """
        Compiles the complete metadata of a survey

        Parameters
        ----------
        survey_folder : str
            Folder path containing contents of a single survey

        Returns
        -------
        [dict]
            List of survey metadata dicts

        """
        survey_meta = []

        if _os.path.isfile(survey_folder):
            survey_folder = _os.path.dirname(survey_folder)

        bag_files = _glob(_os.path.join(survey_folder, '*.bag'))

        for bag in bag_files:
            survey_meta.append(self.__file_meta(bag))

        combine = any([meta['interpolate'] for meta in survey_meta if (len(survey_meta) > 0) and 'interpolate' in meta])
        if combine:
            self.logger.debug(f'combine: {combine}')
            survey_meta = self._parse_groups(survey_meta)

        self.survey_meta = survey_meta

        return survey_meta

    def __file_meta(self, filename: str) -> dict:
        """
        read metadata from file

        Parameters
        ----------
        filename
            file to read

        Returns
        -------
            dict

        """

        try:
            meta_gdal, bag_version = self._parse_bag_gdal(filename)
            meta_xml = self._parse_bag_xml(filename, bag_version=bag_version)
            meta_support = self._known_meta(filename)
            meta_pickle = parse_file_pickle.read_pickle(filename)
            meta_combined = {**meta_pickle, **meta_xml, **meta_gdal, **meta_support}
            meta_final = self._finalize_meta(meta_combined)
            return meta_final
        except ValueError as error:
            _logging.warning(f'Error in reading metadata for {filename}: {error}')
            return {}

    def _parse_groups(self, metadata_list: [dict]) -> [dict]:
        """
        Put the BAGs into groups since some of the old BAGs persist on NCEI.
        This avoids combining data that should not go together, but does not
        solve the problem of adding in data that is not correct (need to
        fix this on NCEI)
        """
        files = [bag_file['from_path'] for bag_file in metadata_list if 'from_path' in bag_file]
        number_groups = {}

        for bag_file in files:
            numbers = _os.path.splitext(bag_file)[0].split('_')[-1].split('of')
            file_meta = [meta_entry for meta_entry in metadata_list if meta_entry['from_path'] == bag_file][0]
            if numbers[-1] in number_groups:
                number_groups[numbers[-1]].append(file_meta)
            else:
                number_groups[numbers[-1]] = [file_meta]

        return_list = []
        for group_number, group_list in number_groups.items():
            if len(group_list) == 1:
                return_list.extend(self._single_surface(group_list))
            else:
                return_list.extend(self._combined_surface(group_list))

        return return_list

    def _zoom_array(self, data_array, zoom_factor, ndv):
        """
        Zoom the provided array with a node centered (not pixel_is_area)
        approach.
        """
        self.logger.debug('Zooming...')
        start = _datetime.datetime.now()
        # build the expanded array
        shape = data_array.shape
        x_dim = zoom_factor * shape[1] - (zoom_factor - 1)
        y_dim = zoom_factor * shape[0] - (zoom_factor - 1)
        elev = _np.zeros((y_dim, x_dim)) + ndv
        elev[::zoom_factor, ::zoom_factor] = data_array

        # build a mask using binary morphology
        bin_arr = _np.zeros(elev.shape)
        idx = _np.nonzero(elev != ndv)
        bin_arr[idx] = 1
        dim = 2 * zoom_factor - 1
        struct = _np.ones((dim, dim), dtype=_np.bool)
        mask = binary_closing(bin_arr, structure=struct)

        # interpolate the gaps in the expanded array
        rel_idx = _np.arange(dim, dtype=_np.int) - int(zoom_factor / 2)
        footprint = _np.tile(rel_idx, (dim, 1))
        x_fp = footprint.flatten()
        y_fp = footprint.T.flatten()
        elev = _interpolate_gaps(elev, mask, y_fp, x_fp, ndv)

        dt = _datetime.datetime.now() - start
        self.logger.debug(f'Zooming completed in {dt} seconds')
        return elev

    def _insert_raster(self, bag_name, comb_bounds, comb_array, comb_res):
        """
        Insert the provided BAG file into the provided array.
        """
        self.logger.debug(f'Inserting {bag_name}')
        # open this file and extract the needed metadata
        bagfile_obj = Open(bag_name, pixel_is_area=False)
        nw, se = bagfile_obj.bounds
        res = bagfile_obj.resolution[0]
        ndv = bagfile_obj.nodata
        # if the data is lower res, expand array and fill in missing nodes
        if res > comb_res:
            if res % comb_res != 0:
                self.logger.warning(f'WARNING!!! Combine process will produce unintended results with {bag_name}')
            zoom_factor = int(res / comb_res)
            elev = self._zoom_array(bagfile_obj.elevation, zoom_factor, ndv)
        else:
            elev = bagfile_obj.elevation
        # determine offset from lower left for this array
        xoff = int(round((nw[0] - comb_bounds[0][0]) / comb_res))
        yoff = int(round((comb_bounds[0][1] - nw[1]) / comb_res))
        # insert the array
        idx = _np.nonzero(elev != ndv)
        yidx = idx[0] + yoff
        xidx = idx[1] + xoff
        comb_array[0, yidx, xidx] = elev[idx]
        del bagfile_obj
        return comb_array

    def _combined_surface(self, metadata_list: [dict]) -> [dict]:
        """
        Build a combine surface from the BAG files referenced in the metadata
        list provided.

        This method works with the data as nodes rather than cells.
        """
        # distill the information about the provided files.
        num_files = len(metadata_list)
        res = []
        files = []
        x = []
        y = []
        proj = []
        have_support_files = False
        for m in metadata_list:
            f = m['from_path']
            files.append(f)
            res.append(m['res'])
            if 'support_files' in m:
                for support_file in m['support_files']:
                    root, ext = _os.path.splitext(support_file)
                    if ext.lower() in ('.tif', '.tiff', '.gpkg'):
                        have_support_files = True
            bagfile_obj = Open(f, pixel_is_area=False)
            nw, se = bagfile_obj.bounds
            x.append([nw[0], se[0]])
            y.append([se[1], nw[1]])
            proj.append(bagfile_obj.wkt)
            ndv = bagfile_obj.nodata
            del bagfile_obj
        x = _np.array(x).flatten()
        y = _np.array(y).flatten()
        res = _np.array(res)
        resx = res[:, 0]
        minres = _np.abs(resx).min()

        # sort the file names by accending resolution
        order = _np.argsort(resx)

        # build the metadata
        combined_surface_metadata = {}
        combined_surface_metadata = metadata_list[order[0]].copy()
        survey_id = _os.path.basename(combined_surface_metadata['from_filename']).split('_')[0]
        combined_surface_metadata['from_filename'] = f"{survey_id}_Xof{num_files}.combined"

        if have_support_files:
            combined_surface_metadata['from_path'] = _os.path.join(self.out_file_location, f"{survey_id}_Xof{num_files}.bag")
            combined_surface_metadata['interpolate'] = True

            # build an array to house the combined dataset
            comb_bounds = ((min(x), max(y)), (max(x), min(y)))
            comb_shape = (int((max(y) - min(y)) / minres + 1), int((max(x) - min(x)) / minres + 1))
            comb_array = _np.zeros((2, comb_shape[0], comb_shape[1])) + ndv

            for n in order[::-1]:
                bag_file = files[n]
                comb_array = self._insert_raster(bag_file, comb_bounds, comb_array, minres)

            # convert the bounds to cell, aka pixel_is_area, and to (nw, se) format
            bounds = ((comb_bounds[0][0] - minres / 2., comb_bounds[0][1] + minres / 2.),
                      (comb_bounds[1][0] + minres / 2., comb_bounds[1][1] - minres / 2.))
            dataset_wkt = proj[order[0]]
            convert_dataset = BagToGDALConverter()
            convert_dataset.components2gdal(comb_array, bounds, res[order[0]], dataset_wkt)
            output_driver = _gdal.GetDriverByName('BAG')
            output_driver.CreateCopy(combined_surface_metadata['from_path'], convert_dataset.dataset)

            combined_surface_metadata['size'] = self._size_finder(combined_surface_metadata['from_path'])
        else:
            self.logger.warning('No support files found. Creating metadata entry but not creating combined surface for interpoaltion.')

        for bag_file in metadata_list:
            bag_file['interpolate'] = False if 'interpolate' in bag_file else False

        metadata_list.append(combined_surface_metadata)

        return metadata_list

    def _single_surface(self, metadata_list: dict) -> [dict]:
        metadata = metadata_list[0]
        if 'interpolate' in metadata and metadata['interpolate']:
            interp_meta = metadata.copy()
            interp_meta['from_filename'] = f"{metadata['from_filename']}.interpolated"
            # interp_meta['from_path'] = _os.path.join(self.out_file_location, f"{metadata['from_filename']}.bag")
            metadata['interpolate'] = False
            return [metadata, interp_meta]
        else:
            return metadata_list


class BagFile:
    """This class serves as the main container for BAG data."""

    def __init__(self, pixel_is_area=True, logger: _logging.Logger = None):
        if logger is None:
            logger = _logging.getLogger('fuse.read')
        self.logger = logger

        self.name = None
        self.nodata = 1000000.0
        self.elevation = None
        self.uncertainty = None
        self.shape = None
        self.pixel_is_area = pixel_is_area
        self.bounds = None
        self.resolution = None
        self.wkt = None
        self.size = None
        self.outfilename = None
        self.version = None

    def open_file(self, filepath: str, method: str):
        """
        Used to read a BAG file using the method determined by input.

        The current methods available are 'gdal' and 'hyo'

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file
        method : str
            The method used to open the file

        """

        if method == 'gdal':
            self._file_gdal(filepath)
        elif method == 'hack':
            self._file_hack(filepath)
        else:
            raise NotImplementedError('Open method not implemented.')

    def from_gdal(self, dataset: _gdal.Dataset):
        """
        Uses an already existing dataset to assign class attributes

        Parameters
        ----------
        dataset : _gdal.Dataset

        """

        try:
            self._known_data(dataset.GetFileList()[0])
        except TypeError:
            # TODO: Find a solution for this blasphemous hack
            pass
        self.elevation = self._gdalreadarray(dataset, 1)
        self.uncertainty = self._gdalreadarray(dataset, 2)
        self.shape = self.elevation.shape
        self.bounds, self.resolution = self._gt2bounds(dataset.GetGeoTransform(), self.shape)
        self.wkt = dataset.GetProjectionRef()
        self.version = dataset.GetMetadata()

    def _file_gdal(self, filepath: str):
        """
        Used to read a BAG file using OSGEO's GDAL module.

        This function reads and populates this object's attributes

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """

        self._known_data(filepath)
        bag_obj = _gdal.Open(filepath)
        self.elevation = self._gdalreadarray(bag_obj, 1)
        self.uncertainty = self._gdalreadarray(bag_obj, 2)
        self.shape = self.elevation.shape
        self.bounds, self.resolution = self._gt2bounds(bag_obj.GetGeoTransform(), self.shape)
        self.wkt = bag_obj.GetProjectionRef()
        self.version = bag_obj.GetMetadata()

        del bag_obj

    def _file_hack(self, filepath: str):
        """
        Used to read a BAG file using pytables and HDF5.

        This function reads and populates this object's attributes

        Object attribute 'pixel_is_area' will shift the bounds outward from the
        values reported in the BAG to the GDAL standard if set to True.

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """

        self._known_data(filepath)

        with _tb.open_file(filepath, mode='r') as bagfile:
            self.elevation = _np.flipud(bagfile.root.BAG_root.elevation.read())
            self.uncertainty = _np.flipud(bagfile.root.BAG_root.uncertainty.read())
            self.shape = self.elevation.shape
            meta_xml = ''.join(str(line, 'utf-8', 'ignore') for line in bagfile.root.BAG_root.metadata.read())
            xml_tree = _et.XML(_re.sub(r'<\?.+\?>', '', meta_xml))

            self.wkt = self._read_wkt_prj(xml_tree)
            if self.wkt is None:
                self.wkt = BAGRawReader().read_metadata(filepath)['from_horiz_datum']

            self.resolution = self._read_res_x_and_y(xml_tree)
            sw, ne = self._read_corners_sw_and_ne(xml_tree)
            sx, sy = sw

            # BAGs are a node-based convention; 1 cell is subtracted to account
            nx = (sx + (self.resolution[0] * (self.shape[1] - 1)))
            ny = (sy + (self.resolution[0] * (self.shape[0] - 1)))

            if self.pixel_is_area:
                # Convert to cell based-convention
                half_cell = 0.5 * self.resolution[0]
                self.bounds = ([sx - half_cell, ny + half_cell], [nx + half_cell, sy - half_cell])
            else:
                self.bounds = ([sx, ny], [nx, sy])

    def _known_data(self, filepath: str):
        """
        Assigns class attributes that are extension agnostic

        :attr:`name`, :attr:`nodata`, and :attr:`size` can be determined
        without using a specific library or method to open file contents, so
        they are assigned by using this function

        Parameters
        ----------
        filepath : str
            The complete file path of the input BAG file

        """

        _fName = _os.path.split(filepath)[-1]
        self.name = _os.path.splitext(_fName)[0]
        self.nodata = 1000000.0
        self.size = self._size_finder(filepath)

    def _nan2ndv(self, arr: _np.array, nodata: float) -> _np.array:
        """
        Reassigns numpy.nan values with a replacement no data value.

        Parameters
        ----------
        arr : numpy.array
            TODO write description
        nodata : float
            The no data value to be assigned to the numpy.array

        Returns
        -------
        type
            numpy.array

        """

        arr[_np.isnan(arr)] = nodata
        return _np.flipud(arr)

    def _npflip(self, arr: _np.array) -> _np.array:
        """
        Performs :func:`numpy.flipud` on the input numpy.array

        Parameters
        ----------
        arr : numpy.array
            TODO write description

        Returns
        -------
        type
            numpy.array

        """

        return _np.flipud(arr)

    def _meta2bounds(self, meta) -> ((float, float), (float, float)):
        """
        Breaks up and assigns the NW and SE corners from the NE and SW corners
        of a :obj:`hyo2.bag.meta` object

        Parameters
        ----------
        meta : hyo2.bag.meta
            TODO write description

        Returns
        -------
        type
            tuple of bounds

        """

        sx, sy = meta.sw
        nx, ny = meta.ne
        return (sx, ny), (nx, sy)

    def _gt2bounds(self, meta, shape: (int, int)) -> (((float, float), (float, float)), float):
        """
        Formats and returns the bounds and resolution

        This function takes a GeoTransform object and array shape and
        calculates the NW and SE corners.

        Object attribute 'pixel_is_area' will shift the bounds inward from the
        gdal geotransform standard to the BAG standard if set to False.

        Parameters
        ----------
        meta : gdal.GetGeoTransform
            TODO write description
        shape : tuple of int
            (y, x) shape of the bag object

        Returns
        -------
        type
            tuple of bounds, resolution

        """

        y, x = shape
        res = (meta[1], meta[5])
        # ulx, uly = _np.round(meta[0]), _np.round(meta[3])
        ulx, uly = meta[0], meta[3]
        lrx = ulx + (x * res[0])
        lry = uly + (y * res[1])
        # res = (_np.round(meta[1], 2), _np.round(meta[5], 2))
        if not self.pixel_is_area:
            # shift from cell bounds to node bounds
            ulx, lry = ulx - res[0] / 2., lry - res[1] / 2.
            lrx, uly = lrx + res[0] / 2., uly + res[1] / 2.
        return ([ulx, uly], [lrx, lry]), res

    def _gdalreadarray(self, bag_obj, band: int) -> _np.array:
        """
        Returns a numpy.array of a GDAL object

        This function uses a :obj:`gdal.Dataset` object and band number to pull
        and read the appropriate array from the object

        Parameters
        ----------
        bag_obj : gdal.Dataset
            TODO write description
        band : int
            raster band number

        Returns
        -------
        type
            numpy.array

        """

        return bag_obj.GetRasterBand(band).ReadAsArray()

    def _size_finder(self, filepath: str) -> int:
        """
        Returns the rounded size of a file in MB as an integer

        Parameters
        ----------
        filepath : str, os.Pathlike
            TODO write description

        Returns
        -------
        int

        """

        return int(_np.round(_os.path.getsize(filepath) / 1000))

    def generate_name(self, outlocation: str, io: bool):
        """
        Assigns the :attr:`outfilename` of the bag object.

        This function uses the input bool ``io`` to determine the naming
        convention of the output file name. This requires an assigned
        :attr:`name` for this object that is not ``None``

        Parameters
        ----------
        outlocation : str, os.Pathlike
            Folder path of the output file
        io : bool
            Boolean determination of interpolated only or full bag data

        """

        if self.name is not None:
            name = f'{self.name}_INTERP_{"ONLY" if io else "FULL"}.bag'
            self.outfilename = _os.path.join(outlocation, name)

    def _read_res_x_and_y(self, xml_tree):
        """ attempts to read resolution along x- and y- axes """

        try:
            ret = xml_tree.xpath('//*/gmd:spatialRepresentationInfo/gmd:MD_Georectified/'
                                 'gmd:axisDimensionProperties/gmd:MD_Dimension/gmd:resolution/gco:Measure',
                                 namespaces=_ns)
        except _et.Error as e:
            self.logger.error(f"unable to read res x and y: {e}")
            return

        if len(ret) == 0:
            try:
                ret = xml_tree.xpath('//*/spatialRepresentationInfo/smXML:MD_Georectified/'
                                     'axisDimensionProperties/smXML:MD_Dimension/resolution/'
                                     'smXML:Measure/smXML:value',
                                     namespaces=_ns2)
            except _et.Error as e:
                self.logger.error(f"unable to read res x and y: {e}")
                return

        try:
            res_x = float(ret[0].text)
            res_y = -float(ret[1].text)
            return res_x, res_y
        except (ValueError, IndexError) as e:
            self.logger.error(f"unable to read res x and y: {e}")
            return

    def _read_corners_sw_and_ne(self, xml_tree):
        """ attempts to read corners SW and NE """
        try:
            ret = xml_tree.xpath('//*/gmd:spatialRepresentationInfo/gmd:MD_Georectified/'
                                 'gmd:cornerPoints/gml:Point/gml:coordinates',
                                 namespaces=_ns)[0].text.split()
        except (_et.Error, IndexError) as e:
            try:
                ret = xml_tree.xpath('//*/spatialRepresentationInfo/smXML:MD_Georectified/'
                                     'cornerPoints/gml:Point/gml:coordinates',
                                     namespaces=_ns2)[0].text.split()
            except (_et.Error, IndexError) as e:
                self.logger.error(f"unable to read corners SW and NE: {e}")
                return
        try:
            sw = [float(c) for c in ret[0].split(',')]
            ne = [float(c) for c in ret[1].split(',')]
            return sw, ne
        except (ValueError, IndexError) as e:
            self.logger.error(f"unable to read corners SW and NE: {e}")
            return

    def _read_wkt_prj(self, xml_tree):
        """ attempts to read the WKT projection string """
        try:
            ret = xml_tree.xpath('//*/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/'
                                 'gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString',
                                 namespaces=_ns)
        except _et.Error as e:
            self.logger.error(f"unable to read the WKT projection string: {e}")
            return
        if len(ret) == 0:
            try:
                ret = xml_tree.xpath('//*/referenceSystemInfo/smXML:MD_CRS',
                                     namespaces=_ns2)
            except _et.Error as e:
                self.logger.error(f"unable to read the WKT projection string: %s: {e}")
                return
            if len(ret) != 0:
                self.logger.warning("unsupported method to describe CRS")
                return

        try:
            wkt_srs = ret[0].text
            return wkt_srs
        except (ValueError, IndexError) as e:
            self.logger.error(f"unable to read the WKT projection string: {e}")
            return


class Open(BagFile):
    """
    Opens a BAG file using :obj:`BagFile`'s :func:`BagFile.open_file` or
    :func:`BagFile.from_gdal` method
    """

    def __init__(self, dataset, pixel_is_area=True):
        BagFile.__init__(self, pixel_is_area)
        if type(dataset) == str:
            self.open_file(dataset, 'hack')
        elif type(dataset) == _gdal.Dataset:
            self.from_gdal(dataset)


def _from_rasterio(bag_file: str) -> _rio.MemoryFile:
    bagfile_obj = Open(bag_file)
    resx, resy = bagfile_obj.resolution
    nw, se = bagfile_obj.bounds
    transform = _rio.Affine(resx, 0, nw[0], 0, resy, nw[1])
    height = bagfile_obj.shape[0]
    width = bagfile_obj.shape[1]
    dtype = bagfile_obj.elevation.dtype

    with _rio.MemoryFile() as memfile:
        with memfile.open(**{'driver': 'GTiff', 'dtype': dtype, 'nodata': 1000000.0,
                             'width': width, 'height': height, 'count': 2, 'crs': _rio.crs.CRS.from_wkt(bagfile_obj.wkt),
                             'transform': transform, 'tiled': False, 'compress': 'lzw', 'interleave': 'band'}) as dataset:
            dataset.write(_np.array([bagfile_obj.elevation, bagfile_obj.uncertainty]))
            del bagfile_obj

        return memfile.open()


class BagToGDALConverter:
    """This class serves as the main container for converting processed bag
    data into a :obj:`gdal.Dataset` object.

    """
    _descriptions = ['Elevation', 'Uncertainty', 'Interpolated']

    def __init__(self, out_verdat: str = None, flip: bool = False):
        """

        Parameters
        ----------
        out_verdat : str
            Output Vertical Coordinate System (ie. 'MLLW')
        """

        self.dataset = None
        self.out_verdat = out_verdat
        self.flip = flip

    def bag2gdal(self, bag: BagFile):
        """
        Converts a :obj:`bag` object into a :obj:`gdal.Dataset`

        Parameters
        ----------
        bag
            BAG file object
        """

        arrays = [bag.elevation, bag.uncertainty]
        bands = len(arrays)
        nw, se = bag.bounds
        nwx, nwy = nw
        scx, scy = se
        y_cols, x_cols = bag.shape
        res_x, res_y = bag.resolution[0], bag.resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols, bands, _gdal.GDT_Float32)
        target_gt = (nwx, res_x, 0, nwy, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        srs = _osr.SpatialReference(wkt=bag.wkt)
        #        if not srs.IsCompound():
        #            srs.SetVertCS(self.out_verdat, self.out_verdat, 2000)
        wkt = srs.ExportToWkt()
        target_ds.SetProjection(wkt)

        for index, item in enumerate(arrays):
            band = target_ds.GetRasterBand(index + 1)
            band.SetDescription(self._descriptions[index])
            band.SetNoDataValue(bag.nodata)

            if self.flip:
                item = _np.flipud(item)

            band.WriteArray(item)
            del band

        self.dataset = target_ds
        del target_ds

    def components2gdal(self, arrays: [_np.array], bounds: ((float, float), (float, float)),
                        resolution: (float, float), prj: str, nodata: float = 1000000.0):
        """
        Converts raw dataset components into a :obj:`gdal.Dataset` object

        Parameters
        ----------
        arrays : list of numpy.array
            Arrays [elevation, uncertainty]
        bounds : tuple of tuple of float
            (NW, SE) corners of the data
        resolution : tuple of float
            (x_res, y_res) of the data
        prj : str
            WKT string of the projection of the data
        nodata : float, optional
            no data value of the arrays, default is 1000000.0


        """
        dims = arrays.shape
        if len(dims) == 3:
            bands, y_cols, x_cols = dims
        elif len(dims) == 2:
            raise ValueError('Only one band in array found.  Missing uncertainty layer?')
        else:
            raise ValueError('Shape of provided array is nonsensical')
        nw, se = bounds
        nwx, nwy = nw
        scx, scy = se
        res_x, res_y = resolution[0], resolution[1]
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_cols, y_cols, bands, _gdal.GDT_Float32)
        target_gt = (nwx, res_x, 0, nwy, 0, res_y)
        target_ds.SetGeoTransform(target_gt)
        srs = _osr.SpatialReference(wkt=prj)
        #        if not srs.IsCompound():
        #            srs.SetVertCS(self.out_verdat, self.out_verdat, 2000)
        wkt = srs.ExportToWkt()
        target_ds.SetProjection(wkt)
        x = 1

        for item in arrays:
            band = target_ds.GetRasterBand(x)
            band.SetDescription(self._descriptions[x - 1])
            band.SetNoDataValue(nodata)

            if self.flip:
                item = _np.flipud(item)

            band.WriteArray(item)
            del band
            x += 1

        self.dataset = target_ds
        del target_ds


@jit(nopython=True)
def _interpolate_gaps(data, mask, y_footprint, x_footprint, ndv):
    """
    Return the provided 2d data with the missing nodes not covered by the mask.
    """
    y_idx, x_idx = _np.nonzero((mask == 1) & (data == ndv))
    numpts = x_idx.shape[0]
    interp_vals = _np.zeros(numpts)
    for n in range(numpts):
        # need to add a check for bounds
        y = y_footprint + y_idx[n]
        x = x_footprint + x_idx[n]
        count = 0
        tmp = 0
        for m in range(x_footprint.shape[0]):
            val = data[y[m], x[m]]
            if val != ndv:
                count = count + 1
                tmp = tmp + val
        if count != 0:
            interp_vals[n] = tmp / count
    for n in range(numpts):
        data[y_idx[n], x_idx[n]] = interp_vals[n]
    return data
