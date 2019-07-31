# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:39:57 2019

@author: Casiano.Koprowski
"""

import lxml.etree as _et

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


def read_res_x_and_y(xml_tree):
    """ attempts to read resolution along x- and y- axes """

    try:
        ret = xml_tree.xpath('//*/gmd:spatialRepresentationInfo/gmd:MD_Georectified/'
                             'gmd:axisDimensionProperties/gmd:MD_Dimension/gmd:resolution/gco:Measure',
                             namespaces=_ns)
    except _et.Error as e:
        print("unable to read res x and y: %s" % e)
        return

    if len(ret) == 0:
        try:
            ret = xml_tree.xpath('//*/spatialRepresentationInfo/smXML:MD_Georectified/'
                                 'axisDimensionProperties/smXML:MD_Dimension/resolution/'
                                 'smXML:Measure/smXML:value',
                                 namespaces=_ns2)
        except _et.Error as e:
            print("unable to read res x and y: %s" % e)
            return

    try:
        res_x = float(ret[0].text)
        res_y = -float(ret[1].text)
        return res_x, res_y
    except (ValueError, IndexError) as e:
        print("unable to read res x and y: %s" % e)
        return


def read_corners_sw_and_ne(xml_tree):
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
            print("unable to read corners SW and NE: %s" % e)
            return
    try:
        sw = [float(c) for c in ret[0].split(',')]
        ne = [float(c) for c in ret[1].split(',')]
        return sw, ne
    except (ValueError, IndexError) as e:
        print("unable to read corners SW and NE: %s" % e)
        return


def read_wkt_prj(xml_tree):
    """ attempts to read the WKT projection string """
    try:
        ret = xml_tree.xpath('//*/gmd:referenceSystemInfo/gmd:MD_ReferenceSystem/'
                             'gmd:referenceSystemIdentifier/gmd:RS_Identifier/gmd:code/gco:CharacterString',
                             namespaces=_ns)
    except _et.Error as e:
        print("unable to read the WKT projection string: %s" % e)
        return
    if len(ret) == 0:
        try:
            ret = xml_tree.xpath('//*/referenceSystemInfo/smXML:MD_CRS',
                                 namespaces=_ns2)
        except _et.Error as e:
            print("unable to read the WKT projection string: %s" % e)
            return
        if len(ret) != 0:
            print("unsupported method to describe CRS")
            return

    try:
        wkt_srs = ret[0].text
        return wkt_srs
    except (ValueError, IndexError) as e:
        print("unable to read the WKT projection string: %s" % e)
        return
