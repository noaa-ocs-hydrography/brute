# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:17:08 2019

@author: Casiano.Koprowski
"""

import configparser
import datetime
import json
import os
import re
import socket
import urllib
import zipfile
from typing import Tuple, List, Union, TextIO, Any

import numpy as np
import requests
from osgeo import osr, ogr, gdal

"""Known global constants"""
prog_loc = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.split(prog_loc)[0]
config = configparser.ConfigParser()
config.read('config.ini')

if config['Downloads']['Downloads'] == '':
    downloads = os.path.join(prog_loc, 'downloads')
else:
    downloads = config['Downloads']['Downloads']

if not os.path.isdir(os.path.join(downloads)):
    os.mkdir(downloads)

area_to_quarter = {'area': 'CHANNELAREAIDPK', 'quarter': 'CHANNELAREAIDFK', 'reach': 'CHANNELAREAIDFK'}
reach_to_quarter = {'reach': 'CHANNELREACHIDPK',  'quarter': 'CHANNELREACHIDFK'}

def create_polygon(coords: List[Tuple[float, float]]) -> ogr.Geometry:
    """
    Creates an ogr.Geometry/wkbLinearRing object from a list of coordinates.

    with considerations from:
    https://gis.stackexchange.com/q/217165

    Parameters
    ----------
    coords :
        A list of [x, y] points to be made into an ogr.Geometry/wkbLinearRing object

    Returns
    -------
    type
        ogr.Geometry/wkbLinearRing object

    """

    ring = ogr.Geometry(ogr.wkbLinearRing)

    for coord in coords:
        ring.AddPoint(coord[0], coord[1], 1)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def create_multipolygon(polys: List[ogr.Geometry]) -> str:
    """
    Creates an ogr.Geometry/wkbMultiPolygon object from a list of
    ogr.Geometry/wkbLinearRing objects.  The ogr.Geometry/wkbMultiPolygon is
    translated and returned as a WTK Multipolygon object.

    with considerations from:
    https://gis.stackexchange.com/q/217165
    and:
    https://pcjericks.github.io/py-gdalogr-cookbook/geometry.html#create-a-multipolygon

    Parameters
    ----------
    polys :
        A list of ogr.Geometry/wkbLinearRing objects
    polys: List[ogr.Geometry] :


    Returns
    -------
    type
        WTK Multipolygon object

    """

    multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)

    for poly in polys:
        multipolygon.AddGeometry(poly)

    return multipolygon.ExportToWkt()


def geometryToShape(coordinates: list):
    """
    Uses a list of coordinate point 'rings' and creates a WTK Multipolygon
    object from them.  This object represents the survey outline.

    eHydro data object geometries are returned as a list of lists/'rings'
    meaning that a survey may have one or many polygons included in it's
    geometry.

    This function takes each 'ring' and determines it's extents and creates a
    ogr.Geometry object for it using :func:`create_polygon`

    The polygons for each 'ring' are then combined into a single WTK
    Multipolygon object using :func:`create_multipolygon`

    The total extent of the geometry and the WTK Multipolygon object are
    returned

    Parameters
    ----------
    coordinates :
        A list of coordinate point 'rings' returned by an eHydro survey query
        in the Geometry attribute


    Returns
    -------
    type
        A WTK Multipolygon object representing the survey outline

    """

    polys = []
    bounds = []

    for ring in coordinates:
        ring = np.array(ring)
        x = ring[:, 0]
        y = ring[:, 1]
        bound = [[np.amin(x), np.amax(y)], [np.amax(x), np.amin(y)]]
        bounds.extend(bound)
        poly = create_polygon(ring)
        polys.append(poly)

    multipoly = create_multipolygon(polys)
    bounds = np.array(bounds)
    xb = bounds[:, 0]
    yb = bounds[:, 1]
    bounds = (np.amin(xb), np.amax(yb)), (np.amax(xb), np.amin(yb))
    return multipoly


def write_geopackage(out_path: str, metadata: dict,
                     spcs: Union[str, int] = 4326):
    """
    Writes out a geopackage containing the bounding geometry of
    of the given query.

    Derived from:
    https://gis.stackexchange.com/a/52708/8104
    via
    https://gis.stackexchange.com/q/217165

    Parameters
    ----------
    out_path : str, os.Pathlike
        String representing the complete file path for the output geopackage
    metadata : dict
        Dictionary holding object metadata
    spcs :
        The ESPG code for the data

    """

    area_name = metadata['CHANNELAREAIDPK']
    area_poly = metadata['geometry']
    metadata_keys = list(metadata.keys())[:-2]
    output_gpkg = os.path.join(out_path, area_name) + ".gpkg"
    output_esri = os.path.join(out_path, area_name) + ".shp"

    # Reference
    if type(spcs) == str:
        proj = osr.SpatialReference(wkt=spcs)
        proj.MorphFromESRI()
    else:
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(spcs)

    # Now convert it to a geopackage with OGR
    driver = ogr.GetDriverByName('GPKG')
    ds = driver.CreateDataSource(output_gpkg)

    layer = ds.CreateLayer(area_name, proj, ogr.wkbMultiPolygon)

    # Add one attribute
    for key in metadata_keys:
        layer.CreateField(ogr.FieldDefn(key, ogr.OFTString))
    defn = layer.GetLayerDefn()

    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    for key in metadata_keys:
        feat.SetField(key, metadata[key])

    # Make a geometry, from wkt object
    geom = ogr.CreateGeometryFromWkt(area_poly)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)

    del layer, feat, geom

    for quarter in metadata['quarters']:
        quarter_name = quarter['CHANNELQUARTERIDPK']
        quarter_poly = quarter['geometry']
        quarter_keys = list(quarter.keys())[:-1]

        quarter_layer = ds.CreateLayer(quarter_name, proj, ogr.wkbMultiPolygon)

        if quarter_layer is not None:
            # Add one attribute
            for key in quarter_keys:
                quarter_layer.CreateField(ogr.FieldDefn(key, ogr.OFTString))
            quarter_defn = quarter_layer.GetLayerDefn()

            # Create a new feature (attribute and geometry)
            quarter_feat = ogr.Feature(quarter_defn)
            for key in quarter_keys:
                quarter_feat.SetField(key, quarter[key])

            # Make a geometry, from wkt object
            quarter_geom = ogr.CreateGeometryFromWkt(quarter_poly)
            quarter_feat.SetGeometry(quarter_geom)
            quarter_layer.CreateFeature(quarter_feat)

            del quarter_layer, quarter_feat, quarter_geom
        else:
            print(f"Failed to create {quarter_name}")

    try:
        gdal.VectorTranslate(output_esri, output_gpkg, format='ESRI Shapefile')
    except ValueError:
        pass

    # Save and close everything
    del ds




def list_chunks(input_list: list, n: int = 1000):
    """
    Returns a list of lists given a maximum length for each list;

    Parameters
    ----------
    input_list :
        list
    n : int, optional
        maximum lengh of the subdivided lists (Default = 1000)

    """
    for i in range(0, len(input_list), n):
        yield input_list[i:i + n]


def objectID_query(query_id: int = 1) -> [int]:
    """
    Holds the Queries for the eHydro REST API, asks for responses, and uses
    the json library to make it readable by the program. Returns the json
    responses for number of surveys and surveys in the response.

    The function uses the requests library to retrieve the API's response(s) and
    uses the json library to make them readable by the program.

    The function returns the json object containing the contents of the query
    and the integer number of surveys contained by the query

    Parameters
    ----------
    query_id : int
        Query type performed ``1`` for ChannelArea ``3`` for ChannelReach (Default == 1)

    Returns
    -------
    [int]
        List of object ids from query

    """

    # The main query parameters that will determine the contents of the response
    if config['Agencies']['Agencies'] != '':
        areas = ''
        agencies = config['Agencies']['Agencies'].split(',')
        if len(agencies) > 1:
            i = 0
            areas += '('
            while i < len(agencies):
                if i == 0:
                    # %27 is ' (Apostrophe); %25 is % (Percent sign)
                    areas += f'UPPER(USACEDISTRICTCODE)+like+%27%25{agencies[0].strip()}%25%27'
                    i += 1
                else:
                    # %27 is ' (Apostrophe); %25 is % (Percent sign)
                    areas += f'+OR+UPPER(USACEDISTRICTCODE)+like+%27%25{agencies[i].strip()}%25%27'
                    i += 1
            areas += ')+'
        else:
            # %27 is ' (Apostrophe); %25 is % (Percent sign)
            areas = f'UPPER(USACEDISTRICTCODE)+like+%27%25{agencies[0].strip()}%25%27'
    else:
        areas = ''

    if areas != '':
        # areas
        where = areas
    else:
        # 1 = 1; %3D is = (Equals sign)
        where = '1%3D1'

    # The query for returning the object IDs
    objIDs = f'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/National_Channel_Framework/FeatureServer/{query_id}/query' + \
             f'?&where={where}&outFields=*&returnGeometry=false&returnIdsOnly=true&outSR=&f=json'

    # print(objIDs)

    # Initial Query execution
    objectIds_request = requests.get(objIDs)

    # Response for survey Object IDs as a json object
    objectIds_json = objectIds_request.json()
    objectIds = objectIds_json['objectIds']

    return objectIds


def geometry_query(objectIds: [int], query_id: int = 1, geometry: bool = True) -> [dict]:
    """
    Uses the json object return of the each queried survey id and the total
    number of surveys included to compile a list of complete returned survey
    data, as provided in the response.

    The function returns the lists of returned survey data as a dictionary
    containing the field_name - value relationships.

    Parameters
    ----------

    Returns
    -------

    """

    x = 0
    rows = []
    attr_list = []
    id_chunks = list(list_chunks([str(objectID) for objectID in objectIds], n=100))

    for chunk in id_chunks:
        id_string = ','.join(chunk)
        query = f'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/National_Channel_Framework/FeatureServer/{query_id}/query' + \
                f'?where=&objectIds={id_string}&outFields=*&returnGeometry={str(geometry).lower()}&outSR=4326&f=json'
        response = requests.get(query)
        page = response.json()

        if id_chunks.index(chunk) == 0:
            fields = page['fields']
            for attr in fields:
                if attr['name'] == 'OBJECTID':
                    pass
                else:
                    attr_list.append(attr['name'])

        object_num = 0
        for objectID in chunk:
            metadata = {}
            for attribute in attr_list:
                try:
                    if attribute.upper() in ("DATEUPLOADED", "DATEEDITED", "DATECREATED", "DATELASTDREDGE"):
                        date = (page['features'][object_num]['attributes'][attribute])

                        if date is None:
                            continue

                        try:
                            date = datetime.datetime.utcfromtimestamp(date / 1000)
                            metadata[attribute] = f'{date:%Y-%m-%d}'
                        except (OSError, TypeError) as e:
                            print(e, date)
                    elif page['features'][object_num]['attributes'][attribute] is not None:
                        metadata[attribute] = str(page['features'][object_num]['attributes'][attribute])
                except KeyError as e:
                    print(e, attribute)
            if geometry:
                try:
                    coords = page['features'][object_num]['geometry']['rings']
                    metadata['geometry'] = geometryToShape(coords)
                except KeyError as e:
                    print(e, 'geometry')

            rows.append(metadata)
            object_num += 1
            x += 1
            print(x, end=' ')
    print('rows complete')
    return rows


def save_polygons(channel_areas: [dict], channel_reaches: [dict], channel_quarters: [dict]):
    """
    Saves the geometry as a Geopackage (.gpkg) for each area and it's associated quarter(s)

    Parameters
    ----------
    channel_areas : [dict]
        list of area info dicts
    channel_quarters : [dict]
        list of quarter info dicts

    """

    a = 0
    for area in channel_areas:
        district_name = area['USACEDISTRICTCODE']
        district_path = os.path.join(downloads, district_name)

        if not os.path.isdir(district_path):
            os.mkdir(district_path)

        area['quarters'] = []
        area_name = area[area_to_quarter['area']]
        area_path = os.path.join(district_path, area_name)

        if not os.path.isdir(area_path):
            os.mkdir(area_path)


        for reach in channel_reaches:
            reach_area = reach[area_to_quarter['reach']]
            reach_name = reach[reach_to_quarter['reach']]

            reach_quarters = []
            if reach_area == area_name:
                for quarter in channel_quarters:
                    quarter_reach = quarter[reach_to_quarter['quarter']]

                    if quarter_reach == reach_name:
                        reach_quarters.append({**reach, **quarter})
                area['quarters'].extend(reach_quarters)

        if 'geometry' in area:
            write_geopackage(area_path, area)
            a += 1

        print(f"{a} of {len(channel_areas)}: {area_name}, completed")


def main():
    channel_ids = objectID_query()
    channel_areas = geometry_query(channel_ids)
    reach_ids = objectID_query(query_id=2)
    channel_reaches = geometry_query(reach_ids, query_id=2, geometry=False)
    quarter_ids = objectID_query(query_id=3)
    channel_quarters = geometry_query(quarter_ids, query_id=3)
    save_polygons(channel_areas, channel_reaches, channel_quarters)

if __name__ == "__main__":
    main()
