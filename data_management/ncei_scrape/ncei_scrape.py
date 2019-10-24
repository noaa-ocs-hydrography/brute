# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 08:51:53 2019

@author: Casiano.Koprowski
"""

import ast
import configparser
import datetime
import gzip
import json
import numpy as np
import os
import pickle
import re
import shutil
import socket
import urllib
from typing import Union, Dict, List

import requests

#from osgeo import gdal
from osgeo import ogr
from osgeo import osr


"""Known global constants"""
progLoc = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.split(progLoc)[0]
config = configparser.ConfigParser()
config.read('config.ini')

if config['Destination']['Folder'] == '':
    downloads = os.path.join(progLoc, 'downloads')
else:
    downloads = config['Destination']['Folder']

if not os.path.isdir(os.path.join(downloads)):
    os.mkdir(downloads)

# https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H12001-H14000/H12001/[BAG, TIFF]

zList = ['xmin', 'ymin', 'xmax', 'ymax']
attributes = {3: ['Name', 'SURVEY_ID', 'CELL_SIZE'],
              0: ['*'],
              1: ['*']}
date_fields = ['DATE_SURVEY_BEGIN', 'DATE_SURVEY_END', 'DATE_MODIFY_DATA',
               'DATE_SURVEY_APPROVAL', 'START_TIME', 'END_TIME']


def wgs84_to_esri(min_x: float, min_y: float, max_x: float, max_y: float) -> dict:
    """
    Given input coordinates (LatLong; WGS84), uses a NCEI webservice to transform
    into ESRI Pseudo-Mercator

    Parameters
    ----------
    min_x
    min_y
    max_x
    max_y


    Returns
    -------
    dict
        {min_x: float,
         min_y: float,
         max_x: float,
         max_y: float,
        }


    """

    coords = {'geometryType': 'esriGeometryPoint', 'geometries': [{'x': min_x, 'y': min_y}, {'x': max_x, 'y': max_y}]}
    query = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/Utilities/Geometry/GeometryServer/project' + \
            f'?inSR=4326&outSR=102100&geometries={coords}&transformation=&transformForward=true&vertical=false&f=json'

    coordsRequest = requests.get(query)
    coordsRequestJSON = coordsRequest.json()
    bounds = []
    z = 0
    for i in range(len(coordsRequestJSON['geometries'])):
        for k, j in coordsRequestJSON['geometries'][i].items():
            bounds.append((zList[z], j))
            z += 1
    bounds = dict(bounds)
    return bounds


def region_bounds(region_file: str) -> [dict]:
    """
    Uses an input geopackage polygon file to determine the geographic extent of
    the region for use by the query

    Parameters
    ----------
    region_file
        The filepath of the input region's geopackage polygon

    Returns
    -------
        list of dict
            A list of the bounding area(s) of a region

    """
    region_vector = ogr.Open(region_file)
    layer = region_vector.GetLayer()

    bounds = {}

    for feature in layer:
        if feature is not None:
            geom = feature.GetGeometryRef()
            points = geom.GetBoundary()
            if points.GetGeometryName() == 'MULTILINESTRING':
                for i in range(0, geom.GetGeometryCount()):
                    inner_geom = geom.GetGeometryRef(i)
                    inner_points = inner_geom.GetBoundary()
                    x_points, y_points = [], []
                    for i in range(0, inner_points.GetPointCount()):
                        # GetPoint returns a tuple not a Geometry
                        point = inner_points.GetPoint(i)
                        x_points.append(round(point[0], 2))
                        y_points.append(round(point[1], 2))
                    bounds[i] = (x_points, y_points)
            else:
                x_points, y_points = [], []
                for i in range(0, points.GetPointCount()):
                    # GetPoint returns a tuple not a Geometry
                    point = points.GetPoint(i)
                    x_points.append(round(point[0], 2))
                    y_points.append(round(point[1], 2))
                bounds[0] = (x_points, y_points)
                break

    del region_vector

    return_bounds = []

    if len(bounds.keys()) > 1:
        for key in bounds.keys():
            x_points, y_points = bounds[key]
            return_bounds.append(wgs84_to_esri(min(x_points), min(y_points), max(x_points), max(y_points)))
    else:
        x_points, y_points = bounds[0]
        return_bounds.append(wgs84_to_esri(min(x_points), min(y_points), max(x_points), max(y_points)))

    return return_bounds


def survey_objectID_query(bounds: [dict], qId=0) -> ([int], int):
    """
    Queries the ObjectIDs for the surveys within the geographic bounds given.

    This function may also take into account certain time periods, if they are
    defined in the configuration file

    Parameters
    ----------
    bounds
        list of bounding boxes
    qId
        (Default = 0)


    Returns
    -------
    list of int
        The ObjectIDs returned by the query
    int
        The number of ObjectIDs returned by the query

    """

    # Renamed bounds because Glen wanted it to be a que, not a queue
    que = bounds

    objectIDs = []
    objectNum = 0

    # Today (ex. '2018-08-08'), unformatted
    today = datetime.datetime.today()
    # Today - 1 (ex. '2018-08-06'), unformatted
    yesterday = today - datetime.timedelta(1)

    data_type, data_format = config_data_type()

    if config['Timeframe']['Start Date'] != '':
        start_parse =  datetime.datetime.strptime(config['Timeframe']['Start Date'], '%Y-%m-%d')
        start = f'{start_parse:%Y-%m-%d}'
    else:
        start = f'{yesterday:%Y-%m-%d}'

    if config['Timeframe']['End Date'] != '':
        end_parse =  datetime.datetime.strptime(config['Timeframe']['End Date'], '%Y-%m-%d')
        start = f'{end_parse:%Y-%m-%d}'
    else:
        end = f'{today:%Y-%m-%d}'

    if config['Timeframe']['Ignore Date'] in (False, 'False') and config['Timeframe']['Date Queries'] != '':
        query_fields = [date_query.strip() for date_query in config['Timeframe']['Date Queries'].split(',')]
        where = ''
        num_fields = len(query_fields)-1
        for field in query_fields:
            where += f'({field}+BETWEEN+timestamp+%27{start}+00:00:01%27+AND+timestamp+%27{end}+11:59:59%27)'
            if len(query_fields) > 0 and query_fields.index(field) < num_fields:
                where += '+OR+'
    else:
        where = '1%3D1'

    for box in que:
        dataList = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/{data_type}/MapServer/{qId}/query?where={where}&text=&objectIds=&time=&geometry={box}&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=true&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
        dataListRequest = requests.get(dataList)
        dataListRequestJSON = dataListRequest.json()
#            print (dataListRequestJSON)
        if dataListRequestJSON['objectIds'] is not None:
            objectIDs.extend(dataListRequestJSON['objectIds'])
            objectNum += len(objectIDs) - 1
            print(objectNum)
        else:
            objectIDs.extend([])

#    print(objectIDs, objectNum)
    return objectIDs, objectNum


def list_chunks(input_list: list, n: int = 1000) -> [list]:
    """
    Yeilds a list of lists given a maximum length for each list;

    Parameters
    ----------
    input_list :
        list
    n : int, optional
        maximum lengh of the subdivided lists (Default = 1000)

    Yeilds
    ------
    list
        List of lists of at most `n` length


    """
    for i in range(0, len(input_list), n):
        yield input_list[i:i + n]


def date_eval(row: dict, stored: dict) -> bool:
    """
    Evaluates date fields between two survey data objects.

    If either of ('DATE_SURVEY_APPROVAL', 'DATE_MODIFY_DATA') are different
    from the stored survey's data, True is returned. False otherwise

    Parameters
    ----------
    row
        The 'new' survey info returned by the query
    stored
        The 'old' survey info stored from a past query

    Returns
    -------
    bool

    """

    date_items = ('DATE_SURVEY_APPROVAL', 'DATE_MODIFY_DATA')

    results = []

    for key in date_items:
        try:
            if key not in row.keys() or key not in stored.keys():
                raise ValueError(f"Key not found for data: {key} in {row['SURVEY_ID']}")
            new, old = datetime.datetime.strptime(row[key], '%Y-%m-%d'), datetime.datetime.strptime(stored[key], '%Y-%m-%d')
            results.append(new > old)
        except ValueError as e:
            results.append(False)
    if True in results:
        return True
    else:
        return False


def create_polygon(coords: [(float, float)]) -> ogr.Geometry:
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


def create_multipolygon(polys: [ogr.Geometry]) -> str:
    """
    Creates an ogr.Geometry/wkbMultiPolygon object from a list of
    ogr.Geometry/wkbLinearRing objects.  The ogr.Geometry/wkbMultiPolygon is
    transelated and returned as a WTK Multipolygon object.

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


def parse_timestamp(timestamp: int) -> datetime.date:
    """
    Made with answers from https://stackoverflow.com/q/17231711

    Parameters
    ----------
    timestamp : int
        ESRI timestamp returned by REST query

    Returns
    -------
    datetime.date
        Parsed date from timestamp

    """
    if timestamp < 0:
        return datetime.datetime(1970, 1, 1) + datetime.timedelta(seconds=timestamp / 1000)
    else:
        return datetime.datetime.utcfromtimestamp(timestamp / 1000)


def survey_compile(objectIDs: list, num: int, history: [dict], poly: str, qId=0, pb=None) -> (list, [dict]):
    """
    Queries and compiles the data from the input objectIDs and checks against
    the history of past queries to mark them for updates

    Parameters
    ----------
    objectIDs
    num
    history
    qId
        (Default = 0)
    pb
        (Default = None)

    Returns
    -------
    (list, [dict])
        Tuple of the attribute fields returned by the query and the list of
        queried data objects

    """

    rows = []
    attr_list = []
    opts = ','.join(attributes[qId])
    if pb is not None:
        pb.SetRange(num)

    data_type, data_format = config_data_type()

    id_chunks = list(list_chunks([str(objectID) for objectID in objectIDs]))

    region_ds = ogr.Open(poly)
    region_layer = region_ds.GetLayer()

    for feature in region_layer:
        if feature is not None:
            region_geom = feature.GetGeometryRef()
            break

    for chunk in id_chunks:
        id_string = ','.join(chunk)
        query = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/{data_type}/MapServer/{qId}/query' + \
                f'?where=&text=&objectIds={id_string}&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields={opts}&returnGeometry=true&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=4326&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
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
            print(object_num, end=' ')

            row = {}
            try:
                for attribute in attr_list:
                    if attribute in date_fields:
                        if page['features'][object_num]['attributes'][attribute] is not None:
                            date = (page['features'][object_num]['attributes'][attribute])
                            try:
                                row[attribute] = f'{parse_timestamp(date):%Y-%m-%d}'
                            except OSError as e:
                                print(f"OSError - {e} for {attribute}: {date}")
                    elif page['features'][object_num]['attributes'][attribute] is not None:
                        row[attribute] = str(page['features'][object_num]['attributes'][attribute])
            except KeyError as e:
                print(e, page)
            try:
                coords = page['features'][object_num]['geometry']['rings']
                row['poly'] = geometryToShape(coords)
            except KeyError as e:
                print(e, 'invalid geometry')

            if 'BAGS_EXIST' in row:
                bags_exist = True if row['BAGS_EXIST'].upper() == 'TRUE' else False
            else:
                bags_exist = False

            if 'poly' in row:
                try:
                    survey_geom = ogr.CreateGeometryFromWkt(row['poly'])
                    try:
                        intersection = region_geom.Intersection(survey_geom)
                        flag = intersection.ExportToWkt()
                    except AttributeError as e:
                        flag = 'GEOMETRYCOLLECTION EMPTY'
                except TypeError as e:
                    flag = 'GEOMETRYCOLLECTION EMPTY'

                if flag != 'GEOMETRYCOLLECTION EMPTY' and data_format == 'BPS' and not bags_exist:
                    rows.append(row)
                elif flag != 'GEOMETRYCOLLECTION EMPTY':
                    rows.append(row)
            else:
                if data_format == 'BPS' and not bags_exist:
                    rows.append(row)
                elif data_format != 'BPS':
                    rows.append(row)
            if pb is not None:
                pb.SetValue(object_num + 1)
            object_num += 1
    for row in rows:
        if len(history) > 0:
            stored = [item for item in history if item['SURVEY_ID'] == row['SURVEY_ID']]
            if len(stored) > 0:
                for survey in stored:
                    new_flag = date_eval(row, survey)
                    if new_flag:
                        row['update'] = new_flag
                        if 'SURVEY_FILES' in stored:
                            row['SURVEY_FILES'] = stored['SURVEY_FILES']
    print('rows complete')
    return attr_list, rows


def write_geopackage(out_path: str, name: str, poly: str,
                     spcs: Union[str, int]):
    """
    Writes out a geopackage containing the bounding geometry of
    of the given query.

    Derived from:
    https://gis.stackexchange.com/a/52708/8104
    via
    https://gis.stackexchange.com/q/217165

    Parameters
    ----------
    out_path :
        String representing the complete file path for the output geopackage
    name :
        String representing the name of the survey; Used to name the layer
    poly :
        The WTK Multipolygon object that holds the survey's bounding data
    spcs :
        The ESPG code for the data

    """
    # Reference
    if type(spcs) == str:
        proj = osr.SpatialReference(wkt=spcs)
        proj.MorphFromESRI()
    else:
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(spcs)

    # Now convert it to a geopackage with OGR
    driver = ogr.GetDriverByName('GPKG')
    ds = driver.CreateDataSource(out_path)
    layer = ds.CreateLayer(name, proj, ogr.wkbMultiPolygon)
    #    layer = ds.CreateLayer(name, None, ogr.wkbMultiPolygon)
    # Add one attribute
    layer.CreateField(ogr.FieldDefn('Survey', ogr.OFTString))
    defn = layer.GetLayerDefn()

    # If there are multiple geometries, put the "for" loop here

    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField('Survey', name)

    # Make a geometry, from wkt object
    geom = ogr.CreateGeometryFromWkt(poly)

    feat.SetGeometry(geom)

    layer.CreateFeature(feat)

    linear_geom = geom.GetLinearGeometry()
    geojson = linear_geom.ExportToJson()

    # Save and close everything
    del ds, layer, feat, geom

    return geojson


def link_grab(source_url: str, extensions: list) -> list:
    """
    Scrapes an FTP site for file download links with the given extentsions

    Parameters
    ----------
    source_url
        The url at which download links are expected
    extensions
        The list of file extensions to search for from the download links
        displayed on the page

    Returns
    -------
    list
        The list of applicable download links found on the page

    """

    hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
           'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
           'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
           'Accept-Encoding': 'none',
           'Accept-Language': 'en-US,en;q=0.8',
           'Connection': 'keep-alive'}
    req = requests.get(source_url, headers=hdr)
    page = req.text
    file_links = []
    for extension in extensions:
        links = [link.strip('"') for link in re.findall(f'".*{extension}"', page)]
        file_links.extend([f'{source_url}/{link}' for link in links if link != ''])
    valid_links = []
    for link in file_links:
        basename, ext = os.path.splitext(link)
        root, basename = os.path.split(basename)
        add_link = True
        if 'combined' in basename.lower() or 'ellipsoid' in basename.lower():
            add_link = False
        elif ext.lower() in ('.gz') and os.path.splitext(basename)[1] not in extensions:
            add_link = False
        if '.bag' in extensions and re.compile(r'[a-z][0-9]{5}\_[a-z]{2}', re.IGNORECASE).search(basename) is None:
            add_link = False
        elif '.tif' in extensions and re.compile(r'[a-z][0-9]{5}(_SSSAB)', re.IGNORECASE).search(basename) is None:
            add_link = False
        if add_link:
            valid_links.append(link)

    return valid_links


def file_downloader(folder: str, download_links: list, saved_files: list) -> list:
    """
    Downloads files if they do not already exist in the input folder location.

    Parameters
    ----------
    folder
        The target folder for downloaded files
    download_links
        The files to be downloaded
    saved_files
        The files previously downloaded

    Returns
    -------
    list
        The list of files downloaded by this function

    """

    saved_links = []
    for link in download_links:
        saved = os.path.join(folder, os.path.split(link)[1])
        dwntime = datetime.datetime.now()
        while True:
            if os.path.exists(saved) and (saved not in saved_files or saved not in saved_links):
                basename, ext = os.path.splitext(saved)
                if ext in ('.gz') and not os.path.exists(basename):
                    unzip = gzip.open(saved, 'rb')
                    save_obj = open(basename, 'wb')
                    shutil.copyfileobj(unzip, save_obj)
                    unzip.close(), save_obj.close()
                    os.remove(saved)
                    saved_links.append(basename)
                    print(link)
                elif ext in ('.gz') and os.path.exists(basename):
                    pass
                else:
                    saved_links.append(saved)
                    print(link)
                break
            elif not os.path.exists(saved):
                try:
                    urllib.request.urlretrieve(link, saved)
                except socket.timeout:
                    urllib.request.urlretrieve(link, saved)
                except urllib.error.HTTPError as e:
                    print(f'e \n{link} {e}')
                    break
                except urllib.error.URLError as e:
                    print(f'e \n{link} {e}')
                    break
            elif datetime.datetime.now() - dwntime > datetime.timedelta(seconds=265):
                print(f'e \n{link} ')
                break
            else:
                print(f'e \n{link}')
                break
    return saved_links


def survey_download(rows: [dict], region: dict) -> [dict]:
    """
    Downloads and stores file information to the input dict

    Parameters
    ----------
    rows
        A list of survey metadata items (dict)
    region
        A dictionary of items that describes a NBS region

    Returns
    -------
    list of dict
        A list of survey metadata items updated with any actions taken by this
        function

    """
    data_type, data_format = config_data_type()

    if data_type == 'nos_hydro_dynamic' and data_format == 'BAG':
        ncei_head = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        branch_tail = r"\NOAA_NCEI_OCS\BAGs\Original"
        ftp_folders = {'BAG': ['.bag', '.gz'], 'TIFF': ['.tif', '.tiff', '.tfw', '.gz']}
    elif data_type == 'nos_hydro_dynamic' and data_format == 'BPS':
        ncei_head = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
        branch_tail = r"\NOAA_NCEI_OCS\BPS\Original"
        ftp_folders = {'GEODAS': ['.a93', 'htm', '.xyz', '.gz']}
    elif data_type == 'multibeam_dynamic':
        ncei_head = 'https://data.ngdc.noaa.gov/ocean/platforms/ships/'
        branch_tail = r"\NOAA_NCEI_MBBDB\MBs\Original"
        ftp_folders = {'multibeam/data/version1/products': ['.xyz', '.gz']}

    for row in rows:
        row['PROCESSING_REGION'] = f"{region['Processing Branch']}_{region['Region']}"

        if rows.index(row) == 0 and config['Destination']['Structure'] == 'NBS':
            branch_path = fr"{row['PROCESSING_REGION']}{branch_tail}"
            download_to = os.path.join(downloads, branch_path)
            if not os.path.isdir(download_to):
                os.makedirs(download_to)
        elif rows.index(row) == 0 and config['Destination']['Structure'] != 'NBS':
            download_to = downloads

        if 'SURVEY_FILES' not in row:
            row['SURVEY_FILES'] = []

        if 'update' in row and len(row['SURVEY_FILES'] != 0):
            for survey_file in row['SURVEY_FILES']:
                try:
                    os.remove(survey_file)
                except (OSError, FileNotFoundError) as e:
                    pass
                row['SURVEY_FILES'].remove(survey_file)
            del row['update']

        survey = row['SURVEY_ID']
        ncei_sub = '/'.join(os.path.splitext(row['DOWNLOAD_URL'])[0].split('/')[-2:])
        if data_type == 'multibeam_dynamic':
            ncei_sub = re.sub(r"_mb", "", ncei_sub)
        download_links = []
        print(f'Survey - {rows.index(row) + 1} of {len(rows)}: {survey}')
        survey_folder = os.path.join(download_to, survey)
        for folder, extensions in ftp_folders.items():
            source_url = f'{ncei_head}/{ncei_sub}/{folder}'
            download_links.extend(link_grab(source_url, extensions))

            if len(download_links) > 0:
                if not os.path.isdir(survey_folder):
                    os.mkdir(survey_folder)
                row['SURVEY_FILES'].extend(file_downloader(survey_folder, download_links, row['SURVEY_FILES']))

        if os.path.isdir(survey_folder):
#            if 'poly' in row:
#                geopackage_name = os.path.join(survey_folder, f'{survey}.gpkg')
#                write_geopackage(geopackage_name, survey, row[poly], 4362)
#                row['SURVEY_FILES'].append(geopackage_name)

            if len(row['SURVEY_FILES']) > 0:
                survey_files = [file_name for file_name in row['SURVEY_FILES'] if os.path.splitext(file_name)[1] in ('.a93', '.bag', '.xyz')]
                for data_file in survey_files:
                    base = os.path.basename(data_file)
                    name, ext = os.path.splitext(base)
                    pickle_name = f'{os.path.join(survey_folder, name)}.pickle'
                    row['SURVEY_FILES'].extend([pickle_name])

                    with open(pickle_name, 'wb') as metafile:
                        pickle.dump(row, metafile)
                        metafile.close()

    return rows


def region_info_json() -> [dict]:
    """
    Returns
    -------
    [dict]
        Returns the dictionary entries for the Branches defined in the config
        file

    """
    region_file = os.path.join(progLoc, config['Reference']['NBS'])
    with open(region_file) as json_file:
        regions = json.load(json_file)
        branches = [branch.strip() for branch in config['Processing Branches']['Branches'].split(',')]
        json_file.close()
        return [region for region in regions if region['Processing Branch'] in branches]


def survey_list() -> [dict]:
    """
    Returns
    -------
    [dict]
        Returns the dictionary entries for the survey query history from the
        .json file defined in the config file

    """

    surveys_file = os.path.join(progLoc, config['Reference']['History'])
    with open(surveys_file) as json_file:
        surveys = json.load(json_file)
        json_file.close()
        return surveys


def info_save(rows: [dict]):
    """
    Saves a list of survey info data (dict) to a .json file

    Parameters
    ----------
    rows
        A list of survey metadata items (dict)

    """
    surveys_file = os.path.join(progLoc, config['Reference']['History'])
    with open(surveys_file, 'w') as json_file:
        json.dump(rows, json_file)
        json_file.close()


def config_data_type() -> (str, str):
    """
    Reads the config value for data type to be queried by the tool

    Returns
    -------
    (str, str)
        query platform, data format

    """

    if config['Data']['Data Type'].upper() in ('', 'BAG'):
        return 'nos_hydro_dynamic', 'BAG'
    elif config['Data']['Data Type'].upper() in ('BPS'):
        return 'nos_hydro_dynamic', 'BPS'
    elif config['Data']['Data Type'].upper() in ('MB', 'MULTIBEAM'):
        return 'multibeam_dynamic', 'MB'
    else:
        raise ValueError(f"Data Type not supported: {config['Data']['Data Type']}")


def main(pb=None):
    """
    TODO write description

    Parameters
    ----------
    pb
        (Default value = None)

    """

    start = datetime.datetime.now()
    print(start)

    if pb is not None:
        pb.Pulse()

    regions = region_info_json()
    survey_history = survey_list()
    data_type, data_format = config_data_type()

    if data_format in ('BAG', 'MB'):
        qId = 0
    elif data_format in ('BPS'):
        qId = 1

    for region in regions:
        print(f"{region['Processing Branch']}_{region['Region']}")
        region_poly = os.path.join(parent_dir, region['Shape'])
        bounds = region_bounds(region_poly)
        objectIDs, bagNum = survey_objectID_query(bounds, qId=qId)
        if bagNum > 0:
            attr_list, rows = survey_compile(objectIDs, bagNum, survey_history, region_poly, qId=qId, pb=pb)
            if len(rows) > 0:
                rows = survey_download(rows, region)
                survey_history.extend(rows)
                info_save(survey_history)
            else:
                print(f'No changes in bag data')
        else:
            print(f'No new items were found within: {region}.')

    end = datetime.datetime.now()
    delta_time = end - start
    print(f'{end}\n{delta_time}')


if __name__ == '__main__':
    main()
