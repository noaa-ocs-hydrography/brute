# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 08:51:53 2019

@author: Casiano.Koprowski
"""


import configparser
import datetime
import gzip
import json
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
progLoc = os.getcwd()
config = configparser.ConfigParser()
config.read('config.ini')


if config['Destination']['Folder'] == '':
    downloads = os.path.join(progLoc, 'downloads')
else:
    downloads = config['Destination']['Folder']

if not os.path.isdir(os.path.join(downloads)):
    os.mkdir(downloads)

ncei_head = 'https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/'
# https://data.ngdc.noaa.gov/platforms/ocean/nos/coast/H12001-H14000/H12001/[BAG, TIFF]

zList = ['xmin', 'ymin', 'xmax', 'ymax']
attributes = {3: ['Name', 'SURVEY_ID', 'CELL_SIZE'],
              0: ['*']}
date_fields = ['DATE_SURVEY_BEGIN', 'DATE_SURVEY_END', 'DATE_MODIFY_DATA',
               'DATE_SURVEY_APPROVAL']


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
        bagList = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/{qId}/query' + \
                  f'?where={where}&text=&objectIds=&time=&geometry={box}&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=true&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
        bagListRequest = requests.get(bagList)
        bagListRequestJSON = bagListRequest.json()
#            print (bagListRequestJSON)
        if bagListRequestJSON['objectIds'] is not None:
            objectIDs.extend(bagListRequestJSON['objectIds'])
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


def survey_compile(objectIDs: list, num: int, history: [dict], qId=0, pb=None) -> (list, [dict]):
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

    id_chunks = list(list_chunks([str(objectID) for objectID in objectIDs]))

    for chunk in id_chunks:
        id_string = ','.join(chunk)
        query = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/{qId}/query' + \
                f'?where=&text=&objectIds={id_string}&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields={opts}&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
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
                                row[attribute] = f'{datetime.datetime.utcfromtimestamp(date / 1000):%Y-%m-%d}'
                            except OSError as e:
                                print(e, date)
                    elif page['features'][object_num]['attributes'][attribute] is not None:
                        row[attribute] = str(page['features'][object_num]['attributes'][attribute])
                rows.append(row)
            except KeyError as e:
                print(e, page)
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
    return file_links


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
                    saved_links.extend([saved, basename])
                else:
                    print(link)
                    saved_links.append(saved)
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

    for row in rows:
        row['PROCESSING_REGION'] = f"{region['Processing Branch']}_{region['Region']}"

        if rows.index(row) == 0 and config['Destination']['Structure'] == 'NBS':
            branch_path = fr"{row['PROCESSING_REGION']}\NOAA_NCEI_OCS\BAGs\Original"
            download_to = os.path.join(downloads, branch_path)
            if not os.path.isdir(download_to):
                os.makedirs(download_to)
        else:
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
        download_links = []
        print(f'Survey - {rows.index(row) + 1} of {len(rows)}: {survey}')

        for folder, extensions in {'BAG': ['.bag', '.bag.gz'],
                                   'TIFF': ['tif', 'tiff', 'tfw']}.items():
            source_url = f'{ncei_head}/{ncei_sub}/{folder}'
            download_links = link_grab(source_url, extensions)

            if len(download_links) > 0:
                survey_folder = os.path.join(download_to, survey)
                if not os.path.isdir(survey_folder):
                    os.mkdir(survey_folder)
                row['SURVEY_FILES'].extend(file_downloader(survey_folder, download_links, row['SURVEY_FILES']))

        pickle_name = f'{os.path.join(survey_folder, survey)}.pickle'
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


def main(pb=None):
    """
    TODO write description

    Parameters
    ----------
    pb
        (Default value = None)

    """

    if pb is not None:
        pb.Pulse()

    regions = region_info_json()
    survey_history = survey_list()

    for region in regions:
        region_poly = os.path.join(progLoc, region['Shape'])
        bounds = region_bounds(region_poly)
        objectIDs, bagNum = survey_objectID_query(bounds, 0)
        if bagNum > 0:
            attr_list, rows = survey_compile(objectIDs, bagNum, survey_history, pb=pb)
            if len(rows) > 0:
                rows = survey_download(rows, region)
                survey_history.extend(rows)
                info_save(survey_history)
            else:
                print(f'No changes in bag data')
        else:
            print(f'No new items were found within: {region}.')


if __name__ == '__main__':
    main()
