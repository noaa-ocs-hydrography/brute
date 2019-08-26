# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 08:51:53 2019

@author: Casiano.Koprowski
"""


import configparser
import csv
import datetime
import gzip
import os
import pickle
import re
import shutil
import socket
import urllib
import zipfile
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


def coordQuery(nx, ny, sx, sy):
    """


    Parameters
    ----------
    nx :

    ny :

    sx :

    sy :


    Returns
    -------

    """

    coords = {'geometryType': 'esriGeometryPoint', 'geometries': [{'x': nx, 'y': sy}, {'x': sx, 'y': ny}]}
    query = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/Utilities/Geometry/GeometryServer/project' + \
            f'?inSR=4326&outSR=102100&geometries={coords}&transformation=&transformForward=true&vertical=false&f=json'

    coordsRequest = requests.get(query)
    coordsRequestJSON = coordsRequest.json()
    #    print (coordsRequestJSON)
    bounds = []
    z = 0
    for i in range(len(coordsRequestJSON['geometries'])):
        for k, j in coordsRequestJSON['geometries'][i].items():
            bounds.append((zList[z], j))
            z += 1
    bounds = dict(bounds)
    #    print (bounds)
    return bounds


def regionQuery(region_file: str):
    """


    Parameters
    ----------
    nx :

    ny :

    sx :

    sy :


    Returns
    -------

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
            coords = {'geometryType': 'esriGeometryPoint', 'geometries': [{'x': min(x_points), 'y': min(y_points)}, {'x': max(x_points), 'y': max(y_points)}]}
            query = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/Utilities/Geometry/GeometryServer/project' + \
                    f'?inSR=4326&outSR=102100&geometries={coords}&transformation=&transformForward=true&vertical=false&f=json'

            coordsRequest = requests.get(query)
            coordsRequestJSON = coordsRequest.json()
            box = []
            z = 0
            for i in range(len(coordsRequestJSON['geometries'])):
                for k, j in coordsRequestJSON['geometries'][i].items():
                    box.append((zList[z], j))
                    z += 1
            return_bounds.append(dict(box))
    else:
        x_points, y_points = bounds[0]

        coords = {'geometryType': 'esriGeometryPoint', 'geometries': [{'x': min(x_points), 'y': min(y_points)}, {'x': max(x_points), 'y': max(y_points)}]}
        query = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/Utilities/Geometry/GeometryServer/project' + \
                f'?inSR=4326&outSR=102100&geometries={coords}&transformation=&transformForward=true&vertical=false&f=json'

        coordsRequest = requests.get(query)
        coordsRequestJSON = coordsRequest.json()
        box = []
        z = 0
        for i in range(len(coordsRequestJSON['geometries'])):
            for k, j in coordsRequestJSON['geometries'][i].items():
                box.append((zList[z], j))
                z += 1
        return_bounds.append(dict(box))
    return return_bounds

#regionQuery(r'R:\Scripts\vlab-nbs\ncei_scrape\MCD_ProcessingBranches\MCD_PBG_WGS84.gpkg')

def bagIDQuery(bounds, qId=3):
    """
    TODO write description

    Parameters
    ----------
    bounds
    qId


    Returns
    -------

    """

    # Renamed bounds because Glen wanted it to be a que, not a queue
    que = bounds

    objectIDs = []
    objectNum = 0

    for box in que:
        bagList = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/{qId}/query' + \
                  f'?where=&text=&objectIds=&time=&geometry={box}&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=true&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
        bagListRequest = requests.get(bagList)
        bagListRequestJSON = bagListRequest.json()
        #    print (bagListRequestJSON)
        objectIDs.extend(bagListRequestJSON['objectIds'])
        objectNum += len(objectIDs) - 1
        print(objectNum)

#    print(objectIDs, objectNum)
    return objectIDs, objectNum


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


def surveyCompile(objectIDs, num, qId=3, pb=None):
    """
    TODO write description

    Parameters
    ----------
    surveyIDs
    num
    qId

    Returns
    -------

    """
    x = 0
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
            print(x, end=' ')

            row = {}
            try:
                for attribute in attr_list:
                    if page['features'][object_num]['attributes'][attribute] is None:
                        row[attribute] = 'null'
                    elif attribute in date_fields:

                        if page['features'][object_num]['attributes'][attribute] is None:
                            row[attribute] = 'null'
                        else:
                            date = (page['features'][object_num]['attributes'][attribute])

                            try:
                                row[attribute] = f'{datetime.datetime.utcfromtimestamp(date / 1000):%Y-%m-%d}'
                            except OSError as e:
                                print(e, date)
                                row[attribute] = 'error'
                    else:
                        row[attribute] = str(page['features'][object_num]['attributes'][attribute])
                rows.append(row)
            except KeyError as e:
                print(e, page)
            if pb is not None:
                pb.SetValue(object_num + 1)
            object_num += 1
    print('rows complete')
    return attr_list, rows


def linkGrab(source_url, extensions):
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


def fileDownloader(folder, download_links, saved_files):
    saved_links = []
    for link in download_links:
        saved = os.path.join(folder, os.path.split(link)[1])
        dwntime = datetime.datetime.now()
        while True:
            if os.path.exists(saved) and saved not in saved_files or saved not in saved_links:
                basename, ext = os.path.splitext(saved)
                if ext not in ('.gz') and not os.path.exists(basename):
                    unzip = gzip.open(saved, 'rb')
                    saved = open(basename, 'wb')
                    shutil.copyfileobj(unzip, saved)
                    unzip.close(), saved.close()
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


def surveyDownload(rows: [dict]) -> [dict]:
    """
    Downloads and stores file information to the input dict

    Parameters
    ----------
    rows : list of dict
        A list of survey metadata items

    """
    for row in rows:
        if rows.index(row) == 0:
            branch_path = fr"{row['PROCESSING_BRANCH']}\NOAA_NCEI_OCS\BAGs\Original"
            branch_folder = os.path.join(downloads, branch_path)
            if not os.path.isdir(branch_folder):
                os.makedirs(branch_folder)

        row['SURVEY_FILES'] = []
        survey = row['SURVEY_ID']
        ncei_sub = '/'.join(os.path.splitext(row['DOWNLOAD_URL'])[0].split('/')[-2:])
        download_links = []
        print(f'Survey - {rows.index(row) + 1} of {len(rows)}: {survey}')

        for folder, extensions in {'BAG': ['.bag', '.bag.gz'],
                                   'TIFF': ['tif', 'tiff', 'tfw']}.items():
            source_url = f'{ncei_head}/{ncei_sub}/{folder}'
            download_links = linkGrab(source_url, extensions)

            if len(download_links) > 0:
                survey_folder = os.path.join(branch_folder, survey)
                if not os.path.isdir(survey_folder):
                    os.mkdir(survey_folder)
                row['SURVEY_FILES'].extend(fileDownloader(survey_folder, download_links, row['SURVEY_FILES']))

        pickle_name = f'{os.path.join(survey_folder, survey)}.pickle'


        with open(pickle_name, 'wb') as metafile:
            pickle.dump(row, metafile)
            row['SURVEY_FILES'].extend([pickle_name])

    return rows


#def csvWriter(attr_list, csvFile, csvLocation, name, pb=None):
#    """
#    TODO write description
#
#    Parameters
#    ----------
#    csvFile
#    csvLocation
#    name
#    pb
#         (Default value = None)
#
#    Returns
#    -------
#
#    """
#    if name == '':
#        num = 0
#        name = f'{datetime.datetime.now():%Y%m%d}_NCEI_Output'
#        while True:
#            if not os.path.exists(f'{name}_{num}.txt'):
#                name = f'{name}_{num}'
#                break
#            else:
#                num += 1
#    name = os.path.join(csvLocation, f'{name}.txt')
#    csvOpen = open(name, 'w', newline='')
#    save = csv.writer(csvOpen, delimiter=',')
#    save.writerow(attr_list)
#    x = 0
#    if pb is not None:
#        pb.SetRange(len(csvFile))
#    for row in csvFile:
#        save.writerow(row)
#        if pb is not None:
#            pb.SetValue(x)
#        x += 1
#    csvOpen.close()


def main(name, nx, sy, sx, ny, qId=0, pb=None):
    """
    TODO write description

    Parameters
    ----------
    name
    nx
    sy
    sx
    ny
    pb
         (Default value = None)

    Returns
    -------

    """

    if pb is not None:
        pb.Pulse()

    if qId == 3:
        noItems = 'BAG Files'
    else:
        noItems = 'Surveys'

    bounds = coordQuery(nx, ny, sx, sy)
    objectIDs, bagNum = bagIDQuery(bounds, qId)
    if bagNum > 0:
        attr_list, rows = surveyCompile(objectIDs, bagNum, qId, pb)
        rows = surveyDownload(rows)
        return rows
#        csvWriter(attr_list, rows, progLoc, name, pb)
    else:
        cardinal_directions = {'North': ny, 'West': sx,
                               'South': sy, 'East': nx}
        return (f'No {noItems} were found within: {cardinal_directions}.')

#main('', -73.5, 40, -74, 40.5, qId=0)
