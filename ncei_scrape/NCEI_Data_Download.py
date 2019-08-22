# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 08:51:53 2019

@author: Casiano.Koprowski
"""


import configparser
import csv
import datetime
import os
import pickle
import re
import shutil
import socket
import urllib
import zipfile
from typing import Union, Dict, List

import requests

from osgeo import gdal
from osgeo import ogr
from osgeo import osr


"""Known global constants"""
progLoc = os.getcwd()
config = configparser.ConfigParser()
config.read('config.ini')


if config['Downloads']['Folder'] == '':
    downloads = os.path.join(progLoc, 'downloads')
else:
    downloads = config['Downloads']['Folder']

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

    bagList = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/{qId}/query' + \
              f'?where=&text=&objectIds=&time=&geometry={bounds}&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=true&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
    bagListRequest = requests.get(bagList)
    bagListRequestJSON = bagListRequest.json()
    #    print (bagListRequestJSON)
    objectIDs = bagListRequestJSON['objectIds']
    if objectIDs is None:
        return [], 0
    else:
        objectNum = len(objectIDs) - 1
        print(objectIDs, objectNum)
        return objectIDs, objectNum


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
    for objectID in objectIDs:
        print(x, end=' ')
        query = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/{qId}/query' + \
                f'?where=&text=&objectIds={objectID}&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields={opts}&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
        response = requests.get(query)
        page = response.json()
        if x == 0:
            fields = page['fields']
            for attr in fields:
                if attr['name'] == 'OBJECTID':
                    pass
                else:
                    attr_list.append(attr['name'])
        row = {}
        try:
            for attribute in attr_list:
                if page['features'][0]['attributes'][attribute] is None:
                    row[attribute] = 'null'
                elif attribute in date_fields:

                    if page['features'][0]['attributes'][attribute] is None:
                        row[attribute] = 'null'
                    else:
                        date = (page['features'][0]['attributes'][attribute])

                        try:
                            row[attribute] = f'{datetime.datetime.utcfromtimestamp(date / 1000):%Y-%m-%d}'
                        except OSError as e:
                            print(e, date)
                            row[attribute] = 'error'
                else:
                    row[attribute] = str(page['features'][0]['attributes'][attribute])
            rows.append(row)
        except KeyError as e:
            print(e, page)
        #            break
        if pb is not None:
            pb.SetValue(x)
        x += 1
#    print(len(rows))
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
        file_links.extend([f'{source_url}{link}' for link in links if link != ''])
    return file_links


def fileDownloader(folder, download_links):
    for link in download_links:
        saved = os.path.join(folder, os.path.split(link)[1])
        dwntime = datetime.datetime.now()
        while True:
            if os.path.exists(saved):
                break
            elif not os.path.exists(saved):
                try:
                    urllib.request.urlretrieve(link, saved)
                except socket.timeout:
                    urllib.request.urlretrieve(link, saved)
                except urllib.error.HTTPError as e:
                    print(f'e \n{link} {e}')
                    download_links.remove(link)
                    break
                except urllib.error.URLError as e:
                    print(f'e \n{link} {e}')
                    download_links.remove(link)
                    break
            elif datetime.datetime.now() - dwntime > datetime.timedelta(seconds=295):
                print(f'e \n{link} ')
                download_links.remove(link)
                break
            else:
                print(f'e \n{link}')
                download_links.remove(link)
                break
    return download_links


def surveyDownload(rows):
    for row in rows:
        survey = row['SURVEY_ID']
        ncei_sub = '/'.join(os.path.splitext(row['DOWNLOAD_URL'])[0].split('/')[-2:])
        download_links = []
        for folder, extensions in {'BAG': ['.bag', '.bag.gz'],
                                   'TIFF': ['tif', 'tiff']}:
            source_url = f'{ncei_head}/{ncei_sub}/{folder}'
            download_links = linkGrab(source_url, extensions)
            if len(download_links) > 0:
                survey_folder = os.path.join(downloads, survey)
                if not os.path.isdir(survey_folder):
                    os.mkdir(survey_folder)
                row['SURVEY_FILES'] = fileDownloader(survey_folder, download_links)

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
