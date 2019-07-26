# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 13:34:27 2019

@author: Casiano.Koprowski
"""

import csv
import datetime
import os

import requests

'''Known global constants'''
# progLoc is the program's own file location / current working directory (cwd)
progLoc = os.getcwd()

zList = ['xmin', 'ymin', 'xmax', 'ymax']
attributes = {3: ['Name', 'SURVEY_ID', 'CELL_SIZE', 'DOWNLOAD_URL'],
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

    nxStr, nyStr, sxStr, syStr = str(nx), str(ny), str(sx), str(sy)
    corner = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/Utilities/Geometry/GeometryServer/project' + \
             f'?inSR=4326&outSR=102100&geometries=%7B"geometryType"+%3A+"esriGeometryPoint"%2C+"geometries"+%3A+%5B%0D%0A+++++%7B%0D%0A+++++++"x"+%3A+{nxStr}%2C%0D%0A+++++++"y"+%3A+{syStr}%0D%0A+++++%7D%2C%7B%0D%0A+++++++"x"+%3A+{sxStr}%2C%0D%0A+++++++"y"+%3A+{nyStr}%0D%0A+++++%7D%0D%0A++%5D%0D%0A%7D&transformation=&transformForward=true&vertical=false&f=json'
    cornerRequest = requests.get(corner)
    cornerRequestJSON = cornerRequest.json()
    #    print (cornerRequestJSON)
    bounds = []
    z = 0
    for i in range(len(cornerRequestJSON['geometries'])):
        for k, j in cornerRequestJSON['geometries'][i].items():
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


def surveyCompile(surveyIDs, num, qId=3, pb=None):
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
    for num in surveyIDs:
        print(x, end=' ')
        bagID = str(num)
        query = f'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/{qId}/query' + \
                f'?where=&text=&objectIds={bagID}&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields={opts}&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
        response = requests.get(query)
        page = response.json()
        if x == 0:
            fields = page['fields']
            for attr in fields:
                if attr['name'] == 'OBJECTID':
                    pass
                else:
                    attr_list.append(attr['name'])
        row = []
        try:
            for attribute in attr_list:
                if page['features'][0]['attributes'][attribute] is None:
                    row.append('null')
                elif attribute in date_fields:

                    if page['features'][0]['attributes'][attribute] is None:
                        row.append('null')
                    else:
                        date = (page['features'][0]['attributes'][attribute])

                        try:
                            date = datetime.datetime.utcfromtimestamp(date / 1000)
                            row.append(str(date.strftime('%Y-%m-%d')))
                        except OSError as e:
                            print(e, date)
                            row.append('error')
                else:
                    row.append(str(page['features'][0]['attributes'][attribute]))
            rows.append(row)
        except KeyError as e:
            print(e, page)
        #            break
        if pb is not None:
            pb.SetValue(x)
        x += 1
    print(len(rows))
    print('rows complete')
    return attr_list, rows


def csvWriter(attr_list, csvFile, csvLocation, name, pb=None):
    """
    TODO write description

    Parameters
    ----------
    csvFile
    csvLocation
    name
    pb
         (Default value = None)

    Returns
    -------

    """
    if name == '':
        num = 0
        name = f'{datetime.datetime.now():%Y%m%d}_NCEI_Output'
        while True:
            if not os.path.exists(f'{name}_{num}.txt'):
                name = f'{name}_{num}'
                break
            else:
                num += 1
    name = os.path.join(csvLocation, f'{name}.txt')
    csvOpen = open(name, 'w', newline='')
    save = csv.writer(csvOpen, delimiter=',')
    save.writerow(attr_list)
    x = 0
    if pb is not None:
        pb.SetRange(len(csvFile))
    for row in csvFile:
        save.writerow(row)
        if pb is not None:
            pb.SetValue(x)
        x += 1
    csvOpen.close()


def main(name, nx, sy, sx, ny, qId=3, pb=None):
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
    print(name, nx, sy, sx, ny)

    if pb is not None:
        pb.Pulse()

    if qId == 3:
        noItems = 'BAG Files'
    else:
        noItems = 'Surveys'

    bounds = coordQuery(nx, ny, sx, sy)
    bagIDs, bagNum = bagIDQuery(bounds, qId)
    if bagNum > 0:
        attr_list, rows = surveyCompile(bagIDs, bagNum, qId, pb)
        csvWriter(attr_list, rows, progLoc, name, pb)
    else:
        cardinal_directions = {'North': ny, 'West': sx, 'South': sy, 'East': nx}
        return (f'No {noItems} were found within: {cardinal_directions}.')
