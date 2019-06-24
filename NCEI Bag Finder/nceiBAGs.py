# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 13:34:27 2019

@author: Casiano.Koprowski
"""

import csv
import os

import requests

'''Known global constants'''
# progLoc is the program's own file location / current working directory (cwd)
progLoc = os.getcwd()

zList = ['xmin', 'ymin', 'xmax', 'ymax']
attributes = ['Name','SURVEY_ID', 'CELL_SIZE', 'DOWNLOAD_URL',]

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
    corner = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/Utilities/Geometry/GeometryServer/project?inSR=4326&outSR=102100&geometries=%7B"geometryType"+%3A+"esriGeometryPoint"%2C+"geometries"+%3A+%5B%0D%0A+++++%7B%0D%0A+++++++"x"+%3A+' + nxStr + '%2C%0D%0A+++++++"y"+%3A+' + syStr + '%0D%0A+++++%7D%2C%7B%0D%0A+++++++"x"+%3A+' + sxStr + '%2C%0D%0A+++++++"y"+%3A+' + nyStr + '%0D%0A+++++%7D%0D%0A++%5D%0D%0A%7D&transformation=&transformForward=true&vertical=false&f=json'
    cornerRequest = requests.get(corner)
    cornerRequestJSON = cornerRequest.json()
#    print (cornerRequestJSON)
    bounds = []
    z = 0
    for i in range(len(cornerRequestJSON['geometries'])):
        for k, j in cornerRequestJSON['geometries'][i].items():
            bounds.append((zList[z],j))
            z += 1
    bounds = dict(bounds)
#    print (bounds)
    return bounds

def bagIDQuery(bounds):
    """
    

    Parameters
    ----------
    bounds :
        

    Returns
    -------

    """
    bounds = str(bounds)
    bagList = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/3/query?where=&text=&objectIds=&time=&geometry=' + bounds + '&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=true&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
    bagListRequest = requests.get(bagList)
    bagListRequestJSON = bagListRequest.json()
#    print (bagListRequestJSON)
    objectIDs = bagListRequestJSON['objectIds']
    objectNum = len(objectIDs) - 1
    print (objectIDs, objectNum)
    return objectIDs, objectNum

def surveyCompile(surveyIDs,num,pb=None):
    """
    

    Parameters
    ----------
    surveyIDs :
        
    num :
        
    pb :
         (Default value = None)

    Returns
    -------

    """
    x = 0
    rows = []
    if pb != None:
        pb.SetRange(num)
    for num in surveyIDs:
        print (x, end=' ')
        bagID = str(num)
        query = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/nos_hydro_dynamic/MapServer/3/query?where=&text=&objectIds=' + bagID + '&time=&geometry=&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=Name,SURVEY_ID,CELL_SIZE,DOWNLOAD_URL&returnGeometry=false&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&f=json'
        response = requests.get(query)
        page = response.json()
        row = []
        try:
            for attribute in attributes:
                if page['features'][0]['attributes'][attribute] == None:
                    row.append('null')
                else:
                    row.append(str(page['features'][0]['attributes'][attribute]))
            rows.append(row)
        except KeyError as e:
            print (e, page)
#            break
        if pb != None:
            pb.SetValue(x)
        x += 1
    print (len(rows))
    print ('rows complete')
    return rows

def csvWriter(csvFile, csvLocation, name, pb=None):
    """
    

    Parameters
    ----------
    csvFile :
        
    csvLocation :
        
    name :
        
    pb :
         (Default value = None)

    Returns
    -------

    """
    name = csvLocation + '\\' + name + '.txt'
    csvOpen = open(name, 'w', newline='')
    save = csv.writer(csvOpen, delimiter = ',')
    save.writerow(attributes)
    x = 0
    if pb != None:
        pb.SetRange(len(csvFile))
    for row in csvFile:
        save.writerow(row)
        if pb != None:
            pb.SetValue(x)
        x += 1
    csvOpen.close()

def main(name,nx,sy,sx,ny,pb=None):
    """
    

    Parameters
    ----------
    name :
        
    nx :
        
    sy :
        
    sx :
        
    ny :
        
    pb :
         (Default value = None)

    Returns
    -------

    """
    print(name,nx,sy,sx,ny)
    if pb != None:
        pb.Pulse()
    bounds = coordQuery(nx,ny,sx,sy)
    bagIDs, bagNum = bagIDQuery(bounds)
    rows = surveyCompile(bagIDs,bagNum,pb)
    csvWriter(rows,progLoc,name,pb)
    return True
