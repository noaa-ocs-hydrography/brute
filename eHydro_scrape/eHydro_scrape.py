# -*- coding: utf-8 -*-
"""
Created on Tue Oct 2 13:44:14 2018

Last Modified: Apr 19 12:12:42 2019

@author: Casiano.Koprowski <casiano.koprowski@noaa.gov>
"""

import configparser
import csv
import datetime
import os
import pickle
import re
import socket
import time
import urllib
import zipfile
from typing import Tuple, List, Union, TextIO, Any

import numpy as np
import requests
from osgeo import osr, ogr

__version__ = '1.2'

"""Known global constants"""
# print (datetime.datetime.now().strftime('%b %d %X %Y'))
# progLoc = '\\eHydro_scrape'
progLoc = os.getcwd()
"""progLoc is the program's own file location / current working directory (cwd)
obtained by :func:`os.getcwd()`"""
full = re.compile(r'FULL.xyz', re.IGNORECASE)
"""regex object for searching zipfile contents for high-res data ending in
``FULL.xyz``
"""
full_a = re.compile(r'_A.xyz', re.IGNORECASE)
"""regex object for searching zipfile contents for high-res data ending in
``_A.xyz``
"""
# location of data items
config = configparser.ConfigParser()
config.read('config.ini')

csvName = 'eHydro_csv.txt'
"""Default name for csv.txt output"""
csvLocation = os.path.join(progLoc, csvName)
# csvLocation = progLoc + '\\' + csvName
"""Default location for :attr:`csvName`"""
logName = 'eHydro_log.txt'
"""Default name for log.txt output"""
logLocation = os.path.join(progLoc, logName)
# logLocation = progLoc + '\\' + logName
"""Default location for :attr:`logName`"""
# """Default location for """
holding = progLoc + '\\downloads\\'
"""Default location for all downloaded data, regardless of the type of query
performed. Data is broken up into folders representing each district (ex.
``\\downloads\\CEMVN``, ``\\downloads\\CENWP``, etc.)
"""
logging = progLoc + '\\logs\\'
"""Default location for individual query logs. These are named like
``YYYYMMDD_0_eHydro_log.txt``
"""
running = progLoc + '\\runs\\'
"""Default location for individual query csv outputs. These are named like
``YYYYMMDD_0_eHydro_csv.txt``
"""
# eHydro survey entry attributes
attributes = ["OBJECTID", "SURVEYJOBIDPK", "SURVEYAGENCY", "CHANNELAREAIDFK",
              "SDSFEATURENAME", "SOURCEPROJECTION", "SOURCEDATALOCATION",
              "SURVEYDATEUPLOADED", "SURVEYDATEEND", "SURVEYDATESTART",
              "SURVEYTYPE", "PROJECTEDAREA", "SOURCEDATAFORMAT",
              "Shape__Area", "Shape__Length"]
"""The specific attributes queried for each survey in :func:`surveyCompile`"""

# check to see if the downloaded data folder exists, will create it if not
if not os.path.exists(holding):
    os.mkdir(holding)
if not os.path.exists(logging):
    os.mkdir(logging)
if not os.path.exists(running):
    os.mkdir(running)


def query() -> Tuple[List[str], int, str]:
    """
    Holds the Queries for the eHydro REST API, asks for responses, and uses
    the json library to make it readable by the program. Returns the json
    responses for number of surveys and surveys in the response.

    -REMOVED- It also saves a prettyprinted version of the response as a text
    file.

    The funtion uses the requests library to retrieve the API's response(s) and
    uses the json library to make them readable by the program.

    The function returns the json object containing the contents of the query
    and the integer number of surveys contained by the query

    :returns: List of survey ids from query, total number of surveys returned by the query, and a string containing the parameters gathered from the config file
    """

    # Today (ex. '2018-08-08'), unformatted
    today = datetime.datetime.today()
    strToday = str(today.strftime('%Y-%m-%d'))
    # Today - 1 (ex. '2018-08-06'), unformatted
    yesterday = today - datetime.timedelta(1)
    strYesterday = str(yesterday.strftime('%Y-%m-%d'))
    if config['Timeframe']['Start Date'] != '':
        start = config['Timeframe']['Start Date']
    else:
        start = strYesterday
    if config['Timeframe']['End Date'] != '':
        end = config['Timeframe']['End Date']
    else:
        end = strToday
    if (config['Agencies']['Only Listed'] == 'yes'
            and config['Agencies']['Agencies'] != ''):
        areas = ''
        agencies = config['Agencies']['Agencies'].split(',')
        if len(agencies) > 1:
            i = 0
            areas += '('
            while i < len(agencies):
                if i == 0:
                    # %27 is ' (Apostrophe); %25 is % (Percent sign)
                    areas += ('UPPER(SURVEYAGENCY)+like+%27%25'
                              + agencies[0].strip()
                              + '%25%27')
                    i += 1
                else:
                    # %27 is ' (Apostrophe); %25 is % (Percent sign)
                    areas += ('+OR+UPPER(SURVEYAGENCY)+like+%27%25'
                              + agencies[i].strip()
                              + '%25%27')
                    i += 1
            areas += ')+'
        else:
            # %27 is ' (Apostrophe); %25 is % (Percent sign)
            areas = ('UPPER(SURVEYAGENCY)+like+%27%25'
                     + agencies[0].strip()
                     + '%25%27')
    else:
        areas = ''

    # The main query parameters that will determine the contents of the response
    # Survey Date Uploaded
    if config['Timeframe']['Ignore Date'] == 'no' and areas != '':
        print('\nStart:', start, '\nEnd:', end)
        # %27 is ' (Apostrophe)
        where = ('SURVEYDATEUPLOADED+>=+%27'
                 + start
                 + 'T00:01:00.000Z%27+AND+SURVEYDATEUPLOADED+<=+%27'
                 + end
                 + 'T11:59:00.000Z%27+AND+'
                 + areas)
    elif config['Timeframe']['Ignore Date'] == 'no' and areas == '':
        print('\nStart:', start, '\nEnd:', end)
        # %27 is ' (Apostrophe)
        where = ('SURVEYDATEUPLOADED+>=+%27'
                 + start
                 + 'T00:01:00.000Z%27+AND+SURVEYDATEUPLOADED+<=+%27'
                 + end
                 + 'T11:59:00.000Z%27')
    else:
        print('\nStart: Ignored', '\nEnd: Ignored')
        start = 'Ignored'
        end = 'Ignored'
        if areas != '':
            # areas
            where = areas
        else:
            # 1 = 1; %3D is = (Equals sign)
            where = '1%3D1'

    # The query for determining how many responses will be returned
    newSurveys = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?where=' + where + '&outFields=*&returnGeometry=false&returnCountOnly=true&outSR=&f=json'

    # The query for returning the object IDs for the given timeframe
    objIDs = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?&where=' + where + '&outFields=*&returnGeometry=false&returnIdsOnly=true&outSR=&f=json'

    print(objIDs, newSurveys)

    # Initial Query execution
    surveyNumRequest = requests.get(newSurveys)
    surveyIDsRequest = requests.get(objIDs)

    # Response for number of surveys used as a json object
    surveyNumJSON = surveyNumRequest.json()
    print(surveyNumJSON)
    newSurveysNum = int(surveyNumJSON['count'])

    # Response for survey Object IDs as a json object
    surveyIDsJSON = surveyIDsRequest.json()
    surveyIDs = surveyIDsJSON['objectIds']
    print(surveyIDs, newSurveysNum)

    if areas == '':
        dist = 'none'
    else:
        dist = config['Agencies']['Agencies']

    paramString = str('\tParameters:\n\t\tStart Date: ' + start
                      + '\n\t\tEnd Date: ' + end
                      + '\n\t\tDistricts: ' + dist
                      + '\n\t\tQuery Only Districts: ' + config['Agencies']['Only Listed']
                      + '\n\t\tKeep All Data: ' + config['Resolutions']['Override'])

    return surveyIDs, newSurveysNum, paramString


def create_polygon(coords: List[Tuple[float, float]]) -> ogr.Geometry:
    """
    Creates an ogr.Geometry/wkbLinearRing object from a list of coordinates.

    with considerations from:
    https://gis.stackexchange.com/q/217165

    :param coords: A list of [x, y] points to be made into an ogr.Geometry/wkbLinearRing object
    :returns: ogr.Geometry/wkbLinearRing object
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
    transelated and returned as a WTK Multipolygon object.

    with considerations from:
    https://gis.stackexchange.com/q/217165
    and:
    https://pcjericks.github.io/py-gdalogr-cookbook/geometry.html#create-a-multipolygon

    :param polys: A list of ogr.Geometry/wkbLinearRing objects
    :returns: WTK Multipolygon object
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

    :param coordinates: A list of coordinate point 'rings' returned by an eHydro survey query in the Geometry attribute
    :returns: A WTK Multipolygon object representing the survey outline
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


def surveyCompile(surveyIDs: list, newSurveysNum: int, pb=None) -> list:
    """
    Uses the json object return of the each queried survey id and the total
    number of surveys included to compile a list of complete returned survey
    data, as provided in the response. The function also takes into account
    that the survey data returns for any date/time are returned as timestamps.
    The function looks for these fields and converts them to datetime objects
    and finally strings.

    The function returns the lists of returned survey data as a list 'rows'.
    The specific :attr:`attributes` for each survey are:

    - OBJECTID.
    - SDSFEATURENAME.
    - SURVEYTYPE.
    - CHANNELAREAIDFK.
    - SURVEYAGENCY.
    - SURVEYDATEUPLOADED.
    - SURVEYDATESTART.
    - SURVEYDATEEND.
    - SOURCEDATALOCATION.
    - SOURCEPROJECTION.
    - SURVEYJOBIDPK.
    - PROJECTEDAREA.
    - SURVEYTYPE.
    - SOURCEDATAFORMAT.
    - Shape__Area.
    - Shape__Length.

    Added to the end of this list but not included in the list for csv export
    is a dictionary of the same information and the survey outline/shape as a
    WTK Multipolygon object. This data is used to writa a metadata pickle
    output and a geopackage in func:`downloadAndCheck`

    :param surveyIDs: List of survey ids usualy generated by :func:`query`
    :param newSurveysNum: Total number of surveys usualy returned by :func:`query`
    :returns: A list compiled of the attributes for every survey in surveyIDs
    """

    x = 0
    rows = []
    if pb != None:
        pb.SetRange(newSurveysNum)
        pb.SetValue(x)
    while x < newSurveysNum:
        print(x, end=' ')
        query = (
                'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?where=OBJECTID+=+'
                + str(surveyIDs[x])
                + '&outFields=*&returnGeometry=true&outSR=4326&f=json')
        response = requests.get(query)
        page = response.json()
        row = []
        metadata = {}
        metadata['version'] = __version__
        for attribute in attributes:
            try:
                if page['features'][0]['attributes'][attribute] == None:
                    row.append('null')
                    metadata[attribute] = 'null'
                elif (attribute == "SURVEYDATEUPLOADED"
                      or attribute == "SURVEYDATEEND"
                      or attribute == "SURVEYDATESTART"):
                    if page['features'][0]['attributes'][attribute] == None:
                        row.append('null')
                        metadata[attribute] = 'null'
                    else:
                        date = (page['features'][0]['attributes'][attribute])
                        date = datetime.datetime.utcfromtimestamp(date / 1000)
                        row.append(str(date.strftime('%Y-%m-%d')))
                        metadata[attribute] = str(date.strftime('%Y-%m-%d'))
                else:
                    row.append(str(page['features'][0]['attributes'][attribute]))
                    metadata[attribute] = str(page['features'][0]['attributes'][attribute])
            except KeyError as e:
                print(e, page)
                row.append('error')
                metadata[attribute] = 'error'
        try:
            coords = page['features'][0]['geometry']['rings']
            metadata['poly'] = geometryToShape(coords)
        except KeyError as e:
            print(e, page)
            metadata['poly'] = 'error'
        row.append(metadata)
        rows.append(row)
        x += 1
        if pb != None:
            pb.SetValue(x)
    print(len(rows))
    print('rows complete')
    return rows


def write_geopackage(out_path: str, name: str, poly: str, spcs: Union[str, int]):
    """
    Writes out a geopackage containing the bounding geometry of
    of the given query.

    Derived from:
    https://gis.stackexchange.com/a/52708/8104
    via
    https://gis.stackexchange.com/q/217165

    :param out_path: String representing the complete file path for the output geopackage
    :param name: String representing the name of the survey; Used to name the layer
    :param poly: The WTK Multipolygon object that holds the survey's bounding data
    :param spcs: The ESPG code for the data
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

    ## If there are multiple geometries, put the "for" loop here

    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField('Survey', name)

    # Make a geometry, from Shapely object
    geom = ogr.CreateGeometryFromWkt(poly)
    feat.SetGeometry(geom)

    layer.CreateFeature(feat)
    feat = geom = None  # destroy these

    # Save and close everything
    ds = layer = feat = geom = None


def contentSearch(contents: List[str]) -> int:
    """
    This funtion takes a list of zipfile contents.

    Using the zipfile contents, it parses the files for any file containing the
    full string '_FULL.xyz'.  If a file name contains this string, it returns a
    Boolean True

    :param contents: A list of file names
    :returns: If search conditions are met, return 1 for True
    """

    x = 0
    for content in contents:
        if full.search(content) or full_a.search(content):
            print('\nvive le resolution', content, end=' ')  # link + '\n')
            x = 1
            return x
        else:
            x = 0


def downloadAndCheck(rows: list, pb=None, to=None) -> Tuple[list, int]:
    """
    This function takes a list of complete survey data as provided by the
    query response `rows`.

    For each survey object in rows, along with the list of arrtibute values, a
    dictionary of the same values is inlcuded as well as a WTK Multipolygon
    object

    For each link provided, it saves the returned data localy. All downloaded
    files are expected to be zip files.  The funtion attempts to open the files
    as zipfile objects:

    1. If it succeeds, the function requests a list of included contents
    of the zipfile objects and passes that list, the download link, and
    location of the local download to contentSearch() to determine if the
    desired data exists. The downloaded content is kept.
    2. If it fails, the link, downloaded
    contents, and survey data in 'rows' are immediately removed.

    Through each of these steps a value is appended to the end of each 'row' or
    survey data object in list 'rows' to indicate the result of the search for
    highest resolution data.  The full range of possible values are:

    1. Yes; the survey contains a FULL.xyz file.
    2. No; the survey does not contain a FULL.xyz file.
    3. BadURL; the survey URL from query response was bad/yeilded no results.
    4. BadZip; the resulting .zip downloaded was corrupt/unable to be opened.

    :param rows: A list compiled of the attributes for every survey in surveyIDs usually generated by :func:`surveyCompile`
    :param pb:
    :param to:
    :returns: The list compiled of the attributes for every survey in surveyIDs with the results of checking for high-res data appended to the end of each entry and the total number of high-res surveys found
    """

    x = len(rows)
    hr = 0
    agencies = config['Agencies']['Agencies']
    if pb != None:
        pb.SetRange(x)
        i = 0
        pb.SetValue(i)
    for row in rows:
        link = row[6]
        surname = row[1]
        agency = row[2]
        meta = row[-1]
        spcs = 4326
        #        spcs = row[5]
        poly = meta['poly']
        name = link.split('/')[-1]
        saved = holding + '/' + agency + '/' + name
        if os.path.exists(holding + '/' + agency):
            pass
        else:
            os.mkdir(holding + '/' + agency)
        saved = os.path.normpath(saved)
        if os.path.exists(saved):
            os.remove(saved)

        if poly != 'error':
            shpfilename = os.path.join(holding + '\\' + agency,
                                       surname + '.gpkg')
            meta['poly_name'] = surname + '.gpkg'
            write_geopackage(shpfilename, surname, poly, spcs)
            sfile = os.path.relpath(surname + '.gpkg')

        metafilename = os.path.join(holding + '\\' + agency,
                                    surname + '.pickle')
        with open(metafilename, 'wb') as metafile:
            pickle.dump(meta, metafile)
        pfile = os.path.relpath(surname + '.pickle')

        print(x, agency, end=' ')
        dwntime = datetime.datetime.now()
        while True:
            if os.path.exists(saved):
                print('x', end=' ')
                break
            else:
                if link != 'null' and link != 'Not in cloud':
                    try:
                        urllib.request.urlretrieve(link, saved)
                    except socket.timeout:
                        urllib.request.urlretrieve(link, saved)
                    except urllib.error.HTTPError as e:
                        print('e \n' + link, e)
                        row.append('No')
                        row.append('BadURL')
                        break
                    except urllib.error.URLError as e:
                        print('e \n' + link, e)
                        row.append('No')
                        row.append('BadURL')
                        break
                elif time.time() - dwntime > 295:
                    print('e \n' + link)
                    row.append('No')
                    row.append('TimeOut')
                    break
                else:
                    print('e \n' + link)
                    row.append('No')
                    row.append('BadURL')
                    break
        if os.path.exists(saved):
            if (config['Resolutions']['Override'] == 'yes'
                    and (agency in agencies or agencies == '')):
                try:
                    zipped = zipfile.ZipFile(saved, mode='a')
                    os.chdir(holding + '/' + agency + '/')
                    zipped.write(pfile)
                    os.remove(pfile)
                    if poly != 'error':
                        zipped.write(sfile)
                        os.remove(sfile)
                    os.chdir(progLoc)
                    contents = zipped.namelist()
                    if contentSearch(contents) != True:
                        print('n', end=' ')
                        zipped.close()
                        row.append('No')
                    else:
                        zipped.close()
                        print('y', end=' ')
                        row.append('Yes')
                        hr += 1
                except zipfile.BadZipfile:
                    os.remove(saved)
                    print('z', end=' ')
                    row.append('BadZip')
                print('o', end=' ')
                if to != None:
                    to.write('\t\t' + agency + '\\' + name + '\n')
                row.append('Yes')
            else:
                try:
                    zipped = zipfile.ZipFile(saved, mode='a')
                    os.chdir(holding + '/' + agency + '/')
                    zipped.write(pfile)
                    os.remove(pfile)
                    if poly != 'error':
                        zipped.write(sfile)
                        os.remove(sfile)
                    os.chdir(progLoc)
                    contents = zipped.namelist()
                    if contentSearch(contents) != True:
                        print('n', end=' ')
                        zipped.close()
                        os.remove(saved)
                        print('r', end=' ')
                        row.append('No')
                    else:
                        zipped.close()
                        print('y', end=' ')
                        row.append('Yes')
                        hr += 1
                        if to != None:
                            to.write('\t\t' + agency + '\\' + name + '\n')
                except zipfile.BadZipfile:
                    os.remove(saved)
                    print('z', end=' ')
                    row.append('BadZip')
                row.append('No')
        x -= 1
        if pb != None:
            i += 1
            pb.SetValue(i)
    if to != None:
        to.write('\n')
    print('\nrow downloads verified')
    return rows, hr


def csvCompare(rows: list, csvFile: List[str], newSurveysNum: int, pb=None) -> Tuple[Union[list, str], str]:
    """
    Takes list 'rows' and list 'csvFile'.  It proceeds to compare each list
    item's contents against each other.  If they match, the relevant list item
    is removed from list 'rows'.  If all items are identical, the function
    returns a string 'No Changes'. If not empty, returns the list 'changes'

    :param rows: The list compiled of the attributes for every survey in surveyIDs. This uses the list generated by :func:`downloadAndCheck` but can also use the original list generated by :func:`surveyCompile`
    :param csvFile: The list of all compiled survey attributes stored in a dedicated .txt file generated by :func:`csvOpen`
    :param newSurveysNum: Total number of surveys usualy returned by :func:`query`
    :param pb:
    :returns: A list compiled of the attributes for every survey in surveyIDs minus the surveys found in the csvFile list or returns a string 'No Changes'
    """

    print(len(rows), end=' ')
    before = str(len(rows))
    if pb != None:
        pb.SetRange(len(rows))
        pb.SetValue(0)
    for line in csvFile:
        x = 0
        y = len(rows)
        while x < y:
            row = rows[x]
            if line[0] == row[0]:
                rows.remove(row)
                y = len(rows)
            x += 1
            if pb != None:
                pb.SetValue(x)
    print(len(rows))
    after = str(len(rows))
    numstring = ('\t\tSurveys in Query: ' + before
                 + '\n\t\tNew Surveys: ' + after)
    if len(rows) != 0:
        return rows, numstring
    else:
        return 'No Changes', numstring


def csvOpen() -> List[str]:
    """
    Uses global variable csvLocation to open eHydro_csv.txt for use.
    Populates a list 'csvFile' with it's contents. Returns list and
    closes file

    :returns: The list of all compiled survey attributes stored in eHydro_csv.txt
    """

    if not os.path.exists(csvLocation):
        create = open(csvLocation, 'w', newline='')
        create.close()
    fileOpened = open(csvLocation, 'r', newline='')
    opened = csv.reader(fileOpened, delimiter=',')
    csvFile = []
    for row in opened:
        csvFile.append(row[:13])
    fileOpened.close()
    return csvFile[1:]


def csvWriter(csvFile: List[str], csvLocation: str, pb=None):
    """
    Uses global variable csvLocation. Opens file at csvLocation for
    writing. Iterates line by line through csvFile and imediatly writes
    to eHydro_csv.txt.

    :param csvFile: A list of survey attributes each containing a row of data to be written to eHydro_csv.txt
    :param csvLocation: Complete file path string for a text file to be created
    :param pb:

    """
    csvOpen = open(csvLocation, 'w', newline='')
    save = csv.writer(csvOpen, delimiter=',')
    if pb != None:
        pb.SetRange(len(csvFile))
        pb.SetValue(0)
        x = 0
    for row in csvFile:
        truncate = row[:12]
        truncate.extend(row[-2:])
        save.writerow(truncate)
        if pb != None:
            x += 1
            pb.SetValue(x)
    csvOpen.close()


def logOpen(logType: Union[str, bool], to=None) -> Tuple[Tuple[TextIO, Any], str]:
    """
    Uses global variable logLocation. Opens file at logLocation
    for appending. Writes text stating when the function was called.
    Returns the file object for future writing.

    :param logType: Boolean expression that determines whether or not the object for the continuous log is passed back or a new log is created for an individual run
    :param to:
    :returns: A text document object representing the output log, and the file path for the output log object
    """

    timestamp = ntime()
    message = timestamp + ' - Program Initiated, Log Opened'
    if logType == 'False' or False:
        fo = open(logLocation, 'a')
        nameLog = logLocation
    elif logType == 'True' or True:
        x = 0
        datestamp = date()
        while True:
            name = datestamp + '_' + str(x) + '_' + logName
            logPath = logging + name
            if os.path.exists(logPath):
                x += 1
            else:
                break
        fo = open(logPath, 'w')
        nameLog = logPath
    fileLog = (fo, to)
    logWriter(fileLog, message)
    return fileLog, nameLog


def logWriter(fileLog: Tuple[TextIO, Any], message: str):
    """
    Takes a file object 'fileLog' and a string 'message'. Writes
    'messege' to 'fileLog'

    :param fileLog: A text document object representing the output log
    :param message: A string of text to be written to the fileLog input
    """

    print(message)
    fl, to = fileLog
    fl.write(message + '\n')
    if to != None:
        to.write(message + '\n')


def logClose(fileLog: Tuple[TextIO, Any]):
    """
    Takes a file object 'fileLog'. Writes text stating when the
    function was called. Closes the file object upon completion

    :param fileLog: A text document object representing the output log
    """

    fo = fileLog[0]
    timestamp = ntime()
    message = timestamp + ' - Program Finished, Log Closed\n'
    logWriter(fileLog, message)
    fo.close()


def ntime() -> str:
    """
    Creates and returns a string 'timestamp' that contains a
    formated current date and time at the time of calling.  Generated by
    :obj:`datetime.datetime.now` and formated using
    :obj:`datetime.datetime.strftime` to reflect 'YYYY-MM-DD Time'

    :returns: A string with the current date and time
    """

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %X')
    return timestamp


def date() -> str:
    """
    Creates and returns a string 'datestamp' that contains a
    formated current date at the time of calling.  Generated by
    :obj:`datetime.datetime.now` and formated using
    :obj:`datetime.datetime.strftime` to reflect 'YYYY-MM-DD'

    :returns: A string with the current date and time
    """

    datestamp = datetime.datetime.now().strftime('%Y-%m-%d')
    return datestamp


def fileTime() -> str:
    """
    Creates and returns a string 'timestamp' that contains a
    formated current date and time at the time of calling.  Generated by
    :obj:`datetime.datetime.now` and formated using
    :obj:`datetime.datetime.strftime` to reflect 'YYYYMMDD'

    :returns: A string with the current date and time
    """

    timestamp = datetime.datetime.now().strftime('%Y%m%d')
    return timestamp


def main(pb=None, to=None):
    """
    Main function.

    :param pb:
    :param to:
    """

    runType = config['Data Checking']['Override']
    logType = config['Output Log']['Log Type']
    fileLog, nameLog = logOpen(logType, to)
    try:
        surveyIDs, newSurveysNum, paramString = query()
        logWriter(fileLog, '\tSurvey IDs queried from eHydro\n' + paramString)
        logWriter(fileLog, '\tCompiling survey objects from Survey IDs')
        rows = surveyCompile(surveyIDs, newSurveysNum, pb)
        logWriter(fileLog, '\tSurvey objects compiled from eHydro')
    except:
        logWriter(fileLog, '\teHydro query failed')
    if runType == 'no':
        try:
            csvFile = csvOpen()
            logWriter(fileLog, '\teHydro_csv.txt opened for reading')
        except:
            logWriter(fileLog, '\teHydro_csv.txt unable to be opened')
        try:
            logWriter(fileLog, '\tComparing query results to eHydro_csv.txt')
            changes, numstring = csvCompare(rows, csvFile, newSurveysNum)
            logWriter(fileLog, numstring)
        except:
            logWriter(fileLog,
                      '\t\tUnable to compare query results to eHydro_csv.txt')
    elif runType == 'yes':
        csvFile = []
        changes = rows
    try:
        logWriter(fileLog, '\tParsing new entries for resolution:')
        attributes.append('Hi-Res?')
        attributes.append('Override?')
        if changes != 'No Changes':
            checked, hiRes = downloadAndCheck(changes, pb, to)
            csvFile.extend(checked)
            if config['Output Log']['Query List'] == 'yes':
                logWriter(fileLog, '\tNew Survey Details:')
                for row in checked:
                    txt = ''
                    for i in [1, 4, 5, 6, -2]:
                        txt = txt + attributes[i] + ' : ' + row[i] + '\n\t\t'
                    logWriter(fileLog, '\t\t' + txt)
            logWriter(fileLog, '\t\tTotal High Resloution Surveys: '
                      + str(hiRes) + '/' + str(len(changes)) + '\n')
        else:
            logWriter(fileLog, '\t\t' + changes)
    except:
        logWriter(fileLog, '\tParsing for resolution failed')
    try:
        csvFile.insert(0, attributes)
        csvSave = csvFile
        if runType == 'no':
            csvPath = csvLocation
            csvWriter(csvSave, csvPath, pb)
        elif runType == 'yes':
            x = 0
            datestamp = date()
            while True:
                name = datestamp + '_' + str(x) + '_' + csvName
                csvPath = running + name
                if os.path.exists(csvPath):
                    x += 1
                else:
                    break
            csvWriter(csvSave, csvPath, pb)
        logWriter(fileLog, '\tAdding results to ' + csvPath)
    except:
        logWriter(fileLog, '\tUnable to add results to ' + csvPath)

    logWriter(fileLog, '\tOutput Log saved as ' + nameLog)
    logClose(fileLog)
    print('log closed')


if __name__ == '__main__':
    """Function call to initiate program"""
    print(__version__)
    main()
