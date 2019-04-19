# -*- coding: utf-8 -*-
"""
Created on Tue Oct 2 13:44:14 2018

Last Modified: Apr 19 12:12:42 2019

@author: Casiano.Koprowski <casiano.koprowski@noaa.gov>
"""

import os
import re
import csv
import socket
import urllib
import zipfile
import datetime
import requests
import configparser

"""Known global constants"""
# print (datetime.datetime.now().strftime('%b %d %X %Y'))
#progLoc = '\\eHydro_scrape'
progLoc = os.getcwd()
"""progLoc is the program's own file location / current working directory (cwd)
obtained by :func:`os.getcwd()`"""
full = re.compile(r'FULL.xyz', re.IGNORECASE)
"""regex object for searching zipfile contents for high-res data"""
full_a = re.compile(r'_A.xyz', re.IGNORECASE)
"""regex object for searching zipfile contents for high-res data"""
# location of data items
config = configparser.ConfigParser()
config.read('config.ini')
#txtName = 'eHydro_txt.txt'
#txtLocation = os.path.join(progLoc, txtName)
csvName = 'eHydro_csv.txt'
"""Default name for csv.txt output"""
csvLocation = os.path.join(progLoc, csvName)
"""Default location for :attr:`csvName`"""
logName = 'eHydro_log.txt'
"""Default name for log.txt output"""
logLocation = os.path.join(progLoc, logName)
"""Default location for :attr:`logName`"""
#"""Default location for """
holding = progLoc + '\\downloads\\'
"""Default location for downloaded data"""
logging = progLoc + '\\logs\\'
"""Default location for individual query logs"""
running = progLoc + '\\runs\\'
"""Default location for individual query csvs"""
# eHydro survey entry attributes
attributes = [ "OBJECTID", "SURVEYJOBIDPK", "SURVEYAGENCY", "CHANNELAREAIDFK",
              "SDSFEATURENAME", "SOURCEPROJECTION", "SOURCEDATALOCATION",
              "SURVEYDATEUPLOADED", "SURVEYDATEEND", "SURVEYDATESTART",
              "SURVEYTYPE", "PROJECTEDAREA"]
"""The specific attributes queried for each survey in :func:`surveyCompile`"""

# check to see if the downloaded data folder exists, will create it if not
if not os.path.exists(holding):
    os.mkdir(holding)
if not os.path.exists(logging):
    os.mkdir(logging)
if not os.path.exists(running):
    os.mkdir(running)


def query():
    """Holds the Queries for the eHydro REST API, asks for responses, and uses
    the json library to make it readable by the program. Returns the json
    responses for number of surveys and surveys in the response.

    -REMOVED- It also saves a prettyprinted version of the response as a text
    file.

    The funtion uses the requests library to retrieve the API's response(s) and
    uses the json library to make them readable by the program.

    The function returns the json object containing the contents of the query
    and the integer number of surveys contained by the query

    Returns
    -------
    surveyIDs : list
        List of survey ids from query
    newSurveysNum : int
        Total number of surveys returned by the query
    paramString : str
        String containing the parameters gathered from the config file

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
    print ('\nStart:', start, '\nEnd:', end)
    if (config['Agencies']['Only Listed'] == 'yes'
        and config['Agencies']['Agencies'] != ''):
        areas = ''
        agencies = config['Agencies']['Agencies'].split(',')

        if len(agencies) > 1:
            i = 0
            areas += '('
            while i < len(agencies):
                if i == 0:
                    areas += ('UPPER(SURVEYAGENCY)%20like%20%27%25'
                              + agencies[0]
                              + '%25%27')
                    i += 1
                else:
                    areas += ('%20OR%20UPPER(SURVEYAGENCY)%20like%20%27%25'
                              + agencies[i]
                              + '%25%27')
                    i += 1
            areas += ')%20'
        else:
            areas = ('UPPER(SURVEYAGENCY)%20like%20%27%25'
                     + agencies[0]
                     + '%25%27')
    else:
        areas = ''

    # The main query parameters that will determine the contents of the response
    # Survey Date Uploaded
    if config ['Timeframe']['Ignore Date'] == 'no' and areas != '':
        where = ('SURVEYDATEUPLOADED%20%3E%3D%20%27'
                 + start
                 + 'T00%3A01%3A00.000Z%27%20AND%20SURVEYDATEUPLOADED%20%3C%3D%20%27'
                 + end
                 + 'T11%3A59%3A00.000Z%27%20AND%20'
                 + areas)
    elif config ['Timeframe']['Ignore Date'] == 'no' and areas == '':
        where = ('SURVEYDATEUPLOADED%20%3E%3D%20%27'
                 + start
                 + 'T00%3A01%3A00.000Z%27%20AND%20SURVEYDATEUPLOADED%20%3C%3D%20%27'
                 + end
                 + 'T11%3A59%3A00.000Z%27')
    else:
        if areas != '':
            where = areas
        else:
            where = '1%3D1'

    # The query for determining how many responses will be returned
    newSurveys = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?where=' + where + '&outFields=*&returnGeometry=false&returnCountOnly=true&outSR=&f=json'

    # The query for returning the object IDs for the given timeframe
    objIDs = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?&where=' + where + '&outFields=*&returnGeometry=false&returnIdsOnly=true&outSR=&f=json'

    print (objIDs, newSurveys)

    # Initial Query execution
    surveyNumRequest = requests.get(newSurveys)
    surveyIDsRequest = requests.get(objIDs)

    # Response for number of surveys used as a json object
    surveyNumJSON = surveyNumRequest.json()
    print (surveyNumJSON)
    newSurveysNum = int(surveyNumJSON['count'])

    # Response for survey Object IDs as a json object
    surveyIDsJSON = surveyIDsRequest.json()
    surveyIDs = surveyIDsJSON['objectIds']
    print (surveyIDs, newSurveysNum)

    if areas == '':
        dist = 'none'
    else:
        dist = config['Agencies']['Agencies']

    paramString = str('\tParameters:\n\t\tStart Date: ' + start
                      + '\n\t\tEnd Date: ' + end
                      + '\n\t\tDistricts: ' + dist
                      + '\n\t\tQuery Only Districts: ' + config['Agencies']['Only Listed']
                      + '\n\t\tKeep All Data: ' + config['Resolutions']['Override'])

    return (surveyIDs, newSurveysNum, paramString)

def surveyCompile(surveyIDs, newSurveysNum):
    """Uses the json object return of the each queried survey id and the total
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

    Parameters
    ----------
    surveyIDs : list
        List of survey ids usualy generated by :func:`query`
    newSurveysNum : int
        Total number of surveys usualy returned by :func:`query`

    Returns
    -------
    rows : list
        A list compiled of the attributes for every survey in surveyIDs

    """
    x = 0
    rows = []
    while x < newSurveysNum:
        print (x, end=' ')
        query = ('https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?where=OBJECTID%20%3D%20'
                 + str(surveyIDs[x])
                 + '&outFields=OBJECTID,SDSFEATURENAME,SURVEYTYPE,CHANNELAREAIDFK,SURVEYAGENCY,SURVEYDATEUPLOADED,SURVEYDATESTART,SURVEYDATEEND,SOURCEDATALOCATION,SOURCEPROJECTION,SURVEYJOBIDPK,PROJECTEDAREA&returnGeometry=false&outSR=&f=json')
        response = requests.get(query)
        page = response.json()
        row = []
        for attribute in attributes:
            try:
                if page['features'][0]['attributes'][attribute] == None:
                    row.append('null')
                elif (attribute == "SURVEYDATEUPLOADED"
                    or attribute == "SURVEYDATEEND"
                    or attribute == "SURVEYDATESTART"):
                        if page['features'][0]['attributes'][attribute] == None:
                            row.append('null')
                        else:
                            date = (page['features'][0]['attributes'][attribute])
                            date = datetime.datetime.utcfromtimestamp(date/1000)
                            row.append(str(date.strftime('%Y-%m-%d')))
                else:
                    row.append(str(page['features'][0]['attributes'][attribute]))
            except KeyError as e:
                print (e, page)
                row.append('error')
        rows.append(row)
        x += 1
    print (len(rows))
    print ('rows complete')
    return rows

def contentSearch(contents):
    """This funtion takes a list of zipfile contents, the download link the
    contents came from and the current local file location where the downloaded
    data resides.

    Using the zipfile contents, it parses the files for any file containing the
    full string '_FULL.xyz'.  If a file name contains this string, it returns a
    Boolean True

    Parameters
    ----------
    contents : list
        A list of file names

    Returns
    -------
    int, bool
        If search conditions are met, return 1 for True

    """
    x = 0
    for content in contents:
        if full.search(content) or full_a.search(content):
            print ('\nvive le resolution', content, end=' ') #link + '\n')
            x = 1
            return x
        else:
            x = 0

def downloadAndCheck(rows):
    """This function takes a list of complete survey data as provided by the
    query response ('rows').

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

    Parameters
    ----------
    rows : list
        A list compiled of the attributes for every survey in surveyIDs usually
        generated by :func:`surveyCompile`

    Returns
    -------
    rows : list
        The list compiled of the attributes for every survey in surveyIDs with
        the results of checking for high-res data appended to the end of each
        entry
    hr : int
        The total number of high-res surveys found

    """
    x = len(rows)
    hr = 0
    agencies = config['Agencies']['Agencies']
    for row in rows:
        link = row[6]
        agency = row[2]
        name = link.split('/')[-1]
        saved = holding + '/' + agency + '/' + name
        if os.path.exists(holding + '/' + agency):
            pass
        else:
            os.mkdir(holding + '/' + agency)
        saved = os.path.normpath(saved)
        print  (x, agency,  end=' ')
        while True:
            if os.path.exists(saved):
                print ('x', end=' ')
                break
            else:
                if link != 'null' and link != 'Not in cloud':
                    try:
                        urllib.request.urlretrieve(link, saved)
                    except socket.timeout:
                        urllib.request.urlretrieve(link, saved)
                    except urllib.error.HTTPError as e:
                        print ('e \n' + link)
                        row.append('No')
                        row.append('BadURL')
                        break
                    except urllib.error.URLError as e:
                        print ('e \n' + link)
                        row.append('No')
                        row.append('BadURL')
                        break
                else:
                    print ('e \n' + link)
                    row.append('No')
                    row.append('BadURL')
                    break
        if os.path.exists(saved):
            if config['Resolutions']['Override'] == 'yes' and (agency in agencies or agencies == ''):
                try:
                    zipped = zipfile.ZipFile(saved)
                    contents = zipped.namelist()
                    if contentSearch(contents) != True:
                        print ('n', end=' ')
                        zipped.close()
                        row.append('No')
                    else:
                        zipped.close()
                        print ('y', end=' ')
                        row.append('Yes')
                        hr += 1
                except zipfile.BadZipfile:
                    os.remove(saved)
                    print ('z', end=' ')
                    row.append('BadZip')
                print ('o', end=' ')
                row.append('Yes')
            else:
                try:
                    zipped = zipfile.ZipFile(saved)
                    contents = zipped.namelist()
                    if contentSearch(contents) != True:
                        print ('n', end=' ')
                        zipped.close()
                        os.remove(saved)
                        print ('r', end=' ')
                        row.append('No')
                    else:
                        zipped.close()
                        print ('y', end=' ')
                        row.append('Yes')
                        hr += 1
                except zipfile.BadZipfile:
                    os.remove(saved)
                    print ('z', end=' ')
                    row.append('BadZip')
                row.append('No')
        x -= 1
    print ('\nrow downloads verified')
    return rows, hr

def csvCompare(rows, csvFile, newSurveysNum):
    """Takes list 'rows' and list 'csvFile'.  It proceeds to compare each list
    item's contents against each other.  If they match, the relevant list item
    is removed from list 'rows'.  If all items are identical, the function
    returns a string 'No Changes'. If not empty, returns the list 'changes'

    Parameters
    ----------
    rows : list
        The list compiled of the attributes for every survey in surveyIDs. This
        uses the list generated by :func:`downloadAndCheck` but can also use
        the original list generated by :func:`surveyCompile`
    csvFile : list
        The list of all compiled survey attributes stored in a dedicated .txt
        file generated by :func:`csvOpen`
    newSurveysNum : int
        Total number of surveys usualy returned by :func:`query`

    Returns
    -------
    list, str
        returns A list compiled of the attributes for every survey in surveyIDs
        minus the surveys found in the csvFile list or returns a string 'No
        Changes'

    """
    print(len(rows), end = ' ')
    for line in csvFile:
        x = 0
        y = len(rows)
        while x < y:
            row = rows[x]
            if line[0] == row[0]:
                rows.remove(row)
                y = len(rows)
            x += 1
    print(len(rows))
    if len(rows) != 0:
        return rows
    else:
        return 'No Changes'

def txtWriter(fileText, txtLocation):
    """String "fileText" is writen to the "txtLocation" save path

    Parameters
    ----------
    fileText : list
        A list of strings each containing a line of text to be written in the
        outputLog
    txtLocation : string
        Complete file path string for a text file to be created

    """
    save = open(txtLocation, 'w')
    for row in fileText:
        save.writelines(row)
    save.close()

def csvOpen():
    """Uses global variable csvLocation to open eHydro_csv.txt for use.
    Populates a list 'csvFile' with it's contents. Returns list and
    closes file

    Returns
    -------
    csvFile : list
        The list of all compiled survey attributes stored in eHydro_csv.txt

    """
    if not os.path.exists(csvLocation):
        create = open(csvLocation, 'w')
        create.close()
    fileOpened = open(csvLocation, 'r', newline='\n')
    opened = csv.reader(fileOpened, delimiter = ',')
    csvFile = []
    for row in opened:
        csvFile.append(row[:13])
    fileOpened.close()
    return csvFile[1:]

def csvWriter(csvFile, csvLocation):
    """Uses global variable csvLocation. Opens file at csvLocation for
    writing. Iterates line by line through csvFile and imediatly writes
    to eHydro_csv.txt.

    Parameters
    ----------
    csvFile : list
        A list of survey attributes each containing a row of data to be written
        to eHydro_csv.txt
    txtLocation : string
        Complete file path string for a text file to be created

    """
    csvOpen = open(csvLocation, 'w')
    save = csv.writer(csvOpen, delimiter = ',')
    for row in csvFile:
        save.writerow(row)
    csvOpen.close()

def logOpen(logType):
    """Uses global variable logLocation. Opens file at logLocation
    for appending. Writes text stating when the function was called.
    Returns the file object for future writing.

    Parameters
    ----------
    logType : str, bool
        Boolean expression that determines whether or not the object for the
        continuous log is passed back or a new log is created for an individual
        run

    Returns
    -------
    fileLog : text file object
        A text document object representing the output log
    nameLog : str
        File path for the output log object

    """
    timestamp = time()
    message = timestamp + ' - Program Initiated, Log Opened'
    if logType == 'False' or False:
        fileLog = open(logLocation, 'a')
        nameLog = logLocation
    elif logType == 'True' or True:
        x = 0
        datestamp = date()
        while True:
            name = datestamp +'_' + str(x) + '_' + logName
            logPath = logging + name
#            print (logPath)
            if os.path.exists(logPath):
                x += 1
            else:
                break
        fileLog = open(logPath, 'w')
        nameLog = logPath
    logWriter(fileLog, message)
    return fileLog, nameLog

def logWriter(fileLog, message):
    """Takes a file object 'fileLog' and a string 'message'. Writes
    'messege' to 'fileLog'

    Parameters
    ----------
    fileLog : text file object
        A text document object representing the output log
    message : string
        A string of text to be written to the fileLog input

    """
    print (message)
    fileLog.write(message + '\n')

def logClose(fileLog):
    """Takes a file object 'fileLog'. Writes text stating when the
    function was called. Closes the file object upon completion

    Parameters
    ----------
    fileLog : text file object
        A text document object representing the output log

    """
    timestamp = time()
    message = timestamp +' - Program Finished, Log Closed\n'
    logWriter(fileLog, message)
    fileLog.close()

def time():
    """Creates and returns a string 'timestamp' that contains a
    formated current date and time at the time of calling.  Generated by
    :obj:`datetime.now` and formated using
    :obj:`datetime.datetime.strftime` to reflect 'YYYY-MM-DD Time'

    Returns
    -------
    timestamp : string
        A string with the current date and time

    """
    datetime.datetime.strftime
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %X')
    return timestamp

def date():
    """Creates and returns a string 'datestamp' that contains a
    formated current date at the time of calling.  Generated by
    :obj:`datetime.datetime.now` and formated using
    :obj:`datetime.datetime.strftime` to reflect 'YYYY-MM-DD'

    Returns
    -------
    timestamp : string
        A string with the current date and time

    """
    datestamp = datetime.datetime.now().strftime('%Y-%m-%d')
    return datestamp

def fileTime():
    """Creates and returns a string 'timestamp' that contains a
    formated current date and time at the time of calling.  Generated by
    :obj:`datetime.datetime.now` and formated using
    :obj:`datetime.datetime.strftime` to reflect 'YYYYMMDD'

    Returns
    -------
    timestamp : string
        A string with the current date and time
    """
    timestamp = datetime.datetime.now().strftime('%Y%m%d')
    return timestamp

def main():
    runType = config['Data Checking']['Override']
    logType = config['Output Log']['Log Type']
    fileLog, nameLog = logOpen(logType)
    try:
        surveyIDs, newSurveysNum, paramString = query()
        rows = surveyCompile(surveyIDs, newSurveysNum)
        logWriter(fileLog, '\tSurveys queried from eHydro\n' + paramString)
    except:
        logWriter(fileLog, '\teHydro query failed')
    if runType == 'no':
        try:
            csvFile = csvOpen()
            txtWriter(csvFile, csvLocation)
            logWriter(fileLog, '\teHydro_csv.txt opened for reading')
        except:
            logWriter(fileLog, '\teHydro_csv.txt unable to be opened')
        try:
            logWriter(fileLog, '\tComparing query results to eHydro_csv.txt')
            changes = csvCompare(rows, csvFile, newSurveysNum)
        except:
            logWriter(fileLog, '\t\tUnable to compare query results to eHydro_csv.txt')
    elif runType == 'yes':
        csvFile = []
        changes = rows
    try:
        logWriter(fileLog, '\tParsing new entries for resolution:')
        attributes.append('Hi-Res?')
        attributes.append('Override?')
        if changes != 'No Changes':
            checked, hiRes = downloadAndCheck(changes)
            csvFile.extend(checked)
            if config['Output Log']['Query List'] == 'yes':
                for row in checked:
                    txt = ''
                    for i in [1,4,5,6,12]:
                        txt = txt + attributes[i] + ' : ' + row[i] + '\n\t\t'
                    logWriter(fileLog, '\t\t' + txt)
            logWriter(fileLog, '\t\tTotal High Resloution Surveys: ' + str(hiRes) + '/' + str(len(changes)) + '\n')
        else:
            logWriter(fileLog, '\t\t' + changes)
    except:
        logWriter(fileLog, '\tParsing for resolution failed')
    try:
        csvFile.insert(0, attributes)
        csvSave = csvFile
        if runType == 'no':
            csvPath = csvLocation
            csvWriter(csvSave, csvPath)
        elif runType == 'yes':
            x = 0
            datestamp = date()
            while True:
                name = datestamp +'_' + str(x) + '_' + csvName
                csvPath = running + name
#                print (csvPath)
                if os.path.exists(csvPath):
                    x += 1
                else:
                    break
            csvWriter(csvSave, csvPath)
        logWriter(fileLog, '\tAdding results to ' + csvPath)
    except:
        logWriter(fileLog, '\tUnable to add results to ' + csvPath)

    logWriter(fileLog, '\tOutput Log saved as ' + nameLog)
    logClose(fileLog)
    print('log closed')

if __name__ == '__main__':
    """Function call to initiate program"""
    main()

