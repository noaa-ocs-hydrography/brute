# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 13:44:14 2018

@author: Casiano.Koprowski
"""

import os
import re
import csv
import json
import requests
import urllib
import zipfile
import datetime
import socket
import configparser

'''Known global constants'''
# progLoc is the program's own file location / current working directory (cwd)
progLoc = os.getcwd()
# regex object for searching zipfile contents
full = re.compile(r'FULL.xyz', re.IGNORECASE)
fullalt = re.compile(r'_A.xyz', re.IGNORECASE)
# location of data items
config = configparser.ConfigParser()
config.read('config.ini')
txtName = 'eHydro_txt.txt'
txtLocation = os.path.join(progLoc, txtName)   
csvName = 'eHydro_csv.txt'
csvLocation = os.path.join(progLoc, csvName)
logName = 'eHydro_log.txt'
logLocation = os.path.join(progLoc, logName)
holding = progLoc + '/downloads/'
# eHydro survey entry attributes
attributes = [ "OBJECTID", "SURVEYJOBIDPK", "SURVEYAGENCY", "CHANNELAREAIDFK",
              "SDSFEATURENAME", "SOURCEPROJECTION", "SURVEYDATEUPLOADED",
              "SOURCEDATALOCATION", "SURVEYDATEEND", "SURVEYDATESTART",
              "SURVEYTYPE", "PROJECTEDAREA"]
# check to see if the downloaded data folder exists, will create it if not
if os.path.exists(holding):
    pass
else:
    os.mkdir(holding)
    
def query():
    '''Holds the Queries for the eHydro REST API, asks for responses, and uses
    the json library to make it readable by the program. Returns the json
    responses for number of surveys and surveys in the response.

    It also saves a prettyprinted version of the response as a text file.

    The funtion uses the requests library to retrieve the API's response(s) and
    uses the json library to make them readable by the program.

    The function returns the json object containing the contents of the query
    and the integer number of surveys contained by the query
    '''
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
            areas += '%20AND%20%20('
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
            areas = ('%20AND%20UPPER(SURVEYAGENCY)%20like%20%27%25'
                     + agencies[0]
                     + '%25%27')
    else:
        areas = ''
        
    # The main query parameters that will determine the contents of the response    
    # Survey Date Uploaded
    if config ['Timeframe']['Ignore Date'] == 'no':
        where = ('SURVEYDATEUPLOADED%20%3E%3D%20%27'
                 + start
                 + 'T04%3A00%3A00.000Z%27%20AND%20SURVEYDATEUPLOADED%20%3C%3D%20%27'
                 + end
                 + 'T04%3A00%3A00.000Z%27'
                 + areas)
    else:
        if areas != '':    
            where = areas
        else:
            where = '1%3D1'
    
    # The query for determining how many responses will be returned
    newSurveys = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?where=' + where + '&outFields=*&returnGeometry=false&returnCountOnly=true&outSR=&f=json'
    
    # The query for returning the object IDs for the given timeframe
    objIDs = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?&where=' + where + '&outFields=*&returnGeometry=false&returnIdsOnly=true&outSR=&f=json'
    
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
    
    return (surveyIDs, newSurveysNum)

def surveyCompile(surveyIDs, newSurveysNum):
    '''Uses the json object return of the query and the total number of surveys
    included to compile a list of complete returned survey data, as provided in
    the response. The function also takes into account that the survey data 
    returns for any date/time are returned as timestamps. The function looks 
    for these fields and converts them to datetime objects and finaly strings.

    The function returns the lists of returned survey data as a list 'rows'.
    '''
    x = 0
    rows = []
    while x < newSurveysNum:
        print (x, end=' ')
        query = ('https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?where=OBJECTID%20%3D%20'
                 + str(surveyIDs[x]) 
                 + '&outFields=OBJECTID,SDSFEATURENAME,SURVEYTYPE,CHANNELAREAIDFK,SURVEYDATEUPLOADED,SURVEYAGENCY,SURVEYDATESTART,SURVEYDATEEND,SOURCEDATALOCATION,SOURCEPROJECTION,SURVEYJOBIDPK,PROJECTEDAREA&returnGeometry=false&outSR=&f=json')
        response = requests.get(query)
        page = response.json()
        row = []
        for attribute in attributes:
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
        rows.append(row)
        x += 1
    print (len(rows))
    print ('rows complete')
    return rows

def contentSearch(contents, link, saved):
    '''This funtion takes a list of zipfile contents, the download link the
    contents came from and the current local file location where the downloaded
    data resides.

    Using the zipfile contents, it parses the files for any file containing the
    full string '_FULL.xyz'.  If a file name contains this string, it returns a
    Boolean True
    '''
    x = 0
    for content in contents:
        if full.search(content):
            print ('\nvive le resolution', content, end=' ') #link + '\n')
            x = 1
            return x
        elif fullalt.search(content):
            print ('\nlong live the resolution', content, end=' ') #link + '\n')
            x = 1
            return x
        else:
            x = 0

def downloadAndCheck(rows):
    '''This function takes a list of complete survey data as provided by the 
    query response ('rows').

    For each link provided, it saves the returned data localy. All downloaded 
    files are expected to be zip files.  The funtion attempts to open the files
    as zipfile objects:
        1) If it succeeds, the function requests a list of included contents
        of the zipfile objects and passes that list, the download link, and
        location of the local download to contentSearch() to determine if the
        desired data exists. The downloaded content is kept.
        2) If it fails, the link, downloaded
        contents, and survey data in 'rows' are immediately removed.

   Through each of these steps a value is appended to the end of each 'row' or
   survey data object in list 'rows' to indicate the result of the search for 
   highest resolution data.  The full range of possible values are:
       1) Yes; the survey contains a FULL.xyz file
       2) No; the survey does not contain a FULL.xyz file
       3) BadURL; the survey URL from query response was bad/yeilded no results
       4) BadZip; the resulting .zip downloaded was corrupt/unable to be opened
    '''
    x = len(rows)
    agencies = config['Agencies']['Agencies']
    for row in rows:
        link = row[7]
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
                try:
                    urllib.request.urlretrieve(link, saved)
                except socket.timeout:
                    urllib.request.urlretrieve(link, saved)
                except urllib.error.HTTPError:
                    print ('e', end=' ')
                    row.append('BadURL')
                    break
        if os.path.exists(saved):
            if config['Resolutions']['Override'] == 'yes' and agency in agencies:
                print ('o', end=' ')
                row.append('Ovrd')
            elif config['Resolutions']['Override'] == 'yes' and agencies == '':
                print ('o', end=' ')
                row.append('Ovrd')
            else:
                try:
                    zipped = zipfile.ZipFile(saved)
                    contents = zipped.namelist()
                    if contentSearch(contents, link, saved) != True:
                        print ('n', end=' ')
                        zipped.close()
                        os.remove(saved)
                        print ('r', end=' ')
                        row.append('No')
                    else:
                        zipped.close()
                        print ('y', end=' ')
                        row.append('Yes')
                except zipfile.BadZipfile:
                    os.remove(saved)
                    print ('z', end=' ')
                    row.append('BadZip')
        x -= 1
    print ('row downloads verified')
    return rows

def csvCompare(rows, csvFile, newSurveysNum):
    '''Takes list 'rows' and list 'csvFile'.  It proceeds to compare each list 
    item's contents against each other.  If they match, the relevant list item 
    is removed from list 'rows'.  If all items are identical, the function 
    returns a string 'No Changes'. If not empty, returns the list 'changes'
    '''
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
    '''String "fileText" is writen to the "txtLocation" save path'''
    save = open(txtLocation, 'w')
    for row in fileText:
        save.writelines(row)
    save.close()

def csvOpen():
    '''Uses global variable csvLocation to open NCEI_csv.txt for use.
    Populates a list 'csvFile' with it's contents. Returns list and
    closes file
    '''
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
    '''Uses global variables txtLocation and csvLocation. Opens file
    at txtLocation for reading and overwrites file at csvLocation for
    writing. Iterates line by line through 'txt' and imediatly writes
    to 'csv'. Closes both opened files.
    '''
    csvOpen = open(csvLocation, 'w')
    save = csv.writer(csvOpen, delimiter = ',')
    for row in csvFile:
        save.writerow(row)
    csvOpen.close()

def logOpen():
    '''Uses global variable logLocation. Opens file at logLocation
    for appending. Writes text stating when the function was called.
    Returns the file object for future writing.
    '''
    timestamp = time()
    fileLog = open(logLocation, 'a')
    message = '\n' + timestamp + ': Program Initiated, Log Opened'
    logWriter(fileLog, message)
    return fileLog

def logWriter(fileLog, message):
    '''Takes a file object 'fileLog' and a string 'message'. Writes
    'messege' to 'fileLog'
    '''
    print (message)
    fileLog.write(message + '\n')

def logClose(fileLog):
    '''Takes a file object 'fileLog'. Writes text stating when the
    function was called. Closes the file object upon completion
    '''
    timestamp = time()
    message = timestamp +': Program Finished, Log Closed'
    logWriter(fileLog, message)
    fileLog.close()

def time():
    '''Creates and returns a string 'timestamp' that contains a
    formated current date and time at the time of calling.
    '''
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %X')
    return timestamp


def main():
    fileLog = logOpen()
    try:
        surveyIDs, newSurveysNum = query()
        rows = surveyCompile(surveyIDs, newSurveysNum)
        logWriter(fileLog, '\tSurveys queried from eHydro')
    except:
        logWriter(fileLog, '\teHydro query failed')
    try:
        csvFile = csvOpen()
        txtWriter(csvFile, txtLocation)
        logWriter(fileLog, '\teHydo_csv.txt opened for reading')
    except:
        logWriter(fileLog, '\teHydo_csv.txt unable to be opened')
    try:
        logWriter(fileLog, '\tComparing query results to eHydo_csv.txt')
        changes = csvCompare(rows, csvFile, newSurveysNum)
    except:
        logWriter(fileLog, '\t\tUnable to compare query results to eHydo_csv.txt')
    try:
        logWriter(fileLog, '\tParsing new entries for resolution:')
        attributes.append("Hi-Res?")
        if changes != 'No Changes':
            checked = downloadAndCheck(changes)
            csvFile.extend(checked)
            for row in checked:
                txt = ''
                for i in [1,4,5,7,12]:
                    txt = txt + attributes[i] + ' : ' + row[i] + '\n\t\t'
                logWriter(fileLog, '\t\t' + txt)
        else:
            logWriter(fileLog, '\t\t' + changes)
    except:
        logWriter(fileLog, '\tParsing for resolution failed')
    try:
        csvFile.insert(0, attributes)
        csvSave = csvFile
        csvWriter(csvSave, csvLocation)
        logWriter(fileLog, '\tAdding results to eHydo_csv.txt')
    except:
        logWriter(fileLog, '\tUnable to add results to eHydo_csv.txt')
    logClose(fileLog)
    print('log closed')

'''Function call to initiate program'''
main()
