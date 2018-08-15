# -*- coding: utf-8 -*-
"""
Created on Fri Jul 06 15:10:17 2018

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

'''Known global constants'''
# progLoc is the program's own file location / current working directory (cwd)
progLoc = os.getcwd()
# regex object for searching zipfile contents
full = re.compile(r'FULL.xyz', re.IGNORECASE)
# location of downloaded data
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
# global list to be populated with the download link and local location of data
resPresent = []


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
    # Today - 1 (ex. '2018-08-06'), unformatted
    yesterday = today - datetime.timedelta(30)
    # Today - 10 (ex. '2018-07-29'), unformatted
    otherday = today - datetime.timedelta(40)

    # Prints of the formated versions of the date used by the query
    print (today.strftime('%Y-%m-%d'),
           yesterday.strftime('%Y-%m-%d'),
           otherday.strftime('%Y-%m-%d'))

    # The main query parameters that will determine the contents of the response
    '''
    # Survey Date End
    where = ('SURVEYDATEEND%20%3C%3D%20timestamp%20%27'
             + str(today.strftime('%Y-%m-%d'))
             + '%2000%3A01%3A00%27%20AND%20SURVEYDATEEND%20%3E%3D%20timestamp%20%27'
             + str(yesterday.strftime('%Y-%m-%d'))
             + '%2000%3A01%3A00%27')
    '''
    # Survey Date Uploaded
#    where = ('SURVEYDATEUPLOADED%20%3C%3D%20timestamp%20%27'
#                     + str(today.strftime('%Y-%m-%d'))
#                     + '%2000%3A01%3A00%27%20AND%20SURVEYDATEUPLOADED%20%3E%3D%20timestamp%20%27'
#                     + str(yesterday.strftime('%Y-%m-%d'))
#                     + '%2000%3A01%3A00%27')
        # Survey Date Uploaded
    where = ('SURVEYDATEUPLOADED%20%3C%3D%20timestamp%20%27'
                     + str(today.strftime('%Y-%m-%d'))
                     + '%2000%3A01%3A00%27%20AND%20SURVEYDATEUPLOADED%20%3E%3D%20timestamp%20%27'
                     + str(yesterday.strftime('%Y-%m-%d'))
                     + '%2000%3A01%3A00%27')

    # The query for determining how many responses will be returned
    newSurveys = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?f=json&where=' + where + '&returnGeometry=false&spatialRel=esriSpatialRelIntersects&outFields=*&outStatistics=%5B%7B%22statisticType%22%3A%22count%22%2C%22onStatisticField%22%3A%22OBJECTID%22%2C%22outStatisticFieldName%22%3A%22value%22%7D%5D'
    # The query for to actually recieve the json formatted responses
    url = 'https://services7.arcgis.com/n1YM8pTrFmm7L4hs/arcgis/rest/services/eHydro_Survey_Data/FeatureServer/0/query?f=json&where=' + where + '&returnGeometry=false&spatialRel=esriSpatialRelIntersects&maxAllowableOffset=19567&outFields=*&orderByFields=SURVEYDATEUPLOADED%20desc%2CSURVEYDATESTART%20desc&outSR=102100&resultOffset=0&resultRecordCount=10000'
    # Query execution
    surveyNumRequest = requests.get(newSurveys)
    test = requests.get(url)
    # Response for number of surveys used as a json object
    surveyNumJSON = surveyNumRequest.json()
    records = json.dumps(surveyNumJSON, indent = 4)
    newSurveysNum = int(surveyNumJSON['features'][0]['attributes']['value'])
    if int(surveyNumJSON['features'][0]['attributes']['value']) <= 1000:
        newSurveysNum = int(surveyNumJSON['features'][0]['attributes']['value'])
    else:
        newSurveysNum = 1000
    print (newSurveysNum)
    # Actuall response of survey details used as a json object
    page = test.json()
    records = json.dumps(page, indent = 4)
    opened = open(progLoc + '/json.txt', 'w')
    opened.writelines(records)
    opened.close()
    # Returns the json object for the query responses, and int for number of.
    return page, newSurveysNum


def surveyCompile(page, newSurveysNum):
    '''Uses the json object return of the query and the total number of surveys
    included to compile two lists:
        1) A list of complete retunred survey data, as provided in the response
        2) A list of only the download links to the survey data, as provided by
        the response
    The function also takes into account that the survey data returns for any
    date/time are returned as timestamps. The function looks for these fields
    and converts them to datetime objects and finaly strings.

    The function returns the lists of returned survey data and data download
    links.
    '''
    x = 0
    rows = []
    links = []
    while x < newSurveysNum:
        print (x, end=' ')
        row = []
        for attribute in attributes:
            if attribute == "SOURCEDATALOCATION":
                links.append(str(page['features'][x]['attributes'][attribute]))
            if page['features'][x]['attributes'][attribute] == None:
                row.append('null')
            elif (attribute == "SURVEYDATEUPLOADED"
                or attribute == "SURVEYDATEEND"
                or attribute == "SURVEYDATESTART"):
                    if page['features'][x]['attributes'][attribute] == None:
                        row.append('null')
                    else:
                        date = (page['features'][x]['attributes'][attribute])
                        date = datetime.datetime.utcfromtimestamp(date/1000)
                        row.append(str(date.strftime('%Y-%m-%d')))
            else:
                row.append(str(page['features'][x]['attributes'][attribute]))
        if len(row) != 0:
            rows.append(row[:13])
        x += 1
    print (len(rows), len(links))
    lines = csvCheck()
    if lines != False:
        for line in lines:
            print (line, end=' ')
            if line in rows:
                rows.remove(line) 
            if line[7] in links:
                links.remove(line[7])
    print (len(rows), len(links))
    print ('rows complete')
    return links, rows

def contentSearch(contents, link, saved):
    '''This funtion takes a list of zipfile contents, the download link the
    contents came from and the current local file location where the downloaded
    data resides.

    Using the zipfile contents, it parses the files for any file containing the
    full string '_FULL.xyz'.  If a file name contains this string, it's
    download link and current local file location are added to a global list of
    like files

    If a file is found, it returns a Boolean True
    '''
    x = 0
    for content in contents:
        if full.search(content):
            if [link, saved] not in resPresent:
                print ('\nviva la resolution', content, end=' ') #link + '\n')
                resPresent.append([link, saved])
            x = 1
            return x
        else:
            x = 0

def downloadAndCheck(links, rows):
    '''This function takes a list of download links for survey contents and the
    list of the complete survey data as provided by the query response ('rows').

    For each link provided, it saves the returned data to the computer. All
    downloaded files are expected to be zip files.  The funtion attempts to
    open the files as zipfile objects:
        1) If it succeeds, the function requests a list of included contents
        of the zipfile objects and passes that list, the download link, and
        location of the local download to contentSearch() to determine if the
        desired data exists. The downloaded content is kept.
        2) If it fails, the link, downloaded
        contents, and survey data in 'rows' are immediately removed.

    The function returns the resulting pruned version of the list of the
    complete survey data as provided by the query response ('rows') only
    populated by data positive results.
    '''
    x = len(links)
    for link in links:
        name = link.split('/')[-1]
        saved = holding + '/' + name
        saved = os.path.normpath(saved)
        while True:
            if os.path.exists(saved):
                print  (x, 'x', end=' ')
                break
            else:
                try:
                    urllib.request.urlretrieve(link, saved)
                except socket.timeout:
                    urllib.request.urlretrieve(link, saved)
                except urllib.error.HTTPError:
                    for row in rows:
                        if row[7] == link:
                            print ('e', end=' ')
                            rows.remove(row)
                            links.remove(link)
                    break
        if os.path.exists(saved):
            try:
                zipped = zipfile.ZipFile(saved)
                contents = zipped.namelist()
                if contentSearch(contents, link, saved) != True:
                    print ('n', end=' ')
                    zipped.close()
                    os.remove(saved)
                    for row in rows:
                        if row[7] == link:
                            print ('r', end=' ')
                            row.append('No')
    
                else:
                    zipped.close()
                    for row in rows:
                        if row[7] == link:
                            print ('y', end=' ')
                            row.append('Yes')
            except zipfile.BadZipfile:
                os.remove(saved)
                for row in rows:
                        if row[7] == link:
                            print ('r', end=' ')
                            rows.remove(row)
        x -= 1

    print ('row downloads verified')
    return rows

def csvCheck():
    '''This looks for the CSV text file of currently held results and opens it
    for reading and converts the contents using the python csv library.  A list
    of rows is then compiled by the function and returned as an array 'lines',

    If the CSV text file of currently held results does not exist, the function
    retunrs a boolean False
    '''
    if os.path.isfile(progLoc + '/eHydro_txt.txt'):
        opened = open(progLoc + '/eHydro_txt.txt', 'r')
        csvRead = csv.reader(opened, delimiter = ',')
        lines = []
        for line in csvRead:
            if len(line) < 0:
                lines.append(line)
        opened.close()
        return lines
    else:
        return False


def csvPopulate(rows):
    '''This funtion takes takes the list of refined query results. It also
    checks to see if a current CSV text file of data results exists.  If not,
    it will create one.  It then opens this file for reading and uses the csv
    library to treat the text file as a csv object. The function then:
        1) Compiles a list of the file's current contents
        2) Writes the results to the csv file if they do not already exist
    The CSV text file is then closed
    '''
    if os.path.isfile(progLoc + '/eHydro_txt.txt'):
        opened = open(progLoc + '/eHydro_txt.txt', 'r')
        csvRead = csv.reader(opened, delimiter = ',')
        lines = []
        for line in csvRead:
            lines.append(line)
        opened.close()
        opened = open(progLoc + '/eHydro_txt.txt', 'a')
        txtFile = csv.writer(opened, delimiter = ',')
        for row in rows:
            if row in lines:
                pass
            else:
                txtFile.writerow(row[:13])
        opened.close()
    else:
        opened = open(progLoc + '/eHydro_txt.txt', 'w')
        txtFile = csv.writer(opened, delimiter = ',')
        attributes.append('FULL.xyz?')
        txtFile.writerow(attributes)
        for row in rows:
            txtFile.writerow(row[:13])
        opened.close()
    print ('csv saved')

def emailWriter(rows):
    name = (str(datetime.datetime.today().strftime('%Y%m%d'))
            + '_eHydroUpdates.txt')
    text = open(name, 'w')

def main():
    page, newSurveysNum = query()
    if newSurveysNum != 0:
        links, rows = surveyCompile(page, newSurveysNum)
        returnedRows = downloadAndCheck(links, rows)
        csvPopulate(returnedRows)
#        emailWriter(returnedRows)
    else:
        print ('no new surveys')

main()