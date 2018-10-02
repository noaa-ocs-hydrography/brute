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