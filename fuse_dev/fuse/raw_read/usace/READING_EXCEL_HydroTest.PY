# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 13:28:16 2018
READING_EXCEL_HydroTest.py
@author: jkinney
"""

import pandas as pd
import os, sys
import csv

inputehydrocsv = r"N:\Individual\Juliet\TestLoad\eHydro_cassiano_ex.csv"
inputxls =r"N:\Individual\Juliet\TestLoad\eHydro_csv2_copy.xlsx"#r"N:\Individual\Juliet\TestLoad\eHydro_test.xlsx"#Dummy multisheet spreadsheet
exportfile2 = r'N:\Individual\Juliet\TestLoad\eHydro_csv2b2.txt'
def tryexcel(inputxls):
    dataframe_excel1d ={}
    sheet_to_df_map_dictionary = {}
    sheet_to_df_map_dictionary_hydro = {}
    xls = pd.ExcelFile(inputxls)
    # list all sheets in the file
    xls.sheet_names

    for sheet_name in xls.sheet_names:
        sheet_to_df_map_dictionary_hydro [sheet_name] = xls.parse(sheet_name)
    for sheet_name in xls.sheet_names:
        df_temp_hydro = sheet_to_df_map_dictionary_hydro[sheet_name]    
        headings = list(df_temp_hydro.columns.values)#OR#headings = list(dataframe_excel1)
    return df_temp_hydro, sheet_to_df_map_dictionary_hydro
    
#temp1_df =sheet_to_df_map_dictionary['fish']#can then call be name of sheet if one wants


def trycsv(inputehydrocsv ):
    ehydro_df = pd.read_csv(inputehydrocsv)
    return ehydro_df
##    file1=csv.DictReader(open(inputehydrocsv), dialect='excel', delimiter=',' )#reads in as ordered dictionary
##    fieldnamelist=file1.fieldnames#fieldnames (keys)
##    single_df= pd.DataFrame.from_dict(file1)

#var1= 'CENAN'#NewYork
#var1_District = 'CENAE'#NewEngland
#var_StateP ="North Carolina"#'New York'
#This values above need to be passed not hard coded here, these values will overwrite the passed values if just left as is.

def searchehydrolist(ehydro_df, var1_District, var_StateP =None, DateOlder=None, DateNewer =None):
    subset_Ehydro = ehydro_df[ehydro_df['SURVEYAGENCY'].str.contains(var1_District)]
    subset_Ehydro_StateP = 'None'
    subset_Ehydro_DateStart_range = 'None'
    subset_Ehydro_DateStart_StateP_range = 'None'
    if var_StateP != None:
        subset_Ehydro_StateP = ehydro_df[ehydro_df['SOURCEPROJECTION'].str.contains(var_StateP )]
    if DateOlder != None and  DateNewer != None:
        try:
            subset_Ehydro_DateStart_range = ehydro_df[(ehydro_df.SURVEYDATESTART > DateOlder) & (ehydro_df.SURVEYDATESTART <DateNewer)]#subset_Ehydro_DateStart_range =ehydro_df[(ehydro_df.SURVEYDATESTART > "2018-08-01")]
            print('exporting Date range subset')
            if subset_Ehydro_StateP  is not 'None':
                print('subset_Ehydro_StateP has a value')
                try:
                    subset_Ehydro_DateStart_StateP_range = subset_Ehydro_StateP[(subset_Ehydro_StateP.SURVEYDATESTART > DateOlder) & (subset_Ehydro_StateP.SURVEYDATESTART <DateNewer)]
                except:
                    print('date after subset problem')
        except:
            print('check date format')
    elif DateOlder != None:
        try:
            subset_Ehydro_DateStart_range = ehydro_df[(ehydro_df.SURVEYDATESTART > DateOlder)]
            if subset_Ehydro_StateP is not 'None':
                subset_Ehydro_DateStart_StateP_range = subset_Ehydro_StateP[(subset_Ehydro_StateP.SURVEYDATESTART > DateOlder)]
        except:
            print('check date format')
    
    elif DateNewer != None:
        try:
            subset_Ehydro_DateStart_range = ehydro_df[(ehydro_df.SURVEYDATESTART <DateNewer)]
            if subset_Ehydro_StateP is not 'None':
                subset_Ehydro_DateStart_StateP_range = subset_Ehydro_StateP[(subset_Ehydro_StateP.SURVEYDATESTART > DateNewer)]
        except:
            print('check date format')   
    #Var3 =ehydro_df[(ehydro_df.SURVEYDATEEND > "2018-08-01")]
    #subset_Ehydro_DateStart_range =ehydro_df[(ehydro_df.SURVEYDATESTART > "2017-08-01") & (ehydro_df.SURVEYDATESTART <"2018-01-0")]#subset_Ehydro_DateStart_range =ehydro_df[(ehydro_df.SURVEYDATESTART > "2018-08-01")]
    surveyarea='SDSFEATURENAME'
    surveyname='SURVEYJOBIDPK'
    channelareaid = 'CHANNELAREAIDFK'
    SUREND='SURVEYDATEEND'
    SURSTA='SURVEYDATESTART'
    return subset_Ehydro, subset_Ehydro_StateP, subset_Ehydro_DateStart_range, subset_Ehydro_DateStart_StateP_range
#example

    
# for pandas version < 0.21.0
#this method works to bring all sheets into a dictionary of dataframes:
#ssdf =pd.read_excel(inputxls, sheetname =None)
#pulls first sheet#dataframe_excel1B= pd.read_excel(inputxls, sheet_name = 1)
# for pandas version >= 0.21.0)
#didn't pull all data frames# dataframe_excel1d= pd.read_excel(inputxls, sheet_name = None)#Default reads all sheets to a df dictionary
#dataframe_excel1B = pd.read_excel(inputxls)#Default reads all sheets to a df dictionary
#pulls sheet of this name
