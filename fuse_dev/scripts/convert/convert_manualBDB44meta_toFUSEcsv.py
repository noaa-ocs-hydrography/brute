# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:31:16 2018
This is supposed to convert information from the old manually entered 
tables used to build the NBS database in CARIS BDB 4.1 to the new FUSE
csv structure

@author: jkinney
"""

import meta2csv as m 
import csv 
import pandas as pd
import sys, os
 
_cols = m._cols#list of attributes
errorfile = "error_44to51.txt"

def usace_convert_manualfile():
    """
    usace_convert_manualfile()
    
    hard coding input table location from old BDB 44 manually entered attributes
    output new csv in FUSE directory
        
    """
    inputtab= r"N:\Data\metadata\NewYork\USACE.txt"
    csvfilename = r"N:\New_Directory_1\NorthEast_PBC\USACE\eHydro_NewYorkDistrict\MLLW\Metadata\ehydrometa.csv"
    file1=inputtab
    meta = _manualkeys(file1,  map_columns, out_map_columns)
    m.write_meta2csv(meta, csvfilename)

def fromNCEI_convert_manualfile():
    """
    fromNCEI_convert_manualfile()
    
    hard coding input table location from old BDB 44 manually entered attributes
    output new csv in FUSE directory
        
    """
    inputtab= r"N:\Data\metadata\NewYork\BAG.txt"
    test=r"N:\New_Directory_1\NorthEast_PBC\NOAA_NCEI_OCS\BAGs\MLLW\Metadata\BAG_fuse_meta_test1.csv"
    #test2 = r"N:\Data\metadata\NewYork\BathyPointStore.txt"
    csvfilename = test
    csvfilename2 =r"N:\New_Directory_1\NorthEast_PBC\NOAA_NCEI_OCS\BAGs\MLLW\Metadata\ncei_files__meta_test.csv"
    meta, out_meta = _OCSmanualkeys(inputtab,  map_columns_OCS, map_columns_OCS2, out_map_columns_OCS)
    m.write_meta2csv(meta, csvfilename)
    m.write_meta2csv(out_meta, csvfilename2)

col_NCEI_template = [    
    'from_filename',
    'from_path',
    'script: to_filename',
    'manual: to_filename',
    'script: start_date',
    'manual: start_date',
    'script: end_date',
    'manual: end_date',
    'script: from_fips',
    'manual: from_fips',
    'script: from_horiz_datum',
    'manual: from_horiz_datum',
    'script: from_horiz_units',
    'manual: from_horiz_units',
    'script: from_horiz_unc',
    'manual: from_horiz_unc',
    'script: to_horiz_datum',
    'manual: to_horiz_datum',
    'script: from_vert_datum',
    'manual: from_vert_datum',
    'script: from_vert_key',
    'manual: from_vert_key',
    'script: from_vert_units',
    'manual: from_vert_units',
    'script: from_vert_unc',
    'manual: from_vert_unc',
    'script: to_vert_datum',
    'manual: to_vert_datum',
    'script: to_vert_units',
    'manual: to_vert_units',
    'script: agency',
    'manual: agency',
    'script: source_indicator',
    'manual: source_indicator',
    'script: source_type',
    'manual: source_type',
    'script: complete_coverage',
    'manual: complete_coverage',
    'script: complete_bathymetry',
    'manual: complete_bathymetry',
    'script: vert_uncert_fixed',
    'manual: vert_uncert_fixed',
    'script: vert_uncert_vari',
    'manual: vert_uncert_vari',
    'script: horiz_uncert',
    'manual: horiz_uncert',
    'script: feat_size',
    'manual: feat_size',
    'script: feat_detect',
    'manual: feat_detect',
    'script: feat_least_depth',
    'manual: feat_least_depth',
    'script: interpolated',
    'manual: interpolated',
    'script_version',
    'reviewed',
    'Last Updated',
    'Notes',
    'ID',
    'File Type',
    'Source',
    'script: LeadLine',
    'manual: LeadLine',
    'script: MBES',
    'manual: MBES',
    'script: SSS',
    'manual: SSS',
    'script: SB',
    'manual: SB',
    'script: Lidar',
    'manual: Lidar',
    'MCD Priority',
    'Changeablity',
    'Calculated Rate',
    'TECSOU',
    ]
map_columns = {'from_filename' : 'FileName',
              'start_date' :  'Startdat1',
              'end_date' : 'Enddat1',
              'source_type' : 's_ftyp',
              'complete_coverage' : 'flcvrg',
              'complete_bathymetry' : 'flbath',
              'feat_size' : 'f_size',
              'feat_detect' : 'f_dtct',
              'feat_least_depth' : 'f_lstd',
              'interpolated' : 'interp',
              'reviewed' : 'ManualReview',
              #'Text Channel or Survey Name'
              'from_vert_datum': 'VDatumTransformation',
              'from_vert_units': 'Source Vertical Units',
              'from_horiz_datum' : 'Source Coordval',              
              'Notes' : 'Note',
              'from_vert_key' : 'VERDAT',
              }

map_columns_OCS = {
              'from_filename' : 'ID',
              'end_date' : 'SUREND',
              'source_type' : 'File Type',#'s_ftyp'
              'complete_coverage' : 'flcvrg',
              'complete_bathymetry' : 'flbath',
              'feat_size' : 'f_size',
              'feat_detect' : 'f_dtct',
              'feat_least_depth' : 'f_lstd',
              'Changeablity' : 'Changeablity',
              'from_vert_datum' : 'VERDAT_override',
              'from_vert_key' : 'VERDAT_override',
              'vert_uncert_fixed' : 'vun_fx',
              'vert_uncert_vari' : 'vun_vb',
              'MBES': 'MBES',
              'LeadLine' : 'LeadLine',
              'SSS' : 'SSS',
              'Lidar' : 'Lidar',
              'SB' : 'SB',       
              }
map_columns_OCS2 = {
              'Calculated Rate': 'Calculated Rate',
              'reviewed' : 'Needs Review',
              'reviewed' : 'ManualReview',
              'ID':'ID',
              'Notes':'Note',
              }

place_holder_OCS={
              'start_date' :  'SURSTA',
              #'from_filename' : 'FileName',
              #'from_vert_datum': 'VDatumTransformation',
              'interpolated' : 'interp',
              'from_vert_key':'VERDAT',
              }
out_map_columns = {map_columns[k] : k for k in map_columns}
                #'from_fips':'man_convert_fips'
                #place holder
                #meta['man_convert_fips']=checkforFIPs(meta['Source Coordval'])

out_map_columns_OCS = {map_columns_OCS[k] : k for k in map_columns_OCS}
place_holder ={'script_version' : 's_scpv',
              'source_indicator' : 'SORIND',
              'vert_uncert_fixed' : 'vun_fx',
              'vert_uncert_vari' : 'vun_vb',
              'to_vert_datum' : 'VERDAT',
              'to_vert_units' : 'DUNITS',
              'horiz_uncert' : 'POSACC',
              'to_horiz_datum' : 'HORDAT',
              'agency' : 'AGENCY',
              'source_indicator' : 'SORIND',
              }

Not_mapping_to_attributeyet = ['MatchtoScriptVariable',
                                'ChannelName', 
                                'Text Channel or Survey Name',
                                'ProjectName',
                                'Survey Number',
                                'Survey Type',
                                'Survey Name',
                                'SoundingDensity_rm (ft)',
                                'Benchmark',
                                'RTK',
                                'TECSOU',
                                'SonarSystem',
                                'SonarManufacturer',
                                'Multibeam Angle',
                                'Vessel Name',
                                'Source Sounding Density 3x3',
                                'Sounding Density 5x5',
                                'SurveyInfo',
                                'Collector for USACE',
                                'Source',
                                '',
                                'MCD Priority', 
                                'csar_1',
                                ]

#Full list of attributes:'MatchtoScriptVariable', 'FileName', 'Note', '', 'ChannelName', 'Text Channel or Survey Name', 'ProjectName', 'Enddat1', 'Startdat1', 'VertDat', 'VDatumTransformation', 'Source Vertical Units', 'Source Coordval', 'Survey Number', 'Survey Type', 'Survey Name', 'SoundingDensity_rm (ft)', 'Benchmark', 'RTK', 'TECSOU', 'SonarSystem', 'SonarManufacturer', 'Multibeam Angle', 'Vessel Name', 'Source Sounding Density 3x3', 'Sounding Density 5x5', 'SurveyInfo', 'Collector for USACE', 'Source', '', 'MCD Priority', 'flcvrg', 'f_dtct', 'f_lstd', 'f_size', 'flbath', 'ManualReview', 'csar_1'
def import_usace_manual(inputtab):
    df1=pd.DataFrame()
    file1= csv.DictReader(open(inputtab), delimiter= '\t')
    fieldnamelist=file1.fieldnames
    print(fieldnamelist)
    for r,row in enumerate(file1):
        #where row is an ordered dictionary
        newOrderedDict=row
        metadata=row
        for key in newOrderedDict.keys():
            df1.loc[r,key]=newOrderedDict[key]

def import_BPS_manual(inputtab):
    df1=pd.DataFrame()
    file1= csv.DictReader(open(inputtab), delimiter= '\t')
    fieldnamelist=file1.fieldnames
    print(fieldnamelist)
    for r,row in enumerate(file1):
        #where row is an ordered dictionary
        newOrderedDict=row
        metadata=row
        for key in newOrderedDict.keys():
            df1.loc[r,key]=newOrderedDict[key]

#meta=_manualkeys(metadata, map_columns, out_map_columns)

def _manualkeys(inputtab, map_columns, out_map_columns ):
    """
    Prepend 'manual: ' to each key in the list of dictionaries such that the
    list goes to the right column when written to a csv.
    """
    meta = csv.DictReader(open(inputtab), delimiter= '\t')
    outputkey =[]
    new_meta = []    
    #inputkeys = map_columns.values()
    inputkeys = out_map_columns.keys()
    outputkey = map_columns.keys()
#    try:
#        for row in meta:
#            print(row)
#            print('testing 2')
#    except csv.Error as e:
#        sys.exit('file {}, line {}: {}'.format(errorfile, meta.line_num, e))
#        print('testing')
    try:
        print('in progress')
        for row in meta:
            print(row)
            try:
                print(row['Text Channel or Survey Name'])
            except:
                print('perhaps value was empty')
            print('inside loop')
            keys = []
            keys = row.keys()
            new_row = {}        
            print(keys)
            print(inputkeys)
            print(outputkey)
            for key in keys:
                print('inside key loop')
                print(key)
                if key in inputkeys:
                    print('this attribute is mapped to the new list')
                    new_row[out_map_columns[key]] = row[key]
                    new_row['manual: ' + out_map_columns[key]] = row[key]
                    print(new_row)
                else:
                    print('this attribute is not mapped to the new file')

            new_meta.append(new_row)
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(errorfile, meta.line_num, e))
        print('weird, not working')
    return new_meta
#m.write_meta2csv(new_meta, csvfilename)

#meta=_manualkeys(metadata, map_columns, out_map_columns)
#meta = _OCSmanualkeys(file1,  map_columns_OCS, map_columns_OCS2, out_map_columns_OCS)
def _OCSmanualkeys(inputtab, map_columns_OCS, map_columns_OCS2, out_map_columns_OCS):
    """
    Prepend 'manual: ' to each key in the list of dictionaries such that the
    list goes to the right column when written to a csv.
    """
    meta = csv.DictReader(open(inputtab), delimiter= '\t')
    outputkey =[]
    new_meta = []
    out_meta = []
    #inputkeys = map_columns.values()
    inputkeys = out_map_columns_OCS.keys()
    outputkey = map_columns_OCS.keys()
    out_map_columns_OCS = {map_columns_OCS[k] : k for k in map_columns_OCS}
    col_NCEI_template.sort()
    try:
        print('in progress')
        for row in meta:
            print(row)
            try:
                print(row['Text Channel or Survey Name'])
            except:
                print('perhaps value was empty')
            print('inside loop')
            keys = []
            keys = row.keys()
            new_row = {}        
            print(keys)
            print(inputkeys)
            print(outputkey)
            for key in keys:
                print('inside key loop')
                print(key)
                if key in inputkeys:
                    print('this attribute is mapped to the new list')
                    new_row[out_map_columns_OCS[key]] = row[key]
                    new_row['manual: ' + out_map_columns_OCS[key]] = row[key]
                    print(new_row)
                elif key in out_map_columns_OCS:
                    print('this attribute is mapped to the new list without manual prefix')
                    new_row[out_map_columns_OCS[key]] = row[key]
                    print(new_row)
                else:
                    print('this attribute is not mapped to the new file')
            if 'from_filename' in new_row:                        
                new_row['from_filename'] = new_row['from_filename'] + '.' + new_row['source_type'].lower()# aka'.bag'#maps from 'source_type' : 'File Type'

            new_meta.append(new_row)
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(errorfile, meta.line_num, e))
        print('weird, not working')
    for row in new_meta:
        keys = []
        keys = row.keys()
        new_row = {}
        for key in keys:
            if key in col_NCEI_template:
                new_row[key] = row[key]
        out_meta.append(new_row)
    return new_meta, out_meta



def main():
    fromNCEI_convert_manualfile()
    print('converting NCEI example file')
    usace_convert_manualfile()
    print('converting USACE example file')
    
if __name__ == '__main__':
    main()   
