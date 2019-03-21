# -*- coding: utf-8 -*-
"""
Edited by Juliet Kinney
extract_ehydro_meta_class_CEMNV.py

Created on Wed Aug  8 10:11:36 2018

@author: grice

Parse the eHydro xyz header using the NY region's files as an example.  The
metadata is returned as a dictionary with keys that match the meta2csv naming
convention where required.
"""
__version__ = 'eHydro_meta2csv 0.0.1'

import os as _os
from glob import glob as _glob
from datetime import datetime as _datetime
import re as _re
import pickle as _pickle
from xml.etree.ElementTree import parse as _parse
try:
    import dateutil.parser as dateparser
except:
    pass

#from BDB.BDB51.preprocessing.fips2wkt import fips2wkt as _fips2wkt

_ussft2m = 0.30480060960121924 # US survey feet to meters

import pandas as pd
##-----------------------------------------------------------------------------
"""
new modules
"""
import eHydro_parsing as eH_p#extract_ehydro_meta_CEMVN_testparse1.py#prior version moved to
#from scripts.E_Hydro_DREG_AOI.DREG_info import DREG_info as D_info#brings in dictionaries

#import DREG_info as D_info#brings in dictionaries

import meta2csv as m2c#from BDB.BDB51.preprocessing import meta2csv as m2c
import csv as _csv
import Extract_Meta_Class as E_M_C#import meta_review_base as E_M_C#
import numpy as _np
##-----------------------------------------------------------------------------

class read_raw_cemvn:
    
    def read_metadata(self, infilename, inputehydrocsv):
        """
        Read all available meta data.
        """
        version='CEMVN'
        self.version = version
        return retrieve_meta_for_Ehydro_out_onefile(infilename, inputehydrocsv)
    
    def read_bathymetry_dat(self, infilename):
        """
        Read the bathymetry.
        """
        # get the dat file for CEMVN
        stub, ext = _os.path.splitext(infilename)
        bathyfilename = stub + '.dat'
        xyz = _np.loadtxt(bathyfilename, delimiter = ' ')
        self.xyz
        return xyz
    
    def read_bathymetry(self, infilename):
        version='CEMVN'
        self.version = version
        #xyz1 = []
        first_instance = eH_p._start_xyz(infilename, version = None)
        if first_instance is not '':    
            xyz = _np.loadtext(infilename, delimeter = ',', skiprows = first_instance)
        else:
            xyz = _np.loadtext(infilename, delimeter = ',')
        return xyz
        #    for line in infile.readlines():
        #        if line == '\n':
        #            continue
        #        elif eH_p._is_xyz_data(line, version):
        #            xyz1.append(line)
        #        else:
        #            break
        #xyz = _np.asarray(xyz1)
        
    def read_pass_meta_tocsv(self, infilename, inputehydrocsv, outputfilename = None):
        if outputfilename == None:
            outputfilename = 'ehydro_meta_dict_to_csv_v1.txt'
            #print("using default in active directory 'ehydro_meta_dict_to_csv_v1.txt'")
            merged_meta, dataframe_nn = self.read_metadata(self, infilename, inputehydrocsv)
            self.meta_dict = merged_meta
        m2c.write_meta2csv_CEMVN(self.meta_dict, outputfilename)
        
    def passdict_tocsv(self, merged_meta, outputfilename = None):
        if outputfilename == None:
            outputfilename = 'ehydro_meta_dict_to_csv_v2.txt'
            #print("using default in active directory 'ehydro_meta_dict_to_csv_v1.txt'")
        m2c.write_meta2csv_CEMVN(merged_meta, outputfilename)               
    
def parse_ehydro_directory(inpath):
    """
    Parse a directory of USACE xyz files from eHydro and return a list of the
    dictionaries containing the meta data from each.
    """
    metalist = []
    flist = _glob(_os.path.join(inpath, '*.xyz'))
    flist.sort()
    for f in flist:
        e_t = Extract_Txt(f)
        meta = e_t.parse_ehydro_xyz(f, meta_source = 'xyz', version='CEMVN', default_meta = '')#parse_ehydro_xyz(infilename, meta_source = 'xyz', version='CEMVN', default_meta = '')#eH_p.parse_ehydro_xyz(f)
        metalist.append(meta)
    return metalist

#
#class extract_xml_meta(filename):
#    def __init__(self,filename):
#        self.filename = filename
#        #methods for xml extraction and attribute mapping
#        
##class dates(datestring):
#    
##class xyz_header(text_object):
#    
#
#class extract_meta(filename):
#    def __init__(self,filename):
#        self.filename = filename
#        #self.district=
#        #self.meta_source
#        #self.version
#        #self.default_meta

def getsubset(highresfolder, ehydrofolder):
    g1 = ehydro_subset_only_ifnotfull_restoo(highresfolder, ehydrofolder)
    return g1
    

def retrieve_meta_for_Ehydro_out_onefile(filename, inputehydrocsv):

    ehydro_df = trycsv(inputehydrocsv)#bring in csv from ehydro website
    #ehydro_df = Rcsv.trycsv(inputehydrocsv)#folders?
    #next if pull the subset of the table in the dataframe related to the list of files passed to it.
    merged_meta = {}
    merge2 = {}
    #merged_meta_rows = {}    
    ehydro_table = Extract_Table(ehydro_df,filename=inputehydrocsv)
    nn = pd.DataFrame()
    nm = pd.DataFrame()
    nnn = pd.DataFrame()
    #g1.sort()
    f = filename
    basename = os.path.basename(f)
    ex_string1 = '*_A.xyz'
    ex_string2 = '*_FULL.xyz'
    ex_string3 = '*_FULL.XYZ'
    ex_string4 = '*_A.XYZ'
    basename = os.path.basename(basename)
    basename = S_f_d.return_surveyid(basename, ex_string1)
    basename = S_f_d.return_surveyid(basename, ex_string2)
    basename = S_f_d.return_surveyid(basename, ex_string3)
    basename = S_f_d.return_surveyid(basename, ex_string4)
    basename = basename.rstrip('.XYZ')  
    basename = basename.rstrip('.xyz')
    meta_from_ehydro, hold_meta2 = ehydro_table.pull_df_by_dict_key_c(basename, searchvalue = None, meta=None, hold_meta2 = None, version = 'ehydro_csv')
    #metadata = ehydro_table.pull_df_by_dict_key(ehydro_df, thisdictionary, basename, searchvalue = None, meta=None, version = None)#other option version = 'casiano_ehydro_csv'
    #ehydro_table.pull_df_by_dict_key(thisdataframe, thisdictionary, basename, searchvalue = None, meta=None, version = None)
    e_t = Extract_Txt(f)
    # xml pull here.
    one_file, v = E_M_C.use_extract_meta(f)#E_M_C.use_extract_meta(test_file_path)
    #since we know its ehydro:
    xmlfilename = one_file.get_xml()
    if os.path.isfile(xmlfilename):
        ###may be useful to add in, right now needs more work for S57 matching      
        ##try:
        ##    e_xml_s57dict = p_usace_xml.extract_s57_dict(xmlfilename)
        ##except:
        ##    print('is this expected path followed?')
        with open(xmlfilename, 'r') as xml_file:
            xml_txt = xml_file.read()
        xmlbasename = os.path.basename(xmlfilename)
        xml_data = p_usace_xml.XML_Meta(xml_txt, filename = xmlbasename)
        if xml_data.version == 'USACE_FGDC':
            meta_xml = xml_data._extract_meta_CEMVN()
        elif xml_data.version == 'ISO-8859-1':
                meta_xml = xml_data._extract_meta_USACE_ISO()
        else:
            meta_xml = xml_data.convert_xml_to_dict2()#some_meta = xml_data.convert_xml_to_dict()
    #relates to Bathy Class get_metadata
        
    meta = e_t.parse_ehydro_xyz(f, meta_source = 'xyz', version='CEMVN', default_meta = '')#parse_ehydro_xyz(infilename, meta_source = 'xyz', version='CEMVN', default_meta = '')#eH_p.parse_ehydro_xyz(f)
    #NEED to see which version can come in #e_xml = Extract_Xml(preloadeddata, version = '', filename = '')
    #NEED to see how to reconcile different possible inputs in QA process.
    list_keys_empty =[]
    combined_row = {}
    subset_row = {}
    subset_no_overlap = {}
    for key in meta:
        if meta[key] == 'unknown' or  meta[key] == '':
            list_keys_empty.append(key)
        else:
            subset_row[key] = meta[key]
            #non blank columns only
            if key in meta_xml:
                if meta[key] == meta_xml[key]:
                    combined_row[key] = meta[key]
                    """
                    only make list within cell if values from different sources are different
                    """
                else:
                    combined_row[key] = meta[key] + ' , ' + meta_xml[key]
            else:
                subset_no_overlap[key] = meta[key]
        #if key in meta_xml:
        #    if list_keys_empty not in meta_xml:
        #        combined_row[key] = meta[key] + ' , ' + meta_xml[key]
        #longer way to say same as above
   #nn_df.to_csv(path_or_buf=r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\df2.txt', encoding='UTF-8', sep ='\t')         
    mydict = meta_xml
    for i0, key1 in enumerate(mydict):#for key1 in mydict:
        key_fromdict=key1
        #use common numberic index i0
        if key_fromdict in combined_row:
            nn.loc[f, key_fromdict] = combined_row[key_fromdict]#[row, column]
        else:
            nn.loc[f, key_fromdict] = mydict[key_fromdict]#[row, column]
        #print(type(mydict[key_fromdict]))
        #nm.loc[i0,key_fromdict] = key_fromdict
    for i0, key1 in enumerate(subset_no_overlap):
        key_fromdict=key1
        nn.loc[f, key_fromdict] = subset_row[key_fromdict]
    for i0, key1 in enumerate(meta_from_ehydro):
        key_fromdict=key1
        nn.loc[f, key_fromdict] = meta_from_ehydro[key_fromdict]
    
    merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row }   #merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row }       
    merged_meta = { **meta, **meta_from_ehydro,**meta_xml }
    #THIS METHOD FOR MERGING DICTIONARIES IS FOR Python 3.5  plus based on PEP 448
    # https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/ 
    #see for more informatio`n
    #merged_dictonary = {**default_values, **override_if_duplicatekeys}# 
    
    #merged_meta = {**e_xml_s57dict, **meta_from_ehydro, **meta}
    #need to fix mapping logic to input to BDB table for xml_meta,
    #most attributes aren't needed at this time. But we can pull them.
    ###m2c.write_meta2csv([merged_meta],meta1csv)
    #m2c.write_meta2csv([merged_meta],metafile1)
    ###m2c.write_meta2csv([merge2],metafile1)
    #m2c.write_meta2csv_full_tab([merged_meta],metafile)#Trying to trouble shoot.
    #write_to_csv([merged_meta],metafile)
    #save pandas dataframe export here
    ###nn.to_csv(path_or_buf=(df_export_to_csv), encoding='UTF-8', sep ='\t')
    #nn.to_csv(path_or_buf=r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\Attempted_combined_df_metafields.txt', encoding='UTF-8', sep ='\t')
    return merged_meta, nn


def retrieve_meta_for_Ehydro_out_df(g1, inputehydrocsv, metafile1, df_export_to_csv):

    ehydro_df = trycsv(inputehydrocsv)#bring in csv from ehydro website
    #ehydro_df = Rcsv.trycsv(inputehydrocsv)#folders?
    #next if pull the subset of the table in the dataframe related to the list of files passed to it.
    merged_meta = {}
    merge2 = {}
    #merged_meta_rows = {}    
    ehydro_table = Extract_Table(ehydro_df,filename=inputehydrocsv)
    nn = pd.DataFrame()
    nm = pd.DataFrame()
    nnn = pd.DataFrame()
    g1.sort()
    ex_string1 = '*_A.xyz'
    ex_string2 = '*_FULL.xyz'
    ex_string3 = '*_FULL.XYZ'
    ex_string4 = '*_A.XYZ'
    for basename in g1:
        f = basename
        basename = os.path.basename(basename)
        basename = S_f_d.return_surveyid(basename, ex_string1)
        basename = S_f_d.return_surveyid(basename, ex_string2)
        basename = S_f_d.return_surveyid(basename, ex_string3)
        basename = S_f_d.return_surveyid(basename, ex_string4)
        basename = basename.rstrip('.xyz')
        basename = basename.rstrip('.XYZ') 
        meta_from_ehydro, hold_meta2 = ehydro_table.pull_df_by_dict_key_c(basename, searchvalue = None, meta=None, hold_meta2 = None, version = 'ehydro_csv')
        #metadata = ehydro_table.pull_df_by_dict_key(ehydro_df, thisdictionary, basename, searchvalue = None, meta=None, version = None)#other option version = 'casiano_ehydro_csv'
        #ehydro_table.pull_df_by_dict_key(thisdataframe, thisdictionary, basename, searchvalue = None, meta=None, version = None)
        e_t = Extract_Txt(f)
        # xml pull here.
        one_file, v = E_M_C.use_extract_meta(f)#E_M_C.use_extract_meta(test_file_path)
        #since we know its ehydro:
        xmlfilename = one_file.get_xml()
        if os.path.isfile(xmlfilename):        
            try:
                e_xml_s57dict = p_usace_xml.extract_s57_dict(xmlfilename)
            except:
                print('is this expected path followed?')
            with open(xmlfilename, 'r') as xml_file:
                xml_txt = xml_file.read()
            xmlbasename = os.path.basename(xmlfilename)
            xml_data = p_usace_xml.XML_Meta(xml_txt, filename = xmlbasename)
            if xml_data.version == 'USACE_FGDC':
                meta_xml = xml_data._extract_meta_CEMVN()
            else:
                meta_xml = xml_data.convert_xml_to_dict2()#some_meta = xml_data.convert_xml_to_dict()
        #relates to Bathy Class get_metadata
            
        meta = e_t.parse_ehydro_xyz(f, meta_source = 'xyz', version='CEMVN', default_meta = '')#parse_ehydro_xyz(infilename, meta_source = 'xyz', version='CEMVN', default_meta = '')#eH_p.parse_ehydro_xyz(f)
        #NEED to see which version can come in #e_xml = Extract_Xml(preloadeddata, version = '', filename = '')
        #NEED to see how to reconcile different possible inputs in QA process.
        list_keys_empty =[]
        combined_row = {}
        subset_row = {}
        subset_no_overlap = {}
        for key in meta:
            if meta[key] == 'unknown' or  meta[key] == '':
                list_keys_empty.append(key)
            else:
                subset_row[key] = meta[key]
                #non blank columns only
                if key in meta_xml:
                    if meta[key] == meta_xml[key]:
                        combined_row[key] = meta[key]
                        """
                        only make list within cell if values from different sources are different
                        """
                    else:
                        combined_row[key] = meta[key] + ' , ' + meta_xml[key]
                else:
                    subset_no_overlap[key] = meta[key]
            #if key in meta_xml:
            #    if list_keys_empty not in meta_xml:
            #        combined_row[key] = meta[key] + ' , ' + meta_xml[key]
            #longer way to say same as above
       #nn_df.to_csv(path_or_buf=r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\df2.txt', encoding='UTF-8', sep ='\t')         
        mydict = meta_xml
        for i0, key1 in enumerate(mydict):#for key1 in mydict:
            key_fromdict=key1
            #use common numberic index i0
            if key_fromdict in combined_row:
                nn.loc[f, key_fromdict] = combined_row[key_fromdict]#[row, column]
            else:
                nn.loc[f, key_fromdict] = mydict[key_fromdict]#[row, column]
            #print(type(mydict[key_fromdict]))
            #nm.loc[i0,key_fromdict] = key_fromdict
        for i0, key1 in enumerate(subset_no_overlap):
            key_fromdict=key1
            nn.loc[f, key_fromdict] = subset_row[key_fromdict]
        for i0, key1 in enumerate(meta_from_ehydro):
            key_fromdict=key1
            nn.loc[f, key_fromdict] = meta_from_ehydro[key_fromdict]
        
        merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row }   #merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row }       
        merged_meta = { **meta, **meta_from_ehydro,**meta_xml }
        #THIS METHOD FOR MERGING DICTIONARIES IS FOR Python 3.5  plus based on PEP 448
        # https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/ 
        #see for more informatio`n
        #merged_dictonary = {**default_values, **override_if_duplicatekeys}# 
        
        #merged_meta = {**e_xml_s57dict, **meta_from_ehydro, **meta}
        #need to fix mapping logic to input to BDB table for xml_meta,
        #most attributes aren't needed at this time. But we can pull them.
        #m2c.write_meta2csv([merged_meta],metafile1)
        m2c.write_meta2csv([merge2],metafile1)
        #m2c.write_meta2csv_full_tab([merged_meta],metafile)#Trying to trouble shoot.
        #write_to_csv([merged_meta],metafile)
    #save pandas dataframe export here
    nn.to_csv(path_or_buf=(df_export_to_csv), encoding='UTF-8', sep ='\t')
    #nn.to_csv(path_or_buf=r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\Attempted_combined_df_metafields.txt', encoding='UTF-8', sep ='\t')
    return merged_meta, nn

def retrieve_meta_for_Ehydro(highresfolder, ehydrofolder, inputehydrocsv, metafile1, metafile, df_export_to_csv):
    g1 = ehydro_subset_only_ifnotfull_restoo(highresfolder, ehydrofolder)
    ehydro_df = trycsv(inputehydrocsv)#bring in csv from ehydro website
    #ehydro_df = Rcsv.trycsv(inputehydrocsv)#folders?
    #next if pull the subset of the table in the dataframe related to the list of files passed to it.
    merged_meta = {}
    merge2 = {}
    #merged_meta_rows = {}    
    ehydro_table = e_meta.Extract_Table(ehydro_df,filename=inputehydrocsv)
    nn = pd.DataFrame()
    nm = pd.DataFrame()
    nnn = pd.DataFrame()
    g1.sort()
    ex_string1 = '*_A.xyz'
    ex_string2 = '*_FULL.xyz'
    ex_string3 = '*_FULL.XYZ'
    ex_string4 = '*_A.XYZ'
    for basename in g1:
        f = basename
        basename = os.path.basename(basename)
        basename = S_f_d.return_surveyid(basename, ex_string1)
        basename = S_f_d.return_surveyid(basename, ex_string2)
        basename = S_f_d.return_surveyid(basename, ex_string3)
        basename = S_f_d.return_surveyid(basename, ex_string4)
        basename = basename.rstrip('.xyz')
        basename = basename.rstrip('.XYZ')        
        meta_from_ehydro, hold_meta2 = ehydro_table.pull_df_by_dict_key_c(basename, searchvalue = None, meta=None, hold_meta2 = None, version = 'ehydro_csv')
        #metadata = ehydro_table.pull_df_by_dict_key(ehydro_df, thisdictionary, basename, searchvalue = None, meta=None, version = None)#other option version = 'casiano_ehydro_csv'
        #ehydro_table.pull_df_by_dict_key(thisdataframe, thisdictionary, basename, searchvalue = None, meta=None, version = None)
        e_t = e_meta.Extract_Txt(f)
        # xml pull here.
        one_file, v = E_M_C.use_extract_meta(f)#E_M_C.use_extract_meta(test_file_path)
        #since we know its ehydro:
        if '_A.xyz' in f:
            xmlfilename = one_file.get_xml_xt('_A.xyz')
        elif '_FULL.xyz' in f:
            xmlfilename = one_file.get_xml_xt('_FULL.xyz')
        elif '_FULL.XYZ' in f:
            xmlfilename = one_file.get_xml_xt('_FULL.XYZ')
        elif '_A.XYZ' in f:
            xmlfilename = one_file.get_xml_xt('_A.XYZ')
        else:
            xmlfilename = one_file.get_xml()
        if os.path.isfile(xmlfilename):        
            try:
                e_xml_s57dict = p_usace_xml.extract_s57_dict(xmlfilename)
            except:
                print('is this expected path followed?')
            with open(xmlfilename, 'r') as xml_file:
                xml_txt = xml_file.read()
            xmlbasename = os.path.basename(xmlfilename)
            xml_data = p_usace_xml.XML_Meta(xml_txt, filename = xmlbasename)
            if xml_data.version == 'USACE_FGDC':
                meta_xml = xml_data._extract_meta_CEMVN()
            elif xml_data.version == 'ISO-8859-1':
                meta_xml = xml_data._extract_meta_USACE_ISO()
            else:
                meta_xml = xml_data.convert_xml_to_dict2()#some_meta = xml_data.convert_xml_to_dict()
        #relates to Bathy Class get_metadata
            
        meta = e_t.parse_ehydro_xyz(f, meta_source = 'xyz', version='CEMVN', default_meta = '')#parse_ehydro_xyz(infilename, meta_source = 'xyz', version='CEMVN', default_meta = '')#eH_p.parse_ehydro_xyz(f)
        #NEED to see which version can come in #e_xml = Extract_Xml(preloadeddata, version = '', filename = '')
        #NEED to see how to reconcile different possible inputs in QA process.
        list_keys_empty =[]
        combined_row = {}
        subset_row = {}
        subset_no_overlap = {}
        for key in meta:
            if meta[key] == 'unknown' or  meta[key] == '':
                list_keys_empty.append(key)
            else:
                subset_row[key] = meta[key]
                #non blank columns only
                if key in meta_xml:
                    if meta[key] == meta_xml[key]:
                        combined_row[key] = meta[key]
                        """
                        only make list within cell if values from different sources are different
                        """
                    else:
                        combined_row[key] = meta[key] + ' , ' + meta_xml[key]
                else:
                    subset_no_overlap[key] = meta[key]
            #if key in meta_xml:
            #    if list_keys_empty not in meta_xml:
            #        combined_row[key] = meta[key] + ' , ' + meta_xml[key]
            #longer way to say same as above
       #nn_df.to_csv(path_or_buf=r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\df2.txt', encoding='UTF-8', sep ='\t')         
        mydict = meta_xml
        for i0, key1 in enumerate(mydict):#for key1 in mydict:
            key_fromdict=key1
            #use common numberic index i0
            if key_fromdict in combined_row:
                nn.loc[f, key_fromdict] = combined_row[key_fromdict]#[row, column]
            else:
                nn.loc[f, key_fromdict] = mydict[key_fromdict]#[row, column]
            #print(type(mydict[key_fromdict]))
            #nm.loc[i0,key_fromdict] = key_fromdict
        for i0, key1 in enumerate(subset_no_overlap):
            key_fromdict=key1
            nn.loc[f, key_fromdict] = subset_row[key_fromdict]
        for i0, key1 in enumerate(meta_from_ehydro):
            key_fromdict=key1
            nn.loc[f, key_fromdict] = meta_from_ehydro[key_fromdict]
        
        merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row }   #merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row }       
        merged_meta = { **meta, **meta_from_ehydro,**meta_xml }
        #THIS METHOD FOR MERGING DICTIONARIES IS FOR Python 3.5  plus based on PEP 448
        # https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/ 
        #see for more informatio`n
        #merged_dictonary = {**default_values, **override_if_duplicatekeys}# 
        
        #merged_meta = {**e_xml_s57dict, **meta_from_ehydro, **meta}
        #need to fix mapping logic to input to BDB table for xml_meta,
        #most attributes aren't needed at this time. But we can pull them.
        #m2c.write_meta2csv([merged_meta],metafile1)
        m2c.write_meta2csv([merge2],metafile1)
        #m2c.write_meta2csv_full_tab([merged_meta],metafile)#Trying to trouble shoot.
        #write_to_csv([merged_meta],metafile)
    #save pandas dataframe export here
    nn.to_csv(path_or_buf=(df_export_to_csv), encoding='UTF-8', sep ='\t')
    #nn.to_csv(path_or_buf=r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\Attempted_combined_df_metafields.txt', encoding='UTF-8', sep ='\t')
    return merged_meta, nn

      
class Extract_Table:
    """
    metadata={}
    ehydro_table =Extract_Table(ehydro_df)  
    """
    def __init__(self, preloadeddata, version = '', filename = ''):

        self.filename = filename
        self.ehydro_pre_download_table=preloadeddata
        #self.filepath = (filenamepath)
        
#    def link_ehydro_table(self, metadata):#ehydro_pre_download_table
#        #dreg_dicts, dreg_d=D_info.DREGi(input1)
#        metadata['SUREND'] = self.ehydro_pre_download_table.SURVEYDATEEND
#        SUREND = self.ehydro_pre_download_table.SURVEYDATEEND
#        #SURSTA= self.ehydro_pre_download_table.SURVEYDATESTART
#        metadata['SOURCEPROJECTION'] = self.ehydro_pre_download_table.SOURCEPROJECTION
#        #FIPS = convert_tofips(SOURCEPROJECTION)
#        #'Survey Type': 'SURVEYTYPE'
#        self.metadata=metadata
#        return SUREND, metadata
    
    def pull_df_by_dict_key(thisdataframe, thisdictionary, basename, searchvalue = None, meta=None, version = None):
        if meta == None:
            meta={}
        if searchvalue == None:
            if version == None:
                searchvalue = 'SURVEYJOBIDPK'
                print('assuming ehydro')
        if version == 'ehydro_csv':
                searchvalue = 'SURVEYJOBIDPK'
        if version == 'casiano_ehydro_csv':
                searchvalue = 'SURVEYJOBIDPK'
        subsetDataFrame = thisdataframe[thisdataframe[searchvalue].str.contains(basename)==True]
        for thesekeys in thisdictionary:
            #meta[thisdictionary[thesekeys]]=thisdataframe[thesekeys]
            meta[thesekeys]=subsetDataFrame[thisdictionary[thesekeys]]
        return meta
    
    
    def pull_df_by_dict_key_c(self, basename, searchvalue = None, meta=None, hold_meta2 = None, version = None):
        """
        self.thisdataframe
        self.thisdictionary
        meta, hold_meta2 = pull_df_by_dict_key_c(self, basename, searchvalue = None, meta=None, hold_meta2 = None, version = None)
        Expansion of:
        meta = pull_df_by_dict_key(thisdataframe, thisdictionary, basename, searchvalue = None, meta=None, version = None)
        """
        v=0
        if meta == None:
            meta={}
        if hold_meta2 == None:
            hold_meta2 = {}
        if searchvalue == None:
            if version == None:
                searchvalue = 'SURVEYJOBIDPK'
                print('assuming ehydro')
        if version == 'ehydro_csv':
                searchvalue = 'SURVEYJOBIDPK'
                thisdataframe = self.ehydro_pre_download_table
                thisdictionary = self.e_hydro_to_output_dict
                thisdictionary2 = self.e_hydro_crosstable_dict
                v=1
        if version == 'casiano_ehydro_csv':
                searchvalue = 'SURVEYJOBIDPK'
                thisdataframe = self.ehydro_pre_download_table
                thisdictionary = self.e_hydro_to_output_dict
                thisdictionary2 = self.e_hydro_casiano_crosstable_dict
                v=2
        #if version == None:
        #    print('version is unknown will return empty dictionaries')
        if v== 1 or 2:
            subsetDataFrame = thisdataframe[thisdataframe[searchvalue].str.contains(basename)==True]
            for thesekeys in thisdictionary:
                """
                output to Glen's table format
                """
                #meta[thisdictionary[thesekeys]]=thisdataframe[thesekeys]
                for values in subsetDataFrame[thisdictionary[thesekeys]]:
                    meta[thesekeys]=values
                if v == 1 or 2:
                    for sourceprojectionstr in subsetDataFrame['SOURCEPROJECTION']:
                        meta['script: from_fips'] = eH_p.convert_tofips( eH_p.SOURCEPROJECTION_dict ,sourceprojectionstr)#SOURCEPROJECTION
            meta['script: start_date'] = self.convert_datefrmt(meta['script: start_date'])
            meta['script: end_date'] = self.convert_datefrmt(meta['script: end_date'])
            """
            converting date string to format for S-57/ MQUAL YYYYMMDD
            """
            if thisdictionary2:
                for thesekeys in thisdictionary2:
                    hold_meta2[thesekeys]=subsetDataFrame[thisdictionary2[thesekeys]]
                    """
                    holding values extended values/ or names for values option in a dictionary
                    """
        return meta, hold_meta2


    def convert_datefrmt(self, date1):
        parsed_date = dateparser.parse(date1)
        ymd_date = parsed_date.strftime('%Y%m%d')#('%Y-%m-%dT%H:%M:%SZ')
        return ymd_date
    
#    def link_ehydro_table(ehydro_pre_download_table, metadata):        
#        #dreg_dicts, dreg_d=D_info.DREGi(input1)
#        metadata['SUREND'] = ehydro_pre_download_table.SURVEYDATEEND
#        SUREND = ehydro_pre_download_table.SURVEYDATEEND
#        #SURSTA= ehydro_pre_download_table.SURVEYDATESTART
#        SOURCEPROJECTION = ehydro_pre_download_table.SOURCEPROJECTION
#        #FIPS = convert_tofips(SOURCEPROJECTION)
#        #'Survey Type': 'SURVEYTYPE'
#        return SUREND, metadata
    
    e_hydro_crosstable_dict ={}
    e_hydro_crosstable_dict = {
            'Survey Type' : 'SURVEYTYPE',
            'script: from_horiz_datum':'SOURCEPROJECTION',#Source Coordval
            'ProjectName': 'SDSFEATURENAME',
            'SUREND':'SURVEYDATEEND',
            'script: end_date':'SURVEYDATEEND',
            #'script: from_fips': 'FIPS',
            'script: agency':'SURVEYAGENCY',
            #'script: agency':'AGENCY',
            'filenamebase':'SURVEYJOBIDPK',
            'EHydro_Channel_ID': 'CHANNELAREAIDFK',
            'script: start_date':'SURVEYDATESTART',
            'SURSTA':'SURVEYDATESTART',
            #'script: Hi-Res':'Hi-Res?',
            'link':'SOURCEDATALOCATION',
            'SURVEYDATEUPLOADED':'SURVEYDATEUPLOADED',}
    
    e_hydro_casiano_crosstable_dict ={}
    e_hydro_casiano_crosstable_dict = {
            'Survey Type' : 'SURVEYTYPE',
            'script: from_horiz_datum':'SOURCEPROJECTION',#Source Coordval
            'ProjectName': 'SDSFEATURENAME',
            'SUREND':'SURVEYDATEEND',
            'script: end_date':'SURVEYDATEEND',
            #'script: from_fips': 'FIPS',
            'script: agency':'SURVEYAGENCY',
            #'script: agency':'AGENCY',
            'filenamebase':'SURVEYJOBIDPK',
            'EHydro_Channel_ID': 'CHANNELAREAIDFK',
            'script: start_date':'SURVEYDATESTART',
            'SURSTA':'SURVEYDATESTART',
            'script: Hi-Res':'Hi-Res?',
            'link':'SOURCEDATALOCATION',
            'SURVEYDATEUPLOADED':'SURVEYDATEUPLOADED',}
    
    e_hydro_to_output_dict ={}
    e_hydro_to_output_dict = {
            'script: from_horiz_datum':'SOURCEPROJECTION',#Source Coordval
            'ProjectName': 'SDSFEATURENAME',
            'script: end_date':'SURVEYDATEEND',
            #'script: from_fips': 'FIPS',
            'script: agency':'SURVEYAGENCY',
            'script: start_date':'SURVEYDATESTART',
            'from_filenamebase':'SURVEYJOBIDPK',
            #'script: Hi-Res':'Hi-Res?',
            }
###---------------------------------------------------------------------------- 
"""
#USING THE CLASS       
metadata={}
ehydro_table =Extract_Table(ehydro_df)  
"""

###---------------------------------------------------------------------------- 
class Extract_Txt(object):
    
    """
    """
    def __init__(self, preloadeddata, version = '', filename = ''):
        self.filename = preloadeddata
        if filename is not "" or None:
            self.filename_1 = filename
            self.errorfile = _os.path.dirname(filename) + 'TEST_extract_ehdyro_meta_class_CEMVN_ErrorFile.txt'
        else:
            self.errorfile = _os.path.dirname(filename) + 'Default_extract_ehdyro_meta_class_CEMVN.txt'
        #self.filepath = (filenamepath)        
    ##NY one
    def parse_ehydro_xyz(self, infilename, meta_source = 'xyz', version= 'CEMVN', default_meta = ''):#need to change version to None
        """
        'CEMVN'
        """
        """
        Parse an USACE eHydro file for the available meta data.
        
        Default metadata (values predetermined for the file but not in the file)
        can be stored at the location defined by 'default_meta' as a pickled
        dicitonary.  If no path is provided the dictionary in the same folder as
        the data but in the file 'default.pkl' will be used.  If this file does
        not exist no default metadata will be loaded.  If the same keyword for
        the metadata exists both in the file metadata and in the default location,
        the file metadata will take precidence.
        """
        name_meta = self.parse_ehydro_filename(infilename)#moved from eH_p
        if meta_source == 'xyz':
            file_meta = self.parse_xyz_header(infilename, version)#moved from eH_p
        #elif meta_source == 'xml':
        #    file_meta = parse_ehydro_xml(infilename)
        default_meta = self.load_default_metadata(infilename, default_meta)#moved from eH_p
        merged_meta = {**default_meta, **name_meta, **file_meta}
        if 'from_horiz_unc' in merged_meta:
            if merged_meta['from_horiz_units'] == 'US Survey Foot':
                val = _ussft2m * float(merged_meta['from_horiz_unc'])
                merged_meta['horiz_uncert'] = val
        if 'from_vert_unc' in merged_meta:
            if merged_meta['from_vert_units'] == 'US Survey Foot':
                val = _ussft2m * float(merged_meta['from_vert_unc'])
                merged_meta['vert_uncert_fixed'] = val
                merged_meta['vert_uncert_vari'] = 0
        sorind = (name_meta['projid'] + '_' + 
                  name_meta['uniqueid'] + '_' + 
                  name_meta['subprojid'] + '_' + 
                  name_meta['start_date'] + '_' + 
                  name_meta['statuscode'])
        merged_meta['source_indicator'] = 'US,US,graph,' + sorind
        merged_meta['script_version'] = __version__
        return merged_meta

    def parse_ehydro_filename(self, infilename):
        """
        Parse the provided infilename for the channel project code, unique id,
        subproject code, survey acquistion start date, the survey code, and
        optional field and return a dictionary of these fields.  The dictionary
        contains the following keys:
            projid : the project id code
            uniqueid : the unique identifier, and can represent items such as
                        separate projects, map scale, or CWIS code.
            subprojid : the subproject within the channel project
            startdate : the start date of acquisition
            statuscode : survey description code, such as for a condition survey
            optional : this is the contents of the condition field
            from_path : this is named to match other scripts downstream
            from_filename : this is also named to match other file downstream
        """
        base = _os.path.basename(infilename)
        name, ext = _os.path.splitext(base)
        splitname = name.split('_')
        
        if len(splitname) >= 4:
            meta = {
                    'from_path' : infilename,
                    'from_filename' : base,
                    'projid' : splitname[0],
                    'uniqueid' : splitname[1],
                    'subprojid' : splitname[2],
                    'start_date' : splitname[3],
                    }
            #'statuscode' : splitname[4],
            if len(splitname) >4:
                meta['statuscode']=splitname[4]
                if len(splitname) > 5:
                    option = splitname[5]
                    if len(splitname) > 6:
                        for n in range(6, len(splitname)):
                            option = option + '_' + splitname[n]
                    meta['optional'] = option
            else:
                meta['statuscode']=''
                
        else:
            print(name + ' appears to have a nonstandard naming convention.')
        return meta
    
    def parse_xyz_header(self, infilename, version=None):
        """
        Parse the xyz file header for meta data and return a dictionary.  The
        key words used to search are
            NOTES
            PROJECT_NAME
            SURVEY_NAME
            DATES_OF_SURVEY
        """
        header = []
        metalist = []
        more_metalist = []
        # get the header
        if version == 'CEMVN':
            with open(infilename, 'r') as infile:
                for line in infile.readlines():
                    if line == '\n':
                        continue
                    elif eH_p._is_header2(line):
                        header.append(line)
                    else:
                        break
            # search the header for lines starting with the key words
            for line in header:
                if line.startswith('NOTES'):
                    metalist.append(eH_p._parse_note(line))
                elif line.startswith('PROJECT_NAME'):#not expecting
                    metalist.append(eH_p._parse_projectname(line))
                elif line.startswith('SURVEY_NAME'):#not expecting
                    metalist.append(eH_p._parse_surveyname(line))
                elif line.startswith('DATES_OF_SURVEY'):#not expecting
                    metalist.append(eH_p._parse_surveydates(line))
                elif line.startswith('SOUNDING_FREQUENCY'):
                    metalist.append(eH_p._parse_sounding_frequency(line))
                elif line.startswith('SURVEY_TYPE'):
                    metalist.append(eH_p._parse_survey_type(line))
                elif line.startswith('SURVEY_CREW'):
                    metalist.append(eH_p._parse_survey_crew(line))
                elif line.startswith('SEA_CONDITION'):
                    metalist.append(eH_p._parse_sea_condition(line))
                elif line.startswith('VESSEL_NAME'):
                    metalist.append(eH_p._parse_vessel_name(line))
                elif line.startswith('LWRP'):
                    metalist.append(eH_p._parse_LWRP_(line))
                elif line.startswith('Gage_Reading'):
                    metalist.append(eH_p._parse_Gage_Reading(line, allcap1 = 1))
                elif line.startswith('GAGE_READING'):
                    metalist.append(eH_p._parse_Gage_Reading(line,allcap1 = 2))
                elif line.startswith('SOUND VELOCITY'):
                    metalist.append(eH_p._parse_sound_velocity(line))    
                elif eH_p._is_RTK(line, version='CEMVN'):
                    more_metalist.append(line)
                    metadata['RTK']='YES'
                    metalist.append(metadata['RTK'])
                    if eH_p._is_RTK_Tide(line, version = 'CEMVN'):
                        metadata['RTK TIDES']='YES'
                        metalist.append(metadata['RTK TIDES'])
                else:
                    more_metalist.append(line)#usually gage offets from MLG and NGVD/NAVD**, tides or river water levels
                    #plus gage name
                    #may include DRAFT, VELOCITY (sound velocity), Index within lines as well
                       
                #elif line.startswith('Ranges'):
                #    metalist.append(_parse_Ranges(line))
            # bring all the dictionaries together
            meta = {}
        #completed over 400 files successfully:
        #got got on     meta = {**meta, **m}
        #TypeError: 'NoneType' object is not a mapping
        try:
            for m in metalist:
                meta = {**meta, **m}
            return meta
        except:
            meta = {}
            errorfile = self.errorfile 
            with open(errorfile,'a') as metafail:
                metafail.write(infilename + '\n')
            return meta
    
    
    def load_default_metadata(self, infilename, default_meta):
        """
        Given the file name for data and a default metadata file (containing a
        picked dictionary), look for the default file.  If that files does not
        exist, look for a file named 'default.pkl' in the same directory as the
        provided file name.
        """
        if len(default_meta) == 0:
            path, infile = _os.path.split(infilename)
            default_meta = _os.path.join(path,'default.pkl')
        if _os.path.exists(default_meta):
            with open(default_meta, 'rb') as metafile:
                meta = _pickle.load(metafile)
        else:
            meta = {}
        return meta
##-----------------------------------------------------------------------------        
class Extract_Xml(object):
    """
    #NEED TO MERGE METHOD WITH parse_usace_xml_1.py
    """
    def __init__(self, preloadeddata, version = '', filename = ''):

        self.filename = filename
        #self.filepath = (filenamepath)

    def parse_ehydro_xml(infilename):
        """
        Parse the eHydro XML file as provided by Wilmington, Charleston, and
        Norfolk Districts.
        """
        xml_meta = eH_p._parse_xml(infilename)
        text_meta = eH_p._parse_xml_text(infilename)
        meta_out = {**xml_meta, **text_meta}
        return meta_out
    
    def _parse_xml(self, infilename):
        """
        Parse the xml portion of the xml file
        """
        xml_meta = {}
        tree = _parse(infilename)
        root = tree.getroot()
        val = root.findtext('.//altdatum')
        if val is not None:
            xml_meta['from_vert_key'] = val
        val = root.findtext('.//altunits')
        if val is not None:
            xml_meta['from_vert_units'] = val
        val = root.findtext('.//abstract')
        if val is not None:
            vals = val.split(' ')
            for n,v in enumerate(vals):
                if v == 'dates':
                    break
            date_str = vals[n+1]
            start, end = date_str.split(',')
            xml_meta['start_date'] = start.replace('-','')
            xml_meta['end_date'] = end.replace('-','')
        return xml_meta
            
    def _parse_xml_text(self, infilename):
        """
        Pase the text portion of the xml file
        """
        txt_meta = {}
        txt_keys = {'Implied_Vertical_Accuracy' : 'from_vert_unc',
                    'Implied_Horizontal_Accuracy' : 'from_horiz_unc',
                    'Horizontal_Zone' : 'from_horiz_datum',
                    'Units' : 'from_horiz_units'}
        keys = txt_keys.keys()
        with open(infilename,'r') as metafile:
            for line in metafile:
                for key in keys:
                    if line.startswith(key):
                        meta_key = txt_keys[key]
                        txt_meta[meta_key] = line
        for key in txt_meta:
            if key == 'from_vert_acc' or key == 'from_horiz_acc':
                line = txt_meta[key]
                val = line.split()[-2]
                txt_meta[key] = val
            elif key == 'from_horiz_datum':
                line = txt_meta[key]
                val = line.split(':')[-1]
                val = val.lstrip(' ')
                val = val.rstrip('\n')
                txt_meta[key] = val
                fips = _re.search('\d\d\d\d', val)
            elif key == 'from_horiz_units':
                line = txt_meta[key]
                val = line.split(':')[-1]
                val = val.lstrip(' ')
                val = val.rstrip('\n')
                txt_meta[key] = val
        if fips is not None:
            txt_meta['fips'] = int(fips.group())
        return txt_meta
##-----------------------------------------------------------------------------        
