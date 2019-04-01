# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:11:04 2019

@author: jkinney
"""
import os

def use_extract_meta_CEMVN(one_file):    
    source = 'USACE'
    source_type = 'ehydro'
    district ='CEMVN'
    versionmeta = district
    one_file_info = Version_Info(one_file)#, filename = '')
    one_file_info.pass_version_Meta(versionmeta, source , source_type, district)
    v = one_file_info.calculate_version_f()
    return one_file_info, v

def use_extract_meta_to_pass(one_file, source, source_type , district, versionmeta = None):    
    if versionmeta == None:
        versionmeta = district
    one_file_info = Version_Info(one_file, filename = '')
    one_file_info = one_file_info.pass_version_Meta(versionmeta, source , source_type, district)
    v = one_file_info.calculate_version_f()
    return one_file_info, v

def pass_metafrom_object(one_file):
    if len(source.USACE_District)>0:
        one_file_info = use_extract_meta_to_pass(one_file, source , source_type, district, versionmeta)
    else:
        one_file_info = use_extract_meta_to_pass(one_file, source , source_type, versionmeta)
    #one_file.versionmeta
    v = one_file_info.calculate_version_f()
    return one_file_info, v
    
    
class Version_Info(object):
    """
    """
    def __init__(self, preloadeddata,  filename = ''):
        if filename == '':
            self.filename = preloadeddata
        elif filename == None:
            self.filename = preloadeddata
        else:
            self.filename = filename                                
        #self.filepath = (filenamepath)
        #if self.filepath == False:
        self.filepath = os.path.dirname(self.filename)
        
    def calculate_version_f(self):
        if self.source:
            v = 'information on version'
        else:
            v = 'unknown'
        #same as in Bath_Data class?
        return v
                                 
    def pass_version_Meta(self, versionmeta = None, source = None , source_type = None, USACE_District = ''):
        """
        based on file paths and naming convention calculate version
        versionmeta option to pass this version information from another source 
        as well.
        """
        if source == None:
            source = ''
        if source_type == None:
            source_type = ''
        if USACE_District == None:
            USACE_District = ''            
        if versionmeta is not None:
            self.source = source 
            self.source_type = source_type
            self.versionmeta = versionmeta
            self.USACE_District_Dict = USACE_District
        else:
            self.source = 'USACE' 
            self.source_type = 'ehydro'
            self.versionmeta = 'CEMVN'
            self.USACE_District_Dict = 'CEMVN'
           
    def get_xml(self):
        if self.source_type =='ehydro':
            filename = self.filename
            basef = filename.rstrip('.xyz')
            basef = filename.rstrip('.XYZ')
            xml_name = basef + '.xml'
            return xml_name            

    def get_xml_xt(self, extension):
        if self.source_type =='ehydro':
            filename = self.filename
            basef = filename.rstrip(extension)
            xml_name = basef + '.xml'
            return xml_name

