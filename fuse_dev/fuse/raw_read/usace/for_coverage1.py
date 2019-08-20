# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:09:19 2019

@author: jkinney
"""

def convert_meta_to_input_coverage(m):
    m['feat_size'] = m['f_size']
    m['script: feat_size'] = m['f_size']
    m['feat_detect'] = m['f_dict'] 
    m['script: feat_detect'] = m['f_dict'] 
    m['script: feat_least_depth'] =m['f_lstd']
    m['feat_least_depth'] =m['f_lstd']
    m['script: complete_coverage'] = m['flcvrg']
    m['complete_coverage'] = m['flcvrg']
    m['script: complete_bathymetry'] = m['flbath']
    m['complete_bathymetry'] = m['flbath']
    return m