# -*- coding: utf-8 -*-
"""
score.py

grice 20190725
V 0.0.1 20190725

Utilities for calculating metrics for scoring the quality of data.
"""
from datetime import datetime as _datetime
import math as _math

def catzoc(metadata: dict) -> int:
    """
    Return an enumeration representing the catzoc assocaited with the provided
    metrics.
    
    The enumeration (catzoc : value) is as follows:
        A1 : 1
        A2 : 2
        B  : 3
        C  : 4
        D  : 5
        U  : 6
        
    The provided metadata is expected to contain the following metadata values:
        
    """
    s = supersession(metadata)
    if s >= 100:
        return 1
    elif s >= 90:
        return 2
    elif s >= 80:
        return 3
    elif s >= 70:
        return 4
    elif s >= 60:
        return 5
    else:
        return 6

def supersession(metadata: dict) -> float:
    """
    Return the superssion score as defined in Wyllie 2017 at US Hydro for the
    catzoc score.
    """
    if "feat_detect" in metadata and 'complete_coverage' in metadata and 'horiz_uncert_fixed' in metadata and 'vert_uncert_fixed' in metadata:
        feat_score = _get_feature_detection(metadata)
        cov_score = _get_coverage(metadata)
        horz_score = _get_horizontal_uncertainty(metadata)
        vert_score = _get_vertical_uncertainty(metadata)
        return min(feat_score, cov_score, horz_score, vert_score)
    else:
        survey_name = metadata['from_filename']
        raise ValueError(f'Metadata is not available to score {survey_name}')

def _get_feature_detection(metadata: dict) -> float:
    """
    Determine the feature detection capability from the ability to detect
    features, detect the least depth, and the size of the feature.
    """
    detected = metadata['feat_detect'].upper() == 'TRUE'
    least_depth = metadata['feat_least_depth'].upper() == 'TRUE'
    size_okay = False
    if 'feat_size' in metadata:
        size = float(metadata['feat_size'])
        if size <= 2: 
            size_okay = True
    if detected and least_depth and size_okay:
        return 100
    else:
        return 80
    
def _get_coverage(metadata: dict) -> float:
    """
    Determine the coverage score and return.
    """
    cov = metadata['complete_coverage'].upper() == 'TRUE'
    if cov:
        return 100
    else:
        return 80
    
def _get_horizontal_uncertainty(metadata: dict) -> float:
    """
    Determine the horizontal uncertainty score and return.
    """
    h_fix = float(metadata['horiz_uncert_fixed'])
    h_var = float(metadata['horiz_uncert_vari'])
    if h_fix <= 5 and h_var <= 0.05:
        s = 100
    elif h_fix <= 20:
        s = 90
    elif h_fix <= 50:
        s = 80
    elif h_fix <= 500: 
        s = 70
    else:
        s = 60
    return s

def _get_vertical_uncertainty(metadata: dict) -> float:
    """
    Determine the vertical uncertainty score and return.
    """
    v_fix = float(metadata['vert_uncert_fixed'])
    v_var = float(metadata['vert_uncert_vari'])
    if v_fix <= 0.5 and v_var <= 0.01:
        s = 100
    elif v_fix <= 1 and v_var <= 0.02:
        s = 90
    elif v_fix <= 2 and v_var <= 0.05: 
        s = 70
    else:
        s = 60
    return s

def decay(metadata: dict, date: _datetime, alpha: float = 0.022 ) -> float:
    """
    Return the decayed supersession_score.
    """
    sd = _datetime.strptime(metadata['end_date'],'%Y%m%d')
    ss = float(metadata['supersession_score'])
    dt = date - sd
    days = dt.days + dt.seconds / (24 * 60 * 60)
    years = days / 365
    ds = ss * _math.exp(-alpha * years)
    return ds