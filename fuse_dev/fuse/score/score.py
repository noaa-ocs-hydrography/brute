# -*- coding: utf-8 -*-
"""
score.py

grice 20190725
V 0.0.1 20190725

Utilities for calculating metrics for scoring the quality of data.
"""
from datetime import datetime
import math

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
    feat_score = _get_feature_detection(metadata)
    cov_score = _get_coverage(metadata)
    horz_score = _get_horizontal_uncertainty(metadata)
    vert_score = _get_vertical_uncertainty
    return min(feat_score, cov_score, horz_score, vert_score)

def _get_feature_detection(metadata: dict) -> float:
    """
    Determine the feature detection capability from the ability to detect
    features, detect the least depth, and the size of the feature.
    """
    detected = metadata['feat_detect'].upper() == 'TRUE'
    least_depth = metadata['feat_least_depth'].upper() == 'TRUE'
    size = float(metadata['feat_size'])
    if size <= 2: 
        size_okay = True
    else:
        size_okay = False
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
    h_fix = float(['horiz_uncert_fixed'])
    h_var = float(['horiz_uncert_vari'])
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

def decay(metadata: dict, date: datetime, alpha: float = 0.022 ) -> float:
    """
    Return the decayed supersession_score.
    """
    survey_date = metadata['end_date']
    ss = metadata['supersession_score']
    dt = date - survey_date
    ds = ss * math.exp(-alpha * dt)
    return ds