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
        A3 : 3
        B  : 4
        C  : 5
        D  : 6
        
    The provided metadata is expected to contain the following metadata values:
        
    """
    pass

def supersession(metadata: dict) -> float:
    """
    Return the superssion score as defined in Wyllie 2017 at US Hydro for the
    catzoc score.
    """
    pass

def decay(metadata: dict, date: datetime, alpha: float = 0.022 ) -> float:
    """
    Return the decayed supersession_score.
    """
    survey_date = metadata['end_date']
    ss = metadata['supersession_score']
    dt = date - survey_date
    ds = ss * math.exp(-alpha * dt)
    return ds