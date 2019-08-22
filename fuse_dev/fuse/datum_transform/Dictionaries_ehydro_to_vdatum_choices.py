# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:23:13 2019

@author: jkinney
"""




e_hydro_to_vdatum_vertical_codes = {
        'lwrp' :'mississippi_lwrp',
        'ngvd29': 'ngvd29',
        'mllw' : 'MLLW',
        'mlw' : 'MLW',
        'mhw' : 'MHW',
        'mhhw' : 'MHHW',
        'lwd' : 'LWD',
        'navd88' : 'navd88',
        }


#NEED check for unusual datums and passing to correct vdatum versions: Probably need geospatial query /polygon/grid to 

#NEED check for unusual local datums and determining if we need 0 vertical changes, as datums is already in chart datum
#NEED to add in manual checks as well

#NEED to track vdatum codes versus S-57 and codes needed for MQUAL/ source diagram outputs

#mll
#
#hrd
#crd

