# -*- coding: utf-8 -*-
"""
csar2array.py

Created on Fri Sep 20 10:29:21 2019

@author: grice

This is a first hack at extracting the information to create BAGs from combined
csars.  The primary assumption is the Contributor layer in the CSAR has the
interpolation flag set as 0 or 1 at the end of the contributor string.
"""

import os
import numpy as np
import caris
import caris.coverage

path = r'N:\NationalBathymetricSource\20190509'
fname = 'LALB_mllw_2m_interp.csar'
#open the file
filepath = os.path.join(path, fname)
csar = caris.open(file_name = filepath)
# get the contribtor layer and set the interpolated data to the no data value
if 'Contributor' in csar.band_info:
    mask = csar.read_np_array(band_name = 'Contributor', area = ((0,0),csar.dims))
    contrib_info = csar.band_info['Contributor']
    # ndv = getattr(contrib_info,'ndv')
    table = getattr(contrib_info,'string_table')
    ndv = len(table)
    for n,t in enumerate(table):
        vals = t.split('|')
        source = bool(int(vals[-1]))
        if not source:
            idx = np.nonzero(mask == n)
            mask[idx] = ndv
# close the csar
csar = None
# save the mask out
root, ext = os.path.splitext(filepath)
maskname = root + '_mask.npy'
np.save(maskname, mask)