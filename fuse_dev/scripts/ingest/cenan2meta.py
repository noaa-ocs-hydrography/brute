"""
cenan2meta.py

grice
20190408
V.0.0.1 20190408

This is a script to demonstrate importing data from the USACE district CENAN
into the metadata file for qualification.
"""

import os
from glob import glob
import fuse.fuse_ehydro as ffe

cenan = ffe.fuse_ehydro('cenan.config') # this config is local for testing
for path in cenan.rawdata_path:
    flist = glob(os.path.join(path,'*.xyz'))
    for f in flist:
        cenan.read(f)
        cenan.process(f)