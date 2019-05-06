"""
cemvn2meta.py
ed jkinney
20190417
based on cenan2meta.py by
grice
20190408
V.0.0.1 20190408

This is a script to demonstrate importing data from the USACE district CEMVN
into the metadata file for qualification.
"""

import os
from glob import glob
import fuse.fuse_ehydro as ffe

cemvn = ffe.fuse_ehydro('cemvn.config') # this config is local for testing
for path in cemvn.rawdata_path:
    flist = glob(os.path.join(path,'*.xyz'))
    for f in flist:
        cemvn.read(f)
        cemvn.process(f)

