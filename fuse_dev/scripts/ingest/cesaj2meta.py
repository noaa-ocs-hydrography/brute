"""
cesaj2meta.py
ed jkinney
20190417
based on cenan2meta.py by
grice
20190408
V.0.0.1 20190408

This is a script to demonstrate importing data from the USACE district CESAJ
into the metadata file for qualification.
"""

import os
from glob import glob

import fuse.fuse_processor as ffp

if __name__ == '__main__':
    cesaj =  ffp.FuseProcessor('cesaj.config')  # this config is local for testing
    for path in cesaj.rawdata_path:
        flist = glob(os.path.join(path, '*.xyz'))
        for f in flist:
            cesaj.read(f)
            cesaj.process(f)
