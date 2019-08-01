"""
ceswg2meta.py
ed jkinney
20190418
based on ceswg2meta.py by
grice
20190408
V.0.0.1 20190408

This is a script to demonstrate importing data from the USACE district CESWG
into the metadata file for qualification.
"""

import os
from glob import glob

import fuse.fuse_ehydro as ffe

if __name__ == '__main__':
    ceswg = ffe.FuseProcessor_eHydro('ceswg.config')  # this config is local for testing
    for path in ceswg.rawdata_path:
        flist = glob(os.path.join(path, '*.xyz'))
        for f in flist:
            ceswg.read(f)
            ceswg.process(f)
