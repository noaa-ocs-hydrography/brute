"""
cespl2meta.py
ed ckoprowski
20190626
based on cenan2meta.py by
grice
20190408
V.0.0.1 20190408

This is a script to demonstrate importing data from the USACE district CESPL
into the metadata file for qualification.
"""

import os
from glob import glob

import fuse.fuse_ehydro as ffe

cespl = ffe.fuse_ehydro('cespl.config')  # this config is local for testing
root = cespl.rawdata_path[0]
top = [os.path.join(root, name) for name in os.listdir(root)]
total = str(len(top))
x = 1
for path in top:
    print('\n\n', str(x), 'of', total, path, sep='')
    flist = glob(os.path.join(path, '*.xyz'))
    for f in flist:
        cespl.read(f)
        cespl.process(f)
    x += 1
