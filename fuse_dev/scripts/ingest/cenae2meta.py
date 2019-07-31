# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:32:03 2019

@author: Casiano.Koprowski
"""

import os
from glob import glob

import fuse.fuse_ehydro as ffe

if __name__ == '__main__':
    cenae = ffe.fuse_ehydro('cenae.config')  # this config is local for testing
    for path in cenae.rawdata_path:
        flist = glob(os.path.join(path, '*.xyz'))
        for f in flist:
            cenae.read(f)
            try:
                cenae.process(f)
            except ValueError as e:
                print(e)