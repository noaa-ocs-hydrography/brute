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

if __name__ == '__main__':
    cenan = ffe.FuseProcessor_eHydro('cenan.config')  # this config is local for testing
    for path in cenan.rawdata_path:
        print(f'Begin working in {path}:')
        c = 1
        flist = glob(os.path.join(path, '*.xyz'))
        for f in flist:
            p,fname = os.path.split(f)
            print(f'{c}:Reading {fname}', end = ', ')
            cenan.read(f)
            try:
                print(f'processing', end = ', ')
                cenan.process(f)
                print(f'done.')
                c += 1
            except ValueError as e:
                print('\n')
                print(e)
                print('\n')
            if c > 3:
                break
