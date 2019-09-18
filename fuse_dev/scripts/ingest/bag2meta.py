# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:09:10 2019

@author: Casiano.Koprowski
"""

import os
from glob import glob

import fuse.fuse_processor as ffp


if __name__ == '__main__':
    bags = ffp.FuseProcessor('bag.config')
    root = r'R:\Scripts\Testing Files\BAGs and SSS Mosaics for Interpolation'
    top = [os.path.join(root, name) for name in os.listdir(root)]
    total = len(top)
    x = 1
    for path in top:
        print(f'\n\n{x}of{total}\n{path}')
        flist = glob(os.path.join(path, '*.bag'))
        for f in flist:
            print(f)
            if 'INTERP' not in f:
                bags.read(f)
        x += 1