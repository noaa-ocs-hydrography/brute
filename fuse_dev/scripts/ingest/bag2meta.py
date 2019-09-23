# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:09:10 2019

@author: Casiano.Koprowski
"""

import os
from glob import glob
import logging as _logging

import datetime

import fuse.fuse_processor as ffp
import fuse.wx_helper.process_title as wx_window


def main():
    start = datetime.datetime.now()
    print(start)
    wx_frame = wx_window.Open_Frame('BAG')
    bags = ffp.FuseProcessor('bag.config')
    c = 1
    root = bags.rawdata_path[0]
    top = [os.path.join(root, name) for name in os.listdir(root)]
#    stop = int(len(top)/2)
    for path in top:
        print(f'Begin working in {path}:')
        flist = glob(os.path.join(path, '*.bag'))
        for f in flist:
            p,fname = os.path.split(f)
            print(f'{c}:Reading {fname}', end = ', ')
            bags.read(f)
            c += 1
            try:
                print(f'processing', end = ', ')
                bags.process(f)
                print(f'done.')
            except Exception as e:
                print('\n')
                print(e)
                bags.logger.log(_logging.DEBUG, e)
                print('\n')

    end = datetime.datetime.now()
    time_delta = end - start
    wx_frame.close()
    print(f'{end}\n{time_delta}')  


if __name__ == '__main__':
    main()