# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:54:19 2019

@author: grice
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
    wx_frame = wx_window.Open_Frame('BPS')
    bps = ffp.FuseProcessor('bps.config')
    path = bps.rawdata_path[0]
    print(f'Begin working in {path}:')
    flist = glob(os.path.join(path, '*.xyz'))
    for c,f in enumerate(flist):
        p,fname = os.path.split(f)
        print(f'{c}:Reading {fname}', end = ', ')
        bps.read(f)
        try:
            print(f'processing', end = ', ')
            bps.process(f)
            print(f'done.')
        except Exception as e:
            print('\n')
            print(e)
            bps.logger.log(_logging.DEBUG, e)
            print('\n')
        if c > 50:
            break

    end = datetime.datetime.now()
    time_delta = end - start
    wx_frame.close()
    print(f'{end}\n{time_delta}')  


if __name__ == '__main__':
    main()