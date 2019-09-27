# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:32:03 2019

@author: Casiano.Koprowski
"""

import os
from glob import glob
import logging as _logging

import datetime

import fuse.fuse_processor as ffp
import fuse.wx_helper.process_title as wx_window

if __name__ == '__main__':
    start = datetime.datetime.now()
    print(start)
    wx_frame = wx_window.Open_Frame('CENAE')
    cenae = ffp.FuseProcessor('cenae.config')  # this config is local for testing
    c = 1
    root = cenae.rawdata_path[0]
    top = [os.path.join(root, name) for name in os.listdir(root)]
    for path in top:
        print(f'{c} - Begin working in {path}:')
        f = cenae.read(path)
        try:
            print(f'processing {f}', end = ', ')
            cenae.process(f)
            print(f'done.')
            c += 1
        except Exception as e:
            print('\n')
            print(e)
            cenae.logger.log(_logging.DEBUG, e)
            print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
    wx_frame.close()
    print(f'{end}\n{time_delta}')