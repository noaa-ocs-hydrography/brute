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


SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    start = datetime.datetime.now()
    print(start)
#    wx_frame = wx_window.Open_Frame('BPS')
    config_list = glob(os.path.join(SCRIPT_DIRECTORY, 'bps*.config'))
    for n,config in enumerate(config_list):
        bps = ffp.FuseProcessor(config)
        root = bps.rawdata_path[0]
        top = [os.path.join(root, name) for name in os.listdir(root)]
        for m,path in enumerate(top):
            print(f'{n}.{m} - Begin working in {path}:')
            paths = bps.read(path)
            for f in paths:
                try:
                    print(f'processing {f} @ {datetime.datetime.now()}', end = ', ')
                    bps.process(f)
                    print(f'done.')
                except Exception as e:
                    print('\n')
                    print(e)
                    bps.logger.log(_logging.DEBUG, e)
                    print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
    # wx_frame.close()
    print(f'{end}\n{time_delta}')