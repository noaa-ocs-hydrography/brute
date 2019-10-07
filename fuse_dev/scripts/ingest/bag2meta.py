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


if __name__ == '__main__':
    start = datetime.datetime.now()
    print(start)
#    wx_frame = wx_window.Open_Frame('BAG')
    config_list = glob('bag*.config')
    for n,config in enumerate(config_list):
        bag = ffp.FuseProcessor(config)
        root = bag.rawdata_path[0]
        top = [os.path.join(root, name) for name in os.listdir(root)]
        for m, path in enumerate(top):
            print(f'{n}.{m} - Begin working in {path}:')
            flist = glob(os.path.join(path, '*.bag'))
            for m, file in enumerate(flist):
                try:
                    print(f'reading {os.path.basename(file)}', end = ', ')
                    f = bag.read(file)
                except Exception as e:
                    print('\n')
                    print(f'read error: {e}')
                    bag.logger.log(_logging.DEBUG, e)
                    print('\n')
                try:
                    print(f'processing @ {datetime.datetime.now()}')
                    bag.process(file)
                    print(f'done.')
                except Exception as e:
                    print('\n')
                    print(f'processing error: {e}')
                    bag.logger.log(_logging.DEBUG, e)
                    print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
#    wx_frame.close()
    print(f'{end}\n{time_delta}')