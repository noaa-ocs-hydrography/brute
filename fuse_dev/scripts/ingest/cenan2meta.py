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
import logging as _logging

import datetime

import fuse.fuse_processor as ffp
import fuse.wx_helper.process_title as wx_window

if __name__ == '__main__':
    start = datetime.datetime.now()
    print(start)
    wx_frame = wx_window.Open_Frame('CENAN')
    cenan = ffp.FuseProcessor('cenan.config')  # this config is local for testing
    c = 1
    root = cenan.rawdata_path[0]
    top = [os.path.join(root, name) for name in os.listdir(root)]
    for path in top:
        print(f'{c} - Begin working in {path}:')
        f = cenan.read(path)
        try:
            print(f'processing {f}', end = ', ')
            cenan.process(f)
            print(f'done.')
            c += 1
        except Exception as e:
            print('\n')
            print(e)
            cenan.logger.log(_logging.DEBUG, e)
            print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
    wx_frame.close()
    print(f'{end}\n{time_delta}')
