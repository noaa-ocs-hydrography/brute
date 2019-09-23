"""
cespl2meta.py
ed ckoprowski
20190626
based on cespl2meta.py by
grice
20190408
V.0.0.1 20190408

This is a script to demonstrate importing data from the USACE district CESPL
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
    wx_frame = wx_window.Open_Frame('CESPL')
    cespl = ffp.FuseProcessor('cespl.config')  # this config is local for testing
    c = 1
    root = cespl.rawdata_path[0]
    top = [os.path.join(root, name) for name in os.listdir(root)]
    for path in top:
        print(f'{c} - Begin working in {path}:')
        f = cespl.read(path)
        try:
            print(f'processing {f}', end = ', ')
            cespl.process(f)
            print(f'done.')
            c += 1
        except Exception as e:
            print('\n')
            print(e)
            cespl.logger.log(_logging.DEBUG, e)
            print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
    wx_frame.close()
    print(f'{end}\n{time_delta}') 
