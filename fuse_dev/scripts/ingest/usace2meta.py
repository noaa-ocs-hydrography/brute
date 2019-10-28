"""
usace2meta.py

grice
20190408
V.0.0.3 20190927

This is a script to demonstrate importing data from the USACE district CENAN
into the metadata file for qualification.

All local ce*.config files are run serially.
"""

import os
from glob import glob
import logging as _logging

import datetime

import fuse.fuse_processor as ffp
# import fuse.wx_helper.process_title as wx_window

if __name__ == '__main__':
    start = datetime.datetime.now()
    print(start)
    #wx_frame = wx_window.Open_Frame('USACE')
    config_list = glob('ce*.config')
    for n,config in enumerate(config_list):
        usace = ffp.FuseProcessor(config)
        root = usace.rawdata_path[0]
        top = [os.path.join(root, name) for name in os.listdir(root)]
        for m,path in enumerate(top):
            print(f'{n}.{m} - Begin working in {path}:')
            paths = usace.read(path)
            for f in paths:
                try:
                    print(f'processing {f} @ {datetime.datetime.now()}', end = ', ')
                    usace.process(f)
                    print(f'done.')
                except Exception as e:
                    print('\n')
                    print(e)
                    usace.logger.log(_logging.DEBUG, e)
                    print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
    # wx_frame.close()
    print(f'{end}\n{time_delta}')
