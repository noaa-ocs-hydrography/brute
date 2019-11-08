# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:09:10 2019

@author: Casiano.Koprowski
"""

import os,sys
from glob import glob
import logging as _logging

import datetime

import fuse.fuse_processor as ffp
import fuse.wx_helper.process_title as wx_window

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    start = datetime.datetime.now()
    print(start)
#    wx_frame = wx_window.Open_Frame('BAG')
    config_list = glob(os.path.join(SCRIPT_DIRECTORY, 'bag*.config'))
    for n,config in enumerate(config_list):
        bag = ffp.FuseProcessor(config)
        root = bag.rawdata_path[0]
        top = [os.path.join(root, name) for name in os.listdir(root)]
        for m, path in enumerate(top):
            print(f'{n}.{m} - Begin working in {path}:')
            try:
                print(f'reading {os.path.basename(path)}', end = ', ')
                fnames = glob(os.path.join(path,'*.bag'))
                files = bag.read(path)
            except Exception as error:                
                print('\n')
                _, _, error_traceback = sys.exc_info()
                msg = f'read error:{error} ({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})'
                print(msg)
                bag.logger.log(_logging.DEBUG, msg)
                print('\n')
            if files is None:
                print('\nNo file ID returned by fuse read.\n')
            else:
                for f in files:
                    try:
                        print(f'processing {f} @ {datetime.datetime.now()}')
                        bag.process(f)
                        print(f'done.')
                    except Exception as error:                
                        print('\n')
                        _, _, error_traceback = sys.exc_info()
                        msg = f'process error:{error} ({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})'
                        print(msg)
                        bag.logger.log(_logging.DEBUG, msg)
                        print('\n')
    end = datetime.datetime.now()
    time_delta = end - start
#    wx_frame.close()
    print(f'{end}\n{time_delta}')