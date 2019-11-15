# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:09:10 2019

@author: Casiano.Koprowski
"""

import logging as _logging
import os
import sys
from datetime import datetime
from glob import glob

from fuse.fuse_processor import FuseProcessor

# import fuse.wx_helper.process_title as wx_window

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    # wx_frame = wx_window.Open_Frame('BAG')
    start_time = datetime.now()
    print(f'starting NOAA BAG processing at {start_time}')
    config_filenames = glob(os.path.join(SCRIPT_DIRECTORY, 'bag_configs', 'pbc_*.config'))

    total_files = 0

    for config_index, config_filename in enumerate(config_filenames):
        bag_processor = FuseProcessor(config_filename)
        config_input_root = bag_processor.rawdata_path[0]

        survey_directories = [os.path.join(config_input_root, filename) for filename in os.listdir(config_input_root)]
        for survey_index, survey_directory in enumerate(survey_directories):
            print(f'config {config_index + 1}/{len(config_filenames)}, ' +
                  f'survey {survey_index + 1}/{len(survey_directories)} - ' +
                  f'Begin working in {survey_directory}:')

            bag_filenames = None
            try:
                bag_filenames = bag_processor.read(survey_directory)
            except Exception as error:
                _, _, error_traceback = sys.exc_info()
                message = f'read error: {error.__class__.__name__} {error} ({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})'
                print(message)
                bag_processor.logger.log(_logging.DEBUG, message)

            if bag_filenames is None or len(bag_filenames) == 0:
                message = '\nNo file ID returned by fuse read.\n'
                print(message)
                bag_processor.logger.log(_logging.DEBUG, message)
            else:
                for bag_index, bag_filename in enumerate(bag_filenames):
                    total_files += len(bag_filenames)

                    try:
                        print(f'file {bag_index}/{len(bag_filenames)}: reading {os.path.basename(bag_filename)}', end=', ')
                        print(f'{datetime.now()}: processing {bag_filename}', end=', ')
                        bag_processor.process(bag_filename)
                        print(f'completed BAG {bag_index + 1} of {len(bag_filenames)} ' +
                              f'of survey {survey_index + 1} of {len(survey_directories)}; ' +
                              f'{(datetime.now() - start_time) / ((bag_index + 1) / (total_files + (len(bag_filenames) * len(survey_directories) * (len(config_filenames) - (config_index + 1)))))} remaining')
                    except Exception as error:
                        _, _, error_traceback = sys.exc_info()
                        message = f'processing error: {error.__class__.__name__} {error} ({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})'
                        print(message)
                        bag_processor.logger.log(_logging.DEBUG, message)

    end_time = datetime.now()
    print(f'completed NOAA BAG processing at {end_time} (took {end_time - start_time})')
    # wx_frame.close()
