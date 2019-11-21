# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:54:19 2019

@author: grice
"""

import os
import sys
from datetime import datetime
from glob import glob

from fuse.fuse_processor import FuseProcessor

# import fuse.wx_helper.process_title as wx_window

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    # wx_frame = wx_window.Open_Frame('USACE')
    start_time = datetime.now()
    print(f'starting BPS processing at {start_time}')
    config_filenames = glob(os.path.join(SCRIPT_DIRECTORY, 'bps_*.config'))

    total_files = 0

    for config_index, config_filename in enumerate(config_filenames):
        bps_processor = FuseProcessor(config_filename)
        config_input_root = bps_processor.rawdata_path[0]

        survey_directories = [os.path.join(config_input_root, filename) for filename in os.listdir(config_input_root)]
        total_files += len(survey_directories)
        for survey_index, survey_directory in enumerate(survey_directories):
            print(f'config {config_index + 1} of {len(config_filenames)}, ' +
                  f'survey {survey_index + 1} of {len(survey_directories)} - ' +
                  f'Begin working in {survey_directory}:')

            xyz_filenames = None
            try:
                xyz_filenames = bps_processor.read(survey_directory)
            except Exception as error:
                _, _, error_traceback = sys.exc_info()
                bps_processor.logger.error(
                    f'read error: {error.__class__.__name__} {error} ({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})')

            if xyz_filenames is None or len(xyz_filenames) == 0:
                bps_processor.logger.warning('\nNo file ID returned by fuse read.\n')
            else:
                for xyz_filename in xyz_filenames:
                    try:
                        print(f'{datetime.now()}: processing {xyz_filename}', end=', ')
                        bps_processor.process(xyz_filename)
                    except Exception as error:
                        _, _, error_traceback = sys.exc_info()
                        bps_processor.logger.error(
                            f'processing error: {error.__class__.__name__} {error} ({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})')
                print(f'completed survey {survey_index + 1} of {len(survey_directories)}; ' +
                      f'{(datetime.now() - start_time) / ((survey_index + 1) / (total_files + (len(survey_directories) * (len(config_filenames) - (config_index + 1)))))} remaining')

    end_time = datetime.now()
    print(f'completed BPS processing at {end_time} (took {end_time - start_time})')
    # wx_frame.close()
