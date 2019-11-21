# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 11:09:10 2019

@author: Casiano.Koprowski
"""

import logging
import os
import sys
from datetime import datetime
from glob import glob

from fuse.fuse_processor import FuseProcessor

# import fuse.wx_helper.process_title as wx_window

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    logger = logging.Logger('usace2meta')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(logging.Formatter('[%(asctime)s] %(name)-30s %(levelname)-8s: %(message)s'))
    logger.addHandler(console)

    # wx_frame = wx_window.Open_Frame('BAG')
    start_time = datetime.now()
    logger.info(f'starting NOAA BAG processing at {start_time}')
    config_filenames = glob(os.path.join(SCRIPT_DIRECTORY, 'bag_*.config'))

    total_files = 0

    for config_index, config_filename in enumerate(config_filenames):
        bag_processor = FuseProcessor(config_filename)
        config_input_root = bag_processor.rawdata_path[0]

        survey_directories = [os.path.join(config_input_root, filename) for filename in os.listdir(config_input_root)]
        for survey_index, survey_directory in enumerate(survey_directories):
            logger.info(f'config {config_index + 1}/{len(config_filenames)}, ' +
                        f'survey {survey_index + 1}/{len(survey_directories)} - ' +
                        f'Begin working in {survey_directory}:')

            bag_filenames = None
            try:
                bag_filenames = bag_processor.read(survey_directory)
            except Exception as error:
                _, _, error_traceback = sys.exc_info()
                logger.error(f'read error: {error.__class__.__name__} {error} ' +
                             f'({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})')

            if bag_filenames is None or len(bag_filenames) == 0:
                logger.warning('\nNo file ID returned by fuse read.\n')
            else:
                for bag_index, bag_filename in enumerate(bag_filenames):
                    total_files += len(bag_filenames)

                    try:
                        logger.info(f'file {bag_index}/{len(bag_filenames)}: reading {os.path.basename(bag_filename)}')
                        logger.info(f'processing {bag_filename}')
                        bag_processor.process(bag_filename)
                        logger.info(f'completed BAG {bag_index + 1} of {len(bag_filenames)} ' +
                                    f'of survey {survey_index + 1} of {len(survey_directories)}; ' +
                                    f'{(datetime.now() - start_time) / ((bag_index + 1) / (total_files + (len(bag_filenames) * len(survey_directories) * (len(config_filenames) - (config_index + 1)))))} remaining')
                    except Exception as error:
                        _, _, error_traceback = sys.exc_info()
                        logger.error(f'processing error: {error.__class__.__name__} {error} ' +
                                     f'({os.path.split(error_traceback.tb_frame.f_code.co_filename)[1]}:{error_traceback.tb_lineno})')

    end_time = datetime.now()
    logger.info(f'completed NOAA BAG processing at {end_time} (took {end_time - start_time})')
    # wx_frame.close()
