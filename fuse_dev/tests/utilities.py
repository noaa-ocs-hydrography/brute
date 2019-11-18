import os
from glob import glob

from fuse.fuse_processor import FuseProcessor

TESTING_DIRECTORY = r'D:\TestingResources'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
NOAA_INPUT_ROOT = os.path.join(INPUT_ROOT, 'NOAA')
USACE_INPUT_ROOT = os.path.join(INPUT_ROOT, 'USACE')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')

USACE_CONFIG_ROOT = os.path.join('data', 'USACE')
NOAA_CONFIG_ROOT = os.path.join('data', 'NOAA')


def process_USACE_points(survey_name: str, interpolation_method: str, output_type: str) -> [str]:
    for filename in glob(os.path.join(OUTPUT_ROOT, f'{survey_name}*')):
        os.remove(filename)

    config_path = os.path.join(USACE_CONFIG_ROOT, f'{survey_name}_{interpolation_method}_{output_type}.config')

    fuse_processor = FuseProcessor(config_path)
    xyz_filenames = fuse_processor.read(fuse_processor.rawdata_path[0])
    return [fuse_processor.process(xyz_filename) for xyz_filename in xyz_filenames]


def process_NOAA_raster(survey_name: str, interpolation_method: str, output_type: str) -> [str]:
    for filename in glob(os.path.join(OUTPUT_ROOT, f'{survey_name}.*')):
        os.remove(filename)

    config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_{interpolation_method}_{output_type}.config')
    input_directory = os.path.join(NOAA_INPUT_ROOT, survey_name)

    output_paths = []
    fuse_processor = FuseProcessor(config_path)
    bag_filenames = fuse_processor.read(fuse_processor.rawdata_path[0])
    return [fuse_processor.process(bag_filename) for bag_filename in bag_filenames]
