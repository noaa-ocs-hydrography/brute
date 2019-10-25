import os
from glob import glob

from fuse.fuse_processor import FuseProcessor

TESTING_DIRECTORY = r'\\OCS-VS-NBS01\nbs\TestingResources'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
NOAA_INPUT_ROOT = os.path.join(INPUT_ROOT, 'NOAA')
USACE_INPUT_ROOT = os.path.join(INPUT_ROOT, 'USACE')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')

USACE_CONFIG_ROOT = os.path.join('data', 'USACE')
NOAA_CONFIG_ROOT = os.path.join('data', 'NOAA')


def process_USACE_points(survey_name: str, interpolation_method: str, output_type: str) -> str:
    for filename in glob(os.path.join(OUTPUT_ROOT, f'{survey_name}*')):
        os.remove(filename)

    config_path = os.path.join(USACE_CONFIG_ROOT, f'{survey_name}_{interpolation_method}_{output_type}.config')
    input_path = os.path.join(USACE_INPUT_ROOT, survey_name, f'{survey_name}.XYZ')

    fuse_processor = FuseProcessor(config_path)
    fuse_processor.read(input_path)
    return fuse_processor.process(input_path)


def process_NOAA_raster(survey_name: str, interpolation_method: str, output_type: str) -> [str]:
    for filename in glob(os.path.join(OUTPUT_ROOT, f'{survey_name}.*')):
        os.remove(filename)

    config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_{interpolation_method}_{output_type}.config')
    input_directory = os.path.join(NOAA_INPUT_ROOT, survey_name)
    bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if name[-4:] == '.bag']

    output_paths = []
    fuse_processor = FuseProcessor(config_path)
    for bag_filename in bag_paths:
        if 'INTERP' not in bag_filename:
            fuse_processor.read(bag_filename)
            output_paths.append(fuse_processor.process(bag_filename))

    return output_paths
