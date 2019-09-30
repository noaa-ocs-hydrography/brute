import os
import unittest

from fuse.fuse_processor import FuseProcessor

TESTING_DIRECTORY = r'\\OCS-VS-NBS02\data\testing'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
NOAA_INPUT_ROOT = os.path.join(INPUT_ROOT, 'NOAA')
USACE_INPUT_ROOT = os.path.join(INPUT_ROOT, 'USACE')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')

USACE_CONFIG_ROOT = os.path.join('data', 'USACE')
NOAA_CONFIG_ROOT = os.path.join('data', 'NOAA')


def process_USACE_points(survey_name: str, interpolation_method: str, output_type: str) -> str:
    config_path = os.path.join(USACE_CONFIG_ROOT, f'{survey_name}_{interpolation_method}_{output_type}.config')
    input_path = os.path.join(USACE_INPUT_ROOT, survey_name, f'{survey_name}.XYZ')
    output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{output_type}')

    if os.path.exists(output_path):
        os.remove(output_path)

    fuse_processor = FuseProcessor(config_path)
    fuse_processor.read(input_path)
    fuse_processor.process(input_path)

    return output_path


def process_NOAA_raster(survey_name: str, interpolation_method: str, output_type: str) -> str:
    config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_{interpolation_method}_{output_type}.config')
    input_directory = os.path.join(NOAA_INPUT_ROOT, survey_name)
    bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if name[-4:] == '.bag']
    output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{output_type}')

    if os.path.exists(output_path):
        os.remove(output_path)

    fuse_processor = FuseProcessor(config_path)
    for bag_filename in bag_paths:
        if 'INTERP' not in bag_filename:
            fuse_processor.read(bag_filename)
            fuse_processor.process(bag_filename)

    return output_path


class TestPointLinear(unittest.TestCase):
    def test_small(self):
        output_path = process_USACE_points('NY_05_RHF_20181227_CS_4787_45X', 'linear', 'bag')
        assert os.path.exists(output_path)


class TestPointKriging(unittest.TestCase):
    def test_small(self):
        output_path = process_USACE_points('NY_05_RHF_20181227_CS_4787_45X', 'kriging', 'bag')
        assert os.path.exists(output_path)

    def test_bounds_error(self):
        output_path = process_USACE_points('BR_01_BRH_20130821_BD_4041_30X', 'kriging', 'bag')
        assert os.path.exists(output_path)

    def test_issue_survey(self):
        output_path = process_USACE_points('BR_01_BRH_20190117_CS_4788_40X', 'kriging', 'bag')
        assert os.path.exists(output_path)


class TestRasterLinear(unittest.TestCase):
    def test_small(self):
        output_path = process_NOAA_raster('H12607', 'linear', 'bag')
        assert os.path.exists(output_path)

    def test_large(self):
        output_path = process_NOAA_raster('H12525', 'linear', 'bag')
        assert os.path.exists(output_path)

    def test_issue_survey(self):
        output_path = process_NOAA_raster('F00521', 'linear', 'bag')
        assert os.path.exists(output_path)


class TestRasterKriging(unittest.TestCase):
    def test_small(self):
        output_path = process_NOAA_raster('H12607', 'kriging', 'bag')
        assert os.path.exists(output_path)


if __name__ == '__main__':
    unittest.main()
