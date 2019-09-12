import os
import unittest

from fuse.fuse_processor import FuseProcessor

TESTING_DIRECTORY = r'\\OCS-VS-NBS02\data\testing'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')

USACE_CONFIG_ROOT = os.path.join('data', 'USACE')
NOAA_CONFIG_ROOT = os.path.join('data', 'NOAA')


class TestPointInterpolator(unittest.TestCase):
    def test_linear(self):
        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'

        file_type = 'bag'
        config_path = os.path.join(USACE_CONFIG_ROOT, f'{survey_name}_linear.config')

        input_path = os.path.join(INPUT_ROOT, survey_name, f'{survey_name}.XYZ')
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        fuse_processor.read(input_path)
        fuse_processor.process(input_path)

        assert os.path.exists(output_path)

    def test_kriging(self):
        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'

        file_type = 'bag'
        config_path = os.path.join(USACE_CONFIG_ROOT, f'{survey_name}_kriging.config')

        input_path = os.path.join(INPUT_ROOT, survey_name, f'{survey_name}.XYZ')
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        fuse_processor.read(input_path)
        fuse_processor.process(input_path)

        assert os.path.exists(output_path)


class TestRasterInterpolator(unittest.TestCase):
    def test_linear(self):
        survey_name = 'H12607'

        file_type = 'bag'
        config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}.config')

        input_directory = os.path.join(INPUT_ROOT, survey_name)
        bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        for bag_filename in bag_paths:
            if 'INTERP' not in bag_filename:
                fuse_processor.read(bag_filename)
                fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)

    def test_kriging(self):
        survey_name = 'H12607'

        file_type = 'bag'
        config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_kriging.config')

        input_directory = os.path.join(INPUT_ROOT, survey_name)
        bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        for bag_filename in bag_paths:
            if 'INTERP' not in bag_filename:
                fuse_processor.read(bag_filename)
                fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)


if __name__ == '__main__':
    unittest.main()
