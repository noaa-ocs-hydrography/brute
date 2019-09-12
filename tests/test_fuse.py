import os
import unittest

import rasterio
from fuse.fuse_processor import FuseProcessor

TESTING_DIRECTORY = r'\\OCS-VS-NBS02\data\testing'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')

USACE_CONFIG_ROOT = os.path.join('data', 'USACE')
NOAA_CONFIG_ROOT = os.path.join('data', 'NOAA')


class TestFuse(unittest.TestCase):
    def test_usace_cespl(self):
        survey_name = 'LA_01_LBC_20151118'

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

    def test_noaa_bag(self):
        survey_name = 'H12525'

        file_type = 'bag'
        config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_linear.config')

        input_directory = os.path.join(INPUT_ROOT, survey_name)
        bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_4m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        for bag_filename in bag_paths:
            if 'INTERP' not in bag_filename:
                fuse_processor.read(bag_filename)
                fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)

    def test_noaa_bag_small(self):
        survey_name = 'H12607'

        file_type = 'bag'
        config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_linear.config')

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

    def test_noaa_bag_stripes(self):
        survey_name = 'H12963'

        file_type = 'bag'
        config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_linear.config')

        input_directory = os.path.join(INPUT_ROOT, survey_name)
        bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_4m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        for bag_filename in bag_paths:
            if 'INTERP' not in bag_filename:
                fuse_processor.read(bag_filename)
                fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)

    def test_noaa_bag_holes(self):
        survey_name = 'H12600'

        file_type = 'bag'
        config_path = os.path.join(NOAA_CONFIG_ROOT, f'{survey_name}_linear.config')

        input_directory = os.path.join(INPUT_ROOT, survey_name)
        bag_paths = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_4m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        fuse_processor = FuseProcessor(config_path)
        for bag_filename in bag_paths:
            if 'INTERP' not in bag_filename:
                fuse_processor.read(bag_filename)
                fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)

    def test_noaa_bag_identical(self):
        survey_1_name = 'H12604'
        survey_2_name = 'H12981'

        file_type = 'bag'
        survey_1_config_path = os.path.join('data', f'{survey_1_name}_linear.config')
        survey_2_config_path = os.path.join('data', f'{survey_2_name}_linear.config')

        survey_1_input_directory = os.path.join(INPUT_ROOT, survey_1_name)
        survey_2_input_directory = os.path.join(INPUT_ROOT, survey_2_name)
        survey_1_bag_paths = [os.path.join(survey_1_input_directory, name) for name in
                              os.listdir(survey_1_input_directory) if '.bag' in name]
        survey_2_bag_paths = [os.path.join(survey_2_input_directory, name) for name in
                              os.listdir(survey_2_input_directory) if '.bag' in name]
        survey_1_output_path = os.path.join(OUTPUT_ROOT, f'{survey_1_name}_4m_interp.{file_type}')
        survey_2_output_path = os.path.join(OUTPUT_ROOT, f'{survey_2_name}_4m_interp.{file_type}')

        if os.path.exists(survey_1_output_path):
            os.remove(survey_1_output_path)

        if os.path.exists(survey_2_output_path):
            os.remove(survey_2_output_path)

        survey_1_fuse_processor = FuseProcessor(survey_1_config_path)
        for bag_filename in survey_1_bag_paths:
            if 'INTERP' not in bag_filename:
                survey_1_fuse_processor.read(bag_filename)
                survey_1_fuse_processor.process(bag_filename)

        survey_2_fuse_processor = FuseProcessor(survey_2_config_path)
        for bag_filename in survey_2_bag_paths:
            if 'INTERP' not in bag_filename:
                survey_2_fuse_processor.read(bag_filename)
                survey_2_fuse_processor.process(bag_filename)

        assert os.path.exists(survey_1_output_path)
        assert os.path.exists(survey_2_output_path)

        with rasterio.open(survey_1_output_path) as survey_1_interpolated, rasterio.open(survey_2_output_path) as survey_2_interpolated:
            assert survey_1_interpolated.read() == survey_2_interpolated.read()


if __name__ == '__main__':
    unittest.main()
