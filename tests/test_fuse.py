import os
import unittest

from fuse.fuse_processor import FuseProcessor

TESTING_DIRECTORY = r'\\OCS-VS-NBS01\nbs\TestingResources'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')


class TestFuse(unittest.TestCase):
    def test_usace_cespl(self):
        config_path = os.path.join('data', 'usace_cespl.config')
        survey_name = 'LA_01_LBC_20151118'
        file_type = 'bag'
        input_path = os.path.join(INPUT_ROOT, survey_name, f'{survey_name}.XYZ')
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        assert os.path.exists(output_path)
        os.remove(output_path)

    def test_bag(self):
        config_path = os.path.join('data', 'bag.config')
        survey_name = 'H12607'
        file_type = 'bag'
        input_directory = os.path.join(INPUT_ROOT, 'H12607 - smol')
        bag_filenames = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        if os.path.exists(output_path):
            os.remove(output_path)

        cenan_fuse_processor = FuseProcessor(config_path)
        for bag_filename in bag_filenames:
            if 'INTERP' not in bag_filename:
                cenan_fuse_processor.read(bag_filename)
                cenan_fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)
        os.remove(output_path)


if __name__ == '__main__':
    unittest.main()
