import os
import pathlib
import unittest

from fuse.fuse_processor import FuseProcessor

DATA_PATH = r'\\OCS-VS-NBS01\nbs\NBS_Data'


class TestFuse(unittest.TestCase):
    def test_cespl(self):
        input_directory = os.path.join(DATA_PATH, r'PBD_Pacific\USACE\eHydro_LosAngeles_CESPL\Original')
        processed_directory = os.path.join(DATA_PATH, r'PBD_Pacific\USACE\eHydro_LosAngeles_CESPL\MLLW\Active')

        if not os.path.exists(input_directory):
            raise EnvironmentError(f'data directory not found: {input_directory}')
        if not os.path.exists(processed_directory):
            pathlib.Path(processed_directory).mkdir(parents=True, exist_ok=True)

        survey_name = 'LA_02_LAC_20150915'
        file_type = 'bag'

        input_path = os.path.join(input_directory, survey_name, f'{survey_name}.XYZ')
        config_path = os.path.join('data', 'cespl.config')
        output_path = os.path.join(processed_directory, f'{survey_name}_5m_interp.{file_type}')

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        assert os.path.exists(output_path)

    def test_bag(self):
        input_directory = r'\\OCS-VS-NBS01\nbs\TestingResources\NOAABAG_related\F00626 - This one has seperated areas of bathy'
        processed_directory = os.path.join(DATA_PATH, r'PBC_Northeast\NOAA_NCEI_OCS\BAGs\MLLW\Active')

        if not os.path.exists(input_directory):
            raise EnvironmentError(f'data directory not found: {input_directory}')
        if not os.path.exists(processed_directory):
            pathlib.Path(processed_directory).mkdir(parents=True, exist_ok=True)

        survey_name = 'F00626'
        file_type = 'bag'

        bag_filenames = [os.path.join(input_directory, name) for name in os.listdir(input_directory) if '.bag' in name]

        config_path = os.path.join('data', 'raster_interpolation.config')
        output_path = os.path.join(processed_directory, f'{survey_name}_5m_interp.{file_type}')

        cenan_fuse_processor = FuseProcessor(config_path)

        for bag_filename in bag_filenames:
            if 'INTERP' not in bag_filename:
                cenan_fuse_processor.read(bag_filename)
                cenan_fuse_processor.process(bag_filename)

        assert os.path.exists(output_path)


if __name__ == '__main__':
    unittest.main()
