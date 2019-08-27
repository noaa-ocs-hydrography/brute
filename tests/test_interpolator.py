import os
import pathlib
import unittest

from fuse.fuse_processor import FuseProcessor
from fuse.interpolator.bag_interpolator import coverage
from fuse.raw_read.noaa import bag

DATA_PATH = r"C:\Data\NBS"


class TestPointInterpolator(unittest.TestCase):
    def test_linear(self):
        input_directory = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'Original')
        processed_directory = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'MLLW', 'Data',
                                           'Active')

        if not os.path.exists(input_directory):
            raise EnvironmentError(f'data directory not found: {input_directory}')
        if not os.path.exists(processed_directory):
            pathlib.Path(processed_directory).mkdir(parents=True, exist_ok=True)

        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'
        file_type = 'bag'

        input_path = os.path.join(input_directory, survey_name, f'{survey_name}.XYZ')
        config_path = os.path.join('data', 'cenan_linear.config')

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        output_path = os.path.join(processed_directory, f'{survey_name}_5m_interp.{file_type}')
        assert os.path.exists(output_path)

    def test_kriging(self):
        input_directory = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'Original')
        processed_directory = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'MLLW', 'Data',
                                           'Active')

        if not os.path.exists(input_directory):
            raise EnvironmentError(f'data directory not found: {input_directory}')
        if not os.path.exists(processed_directory):
            pathlib.Path(processed_directory).mkdir(parents=True, exist_ok=True)

        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'
        file_type = 'bag'

        input_path = os.path.join(input_directory, survey_name, f'{survey_name}.XYZ')
        config_path = os.path.join('data', 'cenan_kriging.config')

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        output_path = os.path.join(processed_directory, f'{survey_name}_5m_interp.{file_type}')
        assert os.path.exists(output_path)


class TestRasterInterpolator(unittest.TestCase):
    def test_align2grid(self):
        bag_testing_directory = os.path.join(DATA_PATH, 'testing', 'bag_interpolator', 'H12607')

        if not os.path.exists(bag_testing_directory):
            raise EnvironmentError(f'test directory not found: {bag_testing_directory}')

        bag_path = os.path.join(bag_testing_directory, 'H12607_MB_4m_MLLW_2of2.bag')
        coverage_list = [os.path.join(bag_testing_directory, 'H12607_SSSAB_1m_600kHz_2of2.tif')]

        bag_dataset = bag.BagFile()
        bag_dataset.open_file(bag_path, 'hack')
        bag_dataset.generate_name(bag_testing_directory, False)
        coverage_dataset = coverage.UnifiedCoverage(coverage_list, bag_dataset.wkt, bag_dataset.name)

        assert coverage_dataset.bounds != bag_dataset.bounds
        assert coverage_dataset.shape != bag_dataset.shape

        coverage_dataset = coverage.align2grid(coverage_dataset, bag_dataset.bounds, bag_dataset.shape,
                                               bag_dataset.resolution, bag_dataset.nodata)

        assert coverage_dataset.bounds == bag_dataset.bounds
        assert coverage_dataset.shape == bag_dataset.shape

    def test_linear(self):
        input_directory = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'Original')
        processed_directory = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'MLLW', 'Data',
                                           'Active')

        if not os.path.exists(input_directory):
            raise EnvironmentError(f'data directory not found: {input_directory}')
        if not os.path.exists(processed_directory):
            pathlib.Path(processed_directory).mkdir(parents=True, exist_ok=True)

        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'
        file_type = 'csar'

        input_path = os.path.join(input_directory, survey_name, f'{survey_name}.XYZ')
        config_path = os.path.join('data', 'cenan_kriging.config')
        output_path = os.path.join(processed_directory, f'{survey_name}_5m_interp.{file_type}')

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        assert os.path.exists(output_path)


if __name__ == '__main__':
    unittest.main()
