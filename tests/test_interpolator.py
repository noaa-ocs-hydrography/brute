import os
import unittest

from fuse.fuse_processor import FuseProcessor
from fuse.interpolator.bag_interpolator import coverage
from fuse.raw_read.noaa import bag

TESTING_DIRECTORY = r'\\OCS-VS-NBS01\nbs\TestingResources'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')


class TestPointInterpolator(unittest.TestCase):
    def test_linear(self):
        config_path = os.path.join('data', 'cenan_linear.config')
        file_type = 'bag'
        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'
        input_path = os.path.join(INPUT_ROOT, survey_name, f'{survey_name}.XYZ')
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        assert os.path.exists(output_path)

    def test_kriging(self):
        config_path = os.path.join('data', 'cenan_kriging.config')
        file_type = 'bag'
        survey_name = 'NY_05_RHF_20181227_CS_4787_45X'
        input_path = os.path.join(INPUT_ROOT, survey_name, f'{survey_name}.XYZ')
        output_path = os.path.join(OUTPUT_ROOT, f'{survey_name}_5m_interp.{file_type}')

        cenan_fuse_processor = FuseProcessor(config_path)
        cenan_fuse_processor.read(input_path)
        cenan_fuse_processor.process(input_path)

        assert os.path.exists(output_path)


class TestRasterInterpolator(unittest.TestCase):
    def test_align2grid(self):
        input_directory = os.path.join(INPUT_ROOT, r'H12607 - smol')
        bag_path = os.path.join(input_directory, 'H12607_MB_4m_MLLW_2of2.bag')
        coverage_list = [os.path.join(input_directory, f'H12607_SSS_{id}_1m.tif') for id in (100, 200)]

        bag_dataset = bag.BagFile()
        bag_dataset.open_file(bag_path, 'hack')
        bag_dataset.generate_name(input_directory, False)
        coverage_dataset = coverage.UnifiedCoverage(coverage_list, bag_dataset.wkt, bag_dataset.name)

        assert coverage_dataset.bounds != bag_dataset.bounds
        assert coverage_dataset.shape != bag_dataset.shape

        coverage_dataset = coverage.align2grid(coverage_dataset, bag_dataset.bounds, bag_dataset.shape,
                                               bag_dataset.resolution, bag_dataset.nodata)

        assert coverage_dataset.bounds == bag_dataset.bounds
        assert coverage_dataset.shape == bag_dataset.shape


if __name__ == '__main__':
    unittest.main()
