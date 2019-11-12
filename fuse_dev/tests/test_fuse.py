import os
import unittest

from tests.utilities import process_USACE_points, process_NOAA_raster


class TestFuse(unittest.TestCase):
    def test_usace_cespl(self):
        output_path = process_USACE_points('LA_01_LBC_20151118', 'linear', 'bag')
        assert os.path.exists(output_path)

    def test_output_csar(self):
        output_path = process_USACE_points('NY_05_RHF_20181227_CS_4787_45X', 'linear', 'csar')
        assert os.path.exists(output_path)

    def test_noaa_bag_stripes(self):
        output_paths = process_NOAA_raster('H12963', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path)

    def test_noaa_bag_holes(self):
        output_paths = process_NOAA_raster('H12600', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path)

    def test_noaa_bag_reprojection(self):
        survey_1_output_paths = process_NOAA_raster('H11250', 'linear', 'bag')
        for output_path in survey_1_output_paths:
            assert os.path.exists(output_path)

        survey_2_output_paths = process_NOAA_raster('H12298', 'linear', 'bag')
        for output_path in survey_2_output_paths:
            assert os.path.exists(output_path)


if __name__ == '__main__':
    unittest.main()
