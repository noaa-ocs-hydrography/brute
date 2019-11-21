import os
import unittest

from tests.utilities import process_USACE_points, process_NOAA_raster


class TestPointLinear(unittest.TestCase):
    def test_small(self):
        output_paths = process_USACE_points('NY_05_RHF_20181227_CS_4787_45X', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_issue_survey(self):
        output_paths = process_USACE_points('BR_01_BRH_20190117_CS_4788_40X', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_holes(self):
        output_paths = process_USACE_points('NB_01_MAI_20170131_CS_4571_30X', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'


# class TestPointKriging(unittest.TestCase):
#     def test_small(self):
#         output_paths = process_USACE_points('NY_05_RHF_20181227_CS_4787_45X', 'kriging', 'bag')
#         for output_path in output_paths:
#             assert os.path.exists(output_path), f'output file not created {output_path}'
#
#     def test_bounds_error(self):
#         output_paths = process_USACE_points('BR_01_BRH_20130821_BD_4041_30X', 'kriging', 'bag')
#         for output_path in output_paths:
#             assert os.path.exists(output_path), f'output file not created {output_path}'
#
#     def test_issue_survey(self):
#         output_paths = process_USACE_points('BR_01_BRH_20190117_CS_4788_40X', 'kriging', 'bag')
#         for output_path in output_paths:
#             assert os.path.exists(output_path), f'output file not created {output_path}'


class TestRasterLinear(unittest.TestCase):
    def test_small(self):
        output_paths = process_NOAA_raster('H12607', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_large(self):
        output_paths = process_NOAA_raster('H12525', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_issue_1(self):
        output_paths = process_NOAA_raster('F00521', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_issue_2(self):
        output_paths = process_NOAA_raster('D00223', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_issue_3(self):
        output_paths = process_NOAA_raster('H11709', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_separated_areas(self):
        output_paths = process_NOAA_raster('F00626', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'

    def test_singlebeam(self):
        output_paths = process_NOAA_raster('F00623', 'linear', 'bag')
        for output_path in output_paths:
            assert os.path.exists(output_path), f'output file not created {output_path}'


# class TestRasterKriging(unittest.TestCase):
#     def test_small(self):
#         output_paths = process_NOAA_raster('H12607', 'kriging', 'bag')
#         for output_path in output_paths:
#             assert os.path.exists(output_path), f'output file not created {output_path}'


if __name__ == '__main__':
    unittest.main()
