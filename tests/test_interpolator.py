import os
import unittest

from fuse.interpolator import interpolator
from fuse.interpolator.bag_interpolator import bag, coverage

DATA_PATH = r"C:\Data\NBS"


class TestBagInterpolator(unittest.TestCase):
    def test_align2grid(self):
        bag_path = os.path.join(DATA_PATH, 'H12607_MB_4m_MLLW_2of2.bag')
        coverage_list = [os.path.join(DATA_PATH, 'H12607_SSSAB_1m_600kHz_2of2.tif')]

        bag_dataset = bag.BagFile()
        bag_dataset.open_file(bag_path, 'hack')
        bag_dataset.generate_name(DATA_PATH, False)
        coverage_dataset = coverage.UnifiedCoverage(coverage_list, bag_dataset.wkt, bag_dataset.name)

        assert coverage_dataset.bounds != bag_dataset.bounds
        assert coverage_dataset.shape != bag_dataset.shape

        coverage_dataset = coverage.align2grid(coverage_dataset, bag_dataset.bounds, bag_dataset.shape,
                                               bag_dataset.resolution, bag_dataset.nodata)

        assert coverage_dataset.bounds == bag_dataset.bounds
        assert coverage_dataset.shape == bag_dataset.shape


class TestPointInterpolator(unittest.TestCase):
    def test_kriging(self):
        input_path = os.path.join(DATA_PATH, 'PBC_Northeast', 'USACE', 'eHydro_NewYork_CENAN', 'Original',
                                  'BR_01_BRH_20190117_CS_4788_40X', 'BR_01_BRH_20190117_CS_4788_40X.XYZ')
        kriging_interpolator = interpolator.Interpolator('point', 'kriging', 500)


if __name__ == '__main__':
    unittest.main()
