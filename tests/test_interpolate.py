import os
import unittest

from fuse.interpolator.bag_interpolator import bag
from fuse.interpolator.bag_interpolator import coverage

DATA_PATH = r"C:\Data\NBS"


class TestBagInterpolator(unittest.TestCase):
    def test_align2grid(self):
        bag_path = os.path.join(DATA_PATH, 'H12607_MB_4m_MLLW_2of2.bag')
        coverage_list = [os.path.join(DATA_PATH, 'H12607_SSSAB_1m_600kHz_2of2.tif')]

        bag_dataset = bag.bag_file()
        bag_dataset.open_file(bag_path, 'hack')
        bag_dataset.generate_name(DATA_PATH, False)
        coverage_dataset = coverage.unified_coverage(coverage_list, bag_dataset.wkt, bag_dataset.name)

        assert coverage_dataset.bounds != bag_dataset.bounds
        assert coverage_dataset.shape != bag_dataset.shape

        coverage_dataset = coverage.align2grid(coverage_dataset, bag_dataset.bounds, bag_dataset.shape,
                                               bag_dataset.resolution, bag_dataset.nodata)

        assert coverage_dataset.bounds == bag_dataset.bounds
        assert coverage_dataset.shape == bag_dataset.shape


if __name__ == '__main__':
    unittest.main()
