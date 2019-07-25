import unittest

from fuse.interpolator.bag_interpolator import bag
from fuse.interpolator.bag_interpolator import coverage


def test_interpolate():
    assert 1 == 1


class TestBagInterpolator(unittest.TestCase):
    def test_align2grid(self):
        bag_path = r"C:\Data\NBS\H12607_MB_4m_MLLW_2of2.bag"
        output_path = r"C:\Data\NBS"
        coverage_list = [r"C:\Data\NBS\H12607_SSSAB_1m_600kHz_2of2.tif"]

        bag_dataset = bag.bag_file()
        bag_dataset.open_file(bag_path, 'hack')
        bag_dataset.generate_name(output_path, False)
        coverage_object = coverage.unified_coverage(coverage_list, bag_dataset.wkt, bag_dataset.name)

        assert coverage_object.bounds != bag_dataset.bounds
        assert coverage_object.shape != bag_dataset.shape

        coverage_object = coverage.align2grid(coverage_object, bag_dataset.bounds, bag_dataset.shape,
                                              bag_dataset.resolution, bag_dataset.nodata)

        assert coverage_object.bounds == bag_dataset.bounds
        assert coverage_object.shape == bag_dataset.shape


if __name__ == '__main__':
    unittest.main()
