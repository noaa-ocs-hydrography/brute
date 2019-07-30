import unittest
from tiling import tiling


class TestTiling(unittest.TestCase):
    def test_tiling(self):
        xy_level = '99'
        tiling.build_gpkg(xy_level, [-180, -90, 180, 90])

        assert tiling.get_tile_name(387382) == '0008AWM'
        assert tiling.get_tile_from_point('99', [-70.5, 43]) == (
            [-70.6640625, 42.890625, -70.3125, 43.2421875], '0008AWM')


if __name__ == '__main__':
    unittest.main()
