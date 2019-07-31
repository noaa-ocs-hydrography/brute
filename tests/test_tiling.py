import unittest
from tiling import tiling


class TestTiling(unittest.TestCase):
    def test_tiling(self):
        xy_level = '44'
        tiling.build_gpkg(xy_level, [-180, -90, 180, 90])

        assert tiling.get_tile_name(361) == '00000A1'
        assert tiling.get_tile_from_point(xy_level, [-70.8749, 43.1889]) == (
            [-70.6640625, 42.890625, -70.3125, 43.2421875], '00000A1')

        assert tiling.get_tile_from_point('88', [70.5, 43]) == ([70.3125, 42.890625, 71.015625, 43.59375], '0000H2P')
        assert tiling.get_shapely('99', [42.890625, 70.3125, 43.59375, 71.015625])[-1] == ['000A0SA', '000A0SB',
                                                                                           '000A1KQ', '000A1KR']


if __name__ == '__main__':
    unittest.main()
