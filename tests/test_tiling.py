import unittest

from tiling import tiling


class TestTiling(unittest.TestCase):
    def test_tiling(self):
        xy_level = '44'
        test_point = (-70.87, 43.18)
        tiling.build_gpkg(xy_level, [-180, -90, 180, 90])

        assert tiling.get_tile_name(361) == '00000A1'
        assert tiling.get_tile_from_point(xy_level, test_point) == ([-78.75, 33.75, -67.5, 45.0], '00000A1')

        assert tiling.get_tile_from_point('33', test_point) == ([-90.0, 22.5, -67.5, 45.0], '000002C')
        assert tiling.get_shapely(xy_level, [-90.0, 22.5, -67.5, 45.0])[-1] == ['0000094', '0000095', '00000A0', '00000A1']

        # os.remove(os.path.join(os.curdir, f'{xy_level}_tesselation.gpkg'))


if __name__ == '__main__':
    unittest.main()
