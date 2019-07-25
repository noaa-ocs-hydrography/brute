import fuse.interpolator.bag_interpolator

def test_interpolate():
    assert 1 == 1


class TestBagInterpolator:
    def test_align2grid(self):
        fuse.interpolator.bag_interpolator.coverage.align2grid()
