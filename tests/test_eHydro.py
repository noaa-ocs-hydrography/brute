import os
import unittest

TESTING_DIRECTORY = r'\\OCS-VS-NBS02\data\testing'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')


class TestEHydro(unittest.TestCase):
    def test_scrape(self):
        raise NotImplementedError()

    def test_move(self):
        raise NotImplementedError()


if __name__ == '__main__':
    unittest.main()
