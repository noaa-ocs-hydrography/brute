import os
import unittest

TESTING_DIRECTORY = r'\\OCS-VS-NBS02\data\testing'
INPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'raw')
OUTPUT_ROOT = os.path.join(TESTING_DIRECTORY, 'output')


class TestRawRead(unittest.TestCase):
    def test_raw_read(self):
        raise NotImplementedError()


if __name__ == '__main__':
    unittest.main()
