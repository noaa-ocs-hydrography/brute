"""
meta_review_base.py

20190130 grice

A base class for reading and writing metadata for review.
"""
from typing import List


class meta_review_base:
    def __init__(self, metafilename: str, metakeys: List[str]):
        """
        Initialize with the filename to use for meta data input / output.

        Parameters
        ----------
        metafilename
        metakeys
        """

        self._metafilename = metafilename
        self._metakeys = metakeys

    def read_meta_record(self, record_key: str):
        """
        Read the meta data record from the meta file.

        :param record_key: 
        """

        pass

    def read_meta_file(self):
        """
        Read all records from the meta file.
        """

        pass

    def write_meta_record(self, record_dictionary: dict):
        """
        Write the provided record to the file.

        :param record_dictionary: 
        """

        pass
