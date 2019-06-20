"""
meta_review_base.py

20190130 grice

A base class for reading and writing metadata for review.
"""

class meta_review_base:
    def __init__(self, metafilename, metakeys):

        """Initialize with the filename to use for meta data input / output."""
        self._metafilename = metafilename
        self._metakeys = metakeys
        
    def read_meta_record(self, record_key):
        """Read the meta data record from the meta file.

        :param record_key: 

        """
        pass
        
    def read_meta_file(self):
        """Read all records from the meta file."""
        pass
        
    def write_meta_record(self, record_dictionary):
        """Write the provided record to the file.

        :param record_dictionary: 

        """
