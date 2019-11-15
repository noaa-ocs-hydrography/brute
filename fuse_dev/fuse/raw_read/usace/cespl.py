# -*- coding: utf-8 -*-
"""
Edited by Juliet Kinney
based on extract_ehydro_meta_class_CESAJ.py

This script takes as an input the filename and path of an USACE e-Hydro .xyz
text file and pull the metadata from the xyz header and file name.
Additionally, it passes the matching .xml file and utilizes
the parse_usace_xml script to pull out the relevant metadata if it
is either FGDC or ISO FGDC USACE metadata format.

update 4/5/19
major update April 2, 2019
update July 12,2019 adding in call to pickle reader
Fri Aug  9 08:50:35 2019; Casiano Edits
"""

from fuse.raw_read.usace.usace import USACERawReader


class CESPLRawReader(USACERawReader):
    def __init__(self):
        super().__init__('CESPL')
