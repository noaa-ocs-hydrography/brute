# -*- coding: utf-8 -*-
"""
fuse_score_ceswg.py

grice
20190725
V.0.0.1 20190725
ed jk

This is a script to demonstrate the posting of the decay score of a file with
USACE data into a CARIS BDB database.
"""

from datetime import datetime

import fuse.fuse_ehydro as ffe

if __name__ == '__main__':
    now = datetime.now()
    poster = ffe.FuseProcessor_eHydro('ceswg.config')  # this config is local for testing
    flist = poster._meta_obj.read_metadata()
    for f in flist:
        if 'from_filename' in f:
            infilename = f['from_filename']
            poster.score(infilename, now)
    poster.disconnect()
