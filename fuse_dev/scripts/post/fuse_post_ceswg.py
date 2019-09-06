# -*- coding: utf-8 -*-
"""
fuse_post_ceswg.py

grice
20190703
V.0.0.1 20190703
ed jk

This is a script to demonstrate the posting of a file with USACE data into a
CARIS BDB database.
"""

import fuse.fuse_ehydro as ffe

if __name__ == '__main__':
    poster = ffe.FuseProcessor_eHydro('ceswg.config')  # this config is local for testing
    flist = poster._meta_obj.read_metadata()
    for f in flist:
        if 'to_filename' in f:
            infilename = f['from_filename']
            poster.post(infilename)
    poster.disconnect()
