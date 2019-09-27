# -*- coding: utf-8 -*-
"""
fuse_post.py

grice
20190703
V.0.0.2 20190927

This is a script to demonstrate the posting of a file with USACE data into a
CARIS BDB database.

All configs in the local directory are run serially.
"""

from glob import glob
import fuse.fuse_ehydro as ffe

if __name__ == '__main__':
    config_list = glob('*.config')
    for config in config_list:
        poster = ffe.FuseProcessor_eHydro('cenan.config')
        flist = poster._meta_obj.read_metadata()
        for f in flist:
            if 'to_filename' in f:
                infilename = f['from_filename']
                poster.post(infilename)
        poster.disconnect()
