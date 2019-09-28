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
import fuse.fuse_processor as ffp

if __name__ == '__main__':
    config_list = glob('*.config')
    for n,config in enumerate(config_list):
        poster = ffp.FuseProcessor(config)
        flist = poster._meta_obj.read_metadata()
        db = poster._config['database_name']
        print(f'working on {config}, posting to {db}')
        for m,f in enumerate(flist):
            if 'to_filename' in f:
                infilename = f['from_filename']
                print(f'{n}.{m}:{infilename}')
                poster.post(infilename)
        poster.disconnect()
