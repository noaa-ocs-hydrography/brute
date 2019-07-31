# -*- coding: utf-8 -*-
"""
fuse_rate.py

grice
20190725
V.0.0.1 20190725

This is a script to demonstrate the posting of the decay score of a file with
USACE data into a CARIS BDB database.
"""

import fuse.fuse_ehydro as ffe
import fuse.score as score
from datetime import datetime

infilename = 'NB_01_MAI_20160916_CS_4514_30X.XYZ'

if __name__ == '__main__':
    poster = ffe.fuse_ehydro('cenan.config')  # this config is local for testing
#    flist = poster._meta_obj.read_meta_file()
#    for f in flist:
#        if 'from_filename' in f:
#            infilename = f['from_filename']
    poster._get_stored_meta(infilename)
    poster._set_log(infilename)
    catzoc = score.catzoc(poster._meta)
    supscr = score.supersession(poster._meta)
    poster._meta['CATZOC'] = catzoc
    poster._meta['supersession_score'] = supscr
    poster._meta_obj.write_meta_record(poster._meta)
    now = datetime.now()
    poster.score(infilename, now)
    poster.disconnect()
