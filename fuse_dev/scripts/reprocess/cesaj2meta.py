"""
cesaj2meta.py

grice
20190408
V.0.0.1 20190408
ed jk

This is a script to demonstrate the processing of data that did not process on 
the initial run.
"""

import fuse.fuse_ehydro as ffe

if __name__ == '__main__':
    processor = ffe.FuseProcessor_eHydro('cesaj.config')  # this config is local for testing
    metalist = processor._meta_obj.read_meta_file()
    flist = []
    for f in metalist:
        flist.append(f['from_filename'])
    for f in metalist:
        if 'to_filename' not in f:
            infilename = f['from_path']
            print(f'Begin working on {infilename}')
            processor.process(infilename)
        if 'interpolate' in f and f['interpolate'].upper() == 'TRUE':
            interp_name = processor._interpolator.gettag(f['from_filename'])
            if interp_name not in flist:
                print(f'Begin working on {infilename}')
                processor.process(infilename)
