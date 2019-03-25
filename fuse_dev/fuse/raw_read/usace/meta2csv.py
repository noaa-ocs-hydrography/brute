# -*- coding: utf-8 -*-
"""
meta2csv

Created on Fri Aug 10 14:58:58 2018

@author: grice

This is a collection of methods for reading and writing metadata from surveys
and putting the data into a csv file.
"""
from pathlib import Path as _Path
from tempfile import NamedTemporaryFile as _NamedTemporaryFile
import shutil as _shutil
import csv as _csv

__version__ = 'meta2csv 0.0.1'

_cols = ['from_filename',
        'from_path',
        'to_filename',
        'start_date',
        'end_date',
        'from_fips',
        'from_horiz_datum',
        'from_horiz_units',
        'from_horiz_unc',
        'to_horiz_datum',
        'from_vert_datum',
        'from_vert_key',
        'from_vert_units',
        'from_vert_unc',
        'to_vert_datum',
        'to_vert_units',
        'agency',
        'source_indicator',
        'source_type',
        'complete_coverage',
        'complete_bathymetry',
        'vert_uncert_fixed',
        'vert_uncert_vari',
        'horiz_uncert',
        'feat_size',
        'feat_detect',
        'feat_least_depth',
        'interpolated',
        'script_version',
        ]


# ordered dict to ensure looping through the keys always gets 'manual' last.
_col_root = {'manual' : 'manual: ',
             'script' : 'script: '}

# this map translates the names used here to the ID used in the database
_field_map = {'from_filename' : 'OBJNAM',
              'start_date' : 'SURSTA',
              'end_date' : 'SUREND',
              'horiz_uncert' : 'POSACC',
              'to_horiz_datum' : 'HORDAT',
              'agency' : 'AGENCY',
              'source_indicator' : 'SORIND',
              'source_type' : 's_ftyp',
              'complete_coverage' : 'flcvrg',
              'complete_bathymetry' : 'flbath',
              'to_vert_datum' : 'VERDAT',
              'to_vert_units' : 'DUNITS',
              'vert_uncert_fixed' : 'vun_fx',
              'vert_uncert_vari' : 'vun_vb',
              'feat_size' : 'f_size',
              'feat_detect' : 'f_dtct',
              'feat_least_depth' : 'f_lstd',
              'interpolated' : 'interp',
              'reviewed' : 'r_name',
              'script_version' : 's_scpv',
              'source_indicator' : 'SORIND',
              }


_s57_cols = ['OBJNAM',
              'SURSTA',
              'SUREND',
              'POSACC',
              'HORDAT',
              'AGENCY',
              'SORIND',
              's_ftyp',
              'flcvrg',
              'flbath',
              'VERDAT',
              'DUNITS',
              'vun_fx',
              'vun_vb',
              'f_size',
              'f_dtct',
              'f_lstd',
              'interp',
              'r_name',
              's_scpv',
              'SORIND'
              ]


def write_meta2csv(meta, csvfilename):
    """
    Open the provided file and add the list of metadata in the provided 
    dictionaries.
    """
    infile = _Path(csvfilename)
    if infile.exists():
        _add_to_csv(meta, csvfilename)
    # just write a new file since there is not one already
    else:
        _write_new_csv(meta,csvfilename)

def write_meta2csv_CEMVN(meta, csvfilename):
    """
    Open the provided file and add the list of metadata in the provided 
    dictionaries.
    """
    infile = _Path(csvfilename)
    if infile.exists():
        _add_to_csv_CEMVN(meta, csvfilename)
    # just write a new file since there is not one already
    else:
        _write_new_csv_CEMVN(meta,csvfilename)
        
def _add_to_csv(meta, csvfilename):
    """
    Add the provided metadata to the provide file.
    """
    orig = []
    # get the names of all the files in the new metadata
    new_meta_files = []
    for m in meta:
        new_meta_files.append(m['from_filename'])
    # update the metadata keys
    meta = _scriptkeys(meta)
    # get all the metadata that is in the file already
    with open(csvfilename, 'r') as csvfile:
        reader = _csv.DictReader(csvfile, fieldnames = _make_col_header())
        for row in reader:
            orig.append(row)
    # move the data into a new temp file
    with _NamedTemporaryFile(mode = 'w', newline='', delete=False) as tempfile:
        writer = _csv.DictWriter(tempfile, 
                                 fieldnames = _make_col_header(),
                                 extrasaction = 'ignore')
        # check to see if the new metadata is the same as an existing file
        for row in orig:
            fname = row['from_filename']
            # if the file is being updated, move over the updated info
            if fname in new_meta_files:
                idx = new_meta_files.index(fname)
                m = meta.pop(idx)
                new_meta_files.pop(idx)
                for key in m.keys():
                    row[key] = m[key]
            writer.writerow(row)
        # append what remains to the file
        for row in meta:
            writer.writerow(row)
    # replace the original file with the temp file
    _shutil.move(tempfile.name, csvfilename)        
        
def _add_to_csv_CEMVN(meta, csvfilename):
    """
    Add the provided metadata to the provide file.
    """
    orig = []
    # get the names of all the files in the new metadata
    new_meta_files = []
    for m in meta:
        new_meta_files.append(m['from_filename'])
    # update the metadata keys
    #meta = _scriptkeys(meta)
    # get all the metadata that is in the file already
    with open(csvfilename, 'r') as csvfile:
        reader = _csv.DictReader(csvfile, fieldnames = _make_col_header())
        for row in reader:
            orig.append(row)
    # move the data into a new temp file
    with _NamedTemporaryFile(mode = 'w', newline='', delete=False) as tempfile:
        writer = _csv.DictWriter(tempfile, 
                                 fieldnames = _make_col_header(),
                                 extrasaction = 'ignore')
        # check to see if the new metadata is the same as an existing file
        for row in orig:
            fname = row['from_filename']
            # if the file is being updated, move over the updated info
            if fname in new_meta_files:
                idx = new_meta_files.index(fname)
                m = meta.pop(idx)
                new_meta_files.pop(idx)
                for key in m.keys():
                    row[key] = m[key]
            writer.writerow(row)
        # append what remains to the file
        for row in meta:
            writer.writerow(row)
    # replace the original file with the temp file
    _shutil.move(tempfile.name, csvfilename)
    
def _write_new_csv(meta, csvfilename):
    """
    Write the provided metadata to a new CSV file.
    """
    meta = _scriptkeys(meta)
    with open(csvfilename, 'w', newline='') as csvfile:
        writer = _csv.DictWriter(csvfile, 
                                 fieldnames = _make_col_header(),
                                 extrasaction = 'ignore')
        writer.writeheader()
        for row in meta:
            writer.writerow(row)
            
def _write_new_csv_CEMVN(meta, csvfilename):
    """
    Write the provided metadata to a new CSV file.
    """
    #meta = _scriptkeys(meta)
    with open(csvfilename, 'w', newline='') as csvfile:
        writer = _csv.DictWriter(csvfile, 
                                 fieldnames = _make_col_header(),
                                 extrasaction = 'ignore')
        writer.writeheader()
        for row in meta:
            writer.writerow(row)
                
def _make_col_header():
    """
    Return the column header names.
    """
    csv_cols = []
    for c in _cols:
        if c is 'from_filename' or c is 'from_path' or c is 'script_version':
            csv_cols.append(c)
        else:
            csv_cols.append(_col_root['script'] + c)
            csv_cols.append(_col_root['manual'] + c)
    csv_cols.append('reviewed')
    csv_cols.append('Last Updated')
    csv_cols.append('Notes')
    return csv_cols

    

def _scriptkeys(meta):
    """
    Prepend 'script: ' to each key in the list of dictionaries such that the
    list goes to the right column when written to a csv.
    """
    new_meta = []
    for row in meta:
        new_row = {}
        keys = row.keys()
        for key in keys:
            if key is 'from_filename' or key is 'from_path':
                new_row[key] = row[key]
            elif key is 'script_version':
                new_row[key] = row[key] + ',' + __version__
            else:
                new_row['script: ' + key] = row[key]
        new_meta.append(new_row)
    return new_meta

def read_simplified_csv(csvfilename):
    """
    Open the provide csv file name, extract the metadata, and combine
    duplicative rows, giving precedence to the manually entered values.
    """
    with open(csvfilename,'r') as csvfile:
        metadata = []
        reader = _csv.DictReader(csvfile)
        # get the row
        for row in reader:
            metadata.append(simplify_row(row))
    return metadata

def simplify_row(row):
    """
    Provided a dictionary representing a row from the csv file, combined the
    manual and scripted values, giving precedence to the manual entries.
    """
    metarow = {}
    # make dictionaries for sorting data into
    for name in _col_root:
        metarow[name] = {}
    metarow['base'] = {}
    # sort each key into the right dictionary
    for key in row:
        # only do stuff with keys that have information
        if len(row[key]) > 0:
            named = False
            for name in _col_root:
                if name in key:
                    named = True
                    val = key.replace(_col_root[name], '')
                    metarow[name][val] = row[key]
            if not named:
                metarow['base'][key] = row[key]
    # combine the dictionaries
    simplerow = {}
    names = [*metarow]
    names.sort(reverse = True)
    for name in names:
        simplerow = {**simplerow, **metarow[name]}
    return simplerow
