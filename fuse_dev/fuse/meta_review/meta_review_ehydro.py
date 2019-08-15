# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:47:59 2019

@author: grice
"""

import csv as _csv
import shutil as _shutil
from pathlib import Path as _Path
from tempfile import NamedTemporaryFile as _NamedTemporaryFile
from typing import List

import fuse.meta_review.meta_review as mrb


class MetaReviewer_eHydro(mrb.MetaReviewer):
    """The ehydro metadata object."""

    # ordered dict to ensure looping through the keys always gets 'manual' last.
    _col_root = {
        'manual': 'manual: ',
        'script': 'script: '
    }

    # this map translates the names used here to the ID used in the database
    _field_map = {
        'from_filename': 'OBJNAM',
        'start_date': 'SURSTA',
        'end_date': 'SUREND',
        'horiz_uncert': 'POSACC',
        'to_horiz_datum': 'HORDAT',
        'agency': 'AGENCY',
        'source_type': 's_ftyp',
        'complete_coverage': 'flcvrg',
        'complete_bathymetry': 'flbath',
        'to_vert_datum': 'VERDAT',
        'to_vert_units': 'DUNITS',
        'vert_uncert_fixed': 'vun_fx',
        'vert_uncert_vari': 'vun_vb',
        'feat_size': 'f_size',
        'feat_detect': 'f_dtct',
        'feat_least_depth': 'f_lstd',
        'interpolated': 'interp',
        'reviewed': 'r_name',
        'script_version': 's_scpv',
        'source_indicator': 'SORIND',
        'catzoc': 'CATZOC',
        'supersession_score': 'supscr',
    }

    _vert_datum = {
        'MLWS': '1',
        'MLLWS': '2',
        'MSL': '3',
        'LLW': '4',
        'MLW': '5',
        'ISLW': '8',
        'MLLW': '12',
        'MHW': '16',
        'MHWS': '17',
        'MHHW': '21',
        'LAT': '23',
        'LOC': '24',
        'IGLD': '25',
        'LLWLT': '27',
        'HHWLT': '28',
        'HAT': '30',
        'Unknown': '701',
        'Other': '703',
        'HRD': '24',  # Adding this for the Hudson River Datum
    }

    _horz_datum = {
        'WGS72': '1',
        'WGS84': '2',
        'WGS_1984': '2',
        'NAD27': '74',
        'NAD83': '75',
        'North_American_Datum_1983': '75',
        'Local': '131',
    }

    def __init__(self, metafile_path: str, meta_keys: List[str]):
        """
        TODO write description

        Parameters
        ----------
        metafile_path
        meta_keys
        """

        super().__init__(metafile_path, meta_keys)
        self._fieldnames = self._make_col_header()

    def _make_col_header(self) -> List[str]:
        """TODO write description"""
        csv_cols = []
        for c in self._metakeys:
            if c in ('from_filename', 'from_path', 'script_version'):
                csv_cols.append(c)
            else:
                csv_cols.append(MetaReviewer_eHydro._col_root['script'] + c)
                csv_cols.append(MetaReviewer_eHydro._col_root['manual'] + c)
        csv_cols.append('reviewed')
        csv_cols.append('Last Updated')
        csv_cols.append('Notes')
        return csv_cols

    def write_meta_record(self, meta: List[dict]):
        """
        Open the provided file and add the list of metadata in the provided
        dictionaries.

        Parameters
        ----------
        meta: Union[List[dict], dict] :
            TODO write description

        Returns
        -------

        """

        infile = _Path(self._metafilename)
        if infile.exists():
            if type(meta) == dict:  # just a single record
                self._add_to_csv([meta])
            elif type(meta) == list:  # this is a list of records
                self._add_to_csv(meta)
            else:
                raise ValueError('Unknown meta data container provided')
        # just write a new file since there is not one already
        else:
            self._write_new_csv(meta)

    def _add_to_csv(self, meta: List[dict]):
        """
        Add the provided metadata to the provide file.

        Parameters
        ----------
        meta: List[dict] :
            TODO write description

        Returns
        -------

        """

        orig = []
        # get the names of all the files in the new metadata
        new_meta_files = []
        for m in meta:
            new_meta_files.append(m['from_filename'])
        # update the metadata keys
        meta = self._scriptkeys(meta)
        # get all the metadata that is in the file already
        with open(self._metafilename, 'r', encoding='utf-8') as csvfile:
            reader = _csv.DictReader(csvfile, fieldnames=self._fieldnames)
            for row in reader:
                orig.append(row)
        # move the data into a new temp file
        with _NamedTemporaryFile(mode='w',
                                 newline='',
                                 encoding='utf-8',
                                 delete=False) as tempfile:
            writer = _csv.DictWriter(tempfile,
                                     fieldnames=self._fieldnames,
                                     extrasaction='ignore')
            # check to see if the new metadata is the same as an existing entry
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
        _shutil.move(tempfile.name, self._metafilename)

    def _write_new_csv(self, meta: List[dict]):
        """
        Write the provided metadata to a new CSV file.

        Parameters
        ----------
        meta: List[dict] :
            TODO write description

        Returns
        -------

        """

        if type(meta) == dict:  # just a single record
            meta = self._scriptkeys([meta])
        elif type(meta) == list:  # this is a list of records
            meta = self._scriptkeys(meta)
        else:
            raise ValueError('Unknown meta data container provided')
        with open(self._metafilename, 'w', newline='', encoding='utf-8') as csvfile:
            writer = _csv.DictWriter(csvfile, fieldnames=self._fieldnames, extrasaction='ignore')
            writer.writeheader()
            for row in meta:
                writer.writerow(row)

    def _scriptkeys(self, meta: List[dict]) -> List[dict]:
        """
        Prepend 'script: ' to each key in the list of dictionaries such that the
        list goes to the right column when written to a csv.

        Parameters
        ----------
        meta: List[dict] :
            TODO write description

        Returns
        -------

        """

        new_meta = []
        for row in meta:
            new_row = {}
            keys = row.keys()
            for key in keys:
                if key in ('from_filename', 'from_path'):
                    new_row[key] = row[key]
                elif key == 'script_version':
                    pass
                    # new_row[key] = f'{row[key]},{__version__}'
                else:
                    new_row[f'script: {key}'] = row[key]
            new_meta.append(new_row)
        return new_meta

    def read_meta_file(self) -> List[dict]:
        """
        Open the provide csv file name, extract the metadata, and combine
        duplicative rows, giving precedence to the manually entered values.

        Parameters
        ----------

        Returns
        -------

        """

        with open(self._metafilename, 'r', encoding='utf-8') as csvfile:
            metadata = []
            reader = _csv.DictReader(csvfile)
            # get the row
            for row in reader:
                metadata.append(self._simplify_row(row))
        return metadata

    def read_meta_record(self, meta_value, meta_key: str = 'from_filename') -> dict:
        """
        Open the provide csv file name, extract the metadata row looking for
        the record name that matches the provided key.  Once the record is
        found the record is "simplified" and returned.

        Parameters
        ----------
        meta_value :
            TODO write description
        meta_key: str :
             (Default value = 'from_filename')

        Returns
        -------

        """

        metadata = {}
        with open(self._metafilename, 'r', encoding='utf-8' ) as csvfile:
            reader = _csv.DictReader(csvfile)
            # get the row
            for row in reader:
                if row[meta_key] == meta_value:
                    metadata = self._simplify_row(row)
        return metadata

    def _simplify_row(self, row: dict) -> dict:
        """
        Provided a dictionary representing a row from the csv file, combined
        the manual and scripted values, giving precedence to the manual
        entries.

        Parameters
        ----------
        row: dict :
            TODO write description

        Returns
        -------

        """

        metarow = {}
        # make dictionaries for sorting data into
        for name in MetaReviewer_eHydro._col_root:
            metarow[name] = {}
        metarow['base'] = {}
        # sort each key (that has information) into the right dictionary
        for key in row:
            # only do stuff with keys that have information
            if len(row[key]) > 0:
                named = False
                for name in MetaReviewer_eHydro._col_root:
                    if name in key:
                        named = True
                        val = key.replace(MetaReviewer_eHydro._col_root[name], '')
                        metarow[name][val] = row[key]
                if not named:
                    metarow['base'][key] = row[key]
        # combine the dictionaries
        simplerow = {}
        names = [*metarow]
        names.sort(reverse=True)
        for name in names:
            simplerow = {**simplerow, **metarow[name]}
        return simplerow

    def row2s57(self, row: dict):
        """
        Convert the expanded column names used in the csv to an S57 name and value.

        Parameters
        ----------
        row: dict :
            TODO write description

        Returns
        -------

        """

        s57row = {}

        # remap the keys
        for key in row:
            if key in MetaReviewer_eHydro._field_map:
                s57row[MetaReviewer_eHydro._field_map[key]] = row[key]
                if row[key] in ('TRUE', 'True'):
                    s57row[MetaReviewer_eHydro._field_map[key]] = 0
                elif row[key] in ('FALSE', 'False'):
                    s57row[MetaReviewer_eHydro._field_map[key]] = 1

        # enforce additional required formating
        if 'VERDAT' in s57row:
            s57row['VERDAT'] = MetaReviewer_eHydro._vert_datum[s57row['VERDAT']]
        elif 'HORDAT' in s57row:
            h = s57row['HORDAT']
            for name in MetaReviewer_eHydro._horz_datum:
                if name in h:
                    s57row['HORDAT'] = MetaReviewer_eHydro._horz_datum[name]
                    break
        elif 'SUREND' in s57row:
            s57row['SORDAT'] = s57row['SUREND']

        return s57row
