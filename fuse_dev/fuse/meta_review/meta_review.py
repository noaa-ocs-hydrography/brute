# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:47:59 2019

@author: grice
"""

import ast as _ast
import csv as _csv
import os
import shutil as _shutil
from tempfile import NamedTemporaryFile as _NamedTemporaryFile


class MetaReviewer:
    """The ehydro metadata object."""

    # ordered dict to ensure looping through the keys always gets 'manual' last.
    _prefixes = {
        'manual': 'manual: ',
        'script': 'script: '
    }

    # this map translates the names used here to the ID used in the database
    _database_keys = {
        'from_filename': 'OBJNAM',
        'start_date': 'SURSTA',
        'end_date': 'SUREND',
        'horiz_uncert': 'POSACC',
        'to_horiz_datum': 'HORDAT',
        'agency': 'AGENCY',
        'source_type': 's_ftyp',
        'complete_coverage': 'flcvrg',
        'complete_bathymetry': 'flbath',
        'to_vert_key': 'VERDAT',
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

    _vertical_datums = {
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
        'HRD': '24',  # Hudson River Datum
    }

    _horizontal_datums = {
        'WGS72': '1',
        'WGS84': '2',
        'WGS_1984': '2',
        'NAD27': '74',
        'NAD83': '75',
        'North_American_Datum_1983': '75',
        'Local': '131',
    }

    def __init__(self, metadata_filename: str, metadata_keys: [str]):
        """
        Create new metadata review object from the given CSV filename and a dictionary of metadata.

        Parameters
        ----------
        metadata_filename
            path to CSV file
        metadata_keys
            dictionary of metadata
        """

        self._metadata_filename = metadata_filename
        self._metadata_keys = metadata_keys
        self._fieldnames = self._csv_header()

    def _csv_header(self) -> [str]:
        """get CSV header as a list of strings"""

        csv_cols = []
        for metadata_key in self._metadata_keys:
            if metadata_key in ('from_filename', 'from_path', 'script_version'):
                csv_cols.append(metadata_key)
            else:
                csv_cols.extend(MetaReviewer._prefixes[prefix] + metadata_key for prefix in ('script', 'manual'))
        csv_cols.extend(('reviewed', 'Last Updated', 'Notes'))
        return csv_cols

    def write_meta_record(self, metadata: [dict]):
        """
        Add the list of metadata to the CSV file.

        Parameters
        ----------
        metadata
            dictionary of metadata
        """

        if os.path.exists(self._metadata_filename):
            # check if only a single record was provided
            if type(metadata) is dict:
                self._add_to_csv([metadata])
            else:
                self._add_to_csv(metadata)
        # just write a new file since there is not one already
        else:
            self._write_new_csv([metadata])

    def _add_to_csv(self, metadata: [dict]):
        """
        Add the provided metadata to the CSV file.

        Parameters
        ----------
        metadata
            list of metadata dictionaries
        """

        orig = []
        # get the names of all the files in the new metadata
        new_meta_files = []
        for m in metadata:
            new_meta_files.append(m['from_filename'])
        # update the metadata keys
        metadata = self._prepend_script_to_keys(metadata)
        # get all the metadata that is in the file already
        with open(self._metadata_filename, 'r', encoding='utf-8') as csvfile:
            reader = _csv.DictReader(csvfile, fieldnames=self._fieldnames)
            for row in reader:
                orig.append(row)
        csvfile.close()
        # move the data into a new temp file
        with _NamedTemporaryFile(mode='w', newline='', encoding='utf-8', delete=False) as tempfile:
            writer = _csv.DictWriter(tempfile, fieldnames=self._fieldnames, extrasaction='ignore')
            # check to see if the new metadata is the same as an existing entry
            for row in orig:
                fname = row['from_filename']
                # if the file is being updated, move over the updated info
                if fname in new_meta_files:
                    idx = new_meta_files.index(fname)
                    m = metadata.pop(idx)
                    new_meta_files.pop(idx)
                    for key in m.keys():
                        row[key] = m[key]
                writer.writerow(row)
            # append what remains to the file
            for row in metadata:
                writer.writerow(row)
            del writer
        tempfile.close()
        # replace the original file with the temp file
        _shutil.move(tempfile.name, self._metadata_filename)

    def _write_new_csv(self, metadata: [dict]):
        """
        Write the provided metadata to a new CSV file.

        Parameters
        ----------
        metadata
            list of metadata dictionaries
        """

        # check if only a single record was provided
        if type(metadata) is dict:
            metadata = self._prepend_script_to_keys([metadata])
        else:
            metadata = self._prepend_script_to_keys(metadata)

        with open(self._metadata_filename, 'w', newline='', encoding='utf-8') as csv_file:
            writer = _csv.DictWriter(csv_file, fieldnames=self._fieldnames, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(metadata)
            csv_file.close()

    def _prepend_script_to_keys(self, metadata: [dict]) -> [dict]:
        """
        Prepend 'script: ' to each key in the list of metadata such that the list goes to the right column when written to a CSV

        Parameters
        ----------
        metadata
            list of metadata dictionaries

        Returns
        -------
            list of metadata dictionaries with `script: ` prepended to keys
        """

        # TODO is `script_version` needed? if so, set the value to f'{row[key]},{__version__}' if provided in the key
        return [{f'script: {key}' if key not in ('from_filename', 'from_path') else key: value for key, value in
                 row.items()} for row in metadata]

    def read_metadata(self) -> [dict]:
        """
        Get a list of combined metadata rows from the CSV.

        Returns
        -------
            list of dictionaries of simplified metadata rows
        """

        with open(self._metadata_filename, 'r', encoding='utf-8') as csv_file:
            return [self._simplify_row(row) for row in _csv.DictReader(csv_file)]

    def read_meta_record(self, search_value: str, meta_key: str = 'from_filename') -> dict:
        """
        Extract the simplified metadata row in the given CSV file where the value of the given key matches the given key.

        Parameters
        ----------
        search_value
            value to find in the CSV file
        meta_key
            key to search

        Returns
        -------
            dictionary of matching row
        """
        meta = {}
        with open(self._metadata_filename, 'r', encoding='utf-8') as csv_file:
            for row in _csv.DictReader(csv_file):
                if row[meta_key] == search_value:
                    meta = self._simplify_row(row)
                    break
            csv_file.close()
        return meta

    def _simplify_row(self, row: dict) -> dict:
        """
        Combine `manual: ` and `script: ` fields in the given CSV row, preferring manual entries.

        Parameters
        ----------
        row
            dicitonary of keys and values from a single row of the metadata CSV

        Returns
        -------
            simplified dictionary
        """

        metarow = {}
        # make dictionaries for sorting data into
        for name in MetaReviewer._prefixes:
            metarow[name] = {}
        metarow['base'] = {}
        # sort each key (that has information) into the right dictionary
        for key in row:
            # only do stuff with keys that have information
            if len(row[key]) > 0:
                named = False
                for name in MetaReviewer._prefixes:
                    if name in key:
                        named = True
                        val = key.replace(MetaReviewer._prefixes[name], '')
                        if val == 'support_files':
                            metarow[name][val] = _ast.literal_eval(row[key])
                        else:
                            metarow[name][val] = row[key]
                if not named:
                    metarow['base'][key] = row[key]
        # combine the dictionaries
        simplified_row = {}
        for name in sorted(list(metarow.keys()), reverse=True):
            simplified_row = {**simplified_row, **metarow[name]}
        return simplified_row

    def csv_to_s57(self, row: dict) -> dict:
        """
        Convert the expanded column names used in the CSV to S57 names / values.

        Parameters
        ----------
        row
            dictionary of CSV row

        Returns
        -------
            dictionary with S57 names / values
        """

        s57_row = {}

        # remap the keys
        for key, value in row.items():
            if key in MetaReviewer._database_keys:
                database_key = MetaReviewer._database_keys[key]

                if value in ('TRUE', 'True'):
                    s57_row[database_key] = 1
                elif value in ('FALSE', 'False'):
                    s57_row[database_key] = 0
                else:
                    s57_row[database_key] = value

        # enforce additional required formating
        if 'VERDAT' in s57_row:
            verdat = s57_row['VERDAT'].upper()
            s57_row['VERDAT'] = MetaReviewer._vertical_datums[verdat]
        elif 'HORDAT' in s57_row:
            h_datum = s57_row['HORDAT']
            for name in MetaReviewer._horizontal_datums:
                if name in h_datum:
                    s57_row['HORDAT'] = MetaReviewer._horizontal_datums[name]
                    break
        elif 'SUREND' in s57_row:
            s57_row['SORDAT'] = s57_row['SUREND']

        return s57_row
