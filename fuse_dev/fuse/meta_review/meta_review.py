# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:47:59 2019

@author: grice
"""

import ast as _ast
import csv as _csv
import os
import shutil as _shutil
from abc import ABC, abstractmethod
from collections import OrderedDict
from tempfile import NamedTemporaryFile as _NamedTemporaryFile
from typing import Union, Any

import psycopg2

DATABASE_CREDENTIALS_FILENAME = r"D:\credentials\postgres_scripting.txt"
POSTGRES_DEFAULT_PORT = '5432'

# this map translates the names used here to the ID used in the database
S57_BASE_TRANSLATIONS = {
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
    'supersession_score': 'supscr'
}
S57_VERTICAL_DATUM_TRANSLATIONS = {
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
    'HRD': '24'  # Hudson River Datum
}
S57_HORIZONTAL_DATUM_TRANSLATIONS = {
    'WGS72': '1',
    'WGS84': '2',
    'WGS_1984': '2',
    'NAD27': '74',
    'NAD83': '75',
    'North_American_Datum_1983': '75',
    'Local': '131'
}

FIELDS_EXCLUDED_FROM_PREFIX = ('from_filename', 'from_path')


class MetadataTable(ABC):
    """abstract class representing a table containing metadata records of bathymetric surveys"""

    # order key iteration in column prefixes
    column_prefixes: OrderedDict
    metadata_fields: [str]

    @property
    def column_names(self) -> [str]:
        columns = []
        for metadata_key in self.metadata_fields:
            if metadata_key in FIELDS_EXCLUDED_FROM_PREFIX or metadata_key == 'script_version':
                columns.append(metadata_key)
            else:
                columns.extend(f'{self.column_prefixes[prefix]}{metadata_key}' for prefix in self.column_prefixes)
        columns.extend(('reviewed', 'last_updated', 'notes'))
        return columns

    @property
    @abstractmethod
    def records(self) -> [dict]:
        raise NotImplementedError

    @abstractmethod
    def records_where(self, where: dict) -> [dict]:
        """
        Query table for matching key-value pairs.

        Parameters
        ----------
        where
            dictionary mapping keys to values, with which to match records

        Returns
        -------
        [dict]
            dictionaries of matching records
        """

        raise NotImplementedError

    @abstractmethod
    def __getitem__(self, primary_key_value: Any) -> dict:
        """
        Query table for the given value of the primary key.

        Parameters
        ----------
        primary_key_value
            value to query from primary key

        Returns
        -------
        dict
            dictionary of matching record
        """

        raise NotImplementedError

    def __setitem__(self, primary_key_value: Any, record: dict):
        """
        Insert the given record into the table with the given primary key value.

        Parameters
        ----------
        primary_key_value
            value of primary key at which to insert record
        record
            dictionary record
        """

        self.insert_records([record])

    @abstractmethod
    def insert_records(self, records: [dict]):
        """
        Insert the list of metadata records into the table.

        Parameters
        ----------
        records
            metadata dictionaries
        """

        raise NotImplementedError

    def _prepend_script_to_keys(self, metadata: [dict]) -> [dict]:
        """
        Prepend 'script' prefix to each applicable key in the list of metadata dictionaries.

        Parameters
        ----------
        metadata
            list of metadata dictionaries

        Returns
        -------
        [dict]
            list of metadata dictionaries with `script: ` prepended to keys
        """

        # TODO is `script_version` needed? if so, set the value to f'{row[key]},{__version__}' if provided in the key
        return [{f'{self.column_prefixes["script"]}{key}' if key not in FIELDS_EXCLUDED_FROM_PREFIX else key: value
                 for key, value in row.items()} for row in metadata]

    def _simplify_record(self, record: dict) -> dict:
        """
        Combine `manual` and `script` fields in the given record, preferring manual entries.

        Parameters
        ----------
        record
            dictionary of keys and values from a single record of the metadata

        Returns
        -------
        dict
            simplified dictionary
        """

        # make dictionaries for sorting data into
        entries_by_prefix = {'base': {}}
        for prefix_name in self.column_prefixes:
            entries_by_prefix[prefix_name] = {}

        # sort each key (that has information) into the dictionary by prefix
        for column_name, value in record.items():
            if value is not None and value != '':
                for prefix_name, prefix in self.column_prefixes.items():
                    if prefix in column_name:
                        metadata_key = column_name.replace(prefix, '')
                        if type(value) is str and metadata_key == 'support_files':
                            entries_by_prefix[prefix_name][metadata_key] = _ast.literal_eval(value)
                        elif type(value) is str and (value.capitalize() == 'True' or value.capitalize() == 'False'):
                            entries_by_prefix[prefix_name][metadata_key] = _ast.literal_eval(value.capitalize())
                        else:
                            entries_by_prefix[prefix_name][metadata_key] = value
                        break
                else:
                    entries_by_prefix['base'][column_name] = value

        # combine the dictionaries, overwriting the prefixed columns according to the order specified in the class attribute
        simplified_row = {}
        for prefix_name in self.column_prefixes:
            simplified_row.update(entries_by_prefix[prefix_name])
        return simplified_row


class MetadataDatabase(MetadataTable):
    """PostGreSQL database table containing metadata records of bathymetric surveys"""

    column_prefixes = OrderedDict([('script', 'script_'), ('manual', 'manual_')])

    postgres_field_types = {
        'from_filename': 'VARCHAR',
        'from_path': 'VARCHAR',
        'to_filename': 'VARCHAR',
        'support_files': 'VARCHAR[]',
        'start_date': 'DATE',
        'end_date': 'DATE',
        'from_horiz_datum': 'VARCHAR',
        'from_horiz_frame': 'VARCHAR',
        'from_horiz_type': 'VARCHAR',
        'from_horiz_units': 'VARCHAR',
        'from_horiz_key': 'VARCHAR',
        'from_vert_datum': 'VARCHAR',
        'from_vert_key': 'VARCHAR',
        'from_vert_units': 'VARCHAR',
        'from_vert_direction': 'VARCHAR',
        'to_horiz_frame': 'VARCHAR',
        'to_horiz_type': 'VARCHAR',
        'to_horiz_units': 'VARCHAR',
        'to_horiz_key': 'VARCHAR',
        'to_vert_datum': 'VARCHAR',
        'to_vert_key': 'VARCHAR',
        'to_vert_units': 'VARCHAR',
        'to_vert_direction': 'VARCHAR',
        'from_horiz_unc': 'REAL',
        'from_horiz_resolution': 'REAL',
        'from_vert_unc': 'REAL',
        'complete_coverage': 'BOOL',
        'bathymetry': 'BOOL',
        'vert_uncert_fixed': 'REAL',
        'vert_uncert_vari': 'REAL',
        'horiz_uncert_fixed': 'REAL',
        'horiz_uncert_vari': 'REAL',
        'to_horiz_resolution': 'REAL',
        'feat_size': 'REAL',
        'feat_detect': 'BOOL',
        'feat_least_depth': 'REAL',
        'catzoc': 'VARCHAR',
        'supersession_score': 'VARCHAR',
        'agency': 'VARCHAR',
        'source_indicator': 'VARCHAR',
        'source_type': 'VARCHAR',
        'interpolated': 'BOOL',
        'posted': 'BOOL',
        'license': 'VARCHAR',
        'logfilename': 'VARCHAR',
        'version_reference': 'VARCHAR',
        'interpolate': 'BOOL',
        'file_size': 'REAL'
    }

    def __init__(self, hostname: str, database: str, table: str, fields: [str], primary_key: str = 'from_filename'):
        """
        Create new metadata review object from the given CSV filename and a dictionary of metadata.

        Parameters
        ----------
        hostname
            hostname of PostGres server
        database
            name of metadata database
        table
            name of metadata table
        fields
            list of metadata entry names
        """

        # parse port from URL
        self.hostname, self.port = split_URL_port(hostname)
        if self.port is None:
            self.port = POSTGRES_DEFAULT_PORT

        self.database_name = database
        self.table_name = table
        self.metadata_fields = fields
        self.primary_key = primary_key

        with open(DATABASE_CREDENTIALS_FILENAME) as database_credentials_file:
            lines = [line.strip() for line in database_credentials_file.readlines()]
            database_username, database_password = lines[:2]

        self.connection = psycopg2.connect(database=self.database_name, user=database_username, password=database_password,
                                           host=self.hostname, port=self.port)

        with self.connection:
            with self.connection.cursor() as cursor:
                if not database_has_table(cursor, self.table_name):
                    data_types = {key: value for key, value in self.postgres_field_types.items() if key in FIELDS_EXCLUDED_FROM_PREFIX}
                    data_types.update({f'{prefix}{key}': value for prefix_name, prefix in self.column_prefixes.items()
                                       for key, value in self.postgres_field_types.items() if key not in FIELDS_EXCLUDED_FROM_PREFIX})
                    data_types.update({key: 'VARCHAR' for key in self.column_names if key not in data_types})

                    schema = ['from_filename VARCHAR PRIMARY KEY'] + [f'{field_name} {data_types[field_name]}'
                                                                      for field_name in self.column_names if field_name != 'from_filename']
                    cursor.execute(f'CREATE TABLE {self.table_name} ({", ".join(schema)});')

    @property
    def records(self) -> [dict]:
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(f'SELECT * FROM {self.table_name}')
                records = cursor.fetchall()

        return [self._simplify_record(dict(zip(self.column_names, record))) for record in records]

    def records_where(self, where: dict) -> [dict]:
        where_clause = ' AND '.join(f'{key} = %s' for key in where.keys())

        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(f'SELECT * FROM {self.table_name} WHERE {where_clause}', [*where.values()])
                records = cursor.fetchall()

        return [self._simplify_record(dict(zip(self.column_names, record))) for record in records]

    def __getitem__(self, primary_key_value: Any) -> dict:
        with self.connection:
            with self.connection.cursor() as cursor:
                cursor.execute(f'SELECT * FROM {self.table_name} WHERE {self.primary_key} = %s', [primary_key_value])
                record = cursor.fetchone()

        return self._simplify_record(dict(zip(self.column_names, record)))

    def insert_records(self, records: [dict]):
        if type(records) is dict:
            records = [records]

        assert all(self.primary_key in record for record in records), f'one or more records does not contain "{self.primary_key}"'

        with self.connection:
            with self.connection.cursor() as cursor:
                for record in records:
                    fields_in_record = [field for field in self.metadata_fields if field in record]
                    columns = [self.column_prefixes["script"] + field if field not in FIELDS_EXCLUDED_FROM_PREFIX else field
                               for field in fields_in_record]
                    values = [record[field] for field in fields_in_record]

                    if table_has_record(cursor, self.table_name, record, self.primary_key):
                        cursor.execute(f'UPDATE {self.table_name} SET ({", ".join(columns)}) = %s;', [tuple(values)])
                    else:
                        cursor.execute(f'INSERT INTO {self.table_name} ({", ".join(columns)}) VALUES %s;', [tuple(values)])


class MetadataFile(MetadataTable):
    """CSV table containing metadata records of bathymetric surveys"""

    column_prefixes = OrderedDict([('script', 'script: '), ('manual', 'manual: ')])

    def __init__(self, filename: str, fields: [str], primary_key: str = 'from_filename'):
        """
        Create new metadata review object from the given CSV filename and a dictionary of metadata.

        Parameters
        ----------
        filename
            path to CSV file
        fields
            dictionary of metadata
        """

        self.filename = filename
        self.metadata_fields = fields
        self.primary_key = primary_key

    @property
    def records(self) -> [dict]:
        with open(self.filename, 'r', encoding='utf-8') as csv_file:
            return [self._simplify_record(row) for row in _csv.DictReader(csv_file)]

    def records_where(self, where: dict) -> [dict]:
        records = []
        with open(self.filename, 'r', encoding='utf-8') as csv_file:
            for row in _csv.DictReader(csv_file):
                if all(row[key] == value for key, value in where):
                    records.append(row)

        return [self._simplify_record(record) for record in records]

    def __getitem__(self, primary_key_value: str) -> dict:
        with open(self.filename, 'r', encoding='utf-8') as csv_file:
            for row in _csv.DictReader(csv_file):
                if row['from_filename'] == primary_key_value:
                    return self._simplify_record(row)
            else:
                return {}

    def __setitem__(self, primary_key_value: str, record: dict):
        super().__setitem__(primary_key_value, record)

    def insert_records(self, records: [dict]):
        assert all(self.primary_key in record for record in records), f'one or more records does not contain "{self.primary_key}"'

        if os.path.exists(self.filename):
            # check if only a single record was provided
            if type(records) is dict:
                records = [records]

            self._add_to_csv(records)
        # just write a new file since there is not one already
        else:
            self._write_new_csv([records])

    def _add_to_csv(self, records: [dict]):
        """
        Add the provided metadata to the CSV file.

        Parameters
        ----------
        records
            list of metadata dictionaries
        """

        orig = []
        # get the names of all the files in the new metadata
        new_meta_files = []
        for record in records:
            new_meta_files.append(record['from_filename'])
        # update the metadata keys
        records = self._prepend_script_to_keys(records)
        # get all the metadata that is in the file already
        with open(self.filename, 'r', encoding='utf-8') as csv_file:
            orig.extend(row for row in _csv.DictReader(csv_file, fieldnames=self.column_names))
        # move the data into a new temp file
        with _NamedTemporaryFile(mode='w', newline='', encoding='utf-8', delete=False) as tempfile:
            writer = _csv.DictWriter(tempfile, fieldnames=self.column_names, extrasaction='ignore')
            # check to see if the new metadata is the same as an existing entry
            for row in orig:
                fname = row['from_filename']
                # if the file is being updated, move over the updated info
                if fname in new_meta_files:
                    idx = new_meta_files.index(fname)
                    record = records.pop(idx)
                    new_meta_files.pop(idx)
                    for key in record.keys():
                        row[key] = record[key]
                writer.writerow(row)
            # append what remains to the file
            for row in records:
                writer.writerow(row)
        # replace the original file with the temp file
        _shutil.move(tempfile.name, self.filename)

    def _write_new_csv(self, records: [dict]):
        """
        Write the provided metadata to a new CSV file.

        Parameters
        ----------
        records
            list of metadata dictionaries
        """

        # check if only a single record was provided
        if type(records) is dict:
            records = self._prepend_script_to_keys([records])
        else:
            records = self._prepend_script_to_keys(records)

        with open(self.filename, 'w', newline='', encoding='utf-8') as csv_file:
            writer = _csv.DictWriter(csv_file, fieldnames=self.column_names, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(records)


def csv_to_s57(row: dict) -> dict:
    """
    Convert the expanded column names used in the CSV to S57 names / values.

    Parameters
    ----------
    row
        dictionary of CSV row

    Returns
    -------
    dict
        dictionary with S57 names / values
    """

    s57_row = {}

    # remap the keys
    for key, value in row.items():
        if key in S57_BASE_TRANSLATIONS:
            key_in_s57 = S57_BASE_TRANSLATIONS[key]

            if value in ('TRUE', 'True'):
                s57_row[key_in_s57] = 1
            elif value in ('FALSE', 'False'):
                s57_row[key_in_s57] = 0
            else:
                s57_row[key_in_s57] = value

    # enforce additional required formatting
    if 'VERDAT' in s57_row:
        s57_row['VERDAT'] = S57_VERTICAL_DATUM_TRANSLATIONS[s57_row['VERDAT'].upper()]
    elif 'HORDAT' in s57_row:
        h_datum = s57_row['HORDAT']
        for name in S57_HORIZONTAL_DATUM_TRANSLATIONS:
            if name in h_datum:
                s57_row['HORDAT'] = S57_HORIZONTAL_DATUM_TRANSLATIONS[name]
                break
    elif 'SUREND' in s57_row:
        s57_row['SORDAT'] = s57_row['SUREND']

    return s57_row


def database_has_table(cursor: psycopg2._psycopg.cursor, table: str) -> bool:
    """
    Whether the given table exists within the given database.

    Parameters
    ----------
    cursor
        psycopg2 cursor
    table
        name of table

    Returns
    -------
    bool
        whether table exists
    """

    cursor.execute(f'SELECT EXISTS(SELECT 1 FROM information_schema.tables WHERE table_name=\'{table}\');')
    return cursor.fetchone()[0]


def table_has_record(cursor: psycopg2._psycopg.cursor, table: str, record: dict, primary_key: str) -> bool:
    """
    Whether the given table exists within the given database.

    Parameters
    ----------
    cursor
        psycopg2 cursor
    table
        name of table
    record
        dictionary record
    primary_key
        name of primary key

    Returns
    -------
    bool
        whether table exists
    """

    # if primary_key is None:
    #     cursor.execute(f'SELECT 1 FROM information_schema.table_constraints ' +
    #                    f'WHERE table_name=\'{table}\' AND constraint_type= \'PRIMARY KEY\';')
    #     primary_key_index = cursor.fetchone()[0] - 1
    #
    #     cursor.execute(f'SELECT * FROM information_schema.columns WHERE table_name=\'{table}\';')
    #     primary_key = cursor.fetchall()[primary_key_index]

    cursor.execute(f'SELECT EXISTS(SELECT 1 FROM {table} WHERE {primary_key}=\'{record[primary_key]}\');')
    return cursor.fetchone()[0]


def split_URL_port(url: str) -> (str, Union[str, None]):
    """
    Split the given URL into host and port, assuming port is appended after a colon.

    Parameters
    ----------
    url
        URL string

    Returns
    ----------
    str, Union[str, None]
        URL and port (if found)
    """

    port = None

    if url.count(':') > 0:
        url = url.split(':')
        if 'http' in url:
            url = ':'.join(url[:2])
            if len(url) > 2:
                port = url[2]
        else:
            url, port = url

    return url, port
