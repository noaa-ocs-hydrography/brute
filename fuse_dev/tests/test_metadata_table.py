import os
import unittest
from tempfile import TemporaryDirectory

import psycopg2
from fuse.meta_review.meta_review import MetadataDatabase, MetadataFile, database_has_table, split_URL_port, FIELD_TYPES, \
    POSTGRES_DEFAULT_PORT

POSTGRES_HOSTNAME_FILENAME = r"D:\credentials\postgres_hostname.txt"
DATABASE_CREDENTIALS_FILENAME = r"D:\credentials\postgres_scripting.txt"


class TestMetadataDatabase(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.database_name = 'metadata'

        with open(POSTGRES_HOSTNAME_FILENAME) as postgres_hostname_file:
            hostname = postgres_hostname_file.readline()

        self.hostname, self.port = split_URL_port(hostname)
        if self.port is None:
            self.port = POSTGRES_DEFAULT_PORT

        with open(DATABASE_CREDENTIALS_FILENAME) as database_credentials_file:
            self.username, self.password = [line.strip() for line in database_credentials_file.readlines()][:2]

        self.connection = psycopg2.connect(database=self.database_name, user=self.username, password=self.password, host=self.hostname,
                                           port=self.port)

    def test_table_creation(self):
        table_name = 'test_create_table'
        fields = list(FIELD_TYPES.keys())

        with self.connection:
            with self.connection.cursor() as cursor:
                if database_has_table(cursor, table_name):
                    cursor.execute(f'DROP TABLE "{table_name}";')

        metadata_table = MetadataDatabase(f'{self.hostname}:{self.port}', 'metadata', table_name, fields, primary_key=fields[0])

        with self.connection:
            with self.connection.cursor() as cursor:
                assert database_has_table(cursor, table_name)
                cursor.execute(f'DROP TABLE "{table_name}";')


class TestMetadataFile(unittest.TestCase):
    def test_file_creation(self):
        fields = list(FIELD_TYPES.keys())

        with TemporaryDirectory() as temporary_directory:
            csv_filename = os.path.join(temporary_directory, 'metadata.csv')
            assert not os.path.exists(csv_filename)
            metadata_table = MetadataFile(csv_filename, fields, primary_key=fields[0])
            assert os.path.exists(csv_filename)


if __name__ == '__main__':
    unittest.main()
