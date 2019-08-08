# -*- coding: utf-8 -*-
"""
_bdb51.py

grice 20190509

These are the tools for interfacting directly with CARIS Bathy DataBASE 5.1
server from within the CARIS conda python environment.
"""

import pickle
import socket
import sys
import time

import caris
import caris.bathy.db as bdb
from get_access import *


class BDB51IO:
    """
    An object for handling CARIS BDB I/O through the Python interface.
    """

    def __init__(self, port, buffer_size):
        """
        Initialize the object.

        Parameters
        ----------
        port
        """

        self.sock = None
        self.port = port
        self.host = 'localhost'
        self.buffer_size = buffer_size
        self.alive = True
        self.connected = False
        self.node_manager = None
        self.database = None
        self._nm = None
        self._db = None
        self._command = []
        self._response = []
        # open the socket
        self._call_home(self.host, self.port)
        # listen for a command
        self.get_commands()

    def _call_home(self, host, port):
        """
        Connect back to the calling environment.
        """
        # set up the socket
        try:
            self.sock = socket.create_connection((host, port))
        except OSError:
            self.sock.close()
            self.sock = None
        # send a message saying we are alive
        if self.sock is None:
            print('could not call home from subprocess')
            sys.exit(1)
        else:
            hello_im_here = self.status({'command': None})
            self.sock.sendall(pickle.dumps(hello_im_here))

    def get_commands(self):
        """
        Listen on the provided socket and wait for something to do.
        """
        while self.alive:
            try:
                data = self.sock.recv(self.buffer_size)
                if data:
                    command = pickle.loads(data)
                    response = self.take_commands(command)
                    print('Response', response, end = '\n\n')
                    data = pickle.dumps(response)
                    # print('pickled data size is {}, buffer size is {}'.format(len(data), self.buffer_size))
                    lensent = self.sock.send(data)
                    # print('sent {} bytes'.format(str(lensent)))
                    if lensent != len(data):
                        print('Failed to send all data. {} sent.').format(lensent)
            except ConnectionResetError:
                print('Remote connection closed.  Assuming die command.')
                self.die({})
        self.sock.close()

    def take_commands(self, command_dict: dict) -> dict:
        """
        Act on commands the provided command dictionary, such as to upload
        data, read and return data, destroy the object and exit the
        environment.

        Parameters
        ----------
        command_dict :
            returns: response
        command_dict: dict :


        Returns
        -------
        type
            response

        """

        # self._command.append(command_dict)
        command = command_dict['command']

        if command == 'connect':
            return self.connect(command_dict)
        elif command == 'status':
            return self.status(command_dict)
        elif command == 'upload':
            return self.upload(command_dict)
        elif command == 'die':
            return self.die(command_dict)

        # self._response.append(response)

    def connect(self, command_dict: dict) -> dict:
        """
        Make a connection to the BDB database using the provided location and
        name. If there is not a predefined password or username one will be
        requested.

        Parameters
        ----------
        command_dict :
            returns: response
        command_dict: dict :


        Returns
        -------
        type
            response

        """

        msg = ''
        self.node_manager = command_dict['node_manager']
        self.database = command_dict['database']

        try:
            if None in (username, password):
                raise ValueError("System environment variable 'nbsscriptuser' and/or 'nbsscriptpass' is not defined")
            else:
                self._nm = bdb.NodeManager(username, password, self.node_manager)
                msg += 'Connected to Node Manager {}\n'.format(self.node_manager)
        except RuntimeError as error:
            msg += str(error)
            command_dict['success'] = False
        except ValueError as error:
            msg += str(error)
            command_dict['success'] = False
        if self._nm is not None:
            try:
                self._db = self._nm.get_database(self.database)
                msg += ', Connected to database {}'.format(self.database)
                self.connected = True
                command_dict['success'] = True
            except RuntimeError as error:
                msg += str(error)
                command_dict['success'] = False

        command_dict['log'] = msg
        return command_dict

    def _check_connection(self):
        """
        Check to see if the database connection is alive.
        """

        pass

    def status(self, command_dict: dict) -> dict:
        """
        Check the status of the object.

        Parameters
        ----------
        command_dict :
            returns: response
        command_dict: dict :


        Returns
        -------
        type
            response

        """

        command_dict['success'] = True
        command_dict['alive'] = self.alive
        command_dict['connected'] = self.connected
        return command_dict

    def upload(self, command_dict: dict) -> dict:
        """
        Send data at the provided location, specifying if the metadata or
        bathymetery should be updated if it already exists.

        Parameters
        ----------
        command_dict :
            returns: response
        command_dict: dict :


        Returns
        -------
        type
            response

        """
        try:
            # what to upload, new or updated data
            action = command_dict['action']
            # the name of the file to get data from
            bathy_path = command_dict['bathy_path']
            with open(command_dict['meta_path'], 'rb') as f:
                metadata = pickle.load(f)
                print(metadata)
            if action == 'new':
                msg = self._upload_new(bathy_path, metadata)
            elif action == 'bathy':
                msg = 'Updating only bathy is not implemented yet'
                # query for the object and replace the bathy
            elif action == 'metadata':
                msg = self._update_metadata(metadata)
            else:
                raise ValueError('Upload action type not understood')

            command_dict['success'] = True
            command_dict['log'] = msg
        except Exception as error:
            command_dict['success'] = False
            command_dict['log'] = str(error)

        return command_dict

    def _upload_new(self, file_path: str, new_metadata: dict) -> str:
        """
        Upload both bathymetry and the metadata.

        Parameters
        ----------
        file_path :

        file_path: str :


        Returns
        -------

        """
        # check to see if file with same object name already exists
        n = -1
        cql = "OBJNAM = '{}'".format(new_metadata['OBJNAM'])
        features = self._db.query('surfac', cql)
        for n,f in enumerate(features):
            pass
        if n > -1:
            up_msg = self._update_metadata(new_metadata)
            objnam = new_metadata['OBJNAM']
            msg = '{} already found in database, updating metadata only...'.format(objnam)
            return msg + up_msg
        else:
            # create a fake feature
            crs = self._db.crs
            fake_coverage = 'POLYGON((0 0,0 1,1 1,1 0,0 0))'
            geom = caris.Geometry(crs, fake_coverage)
            surface = self._db.create_feature('surfac', geom)
            surface['OBJNAM'] = file_path
            surface['srcfil'] = file_path
    
            # get a metadata container to put stuff into
            current_metadata = surface.attributes
            # need to load the metadata dictionary that was put on disk here.
            for key in new_metadata:
                current_metadata[key] = new_metadata[key]
            # commit the feature to the database
            self._db.commit()
    
            # upload coverage
            surface.upload_coverage(file_path)
            return 'Uploaded {} to {}'.format(file_path, self.database)
    
    def _update_metadata(self, new_metadata: dict) -> str:
        """
        Query for the data object and update the metadata with the provided
        dictionary.
        """
        n = -1
        cql = "OBJNAM = '{}'".format(new_metadata['OBJNAM'])
        print(cql)
        features = self._db.query('surfac', cql)
        for n,f in enumerate(features):
            if n > 0:
                msg = "More than one object found in query with provided key!"
                break
            else:
                current_metadata = f.attributes
                for key in new_metadata:
                    current_metadata[key] = new_metadata[key]
                self._db.commit()
                msg = "FID found: {}".format(str(f.id))
        if n == -1:
            msg = "No object with the provided name found."
        return msg

    def die(self, command_dict: dict) -> dict:
        """
        Update the object "alive" flag.

        The command dictionary must contain a key word "action"

        Parameters
        ----------
        command_dict :
            returns: response
        command_dict: dict :


        Returns
        -------
        type
            response

        """
        if 'action' in command_dict:
            time.sleep(int(command_dict['action']))
        self._nm = None
        self._db = None
        command_dict['success'] = True
        command_dict['log'] = 'Stopping I/O with {}'.format(self.database)
        self.alive = False
        return command_dict


def main(port: int, buffer_size: int):
    """
    An event loop waiting for commands.

    Parameters
    ----------
    port :

    port: int :

    buffer_size: int : receive buffersize


    Returns
    -------

    """

    db_io = BDB51IO(port, buffer_size)


if __name__ == '__main__':
    args = sys.argv
    main(int(args[-2]), int(args[-1]))
