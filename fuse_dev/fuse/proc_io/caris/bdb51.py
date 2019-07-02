# -*- coding: utf-8 -*-
"""
BDB51.py

Created on Thu May  9 15:46:49 2019

@author: grice

Classes and methods for working with CARIS Bathy DataBASE 5.1.
"""

import logging
import os
import pickle
import socket
import subprocess
import sys
import threading
from typing import Dict

from fuse.proc_io.caris import helper
from osgeo import gdal


class bdb51:
    """
    A class for doing I/O with the CARIS Bathy DataBASE 5.1 server.

    This class spawns a conda environment containing the CARIS python API and
    communicates with that environment over subprocess pipes.

    The pipes carry pickled python dictionaries.

    General commands that are sent to the CARIS environment include:
        status : Is the system is connected to the server, what was the last
            file uploaded.
        connect : Connect the database.
        upload : Send instruction on surveys to upload and what to upload.
        die : Self destruct the _BDB51 object in the CARIs environment.

    Replies from the CARIS environment can be
        status : connectivity to the database and the last file uploaded.
        log : returning a log message
    """

    def __init__(self, database_loc: str, database_name: str, caris_env_name: str = 'NBS35',
                 host='localhost'):
        """
        Instantiate the object and connect to the database referenced, waiting
        until the database responds.

        Parameters
        ----------
        database_loc
        database_name
        caris_env_name
        """

        self.host = host
        self.port = None
        self.sock = None
        self.conn = None
        self.data = None
        self.alive = False  # is the sub environment talking via TCP
        self.connected = False  # is the sub environment talking to the db
        self._response = None
        self._command = None
        self.msg_id = 0
        self.database_loc = database_loc
        self.database_name = database_name
        self.caris_environment_name = caris_env_name
        self._logger = logging.getLogger('fuse')

        if len(self._logger.handlers) == 0:
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(logging.DEBUG)
            self._logger.addHandler(ch)

        self._thread = threading.Thread(target=self._form_connection)
        self._thread.start()

    def _form_connection(self):
        """
        Open a socket, start the environment and then pass along when CARIS
        env sends a response.

        Parameters
        ----------

        Returns
        -------

        """

        self.sock = socket.socket()
        self.sock.bind((self.host, 0))
        self.port = self.sock.getsockname()[1]
        self.sock.listen(1)
        self._start_env(self.port)

        while True:
            self._conn, addr = self.sock.accept()
            self._logger.log(logging.DEBUG, 'accepted connection from {} at {}'.format(addr, self._conn))
            data = self._conn.recv(1024)

            if len(data) > 0:
                response = pickle.loads(data)

                if response['success']:
                    self.alive = response['alive']
            else:
                print('No response from subprocess received')

            while self.alive:
                ''' 
                TO DO: should we should try to detect if the connection is 
                broken and look for another connection if it is
                '''
                if self._command is not None:
                    try:
                        self._conn.send(pickle.dumps(self._command))
                        self._command = None
                        data = self._conn.recv(1024)
                        if len(data) > 0:
                            response = pickle.loads(data)
                            self._response = response
                    except ConnectionResetError:
                        self.alive = False
                        break
                    except Exception as error:
                        self._logger.log(logging.DEBUG, str(error))
            if not self.alive:
                break

    def _start_env(self, port: int):
        """
        Instantiate an environment containing CARIS BDB51 object for talking
        to the database.

        Parameters
        ----------
        port :

        port: int :


        Returns
        -------

        """

        # set up the path to python in the CARIS environment
        conda_env_name = self.caris_environment_name
        conda_env_path = helper.retrieve_env_path(conda_env_name)
        python_path = os.path.join(conda_env_path, 'python')

        # set the location for running the database i/o object
        start = os.path.realpath(os.path.dirname(__file__))
        db_obj = os.path.join(start, 'wrap_bdb51.py')
        activate_file = helper.retrieve_activate_batch()

        args = ["cmd.exe", "/K", "set pythonpath= &&",  # setup the commandline
                activate_file, conda_env_name, "&&",  # activate the Caris 3.5 virtual environment
                python_path, db_obj, str(port),  # call the script for the object
                ]
        args = ' '.join(args)
        self._logger.log(logging.DEBUG, args)

        try:
            self.db = subprocess.Popen(args, creationflags=subprocess.CREATE_NEW_CONSOLE)
        except:
            err = 'Error executing: {}'.format(args)
            print(err)
            self._logger.log(logging.DEBUG, err)

        try:
            # if len(out) > 0:
            #     msg = out.decode(encoding='UTF-8')
            #     print(msg)
            #     self._logger.log(logging.DEBUG, msg)
            # if len(err) > 0:
            #     msg = err.decode(encoding='UTF-8')
            #     print(msg)
            #     self._logger.log(logging.DEBUG, msg)
            pass
        except Exception as e:
            err = 'Error in handling error output: {}'.format(e)
            print(err)
            self._logger.log(logging.DEBUG, err)

    def connect(self):
        """
        Form and send the connect command to the BDB51 wapper.
        """

        command = {'command': 'connect', 'node_manager': self.database_loc, 'database': self.database_name}
        response = self._set_command(command)

    def status(self):
        """
        Check to see if the subprocess is still communicating and connected to
        the database.

        Parameters
        ----------

        Returns
        -------

        """

        command = {'command': 'status'}
        response = self._set_command(command)

    def upload(self, dataset: gdal.Dataset, instruction: str):
        """
        Send the BDB environment instructions on where data is and what to do
        with it.

        Parameters
        ----------
        dataset :
            param instruction:
        dataset: gdal.Dataset :

        instruction: str :


        Returns
        -------

        """

        pass

    def die(self, delay: int = 0):
        """
        Destroy the BDB51 wrapper object and environment.

        Parameters
        ----------
        delay :
            Default value = 0)
        delay: int :
             (Default value = 0)

        Returns
        -------

        """

        command = {'command': 'die', 'action': int(delay)}
        response = self._set_command(command)

        if response['command'] == 'die' and response['success']:
            self.alive = False

    def _set_command(self, command: Dict[str, str]):
        """
        Set the object command variable and wait for a response from
        the subprocess.
        """

        command['id'] = self.msg_id
        self.msg_id += 1

        if self._command is None and self._response is None:
            self._command = command

            while True:
                if self._response is not None:  # we need a way to check if the connection is alive
                    response = self._response
                    self._logger.log(logging.DEBUG, str(response))

                    if not response['success']:
                        print('{} failed!'.format(response['command']))

                    self._response = None
                    break
        else:
            raise ValueError('command / response state is unexpected')

        return response
