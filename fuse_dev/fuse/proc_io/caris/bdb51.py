# -*- coding: utf-8 -*-
"""
BDB51.py

Created on Thu May  9 15:46:49 2019

@author: grice

Classes and methods for working with CARIS Bathy DataBASE 5.1.
"""

import subprocess
import logging
import os, sys
import pickle
from fuse.proc_io.caris import helper

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
    def __init__(self, database_loc, database_name, caris_env_name = 'NBS35'):
        """
        Instantiate the object and connect to the database referenced, waiting
        until the database responds.
        """
        
        self.caris_environment_name = caris_env_name
        self._logger = logging.getLogger('fuse')
        if len(self._logger.handlers) == 0:
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(logging.DEBUG)
            self._logger.addHandler(ch)
        self._start_env()
        
    def _start_env(self):
        """
        Instantiate an environment containing CARIS BDB51 object for talking
        to the database.
        """
        # set up the path to python in the CARIS environment
        conda_env_name = self.caris_environment_name
        conda_env_path = helper.retrieve_env_path(conda_env_name)
        python_path = os.path.join(conda_env_path, 'python')
        self._port = 65505
        # set the location for running the database i/o object
        start = os.path.realpath(os.path.dirname(__file__))
        db_obj = os.path.join(start, 'wrap_bdb51.py')
        activate_file = helper.retrieve_activate_batch()
        args = ["cmd.exe", "/C", "set pythonpath= &&", # setup the commandline
                activate_file, conda_env_name, "&&",  # activate the Caris 3.5 virtual environment
                python_path, db_obj, str(self._port), # call the script for the object
                ]
        args = ' '.join(args)
        self._logger.log(logging.DEBUG, args)
        try:
            self.db = subprocess.Popen(
                    args,
                    )
        except:
            err = 'Error executing: {}'.foramt(args)
            print(err)
            self._logger.log(logging.DEBUG, err)
        try:
            pass
#            if len(out) > 0:
#                msg = out.decode(encoding='UTF-8')
#                print(msg)
#                self._logger.log(logging.DEBUG, msg)
#            if len(err) > 0:
#                msg = err.decode(encoding='UTF-8')
#                print(msg)
#                self._logger.log(logging.DEBUG, msg)
        except Exception as e:
            err = 'Error in handling error output: {}'.format(e)
            print(err)
            self._logger.log(logging.DEBUG, err)
    
    def connect(self):
        """
        Form and send the connect command to the BDB51 wapper.
        """
        command = {'command':'connect'}
        self._send_command(command)
    
    def status(self):
        """
        Check to see if the subprocess is still communicating and connected to
        the database.
        """
        command = {'command':'status'}
        self._send_command(command)
        
    def upload(self, dataset, instruction):
        """
        Send the BDB environment instructions on where data is and what to do
        with it.
        """
        pass
    
    def die(self, delay = 0):
        """
        Destroy the BDB51 wrapper object and environment.
        """
        command = {'command':'die'}
        if delay == 0:
            command['action'] = 'now'
        else:
            command['action'] = int(delay)
        self._send_command(command)
    
    def _send_command(self, command):
        """
        Send a command to the BDB51 wrapper.
        """
        pc = pickle.dumps(command)
        self.db.stdin.write(pc)
        response = self.db.stdout.read()
        print(response.decode(encoding='UTF-8'))
        
        #self._logger.log(logging.DEBUG, pickle.loads(response))
        #self._logger.log(logging.DEBUG, pickle.loads(err))
        