# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:09:59 2019

@author: Casiano.Koprowski
"""

import os
import sys
import time
import pickle
import socket
import logging
import subprocess

from fuse.proc_io.caris import helper


class bdb51:

    def __init__(self):
        self.host = 'localhost'
        self.port = None
        self.sock = None
        self.conn = None
        self.addr = None
        self.data = None
        self.bmsg = None

    def _start_env(self, port):
        """
        Instantiate an environment containing CARIS BDB51 object for talking
        to the database.
        """
        # set up the path to python in the CARIS environment
        conda_env_path = helper.retrieve_env_path('CARIS35')
        python_path = os.path.join(conda_env_path, 'python')
        # set the location for running the database i/o object
        start = os.path.realpath(os.path.dirname(__file__))
        db_obj = os.path.join(start, 'wrap_bdb51(cpk).py')
        activate_file = helper.retrieve_activate_batch()
        args = ["cmd.exe", "/k", "set pythonpath= &&",  # setup the commandline
                activate_file, conda_env_path, "&&",  # activate the Caris 3.5 virtual environment
                python_path, db_obj, str(port),  # call the script for the object
                ]
        args = ' '.join(args)
        print(args)
        try:
            db = subprocess.Popen(
                    args,
                    creationflags=subprocess.CREATE_NEW_CONSOLE,
                    close_fds=True)
        except:
            err = 'Error executing: {}'.foramt(args)
            print(err)

    def form_connection(self, host='localhost', port=0):
        """
        Open a socket, start the environment and then pass along when CARIS
        env sends a response.
        """
        self.sock = socket.socket()
        self.sock.bind((host, port))
        self.port = self.sock.getsockname()[1]
        self.sock.close()
        for res in socket.getaddrinfo(host, self.port, socket.AF_UNSPEC,
                              socket.SOCK_STREAM, 0, socket.AI_PASSIVE):
            af, socktype, proto, canonname, sa = res
            try:
                self.sock = socket.socket(af, socktype, proto)
            except OSError as msg:
                self.sock = None
                continue
            try:
                self.sock.bind(sa)
                self.sock.listen(1)
            except OSError as msg:
                self.sock.close()
                self.sock = None
                continue
            break
        if self.sock is None:
            print('could not open socket')
            sys.exit(1)

        self._start_env(self.port)
        self.conn, self.addr = self.sock.accept()
        self.bmsg = b'Alive'
        with self.conn:
            while True:
                if len(self.bmsg) > 0:
                    print('Connected by', self.addr)
                    self.data = self.conn.recv(1024)
                    print('Received', repr(self.data))
                    self.conn.send(self.bmsg)
                    self.bmsg = b''
                else:
                    pass

    def domath(self):
        self.data = self._is_alive()
        if self.data == b'Alive':
            self.conn.send(b'domath')
        self.data = self.conn.recv(1024)
        print('Received', repr(self.data))
        if self.data == b'done':
            print('done')

    def die(self):
        self.data = self._is_alive()
        if self.data == b'Alive':
            self.conn.send(b'die')
        self.data = self.conn.recv(1024)
        print('Received', repr(data))
        if self.data == b'done':
            print('dead')

    def _is_alive(self):
#        self.conn, self.addr = self.sock.accept()
        with self.conn:
            while True:
                print('Connected by', self.addr)
                self.data = self.conn.recv(1024)
                print('Received', repr(self.data))
        return self.data

#test = bdb51()
#test.domath()
#test.die()

