# -*- coding: utf-8 -*-
"""
_bdb51.py

grice 20190509

These are the tools for interfacting directly with CARIS Bathy DataBASE 5.1
server from within the CARIS conda python environment.
"""

import os, sys
import pickle
import time
import socket

import caris.bathy.db as bdb
import caris

from get_access import *

class bdb51_io:
    """
    
    """
    def __init__(self, port):
        """
        
        """
        self.port = port
        self.alive = True
        self.connected = False
        self.node_manager = None
        self.database = None
        self._nm = None
        self._db = None
        self._command = []
        self._response = []
        pass
    
    def take_commands(self):
        """
        Call the provided port on the host and ask for something to do.
        """
        conn = socket.create_connection(('',self.port))
        conn.setblocking(True)
        # need to build a packet to send to open the connection
        hello_im_here = self.status({})
        conn.send(hello_im_here)
        while self.alive:
            command = conn.recv()
            response = self.take_commands(command)
            conn.send(response)
        conn.close()
    
    def connect(self, command_dict):
        """
        Make a connection to the BDB database using the provided location and
        name.  If there is not a predefined password or username one will be
        requested.
        """
        msg = ''
        self.node_manager = command_dict['node_manager']
        self.database = command_dict['database']
        try:
            self._nm = bdb.NodeManager(username, password, self.node_manager)
            msg = msg + 'Connected to Node Manager {}'.format(self.node_manager)
        except RuntimeError as e:
            sys.stderr(pickle.dumps(e))
        if self._nm is not None:
            try:
                self._db = self._nm.get_database(self.database)
                msg = msg + ', Connected to database {}'.format(self.database)
                self.connected = True
            except RuntimeError as e:
                sys.stderr(pickle.dumps(e))
        response = {'command':'connect','response':'success', 'log':msg}
        return response
    
    def _check_connection(self):
        """
        Check to see if the database connection is alive.
        """
        pass
    
    def status(self, command_dict):
        """
        Check the status of the object.
        """
        response = {'command':'status','alive':self.alive}
        return response
    
    def upload(self, command_dict):
        """
        Send data at the provided location, specifying if the metadata or
        bathymetery should be updated if it already exists.
        """
        # what to upload, new or updated data
        action = command_dict['action']
        # the name of the file to get data from
        name = command_dict['name']
        if action == 'new':
            msg = self._upload_new(name)
        response = {'command':'upload','response':'success', 'log':msg}
        return response
    
    def _upload_new(self, file_path):
        """
        Upload both bathymetry and the metadata.
        """
        # create a fake feature
        crs = self._db.crs
        fake_coverage = 'POLYGON((0 0,0 1,1 1,1 0,0 0))'
        geom = caris.Geometry(crs, fake_coverage)
        surface = self._db.create_feature('surfac', geom)
        surface['OBJNAM'] = file_path    
        surface['srcfil'] = file_path
        # get a metadata container to put stuff into
#        metadata = surface.attributes
#        # need to load the metadata dictionary that was put on disk here.
#        metafilename = get the name here
#        with pickle.load(metafilename) as new_meta:
#         try:
#             metadata = new_meta
#         except:
#             surface.attribute['OBJNAM']  = 'MetaDataFail'
#             with open('metadata_error_file.txt','a') as metafail:
#                 metafail.write(file_path + '\n')
        # commit the feature to the database
        self._db.commit()
        # upload coverage
        try:
            surface.upload_coverage(file_path)
        except RuntimeError as e:
            sys.stderr(pickle.dumps(e))
        
    def query(self, command_dict):
        """
        Query for and return data from the db.
        """
        pass

    def die(self, command_dict):
        """
        Update the object "alive" flag.
        
        The command dictionary must contain a key word "action"
        """
        action = command_dict['action']
        time.sleep(int(action))
        self._nm = None
        self._db = None
        response = {'response':'log', 'message':'Stopping I/O with {}'.format(self.database)}
        self.alive = False
        return response
    
    def take_commands(self):
        """
        Act on commands the provided command dictionary, such as to upload
        data, read and return data, destroy the object and exit the 
        environment.
        """
        self._command.append(command_dict)
        command = command_dict['command']
        if command == 'connect':
            response = self.connect(command_dict)
        elif command == 'status':
            response = self.status(command_dict)
        elif command == 'upload':
            response = self.upload(command_dict)
        elif command == 'die':
            response = self.die(command_dict)
        self._response.append(response)
        return response
    
def main(port):
    """
    An event loop waiting for commands.
    """
    db_io = bdb51_io(port)
    db_io.request_commands()        

if __name__ == '__main__':
    args = sys.argv
    main(int(args[-1]))