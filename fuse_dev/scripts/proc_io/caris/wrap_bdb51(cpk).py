# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:11:10 2019

@author: Casiano.Koprowski
"""

import socket
import sys
import time


class bdb51_io:
    """

    """

    def __init__(self, port):
        """

        """

        self.sock = None
        self.port = port
        self.host = 'localhost'
        self.alive = True
        self.connected = False
        self.node_manager = None
        self.database = None
        self._nm = None
        self._db = None
        self._command = []
        self._response = []
        self.bmsg = None
        self.value = 0

        self._connect(port, self.host)

    def _connect(self, port, host='localhost'):
        self.sock = None

        for res in socket.getaddrinfo(host, port, socket.AF_UNSPEC, socket.SOCK_STREAM):
            af, socktype, proto, canonname, sa = res

            try:
                self.sock = socket.socket(af, socktype, proto)
            except OSError as msg:
                self.sock = None
                continue

            try:
                self.sock.connect(sa)
                self.connected = True
            except OSError as msg:
                self.sock.close()
                self.sock = None
                continue

            break
        if self.sock is None:
            print('could not open socket')
            sys.exit(1)

        with self.sock:
            self.sock.sendall(b'Alive')
            #            try:

            while True:
                data = self.sock.recv(1024)
                print('Received {}'.format(repr(data)))

                if data:
                    self.bmsg = self._command_data(data)

                print(self.bmsg)
                self.sock.sendall(self.bmsg)

    #            except:
    #                pass

    def _command_data(self, data):
        if data == b'domath':
            bmsg = self._domath()
        if data == b'die':
            bmsg = self._die()
        else:
            bmsg = b'Alive'

        print(bmsg)
        return bmsg

    def _domath(self):
        x = 0

        while True:
            print(x)
            if x > 10:
                break
            x += 1

        return b'done'

    def _die(self):
        x = 10

        while True:
            print(x)
            time.sleep(1)
            if x == 0:
                sys.exit(1)
                break
            x -= 1

        return b'done'


def main(port):
    """
    An event loop waiting for commands.
    """

    db_io = bdb51_io(port)
    print(db_io.alive)


if __name__ == '__main__':
    args = sys.argv
    main(int(args[-1]))
