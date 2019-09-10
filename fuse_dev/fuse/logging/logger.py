# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:21:05 2019

@author: Casiano.Koprowski
"""

import datetime as _dt
import logging as _logging
import os as _os
import sys as _sys

class Log:
    """
    TODO: Write Description

    """

    def __init__(self, log_output_path: str, log_filename: str, name: str = 'fuse', incl_date: bool = True):

        date_time = ''
        if incl_date:
            today = _dt.datetime.now()
            date_time = f'{today:%Y%m%d_%H%M}_'
        self.logger = _logging.getLogger(name)
#        metapath, metafile = _os.path.split(log_output_path)
        logname = _os.path.join(log_output_path, f'{date_time}{log_filename}.log')
        self.name = logname
        # remove handlers that might have existed from previous files
        for h in self.logger.handlers:
            self.logger.removeHandler(h)
         # create file handler for this filename
        fh = _logging.FileHandler(logname)
        formatter = _logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.handler = fh.setFormatter(formatter)
        self.logger.addHandler(fh)

    def close(self):
        """
        Close the object logging file by removing log handler(s)

        """
        for h in self.logger.handlers:
            self.logger.removeHandler(h)
        del self.logger

    def debug(self, message):
        with LoggingContext(self.logger, level=_logging.DEBUG):
            self.logger.debug(message)

    def error(self, message):
        with LoggingContext(self.logger, level=_logging.ERROR):
            self.logger.error(message)

    def info(self, message):
        with LoggingContext(self.logger, level=_logging.INFO):
            self.logger.info(message)

    def warning(self, message):
        with LoggingContext(self.logger, level=_logging.WARNING):
            self.logger.warning(message)


class LoggingProcessor:
    """
    TODO: Write Description

    """

    loggers = {}

    def __init__(self, log_output_path: str, log_filename: str, handler_name: str = 'fuse'):

        self.log_output_path = log_output_path
        self.process_handler = handler_name
        self.loggers[handler_name] = Log(log_output_path, log_filename, name=handler_name)

    def add_logger(self, log_filename: str, log_output_path: str = None, handler_name: str = None, incl_date: bool = False):

        if handler_name is None:
            handler_name = log_filename

        if log_output_path is None:
            log_output_path = self.log_output_path

        self.loggers[handler_name] = Log(log_output_path, log_filename, handler_name, incl_date=incl_date)

    def close_logger(self, log_key):

        self.loggers[log_key].close()
        del self.loggers[log_key]

    def close_all_loggers(self):

        logs = self._get_log_list()
        for log in logs:
            self.loggers[log].close()
            del self.loggers[log]

    def _get_log_list(self):
        return list(self.loggers.keys())

#    def debug(self, log, message):
#        self.loggers[log].debug(message)
#
#    def error(self, log, message):
#        self.loggers[self.process_handler].error(message)
#        self.loggers[log].debug(message)
#
#    def info(self, log, message):
#        self.loggers[log].info(message)
#
#    def warning(self, log, message):
#        self.loggers[self.process_handler].warning(message)
#        self.loggers[log].warning(message)


class LoggingContext():
    """
    https://docs.python.org/3/howto/logging-cookbook.html#using-a-context-manager-for-selective-logging
    """

    def __init__(self, logger, level=None, handler=None, close=True):
        self.logger = logger
        self.level = level
        self.handler = handler
        self.close = close

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)
        if self.handler:
            self.logger.addHandler(self.handler)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        if self.handler:
            self.logger.removeHandler(self.handler)
        if self.handler and self.close:
            self.handler.close()
        # implicit return of None => don't swallow exceptions
