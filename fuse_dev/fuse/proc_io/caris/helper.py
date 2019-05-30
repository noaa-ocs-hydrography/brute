# -*- coding: utf-8 -*-
"""
helper.py

Created on Wed May 15 15:13:07 2019

@author: grice

A set of small utilities for working with CARIS and the CARIS environment.
"""
import os
import sys
# helper function to retrieve the path to the "Scripts" folder in PydroXL
def retrieve_scripts_folder():
    install_prefix = sys.exec_prefix
    folder_path = os.path.realpath(os.path.join(install_prefix, os.pardir, os.pardir, "Scripts"))
    if not os.path.exists(folder_path):
        raise RuntimeError("The Scripts folder does not exist at: %s" % folder_path)
    return folder_path

# helper function to retrieve the path to the "activate.bat" batch file in PydroXL
def retrieve_activate_batch():

    scripts_prefix = retrieve_scripts_folder()
    file_path = os.path.realpath(os.path.join(scripts_prefix, "activate.bat"))
    if not os.path.exists(file_path):
        raise RuntimeError("The activate file does not exist at: %s" % file_path)
    return file_path

def retrieve_env_path(env_name):
    """
    Given a conda environement name, find the environment.
    """
    current_env_loc = os.environ['conda_prefix']
    desired_env_loc = os.path.join(current_env_loc, os.pardir, env_name)
    if os.path.exists(desired_env_loc):
        return desired_env_loc
    else:
        raise RuntimeError('{} environment does not exist in current conda installation'.format(env_name))
        