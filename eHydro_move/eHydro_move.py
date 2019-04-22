# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 17:19:27 2019

@author: Casiano.Koprowski
"""

import os as _os
import shutil as _shutil
import configparser as _cp

"""Known global constants"""
progLoc = _os.getcwd()
"""progLoc is the program's own file location / current working directory (cwd)
obtained by :func:`os.getcwd()`"""
config = _cp.ConfigParser(interpolation = _cp.ExtendedInterpolation())
config.read('config.ini')

#sections = config.sections()
downloads = config['Source']['downloads']
destination = config['Destination']['destination']
folder = config['Destination']['folder']
repo = _os.path.split(progLoc)[0]

def regionPath(root, folder):
    """
    """
    regions = []
    for k, v in dict(config.items('Regions')).items():
        dwnlds = []
        districts = [i.strip() for i in v.split(',')]
        for i in districts:
            dwnlds.append(repo + _os.path.join(downloads, i))
        regions.append((k, dwnlds))
    regions = dict(regions)
    return regions

def fileCollect(path):
    """
    """
    zips = []
    if _os.path.exists(path):
        for root, folders, files in _os.walk(path):
            for file in files:
                zips.append(_os.path.join(root, file))
    if len(zips) > 0:
        return zips
    else:
        zips.append('None')
        return zips

def eHydroZIPs(regions):
    """
    """
    hold = []
    for k, v in regions.items():
        zips = []
        for folder in v:
            zips.extend(fileCollect(folder))
        num = len(zips)
        hold.append((k,(zips,num)))
    hold = dict(hold)
    return hold


def fileMove(regionFiles, destination, folder, text_region=None,
             progressBar=None, text_output=None):
    """
    """
    for k, v in regionFiles.items():
        if text_region != None and progressBar!= None and text_output != None:
            text_region.SetValue(k)
            progressBar.SetRange(v[1])
            progressBar.SetValue(0)
            pbv = 0
        newPath = _os.path.join(destination, k, folder)
        if _os.path.isdir(newPath):
            pass
        else:
            _os.makedirs(newPath)
        for item in v[0]:
            name = _os.path.split(item)[-1]
            newName = _os.path.join(newPath, name)
            if _os.path.exists(item) and not _os.path.exists(newName):
                if text_output != None:
                    text_output.write(name+ '\n')
                _shutil.copy2(item, newName)
            if progressBar!= None:
                pbv += 1
                progressBar.SetValue(pbv)

def _main(text_region=None, progressBar=None, text_output=None):
    """
    """
    if progressBar!= None:
        progressBar.Pulse()
    regions = regionPath(repo, downloads)
    regionFiles = eHydroZIPs(regions)
    fileMove(regionFiles, destination, folder, text_region, progressBar,
             text_output)
