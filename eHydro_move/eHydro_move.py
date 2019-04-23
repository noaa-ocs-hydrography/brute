# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 17:19:27 2019

@author: Casiano.Koprowski
"""

import os as _os
import csv as _csv
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
repo = _os.path.split(progLoc)[0]
method = config.getboolean('Method', 'method')

def regionPath(root, folder):
    """Uses repo path derived from the program's location to find the downloads
    folder of eHydro_scrape.

    This function relies on the associated ``config.ini`` having defined
    regions in the **[Regions]** section.  It uses the region's included
    districts to compile a list of the district download folders to be
    associated with the region.

    A region defined as ``NorthEast = CENAN`` would include an eHydro_scrape
    dowload folder like ``root\\eHydro_scrape\\downloads\\CENAN``

    Parameters
    ----------
    root : str
        A string representing the base path of the repository
    folder : str
        A string representing the downloads folder of eHydro_scrape

    Returns
    -------
    regions : dict
        A dictionary with keys of region names and values of a list of district
        downloads folders for that region

    """
    regions = []
#    for k, v in dict(config.items('Regions')).items():
#        dwnlds = []
#        districts = [i.strip() for i in v.split(',')]
#        for i in districts:
#            dwnlds.append(repo + _os.path.join(downloads, i))
#        regions.append((k, dwnlds))
#    regions = dict(regions)
    fileName = open(config['CSVs']['NBS'], 'r')
    opened = _csv.reader(fileName,delimiter=',')
    temp = []
    for row in opened:
        if len(row) > 0:
            temp.append(row)
    fileName.close()
    for row in temp[1:]:
        regionName = row[1].strip() + '_' + row[0].strip()
        dwnlds = []
        districts = [i.strip() for i in row[2].split(',')]
        for i in districts:
             dwnlds.append(repo + _os.path.join(downloads, i))
        regions.append((regionName, dwnlds))
    regions = dict(regions)
    return regions

def fileCollect(path):
    """Given a folder path, this function will return the complete list of
    files in that folder.

    Because this tool is specifically designed to move eHydro_scrape downloads,
    it is expected that all file paths collected will be ``.zip`` files.  If
    there are non-``.zip`` files in an eHydro_scrape download folder, the
    function will not raise an error, but it will collect the file path for use
    in copying the file.

    Parameters
    ----------
    path : str
        A string representing a folder path

    Returns
    -------
    zips : list
        A list of all files found in the provided file path, returns an empty
        list if none are found

    """
    zips = []
    if _os.path.exists(path):
        for root, folders, files in _os.walk(path):
            for file in files:
                zips.append(_os.path.join(root, file))
    if len(zips) > 0:
        return zips
    else:
        zips.append(None)
        return zips

def eHydroZIPs(regions):
    """Creates a dictionary with keys of region names and values of a list of
        all files downloaded from districts within the respective regions

    For each region provided in the dictionary passed to this function, it will
    use the eHydro_scrape download folders associated with it and compile one
    complete list of all the files associated with the districts within the
    defined region.  This list is then reassociated with the region names as a
    dictionary and returned.

    Parameters
    ----------
    regions : dict
        A dictionary with keys of region names and values of a list of district
        downloads folders for that region

    Returns
    -------
    hold : dict
        A dictionary with keys of region names and values of a list of
        all files downloaded from districts within the respective regions

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


def fileMove(regionFiles, destination, method, text_region=None,
             progressBar=None, text_output=None):
    """Takes a dictionary of keys representing regions and values representing
    a list of all files downloaded from districts within the respective regions
    and attempts to copy the files to the new directory structure.

    Uses paramerters desination and folder along with keys representing region
    names to extrapolate the complete path where the region-associated files
    will be moved to.

    Given a destination path such as ``..\\vlab-nbs\\New_Directory\\``, the
    default folder of ``USACE\\eHydro\\original\\`` and a region defined as
    ``NorthEast = CENAN``, the function would produce a destination folder path
    like ``..\\vlab-nbs\\New_Directory\\northeast\\USACE\\eHydro\\original\\``.
    The eHydro_scrape download files associated for each region would be copied
    to folders like this.  If the folder does not exist, it will be created
    using :func:`os.makedirs` and the folder path produced.

    Parameters
    ----------
    regions : dict
        A dictionary with keys of region names and values of a list of district
        downloads folders for that region
    destination : str
        A string representing the base path of the destination directory
    folder : str
        A string representing the sub-directory where the copied files will be
        placed within a named region
    text_region : wx.TextCtrl, optional
        Optional text field used when run via the included GUI; displays the
        current region being copied
    progressBar : wx.Guage, optional
        Optional guage used when run via the included GUI; displays the copy
        progress of the current region
    text_output: wx.TextCtrl, optional
        Optional text field used when run via the included GUI; displays the
        names of files copied

    """
    fileName = open(config['CSVs']['USACE'], 'r')
    opened = _csv.reader(fileName,delimiter=',')
    district_name = []
    for row in opened:
        if len(row) > 0:
            district_name.append(row[:2])
    fileName.close()
    district_name = dict(district_name[1:])
    for k, v in regionFiles.items():
        if text_region != None:
            text_region.SetValue(k)
        if text_output != None:
            progressBar.SetRange(v[1])
        if progressBar != None:
            progressBar.SetValue(0)
            pbv = 0
        for item in v[0]:
            if item != None:
                splits = _os.path.split(item)
                name = splits[-1]
                district_code = splits[-2].split('\\')[-1]
#                print (k, district_code, item)
                district_abbr = district_code[-3:]
                district_full = district_name[district_abbr]
                eHydro_folder = 'USACE\\eHydro_' + district_full + '_district\\original'
                newerPath = _os.path.join(destination, k, eHydro_folder)
                if _os.path.isdir(newerPath):
                    pass
                else:
                    _os.makedirs(newerPath)
                newName = _os.path.join(newerPath, name)
                if _os.path.exists(item) and not _os.path.exists(newName):
                    if text_output != None:
                        text_output.write(eHydro_folder +'\\'+ name + '\n')
                    if method == False:
                        _shutil.copy2(item, newName)
                    elif method == True:
                        _shutil.move(item, newName)
                if progressBar!= None:
                    pbv += 1
                    progressBar.SetValue(pbv)

def _main(text_region=None, progressBar=None, text_output=None):
    """Main function of the program. Responsible for passing items through the
    intended use of this tool.

    Parameters
    ----------
    text_region : wx.TextCtrl, optional
        Optional text field used when run via the included GUI; displays the
        current region being copied
    progressBar : wx.Guage, optional
        Optional guage used when run via the included GUI; displays the copy
        progress of the current region
    text_output: wx.TextCtrl, optional
        Optional text field used when run via the included GUI; displays the
        names of files copied

    """
    if progressBar!= None:
        progressBar.Pulse()
    regions = regionPath(repo, downloads)
    regionFiles = eHydroZIPs(regions)
    fileMove(regionFiles, destination, method, text_region,
             progressBar, text_output)

if __name__ == '__main__':
    _main()
