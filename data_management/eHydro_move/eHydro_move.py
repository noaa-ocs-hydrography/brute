# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 17:19:27 2019

@author: Casiano.Koprowski
"""

import configparser as _cp
import csv as _csv
import os as _os
import pickle as _pickle
import re as _re
import shutil as _shutil
import zipfile as _zf
from typing import Union

from osgeo import ogr as _ogr
from osgeo import osr as _osr

"""Known global constants"""
progLoc = _os.path.dirname(_os.path.abspath(__file__))
"""progLoc is the program's own file location / current working directory (cwd)
obtained by :func:`os.path.dirname(os.path.abspath(__file__))`"""

zreg = _re.compile(r'.zip', _re.IGNORECASE)
config = _cp.ConfigParser(interpolation=_cp.ExtendedInterpolation())
config.read('config.ini')

# sections = config.sections()
downloads = config['Source']['downloads']
destination = config['Destination']['destination']
if config['USACE Reference']['Path'] != '':
    reference_files = config['USACE Reference']['location']
else:
    reference_files = None
repo = _os.path.dirname(progLoc)
method = config.getboolean('Method', 'method')


def regionPath(repo: str, downloads: str) -> dict:
    """
    Uses repo path derived from the program's location to find the downloads
    folder of eHydro_scrape.

    This function uses the the NBS region definitions file defined in the
    provided ``config.ini`` under the header **[CSVs]** and name 'NBS ='.
    Each region within the definitions file is made up of the region name,
    processing branch, list of USACE districts included, and the relative
    location of a shapefile containing the geographic boundaries of the region

    Each region's information is assigned to a dictionary with the
    [Region]_[ProcessingBranch] as the key and the subseqent districts and
    shapefile location are assigned as the value as a list

    Parameters
    ----------
    repo : str
        A string representing the base path of the repository
    downloads : str
        A string representing the downloads folder of eHydro_scrape

    Returns
    -------
    type
        A dictionary with keys of region names and values of a list of district
        downloads folders and relative path for the shapfile for that region

    """

    regions = []

    default_downloads = repo + r'\eHydro_scrape\downloads'

    with open(config['CSVs']['NBS'], 'r') as fileName:
        temp = [row for row in _csv.reader(fileName, delimiter=',') if len(row) > 0]

    for row in temp[1:]:
        regionName = f'{row[1].strip()}_{row[0].strip()}'
        dwnlds = []
        districts = [i.strip() for i in row[2].split(',')]


        for i in districts:
            if downloads != '' and _os.path.isdir(downloads):
                dwnlds.append(_os.path.join(downloads, i))
            elif downloads == '' and _os.path.isdir(default_downloads):
                dwnlds.append(_os.path.join(default_downloads, i))
            else:
                raise RuntimeError(f'Unable to find location of downloaded files. If the files are not stored in {default_downloads}, please ensure the Source\downloads field is set in the config.ini or that {default_downloads} exists')

        bounds = row[-1]
        regions.append((regionName, [dwnlds, bounds]))

    regions = dict(regions)
    return regions


def open_ogr(path) -> tuple:
    """
    Parameters
    ----------
    path :

    Returns
    -------

    """

    ds = _ogr.Open(path)
    ds_layer = ds.GetLayer()

    for feature in ds_layer:
        if feature is not None:
            geom = feature.GetGeometryRef()
            ds_geom = _ogr.CreateGeometryFromWkt(geom.ExportToWkt())
            break

    ds_proj = ds_layer.GetSpatialRef()

    if ds_proj.IsProjected:
        proj_name = ds_proj.GetAttrValue('projcs')
    else:
        proj_name = ds_proj.GetAttrValue('geogcs')

    return ds_geom, ds_proj  # , proj_name


def write_shapefile(out_shp: str, name: str, geom, spcs: Union[str, int]):
    """
    Writes out a geopackage shapefile containing the bounding geometry of
    of the given query.

    Derived from:
    https://gis.stackexchange.com/a/52708/8104
    via
    https://gis.stackexchange.com/q/217165

    Parameters
    ----------
    out_shp
        String representing the complete file path for the output shapefile
    name
        String representing the name of the survey; Used to name the layer
    geom

    spcs

    """

    # Reference
    if type(spcs) == str:
        proj = _osr.SpatialReference(wkt=spcs)
        proj.MorphFromESRI()
    else:
        proj = _osr.SpatialReference()
        proj.ImportFromEPSG(spcs)

    # Now convert it to a shapefile with OGR
    driver = _ogr.GetDriverByName('GPKG')
    ds = driver.CreateDataSource(out_shp)
    layer = ds.CreateLayer(name, proj, _ogr.wkbPolygon)
    #    layer = ds.CreateLayer(name, None, ogr.wkbMultiPolygon)
    # Add one attribute
    layer.CreateField(_ogr.FieldDefn('id', _ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    # If there are multiple geometries, put the "for" loop here

    # Create a new feature (attribute and geometry)
    feat = _ogr.Feature(defn)
    feat.SetField('id', 123)

    # Make a geometry, from Shapely object
    #    geom = _ogr.CreateGeometryFromWkt(poly)
    feat.SetGeometry(geom)

    layer.CreateFeature(feat)

    # Save and close everything
    del ds, layer, feat, geom


def fileCollect(path: str, bounds: str) -> list:
    """
    Given a folder path, this function will return the complete list of
    files in that folder.

    Because this tool is specifically designed to move eHydro_scrape downloads,
    it is expected that all file paths collected will be ``.zip`` files.  If
    there are non-``.zip`` files in an eHydro_scrape download folder, the
    function will not raise an error, but it will collect the file path for use
    in copying the file.

    Parameters
    ----------
    path :
        A string representing a folder path
    bounds :
        returns: A list of all files found in the provided file path, returns
        an empty list if none are found

    Returns
    -------
    type
        A list of all files found in the provided file path, returns an empty
        list if none are found

    """

    zips = []
    bfile = _os.path.join(repo, bounds)
    bname = _os.path.splitext(bounds.split('\\')[1])[0]
    print(f'\n\n{bname}')
    meta_geom, meta_proj = open_ogr(bfile)

    if _os.path.exists(path):
        for root, folders, files in _os.walk(path):
            for item in files:
                if zreg.search(item):
                    zips.append(_os.path.join(root, item))

    slen = len(zips)
    print(f'downloaded zips: {slen}')
    x = 1
    verified_zips = []
    for zfile in zips:
        root = _os.path.dirname(zfile)
        _os.chdir(root)

        try:
            zipped = _zf.ZipFile(zfile)
            contents = zipped.namelist()
            for name in contents:
                if _re.compile(r'\.gpkg$', _re.IGNORECASE).search(name):
                    path = _os.path.join(root, name)
#                    print(x, path)
                    zipped.extract(name)

                    try:
                        ehyd_geom, ehyd_proj = open_ogr(path)
                        try:
                            intersection = meta_geom.Intersection(ehyd_geom)
                            flag = intersection.ExportToWkt()
                        except AttributeError as e:
                            flag = 'GEOMETRYCOLLECTION EMPTY'
                            print(e, bfile, path, meta_geom, ehyd_geom,
                                  sep='\n')
                    except TypeError as e:
                        print(e, meta_proj, ehyd_proj, sep='\n')
                        flag = 'GEOMETRYCOLLECTION EMPTY'

                    if flag != 'GEOMETRYCOLLECTION EMPTY':
#                        print('They did Intersect')
                        verified_zips.append(zfile)
#                    else:
#                        print('They did not Intersect')

                    _os.remove(path)
        except _zf.BadZipfile:
            print('BadZip')

        zipped.close()
        _os.chdir(progLoc)
        x += 1
    print(f'verified zips: {len(verified_zips)}')

    if len(verified_zips) > 0:
        return verified_zips
    else:
        verified_zips.append(None)
        return verified_zips


def eHydroZIPs(regions: {str: [str]}) -> {str: [str]}:
    """
    Creates a dictionary with keys of region names and values of a list of all
    files downloaded from districts within the respective regions

    For each region provided in the dictionary passed to this function, it will
    use the eHydro_scrape download folders associated with it and compile one
    complete list of all the files associated with the districts within the
    defined region.  This list is then reassociated with the region names as a
    dictionary and returned.

    Parameters
    ----------
    regions :
        A dictionary with keys of region names and values of a list of district
        downloads folders for that region

    Returns
    -------
    type
        A dictionary with keys of region names and values of a list of all
        files downloaded from districts within the respective regions

    """

    hold = []

    for region, downloads in regions.items():
        zips = []

        for meta in downloads[0]:
            zips.extend(fileCollect(meta, downloads[1]))

        num = len(zips)
        hold.append((region, (zips, num)))

    hold = dict(hold)
    return hold


def contentSearch(filenames: [str], extensions: [str]) -> [str]:
    """
    This funtion takes a list of zipfile contents.

    Using the zipfile contents, it parses the files for any file containing the
    full string '.xyz', '.xml', or '.pickle'.  If a file name contains this
    string, it is added to a list of found files.

    Parameters
    ----------
    filenames
        A list of file names
    extensions
        list of extensions to filter by

    Returns
    -------
    [str]
        a list of files that met the correct conditions
    """

    return [filename for filename in filenames if
            any(_re.compile(rf'{extension}$', _re.IGNORECASE).search(filename) is not None for extension in extensions)]


def zipManipulate(path: str, name: str):
    """
    TODO write description

    Parameters
    ----------
    path : str

    name : str


    Returns
    -------

    """

    _os.chdir(path)

    try:
        zipped = _zf.ZipFile(name)
        contents = zipped.namelist()

        for item in contentSearch(contents, ['xyz', 'xml', 'pickle', 'gpkg']):
            zipped.extract(item)

        zipped.close()

        if reference_files is not None:
            for item in _os.listdir(path):
                if 'pickle' in item:
                    with open(item, 'rb') as pickle_file:
                        pickle_metadata = _pickle.load(pickle_file)
                        district = pickle_metadata['SURVEYAGENCY']
                        area = pickle_metadata['CHANNELAREAIDFK']
                        area_file = f"{pickle_metadata['CHANNELAREAIDFK']}.gpkg"
                        ref_outline = _os.path.join(reference_files, district, area, area_file)
                        pickle_metadata['poly_name'] = ref_outline

        _os.remove(name)
    except _zf.BadZipfile:
        _os.remove(name)

    _os.chdir(progLoc)


def fileMove(regionFiles: {str: [str]}, destination: str, method,
             text_region=None, progressBar=None, text_output=None):
    """
    Takes a dictionary of keys representing regions and values representing
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
    regionFiles :
        A dictionary with keys of region names and values of a list of district
        downloads folders for that region
    destination :
        A string representing the base path of the destination directory
    method
    text_region: wx.TextCtrl :
         (Default value = None) Optional text field used when run via the
         included GUI; displays the current region being copied
    progressBar: wx.Guage :
         (Default value = None) Optional guage used when run via the included
         GUI; displays the copy progress of the current region
    text_output: wx.TextCtrl :
         (Default value = None) Optional text field used when run via the
         included GUI; displays the names of files copied

    """

    with open(config['CSVs']['USACE'], 'r') as fileName:
        district_name = dict([row[:2] for row in _csv.reader(fileName, delimiter=',') if len(row) > 0][1:])

    for region, downloads in regionFiles.items():
        if text_region is not None:
            text_region.SetValue(region)

        if progressBar is not None:
            progressBar.SetRange(downloads[1])
            progressBar.SetValue(0)
            pbv = 0

        for item in downloads[0]:
            if item is not None:
                splits = _os.path.split(item)
                name = splits[-1]
                surname = _os.path.splitext(name)[0]
                district_code = splits[-2].split('\\')[-1]
                district_abbr = district_code[-3:]
                district_full = f'{district_name[district_abbr]}_{district_code}'
                eHydro_folder = _os.path.join('USACE', f'eHydro_{district_full}', 'Original')
                newerPath = _os.path.join(destination, region, eHydro_folder, surname)

                if _os.path.isdir(newerPath):
                    pass
                else:
                    _os.makedirs(newerPath)

                newName = _os.path.join(newerPath, name)

                if _os.path.exists(item):
                    if text_output is not None:
                        text_output.write(f'{_os.path.join(eHydro_folder, name)}\n')

                    if method:
                        _shutil.move(item, newName)
                    elif not method:
                        _shutil.copy2(item, newName)

                    zipManipulate(newerPath, newName)

                if progressBar is not None:
                    pbv += 1
                    progressBar.SetValue(pbv)


def _main(text_region=None, progressBar=None, text_output=None):
    """
    Main function of the program. Responsible for passing items through the
    intended use of this tool.

    Parameters
    ----------
    text_region
         (Default value = None) Optional text field used when run via the
         included GUI; displays the current region being copied
    progressBar
         (Default value = None) Optional guage used when run via the included
         GUI; displays the copy progress of the current region
    text_output
         (Default value = None) Optional text field used when run via the
         included GUI; displays the names of files copied
    """

    if progressBar is not None:
        progressBar.Pulse()

    regions = regionPath(repo, downloads)
    regionFiles = eHydroZIPs(regions)
    fileMove(regionFiles, destination, method, text_region, progressBar, text_output)


if __name__ == '__main__':
    _main()
