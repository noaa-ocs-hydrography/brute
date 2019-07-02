# -*- coding: utf-8 -*-
"""
proc_io.py

Created on Thu Feb 14 15:13:52 2019

@author: grice
"""

import logging
import os
import pickle
import subprocess
import sys
from tempfile import TemporaryDirectory as tempdir
from typing import Tuple, List

import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
from fuse.proc_io import caris

__version__ = 'Test'


class proc_io:
    """
    A class to abstract the reading and writing of bathymetry.
    """

    def __init__(self, in_data_type: str, out_data_type: str, work_dir: str = None, z_up: bool = True,
                 nodata: float = 1000000.0, caris_env_name: str = 'CARIS35', overwrite: bool = True):
        """
        Initialize with the data type to be worked.

        Parameters
        ----------
        in_data_type : str
        out_data_type : str
        work_dir : str, optional
            A string representing an optional folder path to store working
            data. Default is ``None`` working data is stored in system temp
            data locations
        z_up : bool, optional
            Boolean determining the vertical orienation of the data. Default is
            ``True``, the data is marked as elevation, depth is positive
        nodata : float, optional
            Default is ``1000000.0``. The no data value for the input data.
            This is the value used by CSAR and BAG files for no data
        caris_env_name : str, optional
            Default is ``NBS35``.  Determines which caris environment is used
        overwrite : bool, optional
            Default is ``True``. If a file with an existing name is input, this
            will determine whether the file is overwritten or kept
        """

        self._in_data_type = in_data_type
        self._out_data_type = out_data_type
        self._z_up = z_up
        self._write_nodata = nodata
        self._caris_environment_name = caris_env_name
        self.overwrite = overwrite

        if work_dir is None:
            self._work_dir = tempdir()
            self._work_dir_name = self._work_dir.name
        else:
            self._work_dir_name = work_dir

        self._logger = logging.getLogger('fuse')

        if len(self._logger.handlers) == 0:
            ch = logging.StreamHandler(sys.stdout)
            ch.setLevel(logging.DEBUG)
            self._logger.addHandler(ch)

        if self._out_data_type == "carisbdb51":
            self._bdb51 = caris.bdb51()

    def write(self, dataset: gdal.Dataset, instruction: str, metadata: dict = None):
        """
        Write the provided data to the predefined data type.

        :param dataset: 
        :param instruction: 
        :param metadata:  (Default value = None)
        """

        self._logger.log(logging.DEBUG, f'Begin {self._out_data_type} write')

        if os.path.exists(instruction) and self.overwrite:
            self._logger.log(logging.DEBUG, f'Overwriting {instruction}')
            os.remove(instruction)
            if self._out_data_type == 'bag':
                caris_xml = instruction + '.aux.xml'
                if os.path.exists(caris_xml):
                    os.remove(caris_xml)

        if self._out_data_type == 'csar':
            self._write_csar(dataset, instruction)
        elif self._out_data_type == 'bag':
            self._write_bag(dataset, instruction, metadata)
        elif self._out_data_type == 'gpkg':
            self._write_points(dataset, instruction)
        elif self._out_data_type == 'carisbdb51':
            self._send2bdb51(dataset, instruction)
        else:
            raise ValueError(f'writer type unknown: {self._out_data_type}')

    def _write_csar(self, dataset: gdal.Dataset, outfilename: str):
        """
        Convert the provided gdal dataset into a csar file.
        
        The data and metadata are saved out to a file and then loaded into the
        wrapper around the csar writer.

        :param dataset: 
        :param outfilename: 
        """

        conda_env_name = self._caris_environment_name

        # put the provided data into the right form for the csar conversion.
        if self._in_data_type == 'gdal':
            dataset = self._set_gdalndv(dataset)
            data, metadata = self._gdal2array(dataset)
            metadata['outfilename'] = outfilename
        elif self._in_data_type == 'point':
            data, metadata = self._point2array(dataset)
            splits = os.path.splitext(outfilename)
            metadata['outfilename'] = f'{splits[0]}_Points{splits[1]}'
            print(metadata['outfilename'], outfilename)
        else:
            raise ValueError(f'input data type unknown: {self._in_data_type}')
        metadata['z_up'] = self._z_up
        conda_env_path = caris.helper.retrieve_env_path(conda_env_name)
        python_path = os.path.join(conda_env_path, 'python')
        # save the provided dataset and metadata to a file
        datafilename = os.path.join(self._work_dir_name, 'rasterdata.npy')
        np.save(datafilename, data)
        metafilename = os.path.join(self._work_dir_name, 'metadata')

        with open(metafilename, 'wb') as metafile:
            pickle.dump(metadata, metafile)

        # set the locations for running the wrap_csar script
        start = os.path.realpath(os.path.dirname(__file__))
        write_csar = os.path.join(start, 'caris', 'wrap_csar.py')
        activate_file = caris.helper.retrieve_activate_batch()

        if os.path.exists(write_csar):
            args = ["cmd.exe", "/K", "set pythonpath= &&",  # setup the commandline
                    activate_file, conda_env_name, "&&",  # activate the Caris 3.5 virtual environment
                    python_path, write_csar,  # call the script
                    f'"{datafilename.replace("&", "^&")}"',  # surface path
                    f'"{metafilename.replace("&", "^&")}"',  # metadata path
                    f'"{self._in_data_type.replace("&", "^&")}"',  # data type
                    ]
            args = ' '.join(args)
            print(args)
            self._logger.log(logging.DEBUG, args)

            try:
                proc = subprocess.Popen(args, creationflags=subprocess.CREATE_NEW_CONSOLE)
            except:
                err = 'Error executing: {}'.foramt(args)
                print(err)
                self._logger.log(logging.DEBUG, err)

            try:
                stdout, stderr = proc.communicate()
                self._logger.log(logging.DEBUG, stdout)
                self._logger.log(logging.DEBUG, stderr)
            except:
                err = 'Error in handling error output'
                print(err)
                self._logger.log(logging.DEBUG, err)

            if not os.path.exists(metadata['outfilename']):
                err = f"Unable to create {metadata['outfilename']}"
                self._logger.log(logging.DEBUG, err)
                raise RuntimeError(err)
        else:
            err = f"Unable to overwrite {metadata['outfilename']}"
            self._logger.log(logging.DEBUG, err)
            raise RuntimeError(err)

    def _write_bag(self, dataset: gdal.Dataset, outfilename: str, metadata=None):
        """
        Convert the provided gdal dataset into a bag file.

        :param dataset: 
        :param outfilename: 
        :param metadata:  (Default value = None)
        """

        if self._in_data_type == 'gdal':
            dataset = self._set_gdalndv(dataset)
        else:
            raise ValueError(f'input data type unknown: {self._in_data_type}')

        if metadata is not None:
            raise NotImplementedError('bag xml metadata write has not been implemented')

        print(dataset.GetGeoTransform())

        # Prepare destination file
        driver = gdal.GetDriverByName("BAG")

        # write and close output raster dataset
        dest = driver.CreateCopy(outfilename, dataset)
        dest = None
        self._logger.log(logging.DEBUG, 'BAG file created')

    def _write_points(self, dataset: gdal.Dataset, outfilename: str):
        """
        Convert the provided gdal dataset into a geopackage file.

        :param dataset: 
        :param outfilename: 
        """

        points, meta = self._point2wkt(dataset)
        crs = meta['crs']
        proj = osr.SpatialReference(wkt=crs)

        splits = os.path.splitext(outfilename)
        outfilename = f'{splits[0]}_Points.gpkg'

        if os.path.exists(outfilename):
            os.remove(outfilename)

        driver = ogr.GetDriverByName('GPKG')
        ds = driver.CreateDataSource(outfilename)
        layer = ds.CreateLayer('Point', proj, ogr.wkbMultiPoint)

        layer.CreateField(ogr.FieldDefn('X', ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn('Y', ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn('Z', ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn('Elevation', ogr.OFTReal))
        layer.CreateField(ogr.FieldDefn('Feature Type', ogr.OFTString))
        defn = layer.GetLayerDefn()

        for point in points:
            # Create a new feature (attribute and geometry)
            feat = ogr.Feature(defn)
            feat.SetField('X', point['x'])
            feat.SetField('Y', point['y'])
            feat.SetField('Z', point['z'])
            feat.SetField('Elevation', point['z'])
            feat.SetField('Feature Type', 'Point')

            # Make a geometry, from Shapely object
            geom = ogr.CreateGeometryFromWkt(point['wkt'])
            feat.SetGeometry(geom)

            layer.CreateFeature(feat)
            feat = geom = None  # destroy these

        # Save and close everything
        ds = layer = feat = geom = None

    def _write_vector(self, dataset: gdal.Dataset, outfilename: str):
        """
        TODO write description

        :param dataset: 
        :param outfilename: 
        """

        splits = os.path.split(outfilename)[1]
        name = os.path.splitext(outfilename)[0]
        outfilename = os.path.join(splits[0], f'{name}_Vector.gpkg')

        proj = dataset.GetProjection()
        proj = osr.SpatialReference(wkt=proj)
        band = dataset.GetRasterBand(1)

        driver = ogr.GetDriverByName('GPKG')
        ds = driver.CreateDataSource(outfilename)
        layer = ds.CreateLayer(name, proj, ogr.wkbMultiPolygon)

        # Add one attribute
        layer.CreateField(ogr.FieldDefn('Survey', ogr.OFTString))
        defn = layer.GetLayerDefn()

        # Create a new feature (attribute and geometry)
        feat = ogr.Feature(defn)
        feat.SetField('Survey', name)

        gdal.Polygonize(band, None, layer, 0, [],
                        callback=None)

        ds = band = None

    def _gdal2array(self, dataset: gdal.Dataset) -> np.array:
        """
        Convert the gdal dataset into a numpy array and a dictionary of
        metadata of the geotransform information and return.
        
        The gdal dataset should have he no data value set appropriately.

        :param dataset:
        :returns: array
        """

        meta = {}

        # get the logisitics for converting the gdal dataset to csar
        gt = dataset.GetGeoTransform()
        print(gt)
        meta['resx'] = gt[1]
        meta['resy'] = gt[5]
        meta['originx'] = gt[0]
        meta['originy'] = gt[3]
        meta['dimx'] = dataset.RasterXSize
        meta['dimy'] = dataset.RasterYSize
        print(meta)
        meta['crs'] = dataset.GetProjection()
        rb = dataset.GetRasterBand(1)  # should this be hardcoded for 1?
        meta['nodata'] = rb.GetNoDataValue()

        # get the gdal data raster
        data = rb.ReadAsArray()
        return data, meta

    def _point2array(self, dataset: gdal.Dataset) -> np.array:
        """
        Convert the gdal dataset into a numpy array and a dictionary of
        metadata of the geotransform information and return.
        
        The gdal dataset should have he no data value set appropriately.

        :param dataset:
        :returns: array
        """

        meta = {}

        lyr = dataset.GetLayerByIndex(0)
        crs = lyr.GetSpatialRef().ExportToWkt()

        # Read out of the GDAL data structure
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()
        data = np.zeros((count, 3))
        for n in np.arange(count):
            f = lyr.GetFeature(n)
            data[n, :] = f.geometry().GetPoint()
        data[:, 2] *= -1  # make these heights rather than depths

        meta['crs'] = crs
        print(meta)

        return data, meta

    def _point2wkt(self, dataset: gdal.Dataset) -> Tuple[List[dict], dict]:
        """
        Convert the gdal dataset into a WKT Points object and a dictionary of
        metadata of the geotransform information and return.
        
        The gdal dataset should have he no data value set appropriately.

        :param dataset:
        :returns: points and metadata
        """

        meta = {}
        points = []
        lyr = dataset.GetLayerByIndex(0)
        crs = lyr.GetSpatialRef().ExportToWkt()
        #        multipoint = ogr.Geometry(ogr.wkbMultiPoint)

        # Read out of the GDAL data structure
        lyr = dataset.GetLayerByIndex(0)
        count = lyr.GetFeatureCount()

        for n in np.arange(count):
            info = {}
            point = ogr.Geometry(ogr.wkbPoint)
            f = lyr.GetFeature(n)
            x, y, z = f.geometry().GetPoint()
            point.AddPoint(x, y, z)
            info['x'] = x
            info['y'] = y
            info['z'] = z
            info['wkt'] = point.ExportToWkt()
            points.append(info)
        #            multipoint.AddGeometry(point)
        #        info['wkt'] = multipoint.ExportToWkt()1
        #        points.append(info)

        meta['crs'] = crs
        print(meta)

        return points, meta

    def _set_gdalndv(self, dataset: gdal.Dataset) -> gdal.Dataset:
        """
        Update the gdal raster object no data value and the raster no data
        values in to corrispond with the object no data value.

        :param dataset: 
        """

        # check the no data value
        rb = dataset.GetRasterBand(1)  # should this be hardcoded for 1?
        ndv = rb.GetNoDataValue()

        if self._write_nodata != ndv:
            data = rb.ReadAsArray()
            data = np.where(data == ndv, self._write_nodata, data)
            dataset.GetRasterBand(1).WriteArray(data)

        return dataset
