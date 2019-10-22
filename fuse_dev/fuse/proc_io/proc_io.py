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

import fiona
import fiona.crs
import numpy as np
from fuse.utilities import vectorize_raster, gdal_to_xyz
from osgeo import gdal, osr
from shapely.geometry import MultiPoint

gdal.UseExceptions()
from fuse.proc_io import caris

__version__ = 'Test'


class ProcIO:
    """ A class to abstract the reading and writing of bathymetry data. """

    def __init__(self, in_data_type: str, out_data_type: str, work_dir: str = None, z_up: bool = True, nodata: float = 1000000.0,
                 caris_env_name: str = 'CARIS35', overwrite: bool = True, db_loc: str = None, db_name: str = None):
        """
        Create a new bathymetric data input / output object, reading and writing in the given input and output formats.

        Parameters
        ----------
        in_data_type : str
            input data type
        out_data_type : str
            output data type (one of 'csar', 'bag', 'gpkg', 'carisbdb51')
        work_dir : str, optional
            folder path to store working data; if none is given, a temporary directory will be created
        z_up : bool, optional
            whether positive Z data is elevation (`True`) or depth (`False`)
        nodata : float, optional
            value in input data representing no data; the default nodata value in CSAR and BAG files is 1000000.0
        caris_env_name : str, optional
            name of the Python environment with access to CARIS (Python <= 3.5)
        overwrite : bool, optional
            whether to overwrite existing files
        db_loc : str, optional
            location of bathymetric database
        db_name : str, optional
            name of bathymetric database
        """

        self._in_data_type = in_data_type
        self._out_data_type = out_data_type
        self._z_up = z_up
        self._nodata = nodata
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
            if db_name is not None and db_loc is not None:
                self._bdb = caris.bdb5(db_loc, db_name)
                self._bdb.status()
                self._bdb.connect()
            else:
                raise ValueError('No database name or location provided')

    def write(self, dataset: gdal.Dataset, filename: str, metadata: dict = None, show_console: bool = False):
        """
        Write the provided data to the predefined data type.

        Parameters
        ----------
        dataset
            GDAL dataset
        filename
            filename to write
        metadata
            dictioanry of metadata
        show_console
            whether to show the console of the CARIS subprocess (if writing to a CSAR file)
        """

        self._logger.log(logging.DEBUG, f'Begin {self._out_data_type} write to {filename}')

        if os.path.exists(filename) and self.overwrite:
            self._logger.log(logging.DEBUG, f'Overwriting {filename}')
            os.remove(filename)
            if self._out_data_type == 'csar':
                csar_data = f'{filename}0'
                if os.path.exists(csar_data):
                    os.remove(csar_data)
            elif self._out_data_type == 'bag':
                bag_xml = f'{filename}.aux.xml'
                if os.path.exists(bag_xml):
                    os.remove(bag_xml)


        if self._out_data_type == 'csar':
            self._write_csar(dataset, filename, show_console=show_console)
        elif self._out_data_type == 'bag':
            self._write_bag(dataset, filename, metadata)
        elif self._out_data_type == 'gpkg':
            self._write_points(dataset, filename)
        elif self._out_data_type == 'carisbdb51':
            self._bdb.upload(dataset, filename, metadata)
        else:
            raise ValueError(f'writer type unknown: {self._out_data_type}')

    def close_connection(self):
        """ Closes connection to the bathymetric database """

        if self._out_data_type == 'carisbdb51':
            self._bdb.die()

    def _write_csar(self, raster: gdal.Dataset, filename: str, show_console: bool = False):
        """
        Write the provided GDAL raster dataset to a CSAR file by calling a Python 3.5 environment with access to CARIS.
        The data and metadata are written to a file, then loaded in the CARIS environment.

        Parameters
        ----------
        raster
            GDAL raster dataset
        filename
            filename to write CSAR file
        show_console
            whether to show the console of the CARIS subprocess
        """

        # put the provided data into the right form for the CSAR conversion
        if self._in_data_type == 'gdal':
            raster = self._set_raster_nodata(raster)
            data, metadata = self._gdal_raster_to_array(raster)
        elif self._in_data_type == 'point':
            data, metadata = self._gdal_points_to_array(raster)
        else:
            raise ValueError(f'input data type unknown: {self._in_data_type}')
        metadata['outfilename'] = filename
        metadata['z_up'] = self._z_up
        conda_env_path = caris.helper.retrieve_env_path(self._caris_environment_name)
        python_path = os.path.join(conda_env_path, 'python')

        # save the provided dataset and metadata to a file
        raster_data_filename = os.path.join(self._work_dir_name, 'rasterdata.npy')
        np.save(raster_data_filename, data)
        metadata_filename = os.path.join(self._work_dir_name, 'metadata')
        with open(metadata_filename, 'wb') as metafile:
            pickle.dump(metadata, metafile)

        # set the locations for running `wrap_csar.py`
        write_csar = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'caris', 'wrap_csar.py')
        activate_file = caris.helper.retrieve_activate_batch()

        if os.path.exists(write_csar):
            argument_string = ' '.join((
                # setup the commandline
                'cmd.exe', '/C', 'set pythonpath= &&',
                # activate the Python 3.5 environment with access to CARIS
                activate_file, self._caris_environment_name, '&&',
                # call the script
                python_path, write_csar,
                # filename of saved raster data
                f'"{raster_data_filename.replace("&", "^&")}"',
                # filename of saved metadata
                f'"{metadata_filename.replace("&", "^&")}"',
                # data type
                f'"{self._in_data_type.replace("&", "^&")}"'
            ))
            self._logger.log(logging.DEBUG, argument_string)

            try:
                caris_process = subprocess.Popen(argument_string, creationflags=subprocess.CREATE_NEW_CONSOLE if show_console else 0)
            except Exception as error:
                self._logger.log(logging.DEBUG, f'Error when executing "{argument_string}"\n{error}')

            try:
                stdout, stderr = caris_process.communicate()
                if stdout is not None:
                    self._logger.log(logging.DEBUG, stdout)
                else:
                    self._logger.log(logging.DEBUG, 'No information returned from csar write process.')
                if stderr is not None:
                    self._logger.log(logging.DEBUG, stderr)
                else:
                    self._logger.log(logging.DEBUG, 'No errors returned from csar write process.')
            except Exception as error:
                self._logger.log(logging.DEBUG, f'Error when logging subprocess error\n{error}')

            if not os.path.exists(metadata['outfilename']):
                error_string = f'Unable to create file {metadata["outfilename"]}'
                self._logger.log(logging.DEBUG, error_string)
                raise RuntimeError(error_string)
        else:
            error_string = f'Unable to overwrite file {metadata["outfilename"]}'
            self._logger.log(logging.DEBUG, error_string)
            raise RuntimeError(error_string)

    def _write_bag(self, raster: gdal.Dataset, filename: str, metadata: dict = None):
        """
        Write the given GDAL raster dataset to a BAG file.

        Parameters
        ----------
        raster
            GDAL raster dataset
        filename
            filename of BAG
        metadata
            dictionary of metadata
        """

        if self._in_data_type == 'gdal':
            raster = self._set_raster_nodata(raster)
        else:
            raise NotImplementedError(f'BAG writing not supported for "{self._in_data_type}"')

        if metadata is not None:
            raise NotImplementedError('Writing XML metadata to a BAG has not yet been implemented in GDAL.')

        # Prepare destination file
        bag_driver = gdal.GetDriverByName("BAG")

        # write and close output raster dataset
        try:
            output_raster = bag_driver.CreateCopy(filename, raster)
            del output_raster
            self._logger.log(logging.DEBUG, 'BAG file created')
        except RuntimeError as error:
            self._logger.log(logging.DEBUG, f'Failed to create bag: {error}\n Attempting to pass correct well-known text using OSR')
            old_crs_wkt = raster.GetProjectionRef()
            self._logger.log(logging.DEBUG, f'Old reference: {old_crs_wkt}')
            spatial_reference = osr.SpatialReference(wkt=old_crs_wkt)
            new_crs_wkt = spatial_reference.ExportToWkt()
            self._logger.log(logging.DEBUG, f'New reference: {new_crs_wkt}')
            raster.SetProjection(new_crs_wkt)
            output_raster = bag_driver.CreateCopy(filename, raster)
            del output_raster
            self._logger.log(logging.DEBUG, f'Created BAG file at {filename}')

    def _write_points(self, points: gdal.Dataset, filename: str, layer: int = 0, output_layer: str = 'Elevation'):
        """
        Write provided GDAL point cloud dataset to a geopackage file.

        Parameters
        ----------
        points
            GDAL point cloud dataset
        filename
            filename to write geopackage
        layer
            index of vector layer containing points
        output_layer
            name of layer to create
        """

        filename = f'{os.path.splitext(filename)[0]}_Points.gpkg'
        if os.path.exists(filename):
            os.remove(filename)

        point_layer = points.GetLayerByIndex(layer)
        crs_wkt = point_layer.GetSpatialRef().ExportToWkt()
        del point_layer

        points = gdal_to_xyz(points, layer)
        projection = fiona.crs.from_string(crs_wkt)

        layer_schema = {
            'geometry': 'Point',
            'properties': {
                'X': 'float',
                'Y': 'float',
                'Z': 'float',
                'Elevation': 'float',
                'Feature Type': 'str'
            }
        }

        # in fiona features are input as dictionary "records"
        point_records = [{
            'geometry': {
                'type': 'Point',
                'coordinates': tuple(point)
            },
            'properties': {
                'X': point[0],
                'Y': point[1],
                'Z': point[2],
                'Elevation': point[2],
                'Feature Type': 'Point'
            }
        } for point in points]

        with fiona.open(filename, 'w', 'GPKG', schema=layer_schema, crs=projection, layer=output_layer) as output_file:
            output_file.writerecords(point_records)

    def _write_vectorized_raster(self, filename: str, raster: gdal.Dataset, crs_wkt: str, band_index: int = 1, layer: str = 'Elevation'):
        """
        Write the given GDAL raster dataset to a GDAL vector dataset (vectorizing to a multipolygon).

        Parameters
        ----------
        filename
            file path to write
        raster
            GDAL raster dataset
        band_index
            raster band (1-indexed)
        layer
            name of output layer
        """

        write_geometry(filename, vectorize_raster(raster, band_index), crs_wkt, name=os.path.basename(filename), layer=layer)

    def _gdal_raster_to_array(self, raster: gdal.Dataset) -> (np.array, dict):
        """
        Extract data and metadata from the given band of a GDAL raster dataset.
        The dataset should have the `nodata` value set appropriately.

        Parameters
        ----------
        raster
            GDAL raster array with `nodata` value set
        band_index
            raster band (1-indexed)

        Returns
        ----------
        numpy.array
            array of raster data and a dictionary of metadata
        """

        geotransform = raster.GetGeoTransform()

        raster_band = raster.GetRasterBand(1)
        nodata = raster_band.GetNoDataValue()
        del raster_band

        metadata = {
            'resx': geotransform[1],
            'resy': geotransform[5],
            'originx': geotransform[0],
            'originy': geotransform[3],
            'dimx': raster.RasterXSize,
            'dimy': raster.RasterYSize,
            'crs': raster.GetProjection(),
            'nodata': nodata
        }

        return raster.ReadAsArray(), metadata

    def _gdal_points_to_array(self, points: gdal.Dataset, layer: int = 0) -> (np.array, dict):
        """
        Extract points and metadata from the given layer of a GDAL point cloud dataset.
        The dataset should have the `nodata` value set appropriately.

        Parameters
        ----------
        points
            GDAL point cloud dataset with `nodata` value set
        layer
            index of vector layer containing points

        Returns
        ----------
        numpy.array
            N x M array of points and a dictionary of metadata
        """

        point_layer = points.GetLayerByIndex(layer)
        metadata = {'crs': point_layer.GetSpatialRef().ExportToWkt()}
        del point_layer

        return gdal_to_xyz(points, layer), metadata

    def _gdal_points_to_wkt(self, points: gdal.Dataset, layer: int = 0) -> (str, dict):
        """
        Extract WKT and metadata from the given layer of a GDAL point cloud dataset.
        The dataset should have the `nodata` value set appropriately.

        Parameters
        ----------
        points
            GDAL point cloud dataset with `nodata` value set
        layer
            index of vector layer containing points

        Returns
        ----------
        str, dict
            well-known text of points and a dictionary of metadata
        """

        points, metadata = self._gdal_points_to_array(points, layer)
        return MultiPoint(points.tolist()).wkt, metadata

    def _set_raster_nodata(self, raster: gdal.Dataset, band: int = 1) -> gdal.Dataset:
        """
        Ensure the `nodata` value of the given GDAL raster dataset is equal to the set `nodata` value.

        Parameters
        ----------
        raster
            GDAL raster array with `nodata` value set
        band
            raster band (1-indexed)

        Returns
        ----------
        gdal.Dataset
            dataset with `nodata` value set
        """

        raster_band = raster.GetRasterBand(band)

        # check the no data value
        nodata = raster_band.GetNoDataValue()
        if self._nodata != nodata:
            band_data = raster_band.ReadAsArray()
            band_data[band_data == nodata] = self._nodata
            raster.GetRasterBand(band).WriteArray(band_data)

        del raster_band
        return raster
