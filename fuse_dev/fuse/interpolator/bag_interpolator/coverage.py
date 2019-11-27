# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:29:51 2019

@author: Casiano.Koprowski
"""

import os as _os
from datetime import datetime as _dt

import numpy as _np
from osgeo import gdal as _gdal, gdal
from osgeo import ogr as _ogr
from osgeo import osr as _osr
from scipy.ndimage.interpolation import zoom as _zoom

_gdal.UseExceptions()


# import matplotlib.pyplot as plt


def _maxValue(arr: _np.array):
    """
    Returns the most used value in the array as an integer

    Takes an input array and finds the most used value in the array, this
    value is used by the program to assume the array's nodata value

    Parameters
    ----------
    arr: _np.array :
        An input array

    Returns
    -------

    """

    nums, counts = _np.unique(arr, return_counts=True)
    index = _np.where(counts == _np.amax(counts))
    return int(nums[index])


class GeoTIFF:
    """
    This class serves as the main container for geotiff data.

    In order to read and assign the data relevant to the interpolation process,
    this class heavily depends on OSGEO's GDAL, OGR, and ORS libraries to open
    and populate the data needed

    Parameters
    ----------

    Returns
    -------

    Attributes
    ----------
    nodata : float
        1000000.0 is the nodata value associated with the BAG format
    elevation : numpy.array
        The elevation layer of the BAG
    uncertainty : numpy.array
        The uncertainty layer of the BAG
    shape : tuple
        The (y, x) dimensions of the layer data
    bounds : touple
        The geographic bounds of the data in order (nw, se)=([sx,ny],[nx,sy])
    wkt : str
        The WKT representation of the data CRS
    """

    def __init__(self):
        self.bounds = None
        self.resolution = None
        self.array = None
        self.shape = None
        self.nodata = None
        self.name = None

    def open_file(self, filename: str):
        """
        TODO write description

        Parameters
        ----------
        filename: str :
            The complete file path of the input coverage file

        Returns
        -------

        """

        _fName = _os.path.split(filename)[-1]
        self.name = _os.path.splitext(_fName)[0]
        _ds = _gdal.Open(filename)
        self.bounds, self.resolution = self._getBounds(_ds)
        self.array, self.shape, self.nodata = self._getArrayData(_ds)
        self.wkt = _ds.GetProjectionRef()
        del _ds

    def _getBounds(self, gdal_obj):
        """
        TODO write description

        Parameters
        ----------
        gdal_obj :
            TODO write description

        Returns
        -------

        """

        ulx, xres, xskew, uly, yskew, yres = gdal_obj.GetGeoTransform()
        check_origin(ulx, uly, self.name)
        lrx = ulx + (gdal_obj.RasterXSize * xres)
        lry = uly + (gdal_obj.RasterYSize * yres)

        if lrx < ulx:
            s = ulx
            ulx = lrx
            lrx = s

        if uly < lry:
            s = lry
            lry = uly
            uly = s

        meta = ([ulx, uly], [lrx, lry])
        return meta, (xres, yres)

    def _getArrayData(self, gdal_obj):
        """
        TODO write description

        Parameters
        ----------
        gdal_obj :
            TODO write description

        Returns
        -------

        """

        band = gdal_obj.GetRasterBand(1)
        array = band.ReadAsArray()
        del band
        maxVal = _maxValue(array)

        if maxVal != 0:
            array = (array < maxVal).astype(_np.int)
        else:
            array = (array > maxVal).astype(_np.int)

        shape = array.shape
        return array, shape, 0


class Geopackage:
    """
    This class serves as the main container for geopackage data.

    In order to read and assign the data relevant to the interpolation process,
    this class heavily depends on OSGEO's GDAL, OGR, and ORS libraries to open
    and populate the data needed

    A geopackage coverage requres a ``to_crs`` WKT spatial reference system in
    order to ensure that the contents of the geopackage can be properly
    rasterized.

    Parameters
    ----------

    Returns
    -------

    Attributes
    ----------
    nodata : float
        1000000.0 is the nodata value associated with the BAG format
    elevation : numpy.array
        The elevation layer of the BAG
    uncertainty : numpy.array
        The uncertainty layer of the BAG
    shape : tuple
        The (y, x) dimensions of the layer data
    bounds : touple
        The geographic bounds of the data in order (nw, se)=([sx,ny],[nx,sy])
    wkt : str
        The WKT representation of the data CRS
    """

    def __init__(self):
        self.bounds = None
        self.resolution = None
        self.array = None
        self.shape = None
        self.nodata = None
        self.name = None

    def open_file(self, filename: str, to_crs: str, pixel_size: int = 1, nodata: int = 255):
        """
        TODO write description

        Parameters
        ----------
        filename: str :
            The complete file path of the input coverage file
        to_crs: str :
            WKT object with destination spatial reference system
        pixel_size: int :
            TODO write description
             (Default value = 1)
        nodata: int :
            TODO write description
             (Default value = 255)

        Returns
        -------

        """

        self.wkt = to_crs
        self.resolution = (pixel_size, -pixel_size)
        self.nodata = 0
        self.name = self._name(filename)
        self.array, self.shape, self.bounds = self._readAndRasterize(filename,
                                                                     to_crs, pixel_size, nodata)

    def _readAndRasterize(self, filename: str, to_crs: str, pixel_size: int = 1, nodata: int = 255):
        """
        TODO write description

        Parameters
        ----------
        filename: str :
            The complete file path of the input coverage file
        to_crs: str :
            WKT object with destination spatial reference system
        pixel_size: int :
            TODO write description
             (Default value = 1)
        nodata: int :
            TODO write description
             (Default value = 255)

        Returns
        -------

        """

        fName = _os.path.split(filename)[-1]
        splits = _os.path.splitext(fName)
        name = splits[0]

        # Create destination Spatial Ref from wkt
        to_srs = _osr.SpatialReference(wkt=to_crs)

        # Open the data source and read in the extent
        source_ds = _ogr.Open(filename)
        source_layer = source_ds.GetLayer()
        source_srs = source_layer.GetSpatialRef()

        for feature in source_layer:
            if feature is not None:
                geom = feature.GetGeometryRef()
                # print(geom.ExportToWkt())
                ds_geom = _ogr.CreateGeometryFromWkt(geom.ExportToWkt())
                # print(source_srs, to_proj, sep='\n')
                coordTrans = _osr.CoordinateTransformation(source_srs, to_srs)
                ds_geom.Transform(coordTrans)
                driver = _ogr.GetDriverByName('Memory')
                ds = driver.CreateDataSource('temp')
                layer = ds.CreateLayer(name, to_srs, _ogr.wkbMultiPolygon)

                # Add one attribute
                layer.CreateField(_ogr.FieldDefn('Survey', _ogr.OFTString))
                defn = layer.GetLayerDefn()

                # Create a new feature (attribute and geometry)
                feat = _ogr.Feature(defn)
                feat.SetField('Survey', name)

                # Make a geometry, from Shapely object
                feat.SetGeometry(ds_geom)

                layer.CreateFeature(feat)
                del feat, geom  # destroy these
                break

        x_min, x_max, y_min, y_max = ds_geom.GetEnvelope()
        bounds = ([x_min, y_max], [x_max, y_min])

        # Create the destination data source
        x_dim = int((x_max - x_min) / pixel_size)
        y_dim = int((y_max - y_min) / pixel_size)
        target_ds = _gdal.GetDriverByName('MEM').Create('', x_dim, y_dim,
                                                        1, _gdal.GDT_Float32)
        gt = (x_min, pixel_size, 0, y_max, 0, -pixel_size)
        check_origin(x_min, y_min, fName)
        target_ds.SetGeoTransform(gt)
        target_ds.SetProjection(to_crs)
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(nodata)

        # Rasterize
        _gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])

        array = band.ReadAsArray()
        shape = array.shape

        del band, source_ds, target_ds

        return array, shape, bounds

    def _name(self, filepath: str):
        """
        TODO write description

        Parameters
        ----------
        filepath: str :
            TODO write description

        Returns
        -------

        """

        fName = _os.path.split(filepath)[-1]
        return _os.path.splitext(fName)[0]


class UnifiedCoverage:
    """TODO write description"""

    def __init__(self, coverage_files, bag_wkt=None, bag_name='Test_Data.bag'):
        self.name = None
        self.wkt = bag_wkt
        self.array = None
        self.shape = None
        self.bounds = None
        self.resolution = (1.0, -1.0)
        self.nodata = 0
        _rasters = self._open_data(coverage_files, self.wkt)
        self._align_and_combine(_rasters, bag_name)
        del _rasters

    def _open_data(self, files: [str], bag_wkt: str):
        """
        TODO write description

        Parameters
        ----------
        files: [str] :
            TODO write description
        bag_wkt: str :
            TODO write description

        Returns
        -------

        """

        print('_open')
        bndRasts = []
        y = 0
        error_string = ''

        for item in files:
            fName = _os.path.split(item)[1]
            ext = _os.path.splitext(fName)[1].lower()

            try:
                if ext in ('.tiff', '.tif'):
                    rast = GeoTIFF()
                    rast.open_file(item)
                    bndRasts.append(rast)
                elif ext in ('.gpkg',):
                    rast = Geopackage()
                    rast.open_file(item, bag_wkt)
                    bndRasts.append(rast)
            except ValueError as e:
                error_string += f'\n{e}'
                pass

            y += 1

        if len(bndRasts) == 0:
            raise RuntimeError(f'No valid coverage in {files}{error_string}')

        return bndRasts

    def _align_and_combine(self, rasters: list, name: str):
        """
        TODO write description

        Parameters
        ----------
        rasters: list:
            TODO write description
        name: str :
            TODO write description

        Returns
        -------

        """

        sized, self.bounds = self._align(rasters)
        self.array, self.name, self.shape = self._combine(sized, name)
        del sized

    def _combine(self, rasters: list, name: str):
        """
        Heavy Influence From:
        "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
        (https://gis.stackexchange.com/a/278965)]

        This function takes an array input of coverage objects, a destination
        path (for ouput saving), and a list the names of the input coverage
        objects.

        Takes the arrays of each coverage object and combines them along the z
        axis of each::

            [[za,zb],[za,zb],[za,zb],
             [za,zb],[za,zb],[za,zb],
             [za,zb],[za,zb],[za,zb]]

        Then takes the mean of the values along the z axis. The reslult is then
        modified to a binary raster by making all values that are not the
        nodata value 1 and the values that are nodata 0.

        This binary raster is saved at the input path and also returned with
        the new full path for the output

        Parameters
        ----------
        rasters: list :
            List of coverage objects
        name: str :
            File path of the relevant bag file
        """

        print('_combine')
        shape = [0, 0]
        maxVal = 0
        coverageList = []
        x = 0

        for raster in rasters:
            # print(raster.shape)
            maxVal = raster.nodata
            # print(maxVal)
            cols, rows = raster.shape

            if x == 0:
                # print('original', shape)

                if cols >= shape[0]:
                    shape[0] = cols

                if rows >= shape[1]:
                    shape[1] = rows

                # print('original', shape)
                coverageList.append([x, raster.array])
            elif x != 0:
                if [cols, rows] != shape:
                    # print('nope')
                    pass
                else:
                    coverageList.append([x, raster.array])

            # print('yup')1
            x += 1

        # plt.figure()
        covDict = dict(coverageList)

        for i, grid in covDict.items():
            # plt.imshow(grid)
            # plt.show()
            array = _np.expand_dims(grid, 2)

            if i == 0:
                allarrays = array
            else:
                allarrays = _np.concatenate((allarrays, array), axis=2)

        meanCoverage = _np.nanmean(allarrays, axis=2)
        # plt.imshow(meanCoverage)
        # plt.show()

        if maxVal != 0:
            meanCoverage = (meanCoverage < maxVal).astype(_np.int)
        else:
            meanCoverage = (meanCoverage > maxVal).astype(_np.int)

        # plt.imshow(meanCoverage)
        # plt.show()

        fullname = _os.path.split(name)[1]
        justname = _os.path.splitext(fullname)[0]
        snum = justname.split('_')[0]
        outputname = f'{snum}_COMBINEDPOLY'

        return meanCoverage, outputname, tuple(shape)

    def _align(self, rasters: list):
        """
        Takes an input of an array of coverage objects. The goal of this
        function is to fit the provided coverage objects to the largest
        combined area of the inputs.


        1. Find NW and SE corners of the first input array and the number of rows and columns in the input.
        2. Find NW and SE corners of subsiquent input arrays and record differences as the NW differences (nxd, nyd) SE differences (sxd, syd).

            - also modify the original NW and/or SE corner values if applicable.
            - also modify the original number or rows and columns if applicable.

        3. Resize the area of each input array to the new total combined size of all input arrays based on the stored number of rows and columns.
        4. Reposition the data within each input array to maintain original geographic placement based on the difference of the input array's original NW corner and the recorded 'largest' NW boundry of the combined input areas.
        5. Return coverage objects with the newly resized/repositioned tifs and the new NW and SE extents of the new combined area.

        Parameters
        ----------
        rasters: list :
            List of coverage objects created by :func:`_open_data`
        """

        print('_align')
        x = 0
        nw, se = [0, 0], [0, 0]
        cols, rows = 0, 0
        sxd, nyd = 0, 0
        nxd, syd = 0, 0

        for grid in rasters:
            bounds = grid.bounds
            ul = bounds[0]
            lr = bounds[-1]

            if x == 0:
                nw = ul
                se = lr
                # print(nw, se)
                x += 1
            else:
                ulx, uly = ul
                lrx, lry = lr
                nwx, nwy = nw
                scx, scy = se

                if ulx <= nwx:
                    nxd = nwx - ulx
                    nwx = ulx
                elif ulx >= nwx:
                    nxd = ulx - nwx

                if uly >= nwy:
                    nyd = uly - nwy
                    nwy = uly
                elif uly <= nwy:
                    nyd = nwy - uly

                if lrx >= scx:
                    sxd = lrx - scx
                    scx = lrx
                elif lrx <= scx:
                    sxd = scx - lrx

                if lry <= scy:
                    syd = lry - scy
                    scy = lry
                elif lry >= scy:
                    syd = scy - lry

                # print(nxd, nyd, sxd, syd)
                # print(sxd, nyd)
                nw = [nwx, nwy]
                se = [scx, scy]
                cols, rows = (int(_np.round(nwy - scy)), int(_np.round(scx - nwx)))

        # print('cols: ' , nwy - scy, '\nrows: ', scx - nwx)
        # print(nw, se)
        sizedCoverage = []
        # print('resize?')

        if nxd != 0 or nyd != 0 or sxd != 0 or syd != 0:
            # plt.figure()
            # print('yes')#, _dt.now())
            # ref = _np.full((cols, rows), 0)
            # print(ref.shape)
            for grid in rasters:
                maxVal = grid.nodata
                bndx, bndy = grid.bounds[0]
                nwx, nwy = nw
                # print('old:', bndx, bndy)
                # print('new', nwx, nwy)
                arr = grid.array
                y, x = arr.shape
                # print(y, rows)
                # print(x, cols)

                if x != cols:
                    exp = int(abs(rows - x))
                    # print(exp)
                    add = _np.full((y, exp), maxVal)
                    arr = _np.column_stack([arr, add])
                    y, x = arr.shape

                if y != rows:
                    exp = int(abs(cols - y))
                    # print(exp)
                    add = _np.full((exp, x), maxVal)
                    arr = _np.vstack([arr, add])

                # print(bef, arr.shape)

                if nwx != bndx:
                    rollx = bndx - nwx
                    # print(rollx)
                    arr = _np.roll(arr, int(rollx), axis=1)

                if nwy != bndy:
                    rolly = nwy - bndy
                    # print(rolly)
                    arr = _np.roll(arr, int(rolly), axis=0)

                grid.array = arr
                # plt.imshow(grid.array)
                #  plt.show()
                grid.bounds = (nw, se)
                sizedCoverage.append(grid)

            bounds = (nw, se)
            # print(bounds)
            # print('done')#, _dt.now())
            return sizedCoverage, bounds
        else:
            # print('same')
            bounds = (nw, se)
            return rasters, bounds


def check_origin(x: float, y: float, filename: str) -> ValueError:
    if 0 in (x, y):
        raise ValueError(f'{filename} - Origin is not correctly georeferenced')


def align2grid(coverage, bounds: ((float, float), (float, float)), shape: (int, int), resolution: (float, float),
               nodata: float):
    """
    Takes an input of two arrays representing bag and tif data. These arrays
    hold information like extent, data, and more. The goal of this function is
    to fit and/or shift the provided tif to the area of the bag.

    Steps:

    1. Determine/confirm the resolution of the tif using a mix of file naming, extent, and dimensions.
    2. Determine the resolution of the bag using file naming conventions.
    3. Calculate the tif/bag resolution to determine the value needed to resize the tif.
    4. Resize the tif so that the resolution matches the bag using.
    5. Convert tif values to either numpy.nan or maxVal for use later.
    6. Store the extents of the bag and tif and calculate the difference.
    7. Create an empty array with dimensions equal to the original extents of the tiff data plus the additional space needed to match the total extents of the BAG data.
    8. Align the content of the tiff as needed in order to match the bag and apply the shifted content to the empty array created in Step 7.
    9. Enforce the size of the bag onto the shifted tiff data stored in the array created and modified in Step(s) 7 and 8.

    Parameters
    ----------
    coverage
        Input coverage data object
    bounds
        The ([nx, ny], [sx, sy]) extents to be applied to the input data
    shape
        The (y, x) shape to to be applied to the input data
    resolution
        The (x, y) resolution to be applied to the input data
    nodata
        The nodata value to be applied to the input array object

    Returns
    -------

    """

    print('align2grid')

    # determine the ratio resolution between the BAG resolution and the coverage resolution
    resolution_ratio = coverage.resolution[0] / resolution[0]

    # resample coverage data to match BAG resolution
    if resolution_ratio == 1:
        resampled_coverage_array = coverage.array
    else:
         resampled_coverage_array = _zoom(coverage.array, zoom=[resolution_ratio, resolution_ratio], order=3,
                                         prefilter=False)
#        resampled_coverage_array = _np.kron(coverage.array, _np.ones((resolution_ratio, resolution_ratio)))

    resampled_coverage_array = resampled_coverage_array.astype('float64')
    resampled_coverage_array[resampled_coverage_array > 0] = _np.nan
    resampled_coverage_array[resampled_coverage_array < 1] = float(nodata)

    # clip resampled coverage data to the bounds of the BAG
    output_array = _np.full(shape, nodata)

    cov_ul, cov_lr = _np.array(coverage.bounds[0]), _np.array(coverage.bounds[1])
    bag_ul, bag_lr = _np.array(bounds[0]), _np.array(bounds[1])

    if bag_ul[0] > cov_lr[0] or bag_lr[0] < cov_ul[0] or bag_lr[1] > cov_ul[1] or bag_ul[1] < cov_lr[1]:
        raise ValueError('bag dataset is outside the bounds of coverage dataset')

    ul_index_delta = _np.floor((bag_ul - cov_ul) / _np.array(resolution)).astype(int)
    lr_index_delta = _np.floor((bag_lr - cov_ul) / _np.array(resolution)).astype(int)

    # indices to be written onto the output array
    output_array_index_slices = [slice(0, None), slice(0, None)]

    # BAG leftmost X is to the left of coverage leftmost X
    if ul_index_delta[0] < 0:
        output_array_index_slices[1] = slice(ul_index_delta[0] * -1, output_array_index_slices[1].stop)
        ul_index_delta[0] = 0

    # BAG topmost Y is above coverage topmost Y
    if ul_index_delta[1] < 0:
        output_array_index_slices[0] = slice(ul_index_delta[1] * -1, output_array_index_slices[0].stop)
        ul_index_delta[1] = 0

    # BAG rightmost X is to the right of coverage rightmost X
    if lr_index_delta[0] > resampled_coverage_array.shape[1]:
        output_array_index_slices[1] = slice(output_array_index_slices[1].start,
                                             resampled_coverage_array.shape[1] - lr_index_delta[0])
        lr_index_delta[0] = resampled_coverage_array.shape[1]

    # BAG bottommost Y is lower than coverage bottommost Y
    if lr_index_delta[1] > resampled_coverage_array.shape[0]:
        output_array_index_slices[0] = slice(output_array_index_slices[0].start,
                                             resampled_coverage_array.shape[0] - lr_index_delta[1])
        lr_index_delta[1] = resampled_coverage_array.shape[0]

    # write the relevant coverage data to a slice of the output array corresponding to the coverage extent
    output_array[output_array_index_slices[0], output_array_index_slices[1]] = resampled_coverage_array[
                                                                               ul_index_delta[1]:lr_index_delta[1],
                                                                               ul_index_delta[0]:lr_index_delta[0]]
    del resampled_coverage_array

    coverage.array = output_array
    coverage.bounds = bounds
    coverage.shape = shape
    coverage.resolution = resolution

    return coverage


def write_raster(coverage, outputpath: str, out_verdat: str = 'MLLW', dtype=_gdal.GDT_UInt32, options: int = 0,
                 color_table: int = 0, nbands: int = 1, nodata: bool = False):
    """
    Directly From:
    "What is the simplest way..." on GIS Stack Exchange [Answer by 'Jon'
    (https://gis.stackexchange.com/a/278965)]

    options=['COMPRESS=LZW']

    Parameters
    ----------
    coverage :
        TODO write description
    outputpath: str :
        TODO write description
    out_verdat: str :
        Default value = 'MLLW')
    dtype :
        Default value = _gdal.GDT_UInt32)
    options: int :
        Default value = 0)
    color_table: int :
        Default value = 0)
    nbands: int :
        Default value = 1)
    nodata: float:
        Default value = False)

    Returns
    -------

    """

    print('write_raster')

    height, width = coverage.shape

    # Prepare destination file
    driver = _gdal.GetDriverByName("GTiff")
    if options != 0:
        dest = driver.Create(outputpath, width, height, nbands, dtype, options)
    else:
        dest = driver.Create(outputpath, width, height, nbands, dtype)

    # Write output raster
    if color_table != 0:
        dest.GetRasterBand(1).SetColorTable(color_table)

    dest.GetRasterBand(1).WriteArray(coverage.array)

    if nodata is not False:
        dest.GetRasterBand(1).SetNoDataValue(nodata)

    # Set transform and projection
    x_res, y_res = coverage.resolution
    x_orig, y_orig = coverage.bounds[0]
    gt = (x_orig, x_res, 0, y_orig, 0, y_res)
    dest.SetGeoTransform(gt)
    dest.SetProjection(coverage.wkt)
    dest.SetVertCS(out_verdat, out_verdat, 2000)

    # Close output raster dataset
    del dest


def coverage2gdal(coverage, flip: bool = False) -> gdal.Dataset:
    """
    TODO write description

    Parameters
    ----------
    coverage :
        TODO write description
    flip: bool :
        TODO write description

    Returns
    -------
    type
        gdal dataset

    """

    proj = coverage.wkt
    height, width = coverage.shape
    sw, ne = coverage.bounds
    scx, scy = sw
    nex, ney = ne
    res_x, res_y = coverage.resolution
    gt = (scx, res_x, 0, scy, 0, res_y)
    coverage_gdal = _gdal.GetDriverByName('MEM').Create('', width, height, 1, _gdal.GDT_Float32)
    coverage_gdal.SetGeoTransform(gt)
    coverage_gdal.SetProjection(proj)

    band = coverage_gdal.GetRasterBand(1)
    band.SetNoDataValue(float(coverage.nodata))

    if flip:
        array = _np.flipud(coverage.array)
        band.WriteArray(array)
    else:
        band.WriteArray(coverage.array)

    return coverage_gdal


def write_vector(coverage, outputpath: str, out_verdat: str = 'MLLW', flip: bool = False):
    """
    TODO write description

    Parameters
    ----------
    coverage :
        TODO write description
    outputpath: str :
        TODO write description
    out_verdat: str :
        TODO write description (Default value = 'MLLW')
    flip: bool :
        TODO write description

    Returns
    -------

    """
    float_resolution = abs(coverage.resolution[0])
    if float_resolution < 1:
        resolution = f"{str(float_resolution)[2:]}cm"
    else:
        resolution = f"{str(int(float_resolution))}m"

    name = f'{coverage.name}_{resolution}.gpkg'
    outfilename = _os.path.join(outputpath, name)

    proj = _osr.SpatialReference(wkt=coverage.wkt)
    proj.SetVertCS(out_verdat, out_verdat, 2000)

    cov_ds = coverage2gdal(coverage, flip=flip)
    band = cov_ds.GetRasterBand(1)

    driver = _ogr.GetDriverByName('GPKG')
    ds = driver.CreateDataSource(outfilename)
    layer = ds.CreateLayer(name, proj, _ogr.wkbPolygon)

    # Add one attribute
    layer.CreateField(_ogr.FieldDefn('Survey', _ogr.OFTString))
    defn = layer.GetLayerDefn()

    # Create a new feature (attribute and geometry)
    feat = _ogr.Feature(defn)
    feat.SetField('Survey', name)

    _gdal.Polygonize(band, band, layer, 0, [], callback=None)

    del cov_ds, band
