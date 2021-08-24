import os
import shutil
from shapely import wkt, wkb
from osgeo import ogr, osr, gdal
import pickle
from functools import partial

import numpy

from data_management.db_connection import connect_with_retries
from fuse_dev.fuse.meta_review import meta_review
from nbs.bruty.utils import get_crs_transformer
from nbs.bruty.nbs_postgres import id_to_scoring, get_nbs_records, nbs_survey_sort, connect_params_from_config, make_contributor_csv
from nbs.bruty.world_raster_database import WorldDatabase
from nbs.bruty.generalize import generalize
from nbs.bruty.raster_attribute_table import make_raster_attr_table

"""
1) Get tile geometries and attributes
2) Expand tile geometry based on closing distance
3) Extract data from all necessary Bruty DBs
4) Combine if more than one DB used
5) Binary closing on raster covering expanded area
6) Interpolation only on closed cells that had no bathymetry
7) Export original extents to raster
8) Create Raster Attribute Table (RAT)
"""

# 1) Get tile geometries and attributes
_debug = True
use_gdal = True
b_glen_interp = True

if not _debug:
    user="postgres"
    password = "buildNBS2018!"
    host = 'OCS-VS-NBS01'
    port = '5434'
    connection = connect_with_retries(database="tile_specifications", user=user, password=password, host=host, port=port)
    cursor = connection.cursor()
    cursor.execute(f'SELECT * FROM bruty_bluetopo')
    records = cursor.fetchall()
    fields = [desc[0] for desc in cursor.description]
    pickle.dump([records, fields], open(r"C:\Data\bruty_databases\\bluetopo.pickle", "wb"))

    table_names = ["pbc_utm18n_mllw", "pbc_utm18n_mllw_sensitive", "pbc_utm18n_mllw_prereview",
                   "pbc_utm19n_mllw", "pbc_utm19n_mllw_sensitive", "pbc_utm19n_mllw_prereview"]
    metadata_fields = []
    metadata_records = []
    for table_name in table_names:
        fields, records = get_nbs_records(table_name, "metadata", user, password, hostname=host, port=port)
        metadata_records.append(records)
        metadata_fields.append(fields)
    pickle.dump([metadata_records, metadata_fields], open(r"C:\Data\bruty_databases\\metadata.pickle", "wb"))
else:
    records, fields = pickle.load(open(r"C:\Data\bruty_databases\\bluetopo.pickle", "rb"))
    metadata_records, metadata_fields = pickle.load(open(r"C:\Data\bruty_databases\\metadata.pickle", "rb"))

all_meta_records = {}
all_simple_records = {}
for n, meta_table_recs in enumerate(metadata_records):
    id_col = metadata_fields[n].index('nbs_id')
    for record in meta_table_recs:
        record_dict = dict(zip(metadata_fields[n], record))
        simple_record = meta_review.MetadataDatabase.simplify_record(record_dict)  # todo - make this a static method
        simple_fuse_record = meta_review.records_to_fusemetadata(simple_record)  # re-casts the db values into other/desired types
        all_simple_records[record[id_col]] = simple_fuse_record
    all_meta_records.update({rec[id_col]: rec for rec in meta_table_recs})
name_index = fields.index("name")
resolution_index = fields.index("resolution")
closing_index = fields.index("closing_distance")
sensitive_index = fields.index('data_sensitive')
reviewed_index = fields.index('data_qualified')
prereview_index = fields.index('data_unqualified')
not_for_nav_index = fields.index('quality_filter')

# print(fields)
# ['name',
#  'active',
#  'resolution',
#  'closing_distance',
#  'data_sensitive',
#  'data_qualified',
#  'data_unqualified',
#  'quality_filter',
#  'contributor_attributes',
#  'geometry']
for r in records:

    if "Tile4_PBC18_8b" in r[0]:
        r18_empty = r
    if "Tile29_PBC18_8" in r[0]:
        r18_square_8m = r
    if "Tile20_PBC18_4" in r[0]:
        r18_square = r
    if "Tile1_PBC18_4" in r[0]:
        r18_complex = r  # north_south and skinny

    if "Tile42_PBC19_4" in r[0]:
        r19_complex = r  # self crossing or two polygons?

    if "Tile39_PBC19_4b" in r[0]:
        r19_square = r

# print(r)
# Out[24]:
# ('Tile7_PBC18_4',
#  None,
#  4.0,
#  100.0,
#  False,
#  True,
#  False,
#  False,
#  None,
#  '01030000000100000009000000CDCCCCCCCC7C52C09A9999999919444066666666668652C09A9999999919444066666666668652C066666666663644409A999999999952C066666666663644409A999999999952C0000000000040444066666666666E52C0000000000040444066666666666E52C0CDCCCCCCCC2C4440CDCCCCCCCC7C52C0CDCCCCCCCC2C4440CDCCCCCCCC7C52C09A99999999194440')

for r in records:
    # USING GDAL
    if use_gdal:
        g = ogr.CreateGeometryFromWkb(bytes.fromhex(r[-1]))

        minx, maxx, miny, maxy = g.GetEnvelope()
        g.GetEnvelope()
        # Out[44]: (-74.4, -73.725, 40.2, 40.5)

        g.GetGeometryRef(0).GetPoints()
        # Out[48]:
        # [(-73.95, 40.2), (-74.1, 40.2), (-74.1, 40.425),  (-74.4, 40.425), (-74.4, 40.5),
        #  (-73.725, 40.5), (-73.725, 40.35), (-73.95, 40.35), (-73.95, 40.2)]
    else:
        # *************************************************************************************************
        # USING SHAPELY
        g = wkb.loads(r[-1], hex=True)

        minx, miny, maxx, maxy = g.bounds
        # Out[49]: (-74.4, 40.2, -73.725, 40.5)

        g.boundary.xy
        # Out[28]:
        # (array('d', [-73.95, -74.1, -74.1, -74.4, -74.4, -73.725, -73.725, -73.95, -73.95]),
        #  array('d', [40.2, 40.2, 40.425, 40.425, 40.5, 40.5, 40.35, 40.35, 40.2]))

    # 2) Expand tile geometry based on closing distance

    closing_dist = r[closing_index]

    target_epsg = 4326  # MCD WGS84
    if 'PBC18' in r[0].upper():
        db_epsg = 26918
        base_db = "pbc_utm18n_mllw"
    else:
        db_epsg = 26919
        base_db = "pbc_utm19n_mllw"

    crs_transform = get_crs_transformer(target_epsg, db_epsg)
    if crs_transform:
        xs, ys = [], []
        for x in (minx, maxx):
            for y in (miny, maxy):
                cx, cy = crs_transform.transform(x, y)
                xs.append(cx)
                ys.append(cy)
        minx = min(xs)
        maxx = max(xs)
        miny = min(ys)
        maxy = max(ys)

    # 3) Extract data from all necessary Bruty DBs
    db_base_path = r"C:\Data\bruty_databases"
    output_base_path = r"C:\Data\bruty_databases\output"
    REVIEWED = ""
    NOT_NAV = "_not_for_navigation"
    PREREVIEW = "_prereview"
    SENSITIVE = "_sensitive"
    databases = []
    tables = []
    for datatype, ext in [(reviewed_index, REVIEWED), (prereview_index, PREREVIEW), (sensitive_index, SENSITIVE)]:
        if r[datatype]:  # should we add the main datatype?  (reviewed/qualified, prereview/unqualified, sensitive)
            databases.append(base_db+ext)
            if r[not_for_nav_index]:  # should we add the "not for navigation" data?
                databases.append(base_db + ext + NOT_NAV)

    # @todo Use db.export_into_raster from each DB?
    #   Would just need to retain the score layers too?
    sorted_recs, names_list, sort_dict = id_to_scoring(metadata_fields, metadata_records, for_navigation_flag=(False, False), never_post_flag=(False, False))
    comp = partial(nbs_survey_sort, sort_dict)
    dataset, dataset_score = None, None
    export_filename = os.path.join(output_base_path, r[name_index]+".tif")
    for db_name in databases:
        db = WorldDatabase.open(os.path.join(db_base_path, db_name))
        if dataset is None:
            dataset, dataset_score = db.make_export_rasters(os.path.join(output_base_path, r[name_index]+" - original.tif"), minx-closing_dist, miny-closing_dist,
                                                            maxx+closing_dist, maxy+closing_dist, r[resolution_index])
            exported_filename = dataset.GetFileList()[0]
        db.export_into_raster(dataset, dataset_score, compare_callback=comp)
    # else:
    #     exported_filename = os.path.join(output_base_path, r[name_index]+" - original.tif")
    # FIXME -- remove this temprary copy and exporting to "-original"
    try:
        os.remove(export_filename)
    except FileNotFoundError:
        pass
    # FIXME - remove this change of nodata value once NBS is fixed to accept NaN in addition to 1000000
    shutil.copyfile(exported_filename, export_filename)
    ds = gdal.Open(export_filename, gdal.GA_Update)
    new_nodata = 1000000
    for b in range(1,4):
        band = ds.GetRasterBand(b)
        data = band.ReadAsArray()
        data[numpy.isnan(data)] = new_nodata
        band.SetNoDataValue(new_nodata)
        band.WriteArray(data)
    del band
    del ds


    # 4) Combine if more than one DB used
    # @todo if we use export_into_raster the combine would already be done

    # 5) Binary closing on raster covering expanded area
    # @todo get this code from xipe
    # ?? xipe_dev.xipe.raster.interpolate_within_closing()  # -- this looks to have a bunch of extra logic and such we don't want, Casiano suggested _close_interpolation_coverage

    if b_glen_interp:
        # @todo -- just call copy of process_csar instead
        # srs = osr.SpatialReference()
        # srs.ImportFromEPSG(target_epsg)
        generalize(export_filename, closing_dist, output_crs=target_epsg)  # call the generalize function which used to be process_csar script
        # FIXME - make raster attr table based on database -- was hack.make_enc_raster_attr_table_modified(raster_filename, contributor_table_filename)
        make_raster_attr_table(export_filename, all_simple_records)  # make a raster attribute table for the generalized dataset
        # create_RAT(dataset)  # make a raster attribute table for the raw dataset

    else:
        # interp_values = fuse.interpolator.bag_interpolator.process.RasterInterpolator()
        interpolated_dataset = raster_interp.process.RasterInterpolator().interpolate(raster, 'linear', buffer = closing_distance)

        # binary_mask = fuse.coverage.coverage._close_interpolation_coverage()
        closed_mask = _close_interpolation_coverage(closing_mask, closing_iterations)
        closed_mask = numpy.logical_and(closed_mask, raster_coverage)

        ops = BufferedImageOps(interpolated_dataset)
        for cnt, (ic, ir, cols, rows, col_buffer_lower, row_buffer_lower, nodata, data) in enumerate(
                ops.iterate_gdal(0, 0)):
            altered_mask = numpy.zeros(data[0].shape, dtype=numpy.bool8)
            tile = closed_mask[ir:ir + rows, ic:ic + cols]
            if data[0].shape != tile.shape:
                altered_mask[ir:, ic:] = tile
            else:
                altered_mask = tile
            if not upsampled:
                masked_interpolation = numpy.where(altered_mask, data[0], nodata)
            else:
                masked_interpolation = data[0]
            ops.write_array(masked_interpolation, interpolated_dataset.GetRasterBand(1))


        # 6) Interpolation only on closed cells that had no bathymetry
        # @todo get this code from xipe

        # 7) Export original extents to raster

        # 8) Create Raster Attribute Table (RAT)
        # @todo get this code from xipe
        # @todo contributor_metadata needs to be revised to work on the global indices

        xipe_v2.hack.make_raster_attributes_modified.make_enc_raster_attr_table_modified

        xipe.raster.attribute_tables.create_RAT(raster_filename, contributor_metadata, fields=None, compute_counts=True)
        data = RasterAttributeTableDataset.from_contributor_metadata(contributor_metadata)
        data.write_to_raster(raster_filename, CONTRIBUTOR_BAND_NAME, fields, compute_counts)


