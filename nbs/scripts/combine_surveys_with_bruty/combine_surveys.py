import os
import sys
import time
import traceback
from datetime import datetime
import subprocess
import pathlib
import shutil
from functools import partial
import logging
import io

from nbs.bruty.raster_data import TiffStorage, LayersEnum
from nbs.bruty.history import DiskHistory, RasterHistory, AccumulationHistory
from nbs.bruty.world_raster_database import WorldDatabase, UTMTileBackend, UTMTileBackendExactRes, \
    LockNotAcquired, Lock, EXCLUSIVE, SHARED, NON_BLOCKING
from nbs.bruty.utils import onerr
from nbs.configs import get_logger, iter_configs, set_stream_logging, log_config, parse_multiple_values
from nbs.bruty.nbs_postgres import id_to_scoring, get_nbs_records, nbs_survey_sort, connect_params_from_config

_debug = False

LOGGER = get_logger('bruty.insert')
CONFIG_SECTION = 'insert'


def process_nbs_database(world_db_path, table_names, database, username, password, hostname='OCS-VS-NBS01', port='5434'):
    all_fields = []
    all_records = []
    for table_name in table_names:
        fields, records = get_nbs_records(table_name, database, username, password, hostname=hostname, port=port)
        all_records.append(records)
        all_fields.append(fields)
    sorted_recs, names_list, sort_dict = id_to_scoring(all_fields, all_records)

    db = WorldDatabase.open(world_db_path)
    comp = partial(nbs_survey_sort, sort_dict)
    print('------------   changing paths !!!!!!!!!!')
    while names_list:
        num_names = len(names_list)
        for i in range(num_names-1, -1, -1):
            (filename, sid, path) = names_list[i]
            if _debug:
                pass
                # if sid not in (13425, 13562, 10035):
                #     continue
                # if i > 2:
                #     break

                path_e = path.lower().replace('\\\\nos.noaa\\ocs\\hsd\\projects\\nbs\\nbs_data\\pbc_northeast_utm19n_mllw',
                                              r'E:\Data\nbs\PBC_Northeast_UTM19N_MLLW')
                path_c = path.lower().replace('\\\\nos.noaa\\ocs\\hsd\\projects\\nbs\\nbs_data\\pbc_northeast_utm19n_mllw',
                                              r'C:\Data\nbs\PBC_Northeast_UTM19N_MLLW')
                copy_data = False
                if copy_data:
                    try:
                        os.makedirs(os.path.dirname(path_c), exist_ok=True)
                        if not path_e.endswith("csar"):
                            pass
                            # shutil.copy(path_e, path_c)
                        else:
                            for mod_fname in (f"{path_e}.elev.tif", f"{path_e}.depth.tif", f"{path_e}.csv.zip"):
                                if os.path.exists(mod_fname):
                                    shutil.copy(mod_fname, "c"+mod_fname[1:])
                                    try:
                                        os.remove(path_c)  # take the csar off disk
                                        os.remove(path_c+"0")
                                    except:
                                        pass
                    except FileNotFoundError:
                        print("File missing", sid, path)
                # convert csar names to exported data, 1 of 3 types
                if path.endswith("csar"):
                    for mod_fname in (f"{path_c}.elev.tif", f"{path_c}.depth.tif", f"{path_c}.csv.zip"):
                        if os.path.exists(mod_fname):
                            path = mod_fname
                else:
                    path = path_c
            if not os.path.exists(path):
                print(path, "didn't exist")
                names_list.pop(i)
                continue
            # # @FIXME is contributor an int or float -- needs to be int 32 and maybe int 64 (or two int 32s)
            print('starting', path)
            print(datetime.now().isoformat(), num_names-i, "of", num_names)
            # FIXME there is the possibility that we load metadata looking for SID=xx while it is being processed.
            #    Then it gets written to disk as we figure out what tiles to lock.
            #    We could check in the insert function again (once locks are obtained) to make sure survey=xx is not in the already processed list.
            sid_in_db = True
            if sid not in db.included_ids:
                db.update_metadata_from_disk()  # if not in the cached list then load from disk in case another process changed it in the interim
                if sid not in db.included_ids:
                    sid_in_db = False

            # @todo fixme - make a being processed list and check that the survey is not already being processed.
            #   This is an issue with the zip files overwriting a file being read
            #   but longer term in not starting to process/read the same file - especially for point data where we read the whole
            #   dataset to figure out the tiles to lock

            if not sid_in_db:
                try:
                    lock = Lock(path)  # this doesn't work with the file lock - just the multiprocessing locks
                    if lock.acquire():
                        if path.endswith(".csv.zip"):
                            csv_path = path[:-4]
                            print(f"Extract CSV {path}")
                            p = subprocess.Popen(f'python -m zipfile -e "{path}" "{os.path.dirname(path)}"')
                            p.wait()
                            if os.path.exists(csv_path):
                                try:
                                    # points are in opposite convention as BAGs and exported CSAR tiffs, so reverse the z component
                                    db.insert_txt_survey(csv_path, dformat=[('x', 'f8'), ('y', 'f8'), ('depth', 'f4'), ('uncertainty', 'f4')],
                                                         override_epsg=db.db.epsg, contrib_id=sid, compare_callback=comp, reverse_z=True)
                                except ValueError:
                                    print("Value Error")
                                    print(traceback.format_exc())

                                try:
                                    os.remove(f'{csv_path}')
                                except FileNotFoundError:
                                    print(f'File NOT Found:  {csv_path}')

                                except PermissionError:
                                    time.sleep(1)
                                    try:
                                        os.remove(f'{csv_path}')
                                    except:
                                        print(f"failed to remove{csv_path}")
                            else:
                                print("\n\nCSV was not extracted from zip\n\n\n")
                        elif path.endswith(".npy"):
                            db.insert_txt_survey(path, dformat=[('x', 'f8'), ('y', 'f8'), ('depth', 'f8'), ('uncertainty', 'f8')],
                                                 override_epsg=db.db.epsg, contrib_id=sid, compare_callback=comp, reverse_z=True)
                        else:
                            try:
                                db.insert_survey(path, override_epsg=db.db.epsg, contrib_id=sid, compare_callback=comp)
                            except ValueError:
                                print("Value Error")
                                print(traceback.format_exc())
                        print('inserted', path)
                        names_list.pop(i)
                    else:
                        # print(f"{path} was locked - probably another process is working on it")
                        raise LockNotAcquired()
                except LockNotAcquired:
                    print('files in use for ', sid, path)
                    print('skipping to next survey')
            else:
                print(f"{sid} already in database")
                names_list.pop(i)



def main():
    if len(sys.argv) > 1:
        use_configs = sys.argv[1:]
    else:
        use_configs = pathlib.Path(__file__).parent.resolve()  # (os.path.dirname(os.path.abspath(__file__))

    warnings = ""
    for config_filename, config_file in iter_configs(use_configs):
        stringio_warnings = set_stream_logging("bruty", file_level=logging.WARNING, remove_other_file_loggers=False)
        LOGGER.info(f'***************************** Start Run  *****************************')
        LOGGER.info(f'reading "{config_filename}"')
        log_config(config_file, LOGGER)

        config = config_file[CONFIG_SECTION if CONFIG_SECTION in config_file else 'DEFAULT']
        db_path = pathlib.Path(config['combined_datapath'])
        if not os.path.exists(db_path.joinpath("wdb_metadata.json")):
            try:
                resx, resy = map(float, config['resolution'].split(','))
            except:
                resx = resy = float(config['resolution'])
            # NAD823 zone 19 = 26919.  WGS84 would be 32619
            epsg = int(config['epsg'])
            # use this to align the database to something else (like caris for testing)
            offset_x = float(config['offset_x']) if 'offset_x' in config else 0
            offset_y = float(config['offset_y']) if 'offset_y' in config else 0
            db = WorldDatabase(
                UTMTileBackendExactRes(resx, resy, epsg, RasterHistory, DiskHistory, TiffStorage, db_path, offset_x=offset_x, offset_y=offset_y))
            del db

        if _debug:
            hostname, port, username, password = None, None, None, None
            tablenames_raw, database = config['tablenames'], config['database']
            tablenames = parse_multiple_values(tablenames_raw)
        else:
            tablenames, database, hostname, port, username, password = connect_params_from_config(config)

        process_nbs_database(db_path, tablenames, database, username, password, hostname, port)

    # data_dir = pathlib.Path("c:\\data\\nbs\\test_data_output")  # avoid putting in the project directory as pycharm then tries to cache everything I think
    # def make_clean_dir(name):
    #     use_dir = data_dir.joinpath(name)
    #     if os.path.exists(use_dir):
    #         shutil.rmtree(use_dir, onerror=onerr)
    #     os.makedirs(use_dir)
    #     return use_dir
    # subdir = r"test_pbc_19_exact_multi_locks"
    # db_path = data_dir.joinpath(subdir)
    # make_clean_dir(subdir)

    # # create logger with 'spam_application'
    # logger = logging.getLogger('process_nbs')
    # logger.setLevel(logging.DEBUG)
    # # create file handler which logs even debug messages
    # fh = logging.FileHandler(db_path.joinpath('process_nbs.log'))
    # fh.setLevel(logging.DEBUG)
    # # create console handler with a higher log level
    # ch = logging.StreamHandler()
    # ch.setLevel(logging.INFO)
    # # create formatter and add it to the handlers
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s')
    # fh.setFormatter(formatter)
    # ch.setFormatter(formatter)
    # # add the handlers to the logger
    # logger.addHandler(fh)
    # logger.addHandler(ch)
    #

    # db_path = make_clean_dir(r"test_pbc_19_db")  # reset the database


if __name__ == '__main__':

    # default_config_name = "default.config"

    # turn prints into logger messages
    orig_print = print
    def print(*args, **kywds):
        f = io.StringIO()
        ky = kywds.copy()
        ky['file'] = f
        orig_print(*args, **ky)  # build the string
        LOGGER.info(f.getvalue()[:-1])  # strip the newline at the end
    main()


# "V:\NBS_Data\PBA_Alaska_UTM03N_Modeling"
# UTMN 03 through 07 folders exist
# /metadata/pba_alaska_utm03n_modeling
# same each utm has a table
# \\nos.noaa\OCS\HSD\Projects\NBS\NBS_Data\PBA_Alaska_UTM03N_Modeling
# v drive literal