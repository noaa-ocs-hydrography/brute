import os
import io
import sys
import pathlib
import subprocess
import logging

from nbs.configs import get_logger, iter_configs, set_stream_logging, log_config, parse_multiple_values
from nbs.bruty.nbs_postgres import get_nbs_records, connect_params_from_config

LOGGER = get_logger('bruty.convert_csar')
CONFIG_SECTION = 'convert_csar'

def convert_csar(carisbatch, epsg, table_names, database, username, password, hostname='OCS-VS-NBS01',
                             port='5434', use_zip=False, dest_path=None):
    """Quick script that converts CSAR data using Caris' carisbatch.exe to convert to bag or xyz points"""

    for table_name in table_names:
        fields, records = get_nbs_records(table_name, database, username, password, hostname=hostname, port=port)

        cnt = 0
        for record in records:
            fname = record[fields.index('manual_to_filename')]
            if fname is None or not fname.strip():
                fname = record[fields.index('script_to_filename')]
            if fname is not None and fname.strip().lower().endswith("csar"):
                if record[fields.index('script_resolution')] or record[fields.index('manual_resolution')]:  # has grid filled out
                    if dest_path is not None:
                        local_fname = fname.lower().replace('\\\\nos.noaa\\OCS\\HSD\\Projects\\NBS\\NBS_Data'.lower(), dest_path)
                    else:
                        local_fname = fname
                    if not os.path.exists(f"{local_fname}.csv.zip") and not os.path.exists(f"{local_fname}.csv") and \
                       not os.path.exists(f"{local_fname}.depth.tif") and not os.path.exists(f"{local_fname}.elev.tif"):
                        cnt += 1
                        print(cnt, local_fname)
                        cmd = f'"{carisbatch}" -r ExportRaster --output-format GEOTIFF --compression LZW --include-band Depth --include-band Uncertainty "{local_fname}" "{local_fname}.depth.tif"'
                        p = subprocess.Popen(cmd)
                        p.wait()
                        if not os.path.exists(f"{local_fname}.depth.tif"):
                            cmd = f'"{carisbatch}" -r ExportRaster --output-format GEOTIFF --compression LZW --include-band Elevation --include-band Uncertainty "{local_fname}" "{local_fname}.elev.tif"'
                            p = subprocess.Popen(cmd)
                            p.wait()
                            if not os.path.exists(f"{local_fname}.elev.tif"):
                                cmd = f'"{carisbatch}" -r exportcoveragetoascii --include-band Depth 3 --include-band Uncertainty 3 --output-crs EPSG:{epsg} --coordinate-format GROUND --coordinate-precision 2 --coordinate-unit m "{local_fname}" "{local_fname}.csv"'
                                p = subprocess.Popen(cmd)
                                p.wait()
                                if os.path.exists(f"{local_fname}.csv"):
                                    if use_zip:
                                        p = subprocess.Popen(f'python -m zipfile -c "{local_fname}.csv.zip" "{local_fname}.csv"')
                                        p.wait()
                                        os.remove(f'{local_fname}.csv')
                                    print("was points")
                                else:
                                    print("failed as points and raster????????????????????")

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
        epsg = pathlib.Path(config['epsg'])

        tablenames, database, hostname, port, username, password = connect_params_from_config(config)
        caris_batch_path = config['carisbatch']
        convert_csar(caris_batch_path, epsg, tablenames, database, username, password, hostname, port)

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