import numpy

from nbs.bruty.raster_data import TiffStorage, LayersEnum
from data_management.db_connection import connect_with_retries
from fuse_dev.fuse.meta_review.meta_review import database_has_table, split_URL_port

_debug = True

def connect_params_from_config(config):
    with open(config['URL_FILENAME']) as hostname_file:
        url = hostname_file.readline()
    hostname, port = split_URL_port(url)

    with open(config['CREDENTIALS_FILENAME']) as database_credentials_file:
        username, password = [line.strip() for line in database_credentials_file][:2]
    tablename, database = config['tablename'], config['database']
    return tablename, database, hostname, port, username, password

def make_serial_column(table_id, table_name, database, username, password, hostname='OCS-VS-NBS01', port='5434', big=False):
    connection = connect_with_retries(database=database, user=username, password=password, host=hostname, port=port)
    cursor = connection.cursor()
    # used admin credentials for this
    # cursor.execute("create table serial_19n_mllw as (select * from pbc_utm19n_mllw)")
    # connection.commit()
    if big:
        serial_str = "bigserial"  # use bigserial if 64bit is needed
        bit_shift = 32  # allows 4 billion survey IDs and 4 billion table IDs
    else:
        serial_str = "serial"  # 32bit ints
        bit_shift = 20  # allows for one million survey IDs and 4096 table IDs
    start_val = table_id << bit_shift
    cursor.execute(f'ALTER TABLE {table_name} ADD column sid {serial_str};')
    cursor.execute(f'ALTER SEQUENCE {table_name}_sid_seq RESTART WITH 10000')
    cursor.execute(f"update {table_name} set sid=sid+{start_val}")
    connection.commit()


def get_nbs_records(table_name, database, username, password, hostname='OCS-VS-NBS01', port='5434'):
    if _debug and hostname is None:
        import pickle
        f = open(fr"C:\data\nbs\{table_name}.pickle", 'rb')
        records = pickle.load(f)
        fields = pickle.load(f)
        ## trim back the records to a few for testing
        # print("Thinning survey records for debugging!!!!")
        # filename_col = fields.index('from_filename')
        # id_col = fields.index('sid')
        # thinned_records = []
        # for rec in records:
        #     if rec[id_col] in (12657, 12203, 12772, 10470, 10390):
        #         thinned_records.append(rec)
        # records = thinned_records
    else:
        connection = connect_with_retries(database=database, user=username, password=password, host=hostname, port=port)
        cursor = connection.cursor()
        cursor.execute(f'SELECT * FROM {table_name}')
        records = cursor.fetchall()
        fields = [desc[0] for desc in cursor.description]
    return fields, records



def id_to_scoring(fields, records):
    # Create a dictionary that converts from the unique database ID to an ordering score
    # Basically the standings of the surveys,
    # First place is the highest decay score with a tie breaker of lowest resolution.  If both are the same they will have the same score
    # Alphabetical will have no duplicate standings (unless there is duplicate names) and first place is A (ascending alphabetical)
    # get the columns that have the important data
    decay_col = fields.index("decay_score")
    script_res_col = fields.index('script_resolution')
    manual_res_col = fields.index('manual_resolution')
    script_point_res_col = fields.index('script_point_spacing')
    manual_point_res_col = fields.index('manual_point_spacing')

    filename_col = fields.index('from_filename')
    path_col = fields.index('script_to_filename')
    manual_path_col = fields.index('manual_to_filename')
    id_col = fields.index('sid')
    rec_list = []
    names_list = []
    # make lists of the dacay/res with survey if and also one for name vs survey id
    for rec in records:
        decay = rec[decay_col]
        sid = rec[id_col]
        if decay is not None:
            res = rec[manual_res_col]
            if res is None:
                res = rec[script_res_col]
                if res is None:
                    res = rec[script_point_res_col]
                    if res is None:
                        res = rec[manual_point_res_col]
                        if res is None:
                            print("missing res on record:", sid, rec[filename_col])
                            continue
            path = rec[manual_path_col]
            # A manual string can be an empty string (not Null) and also protect against it looking empty (just a space " ")
            if path is None or not path.strip():
                path = rec[path_col]
                if path is None or not path.strip():
                    print("skipping missing to_path", sid, rec[filename_col])
                    continue
            rec_list.append((sid, res, decay))
            # Switch to lower case, these were from filenames that I'm not sure are case sensitive
            names_list.append((rec[filename_col].lower(), sid, path))  # sid would be the next thing sorted if the names match
    # sort the names so we can use an integer to use for sorting by name
    names_list.sort()
    # do an ordered 2 key sort on decay then res (lexsort likes them backwards)
    rec_array = numpy.array(rec_list)
    sorted_indices = numpy.lexsort([-rec_array[:, 1], rec_array[:, 2]])  # resolution, decay (flip the res so lowest score and largest res is first)
    sorted_recs = rec_array[sorted_indices]
    sort_val = 0
    prev_res, prev_decay = None, None
    sort_dict = {}
    # set up a dictionary that has the sorted value of the decay followed by resolution
    for n, (sid, res, decay) in enumerate(sorted_recs):
        # don't incremenet when there was a tie, this allows the next sort criteria to be checked
        if res != prev_res or decay != prev_decay:
            sort_val += 1
        prev_res = res
        prev_decay = decay
        sort_dict[sid] = [sort_val]
    # the NBS sort order then uses depth after decay but before alphabetical, so we can't merge the name sort with the decay+res
    # add a second value for the alphabetical naming which is the last resort to maintain constistency of selection
    for n, (filename, sid, path) in enumerate(names_list):
        sort_dict[sid].append(n)
    return sorted_recs, names_list, sort_dict


def nbs_survey_sort(id_to_score, pts, existing_arrays, pts_col_offset=0, existing_col_offset=0):
    return nbs_sort_values(id_to_score, pts[LayersEnum.CONTRIBUTOR + pts_col_offset], pts[LayersEnum.ELEVATION + pts_col_offset],
                           existing_arrays[LayersEnum.CONTRIBUTOR+existing_col_offset], existing_arrays[LayersEnum.ELEVATION+existing_col_offset])


def nbs_sort_values(id_to_score, new_contrib, new_elev, accum_contrib, accum_elev):
    # return arrays that merge_arrays will use for sorting.
    # basically the nbs sort is 4 keys: Decay Score, resolution, depth, alphabetical.
    # Decay and resolution get merged into one array since they are true for all points of the survey while depth varies with position.
    # alphabetical is a final tie breaker to make sure the same contributor is picked in the cases where the first three tie.

    # find all the contributors to look up
    unique_contributors = numpy.unique(accum_contrib[~numpy.isnan(accum_contrib)])
    # make arrays to store the integer scores in
    existing_decay_and_res = accum_contrib.copy()
    existing_alphabetical = accum_contrib.copy()
    # for each unique contributor fill with the associated decay/resolution score and the alphabetical score
    for contrib in unique_contributors:
        existing_decay_and_res[accum_contrib == contrib] = id_to_score[contrib][0]
        existing_alphabetical[accum_contrib == contrib] = id_to_score[contrib][1]
    # @FIXME is contributor an int or float -- needs to be int 32 and maybe int 64 (or two int 32s)
    unique_pts_contributors = numpy.unique(new_contrib[~numpy.isnan(new_contrib)])
    decay_and_res_score = new_contrib.copy()
    alphabetical = new_contrib.copy()
    for contrib in unique_pts_contributors:
        decay_and_res_score[new_contrib == contrib] = id_to_score[contrib][0]
        alphabetical[new_contrib == contrib] = id_to_score[contrib][1]

    return numpy.array((decay_and_res_score, new_elev, alphabetical)), \
           numpy.array((existing_decay_and_res, accum_elev, existing_alphabetical)), \
           (False, False, False)




