"""
scrape_ehydro_file_meta_then_move_to_datumfolder_with_fuse2_saj.py

Brings in ehydro all files list, folder of district data:
Then reads metadata in all supporting ehydro files (.xyz headers, xml, and ehydro table):
this is saved out
and datums are used to move data into a new folder based on vertical datum

jkinney February 2019
update 3/1/19
"""

import os

import BDB.BDB51.preprocessing.Searching_files_and_directories as S_f_d
# import fuse.raw_read.usace.cesaj as e_meta
import fuse.raw_read.usace.cesam as cesam
import pandas as pd

try:
    import fuse.raw_read.usace.parse_usace_xml as p_usace_xml  # import parse_usace_xml_1 as p_usace_xml
except:
    print('check if the extract_xml classes have been merged yet see extract_ehydro_meta_class_CEMVN_.py')
import BDB.BDB51.preprocessing.meta2csv as m2c  # from BDB.BDB51.preprocessing import meta2csv as m2c


# ------------------------------------------------------------------
def district_meta(district):
    """
    

    Parameters
    ----------
    district :
        

    Returns
    -------

    """
    district = district.upper()
    if district == 'CESAM':
        import fuse.raw_read.usace.cesam as e_meta
    elif district == 'CEMVN':
        import fuse.raw_read.usace.cemvn as e_meta
    elif district == 'CESAJ':
        import fuse.raw_read.usace.cesaj as e_meta
    elif district == 'CENAN':
        import fuse.raw_read.usace.cenan as e_meta
    else:
        import fuse.raw_read.usace.cemvn as e_meta
    return e_meta


def ehydro_subset_only_ifnotfull_restoo(highresfolder, ehydrofolder):
    """
    

    Parameters
    ----------
    highresfolder :
        
    ehydrofolder :
        

    Returns
    -------

    """
    excludematchesfrom_dir = highresfolder
    directory_path = ehydrofolder
    g1 = S_f_d.exlude_overlap_list_search_xyz(excludematchesfrom_dir, directory_path, search_string='*.xyz',
                                              ex_string='*_A.xyz', splitter='.', splitter2='_A.')
    """
    search_string and ex_string expect '*.xyz' or some other wildcard string or an exact matching string input
    """
    g1.sort()
    return g1


def testonefile():
    """
    Testing
    FUSE module method
    single file
    export out dictionary

    Parameters
    ----------

    Returns
    -------

    """
    # import fuse.raw_read.usace.extract_ehydro_meta_class_CESAM
    infilename = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\downloads\CESAM\PH_01_PSB_20180710_CS.XYZ"
    rr = cesam.read_raw_cesam()
    merged_meta_test = rr.read_metadata(infilename)
    return merged_meta_test


def retrieve_meta_for_Ehydro_fuse_CESAM(highresfolder, ehydrofolder, metafile):
    """
    take input path for ehydro data
    read in metadata using the FUSE module
    and then write to the csv metadata file using Glen's script

    Parameters
    ----------
    highresfolder :
        
    ehydrofolder :
        
    metafile :
        

    Returns
    -------

    """
    # import fuse.raw_read.usace.extract_ehydro_meta_class_CESAM
    rr = cesam.read_raw_cesam()
    g1 = ehydro_subset_only_ifnotfull_restoo(highresfolder, ehydrofolder)
    merge2 = {}
    g1.sort()
    for basename in g1:
        merge2 = rr.read_metadata(basename)
        m2c.write_meta2csv([merge2], metafile)
    return merge2


def retrieve_meta_for_Ehydro_notable(highresfolder, ehydrofolder, district, metafile, df_export_to_csv):
    """
    

    Parameters
    ----------
    highresfolder :
        
    ehydrofolder :
        
    district :
        
    metafile :
        
    df_export_to_csv :
        

    Returns
    -------

    """
    g1 = ehydro_subset_only_ifnotfull_restoo(highresfolder, ehydrofolder)
    merged_meta = {}
    merge2 = {}
    nn = pd.DataFrame()
    g1.sort()
    ex_string1 = '*_A.xyz'
    ex_string2 = '*_FULL.xyz'
    ex_string3 = '*_FULL.XYZ'
    ex_string4 = '*_A.XYZ'

    for basename in g1:
        f = basename
        basename = os.path.basename(basename)
        basename = S_f_d.return_surveyid(basename, ex_string1)
        basename = S_f_d.return_surveyid(basename, ex_string2)
        basename = S_f_d.return_surveyid(basename, ex_string3)
        basename = S_f_d.return_surveyid(basename, ex_string4)
        basename = basename.rstrip('.xyz')
        basename = basename.rstrip('.XYZ')
        meta_from_ehydro = {}
        e_t = e_meta.Extract_Txt(f)

        # xml pull here. #since we know its ehydro:
        if '_A.xyz' in f:
            xmlfilename = e_meta.get_xml_xt(f, '_A.xyz')
        elif '_FULL.xyz' in f:
            xmlfilename = e_meta.get_xml_xt(f, '_FULL.xyz')
        elif '_FULL.XYZ' in f:
            xmlfilename = e_meta.get_xml_xt(f, '_FULL.XYZ')
        elif '_A.XYZ' in f:
            xmlfilename = e_meta.get_xml_xt(f, '_A.XYZ')
        else:
            xmlfilename = e_meta.get_xml(f)

        if os.path.isfile(xmlfilename):
            with open(xmlfilename, 'r') as xml_file:
                xml_txt = xml_file.read()
            xmlbasename = os.path.basename(xmlfilename)
            xml_data = p_usace_xml.XML_Meta(xml_txt, filename=xmlbasename)
            if xml_data.version == 'USACE_FGDC':
                meta_xml = xml_data._extract_meta_CEMVN()
            elif xml_data.version == 'ISO-8859-1':
                meta_xml = xml_data._extract_meta_USACE_ISO()
            else:
                meta_xml = xml_data.convert_xml_to_dict2()  # some_meta = xml_data.convert_xml_to_dict()
        else:
            meta_xml = {}
        meta = e_t.parse_ehydro_xyz(f, meta_source='xyz', version=district,
                                    default_meta='')  # parse_ehydro_xyz(infilename, meta_source = 'xyz', version='CEMVN', default_meta = '')#eH_p.parse_ehydro_xyz(f)

        list_keys_empty = []
        combined_row = {}
        subset_row = {}
        subset_no_overlap = {}

        for key in meta:
            if meta[key] in ('unknown', ''):
                list_keys_empty.append(key)
            else:
                subset_row[key] = meta[key]
                # non blank columns only
                if key in meta_xml:
                    if meta[key] == meta_xml[key]:
                        combined_row[key] = meta[key]
                        # only make list within cell if values from different sources are different
                    else:
                        combined_row[key] = f'{meta[key]} , {meta_xml[key]}'
                else:
                    subset_no_overlap[key] = meta[key]

        mydict = meta_xml

        for i0, key1 in enumerate(mydict):  # for key1 in mydict:
            key_fromdict = key1
            # use common numberic index i0
            if key_fromdict in combined_row:
                nn.loc[f, key_fromdict] = combined_row[key_fromdict]  # [row, column]
            else:
                nn.loc[f, key_fromdict] = mydict[key_fromdict]  # [row, column]

        for i0, key1 in enumerate(subset_no_overlap):
            key_fromdict = key1
            nn.loc[f, key_fromdict] = subset_row[key_fromdict]

        for i0, key1 in enumerate(meta_from_ehydro):
            key_fromdict = key1
            nn.loc[f, key_fromdict] = meta_from_ehydro[key_fromdict]

        merged_meta = {**meta, **meta_from_ehydro, **meta_xml}
        merge2 = {**subset_row, **meta_from_ehydro, **meta_xml, **combined_row}
        # THIS METHOD FOR MERGING DICTIONARIES IS FOR Python 3.5  plus based on PEP 448
        # https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/ 
        # see for more informatio`n
        # merged_dictonary = {**default_values, **override_if_duplicatekeys}#
        m2c.write_meta2csv([merge2],
                           metafile)  # could also try merged_meta, but merge2 outputs a list within a cell if inputs from text and xml for the same variable are different

    # save pandas dataframe export here
    nn.to_csv(path_or_buf=df_export_to_csv, encoding='UTF-8', sep='\t')
    return merged_meta, nn


def test_cesam_ccom_defaultnames():
    """
    TODO write description
    """

    district = 'CESAM'
    highresfolder = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\downloads\cesam_test"
    # CESAM"
    ehydrofolder = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\downloads\cesam_test"
    # CESAM"#r'N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\tst'
    mv_to_dir = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\test_move"
    # "N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\CESAM"#D:\NBS_Data\PBG_GulfCoast_PR_VI\Mississippi\USACE\E-Hydro\download\files\downloads\CEMVN"
    metafile = os.path.join(mv_to_dir, 'metadata',
                            'ehydro_meta_dict_out.txt')  # r"N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\ehydrometa_sm_set.txt"
    df_export_to_csv = os.path.join(mv_to_dir, 'metadata',
                                    'ehydro_allscript_meta_v1.txt')  # r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\Attempted_combined_df_metafields.txt'
    merged_meta1, merged_dataframe = retrieve_meta_for_Ehydro_notable(highresfolder, ehydrofolder, district, metafile,
                                                                      df_export_to_csv)
    move_to_vert_datum_folder_iter(merged_dataframe, mv_to_dir)


# move to a new folder
def move_to_vert_datum_folder_iter(merged_dataframe, mv_to_dir):
    """
    

    Parameters
    ----------
    merged_dataframe :
        
    mv_to_dir :
        

    Returns
    -------

    """
    Ax = merged_dataframe.axes
    Ax_row_index = Ax[0].tolist()
    Ax_row_index.sort()
    for filenamepath in Ax_row_index:
        # 'from_filenamebase'# is in dataframe
        if merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'MLLW':
            # works
            # move to MLLW folder
            S_f_d.look_for_Ehydro_datum_folders(filenamepath, 'MLLW', mv_to_dir)
            print('moving file')
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'MLG':
            # MLG
            S_f_d.look_for_Ehydro_datum_folders(filenamepath, 'MLG', mv_to_dir)
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'LWRP':
            # LWRP
            S_f_d.look_for_Ehydro_datum_folders(filenamepath, 'LWRP', mv_to_dir)
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'NGVD29':
            # NGVD29
            S_f_d.look_for_Ehydro_datum_folders(filenamepath, 'NGVD29', mv_to_dir)
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'MLW':
            # MLW
            S_f_d.look_for_Ehydro_datum_folders(filenamepath, 'MLW', mv_to_dir)
        else:
            if merged_dataframe.at[filenamepath, 'script: from_vert_key'] != "":
                # Vertical datum unknown
                S_f_d.look_for_Ehydro_datum_folders(filenamepath,
                                                    merged_dataframe.at[filenamepath, 'script: from_vert_key'],
                                                    mv_to_dir)
            else:
                S_f_d.look_for_Ehydro_datum_folders(filenamepath)
    print('moving of files complete!')


def main():
    """ """
    test_cesam_ccom_defaultnames()

    """
    UnboundLocalError: local variable 'AllBPSList' referenced before assignment
    occurs when file folder passed has no .xyz files
    """


if __name__ == '__main__':
    main()
