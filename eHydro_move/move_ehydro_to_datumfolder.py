"""
move_ehydro_to_datumfolder.py

Brings in metadata from ehydro all files list, 
given a path for the folder of district data
datums are used to move data into a new folder based on vertical datum

J. Kinney

"""

import os
import shutil
from glob import glob as glob

import pandas as pd


# ------------------------------------------------------------------------------


def test_from_existing_file():
    """
    take an existing metadata .csv file that has all the metadata
    For internal development and testing at this point.

    Parameters
    ----------

    Returns
    -------

    """
    mv_to_dir = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\CEMVN"
    df_export_to_csv = os.path.join(mv_to_dir, 'metadata',
                                    'ehydro_allscript_meta_v1.txt')  # r'N:\New_Directory_1\GulfCoast\USACE\xyz\MLLW\Metadata\Active\Attempted_combined_df_metafields.txt'
    merged_dataframe = pd.read_csv(df_export_to_csv, sep='\t', index_col=0)
    # metafile = os.path.join(mv_to_dir, 'metadata', 'ehydro_meta_dict_out.txt')
    # alternative
    # merged_dataframe = pd.read_csv(metafile, sep = ",", index_col = 1)
    move_to_vert_datum_folder_iter(merged_dataframe, mv_to_dir)


def test_from_existing_file_meta2csvstyle():
    """take an existing metadata .csv file that comes as output from meta2csv"""
    mv_to_dir = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\CEMVN"
    metafile = mv_to_dir + "\metadata\ehydro_meta_dict_out.txt"
    merged_dataframe = pd.read_csv(metafile, sep=",", index_col=1)
    move_to_vert_datum_folder_iter(merged_dataframe, mv_to_dir)


# ------------------------------------------------------------------------------
def look_for_Ehydro_datum_folders(f, datumfolder='unknown', newpathroot=None):
    """
    look_for_Ehydro_datum_folders(f, datumfolder = 'unknown')
    example:look_for_Ehydro_datum_folders(f, 'MLLW')
    moves file to a new directory corresponding to vertical datum

    Parameters
    ----------
    f :
        
    datumfolder :
         (Default value = 'unknown')
    newpathroot :
         (Default value = None)

    Returns
    -------

    """
    if datumfolder == 'unknown' or datumfolder == '':
        d_folder = '\\unknown'
    else:
        d_folder = f'\\{datumfolder}'

    if os.path.exists(f):
        filepath = os.path.dirname(f)
        basename = os.path.basename(f)
        ex_string1 = '*_A.xyz'
        ex_string2 = '*_FULL.xyz'
        ex_string3 = '*_FULL.XYZ'
        ex_string4 = '*_A.XYZ'
        basename = return_surveyid(basename, ex_string1)
        basename = return_surveyid(basename, ex_string2)
        basename = return_surveyid(basename, ex_string3)
        basename = return_surveyid(basename, ex_string4)
        basename = basename.rstrip('.xyz')
        basename = basename.rstrip('.XYZ')
        survey_files = glob(os.path.join(filepath, f'{basename}*'))
        print(survey_files)

        if newpathroot is not None:
            if os.path.exists(newpathroot):
                filedatumpath = os.path.join(newpathroot, d_folder, 'Original')
                # redirect to a new folder path

                if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                    # move the file here
                    for ff in survey_files:
                        shutil.move(ff, filedatumpath)  #
                        print(f'moved to {filedatumpath}')
                    print('did it work?')
                else:
                    # make directory #os.mkdir(r'N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\the_vertical_datum_folder')
                    os.makedirs(filedatumpath)  # os.mkdir(filedatumpath)# add one leaf of folder structure
                    # use os.makedirs(filedatumpath) if one is adding multiple new folder layers

                    if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                        for ff in survey_files:
                            shutil.move(ff, filedatumpath)  # shutil.move(f,filedatumpath)
                            print(f'moved {ff} to datum folder: {datumfolder}')
        elif os.path.dirname(f).find('USACE\ehydro') >= 0 or os.path.dirname(f).find(
                'USACE\E-Hydro') >= 0 or os.path.dirname(f).find('USACE\eHydro' >= 0):
            filepath = os.path.dirname(f)
            filedatumpath = os.path.join(filepath, d_folder, 'Original')  # '\MLLW'
            if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                # move the file here
                for ff in survey_files:
                    shutil.move(ff, filedatumpath)  #
                print(f'moved to {filedatumpath}')
            else:
                # make directory #example os.mkdir(r'N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\the_vertical_datum_folder')
                os.makedirs(filedatumpath)  # os.mkdir(filedatumpath)# add one leaf of folder structure
                if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                    for ff in survey_files:
                        shutil.move(ff, filedatumpath)  # shutil.move(f,filedatumpath)
                        print(f'moved {ff} to datum folder: {datumfolder}')


# ------------------------------------------------------------------------------

def return_surveyid(filenamepath, ex_string):
    """
    

    Parameters
    ----------
    filenamepath :
        
    ex_string :
        

    Returns
    -------

    """

    basename = os.path.basename(filenamepath)
    surveybasename = basename.rstrip(ex_string)
    return surveybasename


# ------------------------------------------------------------------------------
def move_to_vert_datum_folder_iter(merged_dataframe, mv_to_dir):
    """
    move all files with matching extension to new folder based on the datum

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
            # move to MLLW folder
            look_for_Ehydro_datum_folders(filenamepath, 'MLLW', mv_to_dir)
            print('moving file')
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'MLG':
            # MLG
            look_for_Ehydro_datum_folders(filenamepath, 'MLG', mv_to_dir)
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'LWRP':
            # LWRP
            look_for_Ehydro_datum_folders(filenamepath, 'LWRP', mv_to_dir)
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'NGVD29':
            # NGVD29
            look_for_Ehydro_datum_folders(filenamepath, 'NGVD29', mv_to_dir)
        elif merged_dataframe.at[filenamepath, 'script: from_vert_key'] == 'MLW':
            # MLW
            look_for_Ehydro_datum_folders(filenamepath, 'MLW', mv_to_dir)
        else:
            if merged_dataframe.at[filenamepath, 'script: from_vert_key'] != "":
                # Vertical datum unknown
                look_for_Ehydro_datum_folders(filenamepath, merged_dataframe.at[filenamepath, 'script: from_vert_key'],
                                              mv_to_dir)
            else:
                look_for_Ehydro_datum_folders(filenamepath)
    print('moving of files complete!')


def main():
    """ """
    # test_from_existing_file()
    test_from_existing_file_meta2csvstyle()


if __name__ == '__main__':
    main()
