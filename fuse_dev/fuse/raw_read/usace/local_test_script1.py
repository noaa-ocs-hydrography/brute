# -*- coding: utf-8 -*-
"""
local_test_script1.py

Created on Wed Apr  3 16:05:18 2019

@author: jkinney
"""
#import extrac_ehydro_meta_class_CEMVN
import cesam
import cesaj
import ceswg
import cemvn
infilename = r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESWG\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.XYZ"

infilename=r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\test_cesaj\CH_01_CAH_20181025_CS_2018_270_01_A.XYZ"
infilename=r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CEMVN\SW_13_SWP_20190627_CS.xyz"
filename= r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESAJ\IW_11_DAD_20190613_CS_2019_184_01.XYZ"#"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESAJ\IW_11_DAD_20190613_CS_2019_184_01_A.XYZ"
infilename = r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESWG\BT_02_ENT_20190607_CS_46_MLLW_PRE_BD_06_3300_6540.XYZ"

#testing why vdatum output zone error reading from vdatum log
infilename = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\CEMVN\pickle_test\MLLW\Original\test_oddfile1\SW_09_SWP_20190604_CS_FORUM.xyz"
#testing SPCS flag pass CESAM
infilename=r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\CESAM\pickle_test\MLLW\Original\BL_03_SKB_20181023_CS.XYZ"

#CESAJ# r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESAJ\TH_04_TPH_20181025_CS_2019_044_02_A.XYZ"
#r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESAJ\SJ_01_SJH_20190506_CS_2019_146_01.XYZ"
#filewith pickle
#
#infilename = r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\CESAJ\previous_mv\MLLW\Original\PB_01_PBH_20170913_CS_2017_201.XYZ"
#r"N:\New_Directory_1\GulfCoast\USACE\ehydro\EasternGulf\downloads\CESAM\BW_04_HLT_20170320_CS_OLD_LOCK_14_BAR.XYZ"
#r"C:\Users\jkinney\Desktop\testing_etc\AW_01_NSJ_20000501_CS_1999_314_01.XYZ"#'N:\New_Directory_1\GulfCoast\USACE\ehydro\CEMVN\LWRP\Original\MR_48_FRV_20190125_CS.xyz'
rr = CESAJRawReader()#read_raw()
rr = CEMVNRawReader()
rr = CESAMRawReader()
rr= CESWGRawReader()

#rr=cesam.read_raw()
#rr =read_raw_cemvn()
merged_metaECSAJ_A_1 =rr.read_metadata(infilename)
merged_meta_CESWG1e =rr.read_metadata(infilename)
merged_mCEMVN_1 =rr.read_metadata(infilename)
#merged_meta57, meta_df57 =rr.read_metadata(infilename)
merged_mergedCESAM_1b =rr.read_metadata(infilename)

xyz=rr.read_bathymetry(infilename)
numpyout=r"N:\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt"
import numpy as _np
_np.savetxt(numpyout, xyz, delimiter=',')
 
 
#command line
"""
java -jar vdatum.jar ihorz:NAD83:spc:us_ft:102 ivert:ngvd29:us_ft:height ohorz:NAD83:utm:m: overt:mllw:m:height -file:txt:comma,0,1,2,skip0:N:\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt;N:\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\temp

#lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:NAD83:spc:us_ft:102 ivert:ngvd29:us_ft:height ohorz:NAD83:utm:m: overt:mllw:m:height -file:txt:comma,0,1,2,skip0:N:\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt;N:\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\temp_vdatum4
lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:NAD83:spc:us_ft:102 ivert:ngvd29:us_ft:height ohorz:NAD83:utm:m: overt:mllw:m:height -file:txt:comma,0,1,2,skip0:N:\\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt;N:\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\temp_vdatum4

lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:NAD83:spc:us_ft:102 ivert:ngvd29:us_ft:height ohorz:NAD83:utm:m: overt:mllw:m:height -file:txt:comma,0,1,2,skip0:N:\\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt;N:/\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\temp_vdatum4 region:3
lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:nad83:spc:us_ft:102 ivert:ngvd29:us_ft:height ohorz:NAD83:utm:m: overt:mllw:m:height -file:txt:comma,0,1,2,skip0:N:\\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt;N:/\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\temp_vdatum4 region:3

lib\java_home\openjdk-11.0.2\bin\java -jar vdatum.jar ihorz:NAD83_2011:spc:us_ft:102 ivert:ngvd29:us_ft:height ohorz:NAD83_2011:utm:m: overt:mllw:m:height -file:txt:comma,0,1,2,skip0:N:\\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\BH_12_BHA_20181022_CS_37_MLLW_2849_1000.txt;N:/\New_Directory_1\GulfCoast\USACE\ehydro\Galveston\misc_test\temp_vdatum4 region:3

"""


#tabdelimited file
infilename=r"N:\New_Directory_1\GulfCoast\USACE\temp_ehydro\updated_downloadwithpickle\downloads\CESAJ\MH_01_MIH_20180730_CS_2018_219_04_A.XYZ"

keys1 = ['from_filename',
                  'start_date',
                  #'end_date',
                  #'horiz_uncert',
                  #'agency',
                  'source_indicator',
                  'from_vert_datum',
                  'from_vert_key',
                  'from_vert_units',
                  'from_horiz_units',
                  'from_fips',
                  'source_indicator',
                  # from_vert_unc
                  'script_version' ,
                  ]

for key in keys1:
    print(merged_meta_test[key])
 
To_inputs = []                 
                  #'to_vert_datum',
                  #'to_vert_units'    
                  #'to_horiz_datum',    
later_step = [       
                  'source_type',
                  'complete_coverage',
                  'complete_bathymetry',
                  'vert_uncert_fixed' ,
                  'vert_uncert_vari',
                  'feat_size',
                  'feat_detect',
                  'feat_least_depth' ,
                  'interpolated',
                  'reviewed',
                  ]

_field_map = {'from_filename' : 'OBJNAM',
                  'start_date' : 'SURSTA',
                  'end_date' : 'SUREND',
                  'horiz_uncert' : 'POSACC',
                  'to_horiz_datum' : 'HORDAT',
                  'agency' : 'AGENCY',
                  'source_indicator' : 'SORIND',
                  'source_type' : 's_ftyp',
                  'complete_coverage' : 'flcvrg',
                  'complete_bathymetry' : 'flbath',
                  'to_vert_datum' : 'VERDAT',
                  'to_vert_units' : 'DUNITS',
                  'vert_uncert_fixed' : 'vun_fx',
                  'vert_uncert_vari' : 'vun_vb',
                  'feat_size' : 'f_size',
                  'feat_detect' : 'f_dtct',
                  'feat_least_depth' : 'f_lstd',
                  'interpolated' : 'interp',
                  'reviewed' : 'r_name',
                  'script_version' : 's_scpv',
                  'source_indicator' : 'SORIND',}
