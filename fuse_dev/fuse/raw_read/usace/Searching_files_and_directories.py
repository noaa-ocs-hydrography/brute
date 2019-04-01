# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 13:09:55 2019
Searching_files_and_directories.py
@author: jkinney
"""
import os  
from glob import glob as glob
import shutil

def check_afile_against_db(self, db, spreadsheet1):
    if self.version == 'BAG':
        #do something
        filepathname=self.filename
        bag_not_in_db, BOID = check_against_bdb_1file(filepathname, db, spreadsheet1)
        #if bag_not_in_db is not []:
        #   do something
        """
        returns value if, list of files that have yet to be loaded into BDB
        should just be 1 filename, but may be multiple BOID, since it is looking for just 1 file match
        """
        #return bag_not_in_db, BOID#
        pass

def check_against_db(Bathy_filelist, db, spreadsheet1, version = None):
    if version == None:
        version == 'BAG'
    if version == 'BAG':
        bag_not_in_db = check_against_bdb(Bathy_filelist, db, spreadsheet1)
        """
        returns list of files that have yet to be loaded into BDB
        """
        #bag_list=Bathy_filelist 
    return bag_not_in_db


def check_against_table(Bathy_filelist, table1, version ):
    #matches
    #not_a_match
    pass

def check_against_csv(Bathy_filelist, csv_table1, version ):
    pass

#------------------------------------------------------------------------------
def look_for_Ehydro_datum_folders(f, datumfolder = 'unknown', newpathroot = None):
    """
    look_for_Ehydro_datum_folders(f, datumfolder = 'unknown')
    example:look_for_Ehydro_datum_folders(f, 'MLLW')
    moves file to a new directory corresponding to vertical datum
    
    """

    if datumfolder == 'unknown' or datumfolder == '':
        d_folder = '\\unknown'
    else:
        d_folder = '\\' + str(datumfolder)
            
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
        basename=basename.rstrip('.xyz')
        basename = basename.rstrip('.XYZ')  
        survey_files = glob(os.path.join(filepath, basename + '*'))
        print(survey_files)
        if newpathroot is not None:
            if os.path.exists(newpathroot):
                filedatumpath = newpathroot + d_folder + '\\Original'
                #redirect to a new folder path
                if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                    #move the file here
                    for ff in survey_files:
                        shutil.move(ff,filedatumpath)#shutil.move(f,filedatumpath)
                        print('moved to ' + filedatumpath)
                    print('did it work?')
                else:
                    #make directory
                    #os.mkdir(r'N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\the_vertical_datum_folder')
                    os.makedirs(filedatumpath)#os.mkdir(filedatumpath)# add one leaf of folder structure
                    #use os.makedirs(filedatumpath) if one is adding multiple new folder layers
                    if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                        for ff in survey_files:
                            shutil.move(ff,filedatumpath)#shutil.move(f,filedatumpath)
                            print('moved ' + str(ff) + ' to datum folder: ' + str(datumfolder))
        elif os.path.dirname(f).find('USACE\ehydro') >= 0 or os.path.dirname(f).find('USACE\E-Hydro') >= 0  or os.path.dirname(f).find('USACE\eHydro' >= 0 ):
            filepath = os.path.dirname(f)
            filedatumpath = filepath + d_folder + '\\Original'#'\MLLW'
            if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                #move the file here
                for ff in survey_files:
                        shutil.move(ff,filedatumpath)#shutil.move(f,filedatumpath)
                print('moved to ' + filedatumpath)
            else:
                #make directory
                #os.mkdir(r'N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\the_vertical_datum_folder')
                os.makedirs(filedatumpath)#os.mkdir(filedatumpath)# add one leaf of folder structure
                if os.path.isdir(filedatumpath) and os.path.exists(filedatumpath):
                    for ff in survey_files:
                        shutil.move(ff,filedatumpath)#shutil.move(f,filedatumpath)
                        print('moved ' + str(ff) + ' to datum folder: ' + str(datumfolder))
            #test_file=r"N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\MD_56_NO3_20161110\MD_56_NO3_20161110.gdb"
            #test_folder_tomoveto = r'N:\New_Directory_1\GulfCoast\Mississippi\USACE\ehydro\the_vertical_datum_folder'
            #shutil.move(test_file, test_folder_tomoveto)
#------------------------------------------------------------------------------    
            
def bps_bag_search(bagdir, directory_path, search_string):
    """
    bps_bag_search(db, bagdir, directory_path):
        looking at grids and points in .csar files vs. bags
    returns g1, g2 # g3    
    #list of surveys in the case g1 should be points, and g2 the first grid
    This script does a search for all bag survey numbers and then examines
    the list of all BPS survey numbers based on the input files directories
    as passed to the script.  A list of BPS files that exclude surveys 
    existing as .bags to be pulled in is then passed to the next script
    """
    if search_string == None:
        search_string = '*_Active*_NAD83.csar'
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        ### get a list of all the csar's in the directory
        BPSfiles = glob(os.path.join(directory_path, search_string))#   BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_*_NAD83.csar'))#     BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_5m_UTM18N_NAD83.csar'))
        #looking for ACTIVE .csars        
        BPSfiles.sort()
        numbpsfiles = len(BPSfiles)
        # initialize a counter for tracking the successful files
        success_count = 0
        SurveyListB = []
        AllBPSList = []
        for n,BPSfiles1 in enumerate(BPSfiles):
            File_NameBPS = os.path.basename(BPSfiles1)#doesn't make sense later just for testing
            bname1= File_NameBPS.split('_')#use '.' if its just the .csv from the BPS store
            Sname1 = bname1[0]# first column that should equal survey name
            AllBPSList.append(File_NameBPS)
            SurveyListB.append(Sname1)
        #BPSfiles 
    SetBPS = set(AllBPSList)
    if os.path.exists(bagdir) and os.path.isdir(bagdir):
    ### get a list of all the BAGs in the directory for BAGS
        Bagfiles = glob(os.path.join(bagdir,'*.bag'))
        Bagfiles.sort()
        numfiles = len(Bagfiles)
        success_count = 0
        SurveyList = []
        #AllBagList = []       
        #B = []
        for n,Bagfilelist1 in enumerate(Bagfiles):
            File_NameBag=os.path.basename(Bagfilelist1)                
            bname= File_NameBag.split('_')#use '.' if its just the .csv from the BPS store instead of CSARs that we are searching
            Sname = bname[0]# first column that should equal survey name
            SurveyList.append(Sname)
            #AllBagList.append(File_NameBag)            
    SetBPSnames = set(SurveyListB)
    SetBAGsurveynames = set(SurveyList)
    len(SetBPSnames)
    len(SetBAGsurveynames)
    Setnot1 = SetBPSnames.difference(SetBAGsurveynames)
    len(Setnot1)
    SubsetNotBags = list(Setnot1)
    ListBPSnotbags = []    
    for S1 in SubsetNotBags:
        S2 = os.path.join(directory_path, S1)
        S2 = (os.path.join(directory_path, S1) + '*.csar')
        BPSfilesexcludesbags = glob(S2)
        ListBPSnotbags.append(BPSfilesexcludesbags)    
    g = [item for item in ListBPSnotbags if item != []]#remove empty lines in list
    g1 = [item[0] for item in g]#just points
    try:
        g2 = [item[1] for item in g]#first grid size
    except:
        print('not a second item')
        g2 = []
    try:
        g3 = [item[2] for item in g]#second grid size         
    except:
        print('not a third item')
    #will probaby want to adjust which files we want to pull in exactly
    print('outputting a list')
    #print(g1)
    return  g1, g2 # g3    

def bps_bag_search_xyz( bagdir, directory_path, search_string = None):
    """
    g1 = bps(bagdir, directory_path, search_string)
    
    bps_bag_search(db, bagdir, directory_path):
    returns g1, g2 # g3    
    #list of surveys in the case g1 should be points, and g2 the first grid
    This script does a search for all bag survey numbers and then examines
    the list of all BPS survey numbers based on the input files directories
    as passed to the script.  A list of BPS files that exclude surveys 
    existing as .bags to be pulled in is then passed to the next script
    """
    if search_string == None:
        search_string = '*.xyz'
        search_string_ext = search_string
    elif search_string == '*.xyz':
        search_string_ext = search_string
        
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        ### get a list of all the csar's in the directory
        BPSfiles = glob(os.path.join(directory_path, search_string))#   BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_*_NAD83.csar'))#     BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_5m_UTM18N_NAD83.csar'))
        #looking for ACTIVE .csars        
        BPSfiles.sort()
        numbpsfiles = len(BPSfiles)
        # initialize a counter for tracking the successful files
        success_count = 0
        SurveyListB = []
        AllBPSList = []
        for n,BPSfiles1 in enumerate(BPSfiles):
            File_NameBPS = os.path.basename(BPSfiles1)#doesn't make sense later just for testing
            bname1= File_NameBPS.split('.')#use '.' if its just the .csv from the BPS store
            Sname1 = bname1[0]# first column that should equal survey name
            AllBPSList.append(File_NameBPS)
            SurveyListB.append(Sname1)
        #BPSfiles 
    SetBPS = set(AllBPSList)
    if os.path.exists(bagdir) and os.path.isdir(bagdir):
    ### get a list of all the BAGs in the directory for BAGS
        Bagfiles = glob(os.path.join(bagdir,'*.bag'))
        Bagfiles.sort()
        numfiles = len(Bagfiles)
        success_count = 0
        SurveyList = []
        #AllBagList = []       
        #B = []
        for n,Bagfilelist1 in enumerate(Bagfiles):
            File_NameBag=os.path.basename(Bagfilelist1)                
            bname= File_NameBag.split('_')#use '.' if its just the .csv from the BPS store instead of CSARs that we are searching
            Sname = bname[0]# first column that should equal survey name
            SurveyList.append(Sname)
            #AllBagList.append(File_NameBag)            
    SetBPSnames = set(SurveyListB)
    SetBAGsurveynames = set(SurveyList)
    len(SetBPSnames)
    len(SetBAGsurveynames)
    Setnot1 = SetBPSnames.difference(SetBAGsurveynames)
    len(Setnot1)
    SubsetNotBags = list(Setnot1)
    ListBPSnotbags = []    
    for S1 in SubsetNotBags:
        S2 = os.path.join(directory_path, S1)
        S2 = (os.path.join(directory_path, S1) + search_string_ext)
        BPSfilesexcludesbags = glob(S2)
        ListBPSnotbags.append(BPSfilesexcludesbags)    
    g = [item for item in ListBPSnotbags if item != []]#remove empty lines in list
    g1 = [item[0] for item in g]#just points
#    try:
#        g2 = [item[1] for item in g]#first grid size
#    except:
#        print('not a second item')
#        g2 = []
#    try:
#        g3 = [item[2] for item in g]#second grid size         
#    except:
#        print('not a third item')
    #will probaby want to adjust which files we want to pull in exactly
    #designed for multiple types of .csars or other files with the same start.
    print('outputting a list')
    #print(g1)
    return  g1#, g2 # g3

            
#def bps_bag_search(bagdir, directory_path, search_string):
#    """
#    bps_bag_search(db, bagdir, directory_path):
#    returns g1, g2 # g3    
#    #list of surveys in the case g1 should be points, and g2 the first grid
#    This script does a search for all bag survey numbers and then examines
#    the list of all BPS survey numbers based on the input files directories
#    as passed to the script.  A list of BPS files that exclude surveys 
#    existing as .bags to be pulled in is then passed to the next script
#    """
#    if search_string == None:
#        search_string = '*_Active*_NAD83.csar'
#    if os.path.exists(directory_path) and os.path.isdir(directory_path):
#        ### get a list of all the csar's in the directory
#        BPSfiles = glob(os.path.join(directory_path, search_string))#   BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_*_NAD83.csar'))#     BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_5m_UTM18N_NAD83.csar'))
#        #looking for ACTIVE .csars        
#        BPSfiles.sort()
#        numbpsfiles = len(BPSfiles)
#        # initialize a counter for tracking the successful files
#        success_count = 0
#        SurveyListB = []
#        AllBPSList = []
#        for n,BPSfiles1 in enumerate(BPSfiles):
#            File_NameBPS = os.path.basename(BPSfiles1)#doesn't make sense later just for testing
#            bname1= File_NameBPS.split('_')#use '.' if its just the .csv from the BPS store
#            Sname1 = bname1[0]# first column that should equal survey name
#            AllBPSList.append(File_NameBPS)
#            SurveyListB.append(Sname1)
#        #BPSfiles 
#    SetBPS = set(AllBPSList)
#    if os.path.exists(bagdir) and os.path.isdir(bagdir):
#    ### get a list of all the BAGs in the directory for BAGS
#        Bagfiles = glob(os.path.join(bagdir,'*.bag'))
#        Bagfiles.sort()
#        numfiles = len(Bagfiles)
#        success_count = 0
#        SurveyList = []
#        #AllBagList = []       
#        #B = []
#        for n,Bagfilelist1 in enumerate(Bagfiles):
#            File_NameBag=os.path.basename(Bagfilelist1)                
#            bname= File_NameBag.split('_')#use '.' if its just the .csv from the BPS store instead of CSARs that we are searching
#            Sname = bname[0]# first column that should equal survey name
#            SurveyList.append(Sname)
#            #AllBagList.append(File_NameBag)            
#    SetBPSnames = set(SurveyListB)
#    SetBAGsurveynames = set(SurveyList)
#    len(SetBPSnames)
#    len(SetBAGsurveynames)
#    Setnot1 = SetBPSnames.difference(SetBAGsurveynames)
#    len(Setnot1)
#    SubsetNotBags = list(Setnot1)
#    ListBPSnotbags = []    
#    for S1 in SubsetNotBags:
#        S2 = os.path.join(directory_path, S1)
#        S2 = (os.path.join(directory_path, S1) + '*.csar')
#        BPSfilesexcludesbags = glob(S2)
#        ListBPSnotbags.append(BPSfilesexcludesbags)    
#    g = [item for item in ListBPSnotbags if item != []]#remove empty lines in list
#    g1 = [item[0] for item in g]#just points
#    try:
#        g2 = [item[1] for item in g]#first grid size
#    except:
#        print('not a second item')
#        g2 = []
#    try:
#        g3 = [item[2] for item in g]#second grid size         
#    except:
#        print('not a third item')
#    #will probaby want to adjust which files we want to pull in exactly
#    print('outputting a list')
#    #print(g1)
#    return  g1, g2 # g3    

def exlude_overlap_list_search_xyz( excludematchesfrom_dir, directory_path, search_string = None, ex_string = None, splitter = None, splitter2 = None):
    """
    g1 = exlude_overlap_list_search_xyz(xcludematchesfrom_dir, directory_path, search_string = None, ex_string = None, excludematchesfrom_dir, directory_path, search_string = None, ex_string = None)
    
    Set up to look for name matches before end extension "."
    
       returns g1, g2 # g3    #option if one expects multiple matches with slightly different names.
       
    #list of surveys in the case g1 should be points, and g2 the first grid
    This script does a search for all bag survey numbers and then examines
    the list of all BPS survey numbers based on the input files directories
    as passed to the script.  A list of BPS files that exclude surveys 
    existing as .bags to be pulled in is then passed to the next script
    
    example:
        (location of EHydro xyz files, and _A.xyz higher resolution files)
    set exludematchesfrom_dir = 
    directory_path = 
        files could be in the same directory and one would look to find only those xyz files not in the 'full or high resolution list'
        
    g1 = exclude_overlap_list_search_xyz(excludematchesfrom_dir, directory_path, search_string = .xyz, ex_string = _A.xyz, splitter = '.', splitter2 = '_A.')
    """
    if search_string == None:
        search_string = '*.xyz'
        search_string_ext = search_string
    elif search_string == '*.xyz':
        search_string_ext = search_string
    if ex_string == None:
        ex_string = '*.xyz'#exclude matches from this directories list with this extension
    if splitter == None:
        splitter == '.'
        splitter1 = splitter
    elif type(splitter) == list:
        splitter1 = splitter[0]
        if len(splitter)>= 2:
            splitter2 = splitter[1]
    elif type(splitter) == str:
        splitter1 = splitter        
    else:
        print('problem with splitter type')
    if splitter2 == None:
        splitter2 = splitter
            
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        ### get a list of all the csar's in the directory
        bathy_files = glob(os.path.join(directory_path, search_string))#   bathy_files = glob(os.path.join(directory_path,'*_SHOAL_*_NAD83.csar'))#     BPSfiles = glob(os.path.join(directory_path,'*_SHOAL_5m_UTM18N_NAD83.csar'))
        #looking for ACTIVE .csars        
        bathy_files.sort()
        numbathy_files = len(bathy_files)
        # initialize a counter for tracking the successful files
        success_count = 0
        SurveyListB = []
        AllBPSList = []
        for n,bathy_files1 in enumerate(bathy_files):
            File_NameBPS = os.path.basename(bathy_files1)#doesn't make sense later just for testing
            bname1= File_NameBPS.split(splitter1)#use '.' if its just the .csv from the BPS store
            Sname1 = bname1[0]# first column that should equal survey name
            AllBPSList.append(File_NameBPS)
            SurveyListB.append(Sname1)
        #bathy_files 
    SetBPS = set(AllBPSList)
    if os.path.exists(excludematchesfrom_dir) and os.path.isdir(excludematchesfrom_dir):
    ### get a list of all the BAGs in the directory for BAGS
        Bagfiles = glob(os.path.join(excludematchesfrom_dir, ex_string ))
        Bagfiles.sort()
        numfiles = len(Bagfiles)
        success_count = 0
        SurveyList = []
        #AllBagList = []       
        #B = []
        for n,Bagfilelist1 in enumerate(Bagfiles):
            File_NameBag=os.path.basename(Bagfilelist1)                
            bname= File_NameBag.split(splitter2)#use '.' if its just the .csv from the BPS store instead of CSARs that we are searching
            Sname = bname[0]# first column that should equal survey name
            SurveyList.append(Sname)
            #AllBagList.append(File_NameBag)            
    SetBPSnames = set(SurveyListB)
    SetBAGsurveynames = set(SurveyList)
    len(SetBPSnames)
    len(SetBAGsurveynames)
    Setnot1 = SetBPSnames.difference(SetBAGsurveynames)
    len(Setnot1)
    SubsetNotBags = list(Setnot1)
    ListBPSnotbags = []    
    for S1 in SubsetNotBags:
        S2 = os.path.join(directory_path, S1)
        S2 = (os.path.join(directory_path, S1) + search_string_ext)
        bathy_filesexcludesbags = glob(S2)
        ListBPSnotbags.append(bathy_filesexcludesbags)    
    g = [item for item in ListBPSnotbags if item != []]#remove empty lines in list
    g1 = [item[0] for item in g]#just points
    try:
        g2 = [item[1] for item in g]#first grid size
    except:
        print('not a second item')
        g2 = []
    try:
        g3 = [item[2] for item in g]#second grid size         
    except:
        print('not a third item')
    #will probaby want to adjust which files we want to pull in exactly
    #designed for multiple types of .csars or other files with the same start.
    print('outputting a list')
    #print(g1)
    return  g1 #, surveylist, AllBPSList#, g2 # g3            
#------------------------------------------------------------------------------

def return_surveyid(filenamepath, ex_string):
    
    basename = os.path.basename(filenamepath)
    surveybasename = basename.rstrip(ex_string)
    return surveybasename
    
#------------------------------------------------------------------------------
import os    

def check_against_bdb(bag_list, db, spreadsheet1):#Bathy_filelist
    
    BOID = []
    bag_not_in_db = []
    bag_list.sort()    
    for b in bag_list:
        bag_name = os.path.basename(b)
        File_Name1 = bag_name.lower()
        #TRY CQL QUERY APPROACH        
        try:
            cql_1 = "CoverageSource = '" + File_Name1 +"'"#method works #Did not work :"'CoverageSource' = " + File_Name1
            feature_listA = db.query('surfac', cql_1)
            for Bagfeature in feature_listA:
                #print(Bagfeature.id)#print feature id
                BOID = Bagfeature.id
                # #condition = bdb.QueryCondition()#condition = QueryCondition()#changed for 5.1
        except:
            print ('test condition issue')

        if len(BOID) == 0:
            bag_not_in_db.append(b)
    return bag_not_in_db


def check_against_bdb_1file(Bathy_list1, db, spreadsheet1):
    bag_not_in_db = []
    BOID = []
    BOIDLIST = []        
    File_NameBPS=os.path.basename(Bathy_list1)##names of BPS files from folder
    #print(File_NameBPS)
    File_Name1 = File_NameBPS.lower()#string needs to be all lower case in order to match to coverage source
    print(File_Name1)
    print(File_Name1.lower())    
    # convert the path names to a filename list
    cql_1 = "CoverageSource = '" + File_Name1 +"'"#method works #Did not work :"'CoverageSource' = " + File_Name1
    feature_listA = db.query('surfac', cql_1)
    for Bagfeature in feature_listA:#in case the file has multiple entries in the database (this seems to happen mor than we'd want at this point, especiall if load interrupted )
        #print(Bagfeature.id)#print feature id
        BOID = Bagfeature.id
        if len(BOID)== 0:
            bag_not_in_db.append(Bagfeature)
            BOIDLIST.append(BOID)
    return bag_not_in_db, BOID        
        
#------------------------------------------------------------------------------    
#https://geospatial-usace.opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0        
#https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.csv
#https://opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1_0.geojson

