build_caris_env.bat must be run with conda in the path.  For example, cd to the 
C:\PydroXL_19\Scripts directory and run the bat from from that location.  A conda
environment will be created for that conda installation.

To add the CARIS python API place the accompanying nbs.pth file to this file
into the PydroXL_19\envs\CARIS35\Lib\site-packages directory before activating the
environment.

G.Rice 2018-06-12

Updated the bat file to remove all libraries except numpy.  G.Rice 2019-05-03
Changed the name of the bat file to build_caris_env from nbspython.  G.Rice 2019-05-07  