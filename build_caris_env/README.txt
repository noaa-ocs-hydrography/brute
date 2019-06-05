build_caris_env.bat must be run with conda in the path with the CARIS35.yml file.  
For example, cd to the C:\PydroXL_19\Scripts directory and run the bat from from 
that location.  A conda environment will be created for that conda installation.

In addition, the local nbs.pth file will be copied to the conda installation 
environment in site-packages.  This file must reference the location to find the
caris python module.

G.Rice 2018-06-12

Updated the bat file to remove all libraries except numpy.  G.Rice 2019-05-03
Changed the name of the bat file to build_caris_env from nbspython.  G.Rice 2019-05-07
Updated the bat file to reference the explicit conda environment variables of the current environment
build held in the newly included CARIS35.yml file. C.Koprowski 2019-05-14
Updated to include the information concerning the movement of the pth file. G.Rice 2019-06-05