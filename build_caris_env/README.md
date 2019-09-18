### Building the CARIS35 Anaconda environment 

`build_caris_env.bat` must be run with access to both `CARIS35.yml` and `conda.exe`; if you have Pydro installed, follow these steps:
1. edit the text file `nbs.pth` so that it points to your CARIS install on the first line, and to this repository on the second line
2. open Command Prompt (`cmd.exe`)
3. change directories to the Pydro `Scripts` folder (for instance `cd C:\PydroXL_19\Scripts`)
4. run `build_caris_env.bat` by calling its absolute path (for instance `C:\national-bathymetric-source\build_caris_env\build_caris_env.bat`)

An Anaconda environment `CARIS35` will be created. To activate this environment, run `conda activate CARIS35` with access to `conda.exe` (for instance from the Pydro `Scripts` folder or the Anaconda Prompt).

In addition, `nbs.pth` will be copied to the new environment in `site-packages`, so that the CARIS wrapper can see your CARIS installation.

G.Rice 2018-06-12 -         created README
G.Rice 2019-05-03 -         Updated the bat file to remove all libraries except numpy.
G.Rice 2019-05-07 -         Changed the name of the bat file to build_caris_env from nbspython.
C.Koprowski 2019-05-14 -    Updated the bat file to reference the explicit conda environment variables of the current environment build held in the newly included CARIS35.yml file.
G.Rice 2019-06-05 -         Updated to include the information concerning the movement of the pth file.