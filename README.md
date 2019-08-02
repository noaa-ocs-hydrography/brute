### Building the NBS Anaconda environment 

`build_conda_env.bat` must be run with access to both `environment.yml` and `conda.exe`; if you have Pydro installed, follow these steps:
1. open Command Prompt (`cmd.exe`)
2. change directories to the Pydro `Scripts` folder (or instance `cd C:\PydroXL_19\Scripts`)
3. run `build_conda_env.bat` by calling its absolute path (for instance `C:\national-bathymetric-source\build_conda_env.bat`)

An Anaconda environment `NBS` will be created. To activate this environment, run `conda activate NBS` with access to `conda.exe` (for instance from the Pydro `Scripts` folder or the Anaconda Prompt).