conda create -q -y --name CARIS35 python=3.5 numpy
xcopy %~dp0nbs.pth %cd%\..\envs\CARIS35\Lib\site-packages\