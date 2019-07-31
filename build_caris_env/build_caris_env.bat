conda create --name CARIS35 --file %~dp0CARIS35.yml
xcopy %~dp0nbs.pth %cd%\..\envs\CARIS35\Lib\site-packages\