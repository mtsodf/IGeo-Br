Rem cls
@echo off

set folder="C:\Users\lcsou\OneDrive\Documents\Development\GeotechThesis\Casos\Tests"
Geotech.exe %folder% "OneElementTest.geo"

Rem set folder="C:\Users\lcsou\OneDrive\Documents\Development\GeotechThesis\Casos\OffshoreGeomech"
Rem Geotech.exe %folder% "ReservoirCompaction.geo"

set folder="C:\Users\lcsou\OneDrive\Documents\Development\GeotechThesis\Casos\OffshoreGeomech\SlopeStability"
Geotech.exe %folder% "OffShoreSlope.geo"
