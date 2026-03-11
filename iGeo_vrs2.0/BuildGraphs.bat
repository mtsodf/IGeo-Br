@echo on
Rem Compiling Fortran code---------------------------------------------------------------

del BuildGaphs.exe

Rem ifort.exe ^
Rem /O3 /Qm64 /Qpar /Qopenmp /Qmkl:parallel ^
Rem InjeSim_BuildGaphs.f90 -o BuildGaphs.exe

ifort.exe ^
/O3 /Qm64 /Qpar /Qopenmp /Qmkl:parallel ^
InjeSim_BuildGaphs2.f90 -o BuildGaphs.exe

.\BuildGaphs