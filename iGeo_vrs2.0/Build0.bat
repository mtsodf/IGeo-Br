cls
del *.exe *.mod *.obj *.log
Rem del *.exe

@echo on
Rem Compiling Fortran code---------------------------------------------------------------
ifort.exe ^
/Qm64 /Qmkl:parallel ^
SolidMechanicsLib.f90 ^
MainConstants.f90 ^
IO.f90 ^
NumericalLib.f90 ^
MeshAuxiliares.f90 ^
MeshIO.f90 ^
NumericalUtilities.f90 ^
FEMlib.f90 ^
ProblemDefinition.f90 ^
NearWellboreModel.f90 ^
StressLoad.f90 ^
ConstitutiveGeoModels.f90 ^
DiscretizedModel.f90 ^
WellModelBuilder.f90 ^
GravityLoad.f90 ^
PorepressureLoad.f90 ^
UniformLoad.f90 ^
ThermalLoad.f90 ^
InternalForceCalculator.f90 ^
ElementStiffness.f90  ^
StiffnessMatrixCSR.f90 ^
StructuralModel.f90 ^
InterfaceMechanics.f90 ^
HydraulicFractureModel.f90 ^
ParallelRoutines.f90  ^
Damage.f90 ^
ViscoplasticExplicitySolver.f90 ^
PCG_Flow.f90 ^
FluidFlowSolver.f90 ^
Thermal.f90 ^
Report.f90 ^
VTK.f90 ^
VTKgeneral.f90 ^
CreepMod.f90 ^
WF_Excavation.f90 ^
WF_WaterInjection.f90 ^
WellStability.f90 ^
DynamicSolver.f90 ^
WF_MEDSolver.f90 ^
GeomechWell_main_these.f90 ^
iGeo.f90 -o iGeo.exe

RunTests.bat