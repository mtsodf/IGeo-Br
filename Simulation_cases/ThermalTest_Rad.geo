
+++++++++++++++++++++++++++++++++++++++++++++
Geomechanics Modeling
Eng. Luis Carlos de Sousa Jr
+++++++++++++++++++++++++++++++++++++++++++++

RAD MODEL
Fluid flow with thermal solution

#WORKPATH
'C:\Users\cxnx\Desktop\Artigos\2Dfrac\Models'
(ALL FILES SHOULD BE IN THIS PATH)

*************************************************************************
* MODELS   																*
*************************************************************************

COMPLETION DATA----------------------------------------------------------

c #PREFRACTURE
1				    														<- number of initial fracture completions.
1  11    5.0   0.001                                        (lower sand)	<- K1,K2, Xlength, aperture

c #PERFORATIONDEPTH
0.10

FLUID DATA---------------------------------------------------------------

#FLUIDINITIALIZATION                                                    
-5544.28     56464.73                                                        <- OWdepth(m)        pressure(kPa)
10.10085     10.10085                                                        <- oilGrad(kPa/m)    waterGrad(kPa/m)

#WATERPHASEPROPS                                                        
1.0                                                                         <- water_bwI
56464.73																	<- water_pref		kPa
1.028																		<- density_pref
0.0																		    <- water_cw        	1/kPa
1.0																			<- water_viscI     	cP

THERMAL MODEL------------------------------------------------------------

#InitialTemperature
-5552.  100.0  0.0333														<- Z, Tinj(oC), Grad(oC/m)

#ThermalModel
1				  															<- number of set of properties.
1  10	1.0E-5  2.5  830.0   4100.0   0.0									<- K1,K2,ThermalExpCoeff(1/oC),ThermalConductivity(W/mK),RockSpecificHeatCapacity(J/K.Kg),FluidSpecificHeatCapacity(J/Kg.K),rock-fluid Heat Transfer (W/m2K)

*************************************************************************
* SOLVERS																*
*************************************************************************

#NTHREADS
6      

#FLOWFLAG
#THERMALFLAG                               

#PerkinsGonzales                               

*************************************************************************
* WORKFLOW: HYDRAULIC FRACTURE MODEL									*
*************************************************************************

#WF_WaterInj

#WELLPRESSURECONSTRAINT
-5553.0  -5559.0															<- z1,z2
>																			<- time (min); value (kPa ou m³/d); Zref(M)
10	    70400	-5463
100		70400	-5463
1000  	70400	-5463
1000.001       0       0
1000.01        0       0
1000.1         0       0
<

c #WELLRATECONSTRAINT
-5553.0  -5559.0															<- z1,z2
>																			<- time (min); value (m³/d); Zref(M)
10	    1000.	-5463
100		1000.	-5463
1000	1000.	-5463
<

#WaterInjectionTemperature
-5552.  30.0  0.00															<- Z, Tinj(oC), Grad(oC/m)

*************************************************************************
* LAYERS																*
*************************************************************************

-------------------------------------------------------------------------
Name Top  POR   Kh  Kv   Eh  Ev   niH   niV  Kg   Sv   Shmax Shmin
      m         mD  mD	 kPa kPa             GPa  kPa   kPa   kPa 
-------------------------------------------------------------------------
#NWB_STRUCTUREDMESH                                                                
>                                                                                  
ThermalTest																	<- nome do poço                            
RAD																			<- mesh squema: BOX or RAD
20  10   10   1.20															<- Nx,Ny,Nz,ratio
0.10   10.00																<- Lx  Ly
#TABLEofProperties
1	-5551.  	0.10	    0.0	     0.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
2	-5552.  	0.10      0.0	     0.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
3	-5553.  	0.10   1000.0	   100.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
4	-5554.  	0.10   1000.0	   100.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
5	-5555.  	0.10    100.0	    10.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
6	-5556.  	0.10      0.0	     0.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
7	-5557.  	0.10      0.0	     0.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
8	-5558.  	0.10    500.0	    50.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
9	-5559.  	0.10   1000.0	   100.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
10	-5560.  	0.10      0.0	     0.0		7.91E+06	7.91E+06	0.20	0.20	38	133714	127028	72657
<
