 
+++++++++++++++++++++++++++++++++++++++++++++
Geomechanics Modeling
Eng. Luis Carlos de Sousa Jr
+++++++++++++++++++++++++++++++++++++++++++++

2D BOX MODEL
Flow model

#WORKPATH
'C:\Users\cxnx\Desktop\Artigos\2Dfrac\Models'
(ALL FILES SHOULD BE IN THIS PATH)

#NTHREADS
3      

#FLOWFLAG

#CRITICALTIMESTEPFACTOR
100

*************************************************************************
* PROBLEMA DATA   									*
*************************************************************************

Fluid Initialization
{OWdepth(m). Pressure(kPa)}
{oilGrad(kPa/m)} {waterGrad(kPa/m)}
#FLUIDINITIALIZATION
-500.0	50.0E2
0.0         0.0

#PREFRACTURE
1					      <- number of initial fracture completions
1  2  5.0  0.010	    		      <- K1, K2, Xlength, aperture

*************************************************************************
* WORKFLOW: HYDRAULIC FRACTURE MODEL						*
*************************************************************************

#WF_WaterInj

#WELLRATECONSTRAINT
-500.0  -501.0				<- z1,z2
>
1.667E-06   60.0  -500.0		<- time (min) value (m³/d); Zref(m)
6.667E-06   60.0  -500.0	
4.000E-05   60.0  -500.0	
3.200E-04   60.0  -500.0	
3.200e-03   60.0  -500.0
3.840e-02   60.0  -500.0	
<  

*************************************************************************
* PROPERTIES										*
*************************************************************************
																		 
#WATERPHASEPROPS
1.0						<- water_bwI
50.0E2					<- water_pref	kPa
1.028						<- density_pref
0.5E-6					<- water_cw       1/kPa.
1.0						<- water_viscI    cP

-------------------------------------------------------------------------
Name Top  POR   Kh  Kv   Eh	Ev   niH   niV  Kg   Sv   Shmax Shmin
      m         mD  mD	 kPa	kPa             GPa  kPa   kPa   kPa 
-------------------------------------------------------------------------
#NWB_STRUCTUREDMESH                                                                
>                                                                                  
LinearFlow        			<- nome do poço                            
BOX						<- mesh squema: BOX or RAD
5  25  2  1.2				<- Nx,Ny,Nz,ratio
5.0  20.0					<- Lx, Ly (normal and tangential to frature propagation plane).
#TABLEofProperties
layer1  -500.0  0.30  1000.0	1000.0	0.0  0.0  0.30  0.30  38.0  100.0E2  90.0E2  75.0E2
layer1  -500.5  0.30  1000.0	1000.0	0.0  0.0  0.30  0.30  38.0  100.0E2  90.0E2  60.0E2
<
