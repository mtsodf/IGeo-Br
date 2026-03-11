
+++++++++++++++++++++++++++++++++++++++++++++
Geomechanics Modeling
Eng. Luis Carlos de Sousa Jr
+++++++++++++++++++++++++++++++++++++++++++++

RADIAL FINER MODEL.
Implicit solver
Viscoplastic Model

*************************************************************************
* MAIN DATA.															*
*************************************************************************

#WORKPATH
'C:\Users\cxnx\Documents\Doutorado\Codes\WTHM\Cases'
(ALL FILES SHOULD BE IN THIS PATH)

CONTINUOUS DAMAGE MECHANICS----------------------------------------------

#CDM_FRACTURE
3
1     3  100.0  1.0E-1  455.0    (cap shale)								<- K1   K2  T0(kPa)  Max.YYdilation   Gc(Joule)
4    17   10.0  1.0E-1  455.0    (upper sand)
18   20  100.0  1.0E-1  455.0    (intra shale)

PERMEABILITY EVOLUTION---------------------------------------------------
Perm enhancement with shear and plastic strain.
C #PERMENHANCEMENT
C 2
C 41   54  10.0  50.0            (upper sand)    							<- K1   K2  DeltaPerm(mD)/0.0001plasticDilation  DeltaPerm(mD)/0.0001plasticShear
C 58   76  10.0  50.0            (lower sand) 
 
COMPLETION DATA----------------------------------------------------------

#PREFRACTURE
1																			<- number of initial fracture completions.
4  17  1.00  0.001	             (upper sand)   							<- K1,K2, Xlength, aperture

FLUID DATA---------------------------------------------------------------

#FLUIDINITIALIZATION
-2300.70 	23202.55					  									<-	OWdepth(m)				Pressure(kPa)
10.085		  10.085					  									<-	oilGrad(kPa/m)		waterGrad(kPa/m)
							
#WATERPHASEPROPS							
1.0																			<- water_bwI
23202.55																	<- water_pref		kPa
1.028																		<- density_pref
0.5E-6																		<- water_cw        	1/kPa
1.0																			<- water_viscI     	cP

THERMAL MODEL------------------------------------------------------------

#InitialTemperature
-2485.3  50.0  0.0333														<- Z, Tinj(oC/m), Grad(oC/m)
							
#ThermalModel							
1				  															<- number of set of properties.
1  20	1.0E-5  1.7  830.0   4100.0   1.0E1									<- K1,K2,ThermalExpCoeff(1/oC),ThermalConductivity(W/mK),RockSpecificHeatCapacity(J/K.Kg),FluidSpecificHeatCapacity(J/Kg.K),rock-fluid Heat Transfer (W/m2K)

-------------------------------------------------------------------------

*************************************************************************
* SOLVERS																*
*************************************************************************

#NTHREADS
3       

#FLOWFLAG
#STRUCTURALFLAG  
C #THERMALFLAG

C #PerkinsGonzales                               

*************************************************************************
* WORKFLOW: HYDRAULIC FRACTURE MODEL									*
*************************************************************************

#WF_WaterInj

#WELLRATECONSTRAINT
-2485.99  -2494.90															<- z1,z2
>																			<- time (min) value (m³/d); Zref(m)
0.1		  10	-2485.3
60		1000	-2485.3
120		2000	-2485.3
240		3000	-2485.3
240.01	   0	-2485.3
240.10	   0	-2485.3
250  	   0	-2485.3
360  	   0	-2485.3
<

#WaterInjectionTemperature
-2485.3  20.0  0.00															<- Tinj(oC/m), Z, Grad(oC/m)

*************************************************************************
* LAYERS																*
*************************************************************************

"Skin" permeability
1mD = 9.8692×10−13 m²
#EDITprops
1
'KH'     9.8692E-14   1  20  1  1  4  17						    		<- prop.name   value   I1  I2  J1  J2  K1  K2

c #ELASTOPLASTIC_MODEL
3					                                                      	<- number of regions
region 1 (cap shale)*****************************************************
2																		   	<- 1- Damage (EXPLICITY); 2 - Mohr-Coulomb; 3 - CAM-CLAY (modified)
1    3															           	<- K1  K2
2.0E3  20.0  5.0                        	                               	<- c(KPa)  c(degrees) dilation(degrees)
region 2 (sandstone)*****************************************************
2																		    <- 1- Damage (EXPLICITY); 2 - Mohr-Coulomb; 3 - CAM-CLAY (modified)
4   17															            <- K1  K2
0.0  30.0  10.0             	                                           	<- c(KPa)  c(degrees) dilation(degrees)
region 5 (underburden)***************************************************
2																		    <- 1- Damage (EXPLICITY); 2 - Mohr-Coulomb; 3 - CAM-CLAY (modified)
18  20    													                <- K1  K2
2.0E3  20.0  5.0             	                                           	<- c(KPa)  c(degrees) dilation(degrees)
*************************************************************************


-------------------------------------------------------------------------
Name Top  POR   Kh  Kv   Eh		Ev	  niH	 niV  Kg   Sv   Shmax Shmin
      m         mD	mD	 kPa	kPa               GPa  kPa   kPa   kPa 
-------------------------------------------------------------------------
#NWB_STRUCTUREDMESH                                                                
>                                                                                  
MLS3B_test				 		  											<- nome do poço.                             
BOX													                 		<- mesh squema: BOX or RAD
150  20  20  1.30    				                 						<- Nrad,Ntheta,Nz,ratio
150.0  150.0								                 				<- LX  LY
#TABLEofProperties
1	-2483.93	0.35	   0	  0	10600000	10600000	0.384	0.384	38	40616	40210	34772
2	-2484.61	0.35	   0	  0	8640000		 8640000	0.354	0.354	38	40631	40225	34784
3	-2485.30	0.35	   0	  0	6400000		 6400000	0.345	0.345	38	40646	40239	34796
4	-2485.99	0.30	1000	100	1110000		 1110000	0.342	0.342	38	40661	40254	34808
5	-2486.67	0.30	1000	100	1170000		 1170000	0.335	0.335	38	40675	40269	34819
6	-2487.36	0.30	1000	100	1210000		 1210000	0.320	0.320	38	40690	40283	34831
7	-2488.04	0.30	1000	100	1300000		 1300000	0.302	0.302	38	40705	40298	34843
8	-2488.73	0.30	1000	100	1280000		 1280000	0.296	0.296	38	40720	40313	34855
9	-2489.42	0.30	1000	100	1230000		 1230000	0.303	0.303	38	40735	40327	34867
10	-2490.10	0.30	1000	100	1250000		 1250000	0.311	0.311	38	40749	40342	34879
11	-2490.79	0.30	1000	100	1330000		 1330000	0.306	0.306	38	40764	40357	34891
12	-2491.47	0.30	1000	100	1350000		 1350000	0.311	0.311	38	40779	40371	34902
13	-2492.16	0.30	1000	100	1280000		 1280000	0.314	0.314	38	40794	40386	34914
14	-2492.84	0.30	1000	100	1260000		 1260000	0.293	0.293	38	40809	40400	34926
15	-2493.53	0.30	1000	100	1270000		 1270000	0.287	0.287	38	40823	40415	34938
16	-2494.22	0.30	1000	100	1300000		 1300000	0.302	0.302	38	40838	40430	34950
17	-2494.90	0.30	1000	100	1380000		 1380000	0.342	0.342	38	40853	40444	34962
18	-2495.59	0.35	   0	  0	6400000		 6400000	0.345	0.345	38	40868	40459	34973
19	-2496.27	0.35	   0	  0	6400000		 6400000	0.345	0.345	38	40883	40474	34985
20	-2496.96	0.35	   0	  0	6400000		 6400000	0.345	0.345	38	40897	40488	34997
<
