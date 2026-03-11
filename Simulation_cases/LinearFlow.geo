+++++++++++++++++++++++++++++++++++++++++++++++++
IGeo - Injection Geomechanics Modeling
Developed under D.Sc. at Geothecnics program
Departament of Civil and Enviromental Engineering
University of Brasilia
Luis C de Sousa Jr
Adviser: Manoel P C Neto
vrs.1.0: september of 2018
+++++++++++++++++++++++++++++++++++++++++++++++++


******************************************************************************
									CASE DESCREPTION

2D model for linear or piston flow with constant water rate.

Basic data:
	reservoir top					: -500.0 m
	reservoir thickness		: 1.0 m
	horiozontal length			: 5 m x 20 m
	porosity						: 0.30
	permeability					: 9.869233E-13 m2 (1000.0 mD)
	fluid compressibility		: 0.50E-6 1/kPa @ 50.0E2	kPa
	fluid viscosity				: 1.0E-3 Pa.s (1.0 cP) @ 50.0E2 kPa
	injection rate					: 5.787E-3 m3/s (500.0 m3/d)

******************************************************************************
									FLUID DATA

Fluid initial pressure.

#FLUIDINITIALIZATION
-500.0	50.0E2					[ z depth(m), pressure (kPa_ ]
0.0         0.0						[ pressure gradient above (kPa/m), pressure gradient bellow (kPa/m) ]


Fluid properties.

#WATERPHASEPROPS  <<< change for #COMPRESSIBLEFLUIDPROPS
1.0									[fluid formation factor at p_ref]
50.0E2								[reference pressure kPa - p_ref] <<< move to the last line
1.028								[density at p_ref]
0.5E-6								[compressibility 1/kPa at p_ref]
1.0									[viscosity cP at p_ref]


******************************************************************************
									SOLVER PARAMETERS

Fluid flow solver.
#FLOWFLAG

OpenMP(R) threads for parallel processing.
#NTHREADS
2

Minimum time step factor in relation to provided time step.
#INITIALSUBSTEPFACTOR
0.10

******************************************************************************
									INJECTION DATA

Fracture seed

#PREFRACTURE
1					     	[number of initial fracture completions]
1  2  5.0  0.010	   				[K1, K2, fracture lemgth, fracture aperture (m)]			<<< fracture aperture


Injection data

Injection workflow
#WF_WaterInj		<<< #WF_FluidInjection

Injection rates

#WELLRATECONSTRAINT
-500.0  -501.0						[injection interval: z1(m), z2(m)]
>
1.667E-06   60.0  -500.0			[time (min), rate (m³/d), --]
6.667E-06   60.0  -500.0	
4.000E-05   60.0  -500.0	
3.200E-04   60.0  -500.0	
3.200e-03   60.0  -500.0
3.840e-02   60.0  -500.0	
<  

******************************************************************************
									LAYERS AND ROCK PROPERTIES

Table explanation:
	Each linear reffers to a k-layer, from the top (k=1) to bottom (k=NK)
	Name		: layer tag name
	Top			: z depth (m)
	POR			: rock's porosity
	Kh			: horizontal permeability (mD)
	Kv			: vertical permeability (mD)
	E				: Young Modulus (kPa)		<<< trocar por um young só
	ni				: Poisson Voefficient			<<< trocar por um único poisson
	Kg			: grain Bulk Modulus (GPa) <<< trocar por kPa
	Sv			: vertical stress (kPa)
	Sh_max	: maximum vertical stress (kPa)
	Sh_min 	: minimum vertical stress (kPa)

#NWB_STRUCTUREDMESH               		<<< mudar para #FORMATIONS_GRID>  
>                                                   
LinearFlow        	<- nome do poço			<<< remover: o nome é o nome do arquivo.   
5  25  2  1.2						[grid discretization: NI or NX, NJ or NY, NK or NZ, ratio for J or Y]
5.0  20.0							[horizontal lengths: LX, LY]
#TABLEofProperties
layer1  -500.0  0.30  1000.0  1000.0  0.00  0.00  38.0  100.0E2  90.0E2  75.0E2
layer2  -500.5  0.30  1000.0  1000.0  0.00  0.00  38.0  100.0E2  90.0E2  60.0E2
<
