!*****************************************************************************
!COMPUTATIONAL GEOMECHANICS
!Eng. Luis Carlos de Sousa Junior
!Computational Geomechanics
!*****************************************************************************

  MODULE ProblemDefinition

    !***************************************************************************
    !DATA STRUCTURE.
    !***************************************************************************

    TYPE DamageShear
      
      INTEGER:: K1,K2
      REAL(8):: LimitPlasticDeviatoricStrain
      REAL(8):: ResidualDamage
      
    END TYPE
    
    TYPE DamageLiquefaction
      
      INTEGER:: K1,K2
      REAL(8):: MaxMeanEffectiveStress                    !Liquefaction limit.
      REAL(8):: ShearDamage,ElasticDamage
      
    END TYPE
    
    !BASIC WELL DATA.
    TYPE WellboreData

      CHARACTER(LEN=40):: WellName

      !Geometria e malha
      CHARACTER(LEN=3):: MeshType                       !BOX or RAD (radial).
      INTEGER:: NTHETA,NRAD,NZ                          !structured grid.
      REAL(8):: Ratio                                   !ratio at RAD direction.
      REAL(8):: Radius(2)                               !raio interno e externo do modelo.
      REAL(8),ALLOCATABLE:: LayersTops(:)               !topo das camadas.
      CHARACTER(LEN=12),ALLOCATABLE:: LayersNames(:)    !nome das camadas.

      !Layer Properties
      REAL(8),ALLOCATABLE:: EffectivePerms(:,:)         !permeabilidades (H e V) efetiva do fluido principal do po�o.
      REAL(8),ALLOCATABLE:: Porosity(:)                 !porosidade.
      REAL(8),ALLOCATABLE:: YoungModulus(:,:)           !M�dulo de Young: (1,:) vertical, (2,:) horizontal.
      REAL(8),ALLOCATABLE:: Poisson(:,:)                !Raz�o de Poisson: (1,:) vertical, (2,:) horizontal.
      REAL(8),ALLOCATABLE:: Kg(:)                       !M�dulo Bulk dos gr�os.

      !Initialization data
      REAL(8),ALLOCATABLE:: VerticalStress(:)
      REAL(8),ALLOCATABLE:: MaxHorizontalStress(:)
      REAL(8),ALLOCATABLE:: MinHorizontalStress(:)
      !OR
      REAL(8),ALLOCATABLE:: StressTensor(:,:)
      REAL(8),ALLOCATABLE:: Porepressure(:)

      !Mud Weight
      REAL(8):: MudWeight

      !Initial Saturation
      REAL(8),ALLOCATABLE:: Sat(:,:)          !Permeabilidade relativa de cada layer: Swi,Krnw@Swi,Snwr,Krw@Snwr

      !BOUNDARY CONDITION
      INTEGER:: Kint1,Kint2,Kext1,Kext2
      REAL(8):: gradInt,gradExt

      !***************************************************************************
      !CONSTITUTIVE MODEL.
      !***************************************************************************
      INTEGER,ALLOCATABLE:: ConsModel(:,:)              !(1,:) - model number.
                                                        !(2,:) - K1
                                                        !(3,:) - K2
      REAL(8),ALLOCATABLE:: ConsModelParameters(:,:)    !(:,i) - Cons. model parameters.
        
      INTEGER,ALLOCATABLE:: ListOfUndrainedLayers(:)
      REAL(8),ALLOCATABLE:: UndrainedLayerParameters(:,:)
        
      INTEGER,ALLOCATABLE:: ListOfLiquefactionLayers(:)
      REAL(8),ALLOCATABLE:: LiquefactionLayerParameters(:,:)
        
      !***************************************************************************
      !PRE-FRACTURE DATA.
      !***************************************************************************
      LOGICAL:: PreFractureExists
      REAL(8),ALLOCATABLE:: PreFractureLengh(:)       !<-length for each k layer.
      REAL(8),ALLOCATABLE:: PreFractureAperture(:)    !<-aperture for each k_layer.
      
      !INTERFACES.
      INTEGER,ALLOCATABLE:: FractureInterval(:,:)     !(I1-I2-J1-J2-Z1-Z2,counter)
      REAL(8),ALLOCATABLE:: TractionStrengthRatio(:)
        
      !***************************************************************************
      !INJECTION DATA.
      !***************************************************************************
      LOGICAL:: InjectionRate,InjectionPressure
      REAL(8):: InjX,InjZ1,InjZ2
      REAL(8):: Rate
      REAL(8):: Pinj,Zinj
      REAL(8),ALLOCATABLE:: InjPressureTable(:,:),InjRateTable(:,:)
        
      !***************************************************************************
      !WELL BORE DATA.. ( #DirectionalWell )
      !***************************************************************************
      LOGICAL:: HorWellFlag                             !BOX MODEL AND HORIZONTAL WELL (X DIR)
      REAL(8):: WellZ,WellRadius,WellLength,WellPermMult
      LOGICAL,ALLOCATABLE:: WellBoreWallFlag(:)
      
      REAL(8):: Z_phase1                                !9 5/8 casing shoe.
      LOGICAL:: VerticaDisplacementRestriction
      
      !PRODUCTION CASE COMPLETION.
      ! #ProductionCompletion
      CHARACTER(LEN=10):: ProductionCompletionType			!<- OPEN , CASED
			REAL(8):: ProductionCompletionPerfLength          !<- PERFORATION LENGTH (m)
			INTEGER:: ProductionCompletionDistribution			  !<- 1-ONE, 2-TWO or 0-DISTRIBUTED
      
      !DRILLING DATA.
      REAL(8):: OpenIntervalZ(2)
      INTEGER:: DrillingSteps
      REAL(8),ALLOCATABLE:: MudWeightHistory(:)
        
    END TYPE WellboreData

    TYPE FracData
      
      INTEGER:: K1,K2
      REAL(8):: PermIncrementByPlasticDilation      !mD/0.001
      REAL(8):: PermIncrementByPlasticShear         !mD/0.001
      
    END TYPE
    
    !TYPE TractionDamage
    !  
    !  INTEGER:: K1,K2
    !  REAL(8):: T0                                  !Original Traction strength.
    !  REAL(8):: MaxPlasticDilation                  !zero traction strength for higher platic dilation.
    !  REAL(8):: Gc                                  !critical release energy.
    !  
    !END TYPE
    
    TYPE CDM_DATA
      
      INTEGER:: K1,K2
      REAL(8):: ST
      REAL(8):: EyyMax
      REAL(8):: Gc
      
    END TYPE

    TYPE TemperatureVariation
        INTEGER:: Step
        INTEGER:: I1,I2,J1,J2,K1,K2
        REAL(8):: DeltaT
    END TYPE
    
    TYPE InitialStateVariations
        INTEGER:: Step
        INTEGER:: I1,I2,J1,J2,K1,K2
        REAL(8):: DeltaP,DeltaT
    END TYPE
    
    TYPE OutflowModel
        INTEGER:: I1,I2,J1,J2,K1,K2
        INTEGER:: DIR                   !1 OR 2 (X OR Y)
        REAL(8):: J !(m3/d / Pa)
        REAL(8):: Pest,Zref_Pest        !static pressure at a given Z depth.
    END TYPE
    
    TYPE WellPipe
        REAL(8):: Diameter
        INTEGER,ALLOCATABLE:: WellPath(:,:)         !(:,# pipe node) - (4) DIR,i,j,k
                                                    !DIR: 1:X, 2:Y, 3:Z
    END TYPE
    
    TYPE TimeSerie
      CHARACTER(LEN=350):: Name
      REAL(8),ALLOCATABLE:: Table(:,:)    !time, magnitude
    END TYPE

    TYPE TriangularSurface
      !Formed with the dot product between arrays (pt1-pt2) and (pt1-pt3)
      CHARACTER(LEN=350):: Name
      REAL(8):: X1,Y1,Z1,&
                X2,Y2,Z2,&
                X3,Y3,Z3
    END TYPE

    TYPE DistributedForce
      CHARACTER(LEN=350):: Name
      CHARACTER(LEN=350):: TimeSerieName
      CHARACTER(LEN=350):: TriangularSurfaceName
      REAL(8):: Amplitude     !To be multiplied by the magnitude in time series
    END TYPE

    TYPE PD

      CHARACTER(LEN=320):: BaseFile, &   !base file.
                           Path          !work path.
                             
      CHARACTER(LEN=3):: MeshType        !BOX or RAD (radial).
      
      INTEGER:: nthreads
      
      LOGICAL:: Report
        
      !PHYSICS FLAG-------------------------------------------------------------
      LOGICAL:: Flow_flag, Structural_flag, Thermal_flag

      !POREPRESSURE INITIALIZATION----------------------------------------------
      REAL(8):: OWdepth, OWpressure, &   !OW contact.
                wGrad, nwGrad            !weating and non-weating fluid.

      !FLUID PROPERTIES.
      REAL(8):: water_bwI, water_pref,density_pref, water_cw, water_viscI,&     !water props.
                oil_bwI, oil_viscI                                              !oil props.
        
      CHARACTER(LEN=40):: ModelName
        
      !HYDRAULIC FRACTURE MODEL (continuous damage model)-----------------------
      REAL(8),ALLOCATABLE:: FracDamage_Xlength(:),&
                            FracDamage_Aperture(:),&
                            FracDamage_Kn(:)
      INTEGER,ALLOCATABLE:: FracDamage_Kinterval(:,:)
        
      !TRANSIENT FLOW DATA------------------------------------------------------       
      REAL(8):: DeltaTime, TotalTime
      INTEGER:: nSteps
      REAL(8),ALLOCATABLE:: MaxRate(:)
      REAL(8):: CriticalTimeStepFactor
      
      TYPE(InitialStateVariations),ALLOCATABLE:: InitialVariations(:)
                                                    !INTEGER:: Step
                                                    !INTEGER:: I1,I2,J1,J2,K1,K2
                                                    !REAL(8):: DeltaP,DeltaT
        
      !CONSTRUCTION DATA (STRESS INITIALIZATION)--------------------------------
      REAL(8):: GeothermalGrad
      INTEGER,ALLOCATABLE:: LayersInitialProps(:,:)   !K1,K2: layers interval for each line on InitialProps.
      REAL(8),ALLOCATABLE:: InitialProps(:,:)         !densit, porosity, Young (kPa), Poisson, Conductivity (oC/m), kh(mD), kv(mD).
      REAL(8),ALLOCATABLE:: ConstrutionSteps(:,:)     !surface temperature (oC), YY strain, XX strain, +1 for drained OR -1 for undrained.
        
      !OUTPUT CONTROL-----------------------------------------------------------
      INTEGER:: StructuralSolverPrintFreq
      INTEGER:: ElementMonitor                        !IX,IY,IZ
      INTEGER:: NodeMonitor
        
      !FLUID B.C.-------------------------------------------------------------
      LOGICAL:: InfiniteFluidBC
      
      !OUTFLOW (FAR FIELD FLOW)
      TYPE(OutflowModel),ALLOCATABLE:: FarFieldFlow(:)
        
      !PROPERTIES EDITING-----------------------------------------------------
      REAL(8),ALLOCATABLE:: EditProps_values(:)
      INTEGER,ALLOCATABLE:: EditProps_IJKinterval(:,:)
      CHARACTER(LEN=40),ALLOCATABLE:: EditProps_PropsNames(:)
        
      TYPE (WellboreData):: WellData
      
      !FRAC-------------------------------------------------------------------
      TYPE (FracData), ALLOCATABLE:: Frac(:)
      !TYPE (TractionDamage),ALLOCATABLE:: TracDamage(:)
      TYPE (CDM_DATA),ALLOCATABLE:: CDM(:)
      LOGICAL:: FracClosure = .FALSE.
      LOGICAL:: FracUnbalVolume = .FALSE.
      
      !DAMAGE-----------------------------------------------------------------
      TYPE (DamageShear),ALLOCATABLE:: ShearStrainDamage(:)
      TYPE (DamageLiquefaction),ALLOCATABLE:: Liquefaction(:)
      
      !THERMAL----------------------------------------------------------------
      LOGICAL:: PerkinsGonzales
      REAL(8),ALLOCATABLE:: ThermalProps(:,:)
      INTEGER,ALLOCATABLE:: ThermalPropsKinterval(:,:)
      REAL(8):: TopTemperatureZ, TopTemperature, TemperatureGrad
      REAL(8):: Tinj,TinjZ, TinjGrad
      TYPE(TemperatureVariation),ALLOCATABLE:: ExplicitTemperatureChange(:)
                                                    !INTEGER:: Step
                                                    !INTEGER:: I1,I2,J1,J2,K1,K2
                                                    !REAL(8):: DeltaT
      
      !LARGE DISPLACEMENTS----------------------------------------------------
      !#UPDATEGEOMETRY
      LOGICAL:: UpdateGeometry
      
      !LONG TER SIMULATION---------------------------------------------------
      LOGICAL:: LongTermFlag
      INTEGER:: LongTermNumStep
      REAL(8):: LongTermDt
      CHARACTER(LEN=50),ALLOCATABLE:: TimeSteps(:)                      !<- time steps in minutes
      
      !SOIL MODEL------------------------------------------------------------
      CHARACTER(LEN=50),ALLOCATABLE:: CreepTable(:,:)                        !<- specific table for creep models. Contains region and creep parameters.
      
      !SOIL MODEL------------------------------------------------------------
      !TYPE(SoilModel) :: SoilData
      
      !ELASTICITY: SPECIAL MODELS--------------------------------------------
      LOGICAL:: TrasverseIsotropicModel
      REAL(8),ALLOCATABLE:: Dip(:),Strike(:)
      
      !FORMATION DAMAGE------------------------------------------------------
      CHARACTER(LEN=50),ALLOCATABLE:: FluidPermDamage(:,:)                      !K1, K2, time(min), perm.damage(adm)
      
      !WELL PIPE-------------------------------------------------------------
      LOGICAL:: WellPipeModel
      TYPE (WellPipe):: Pipe                                                    !Pipe object
      
      !MECS------------------------------------------------------------------
      TYPE(TimeSerie),ALLOCATABLE:: TimeSeries(:)
      TYPE(TriangularSurface),ALLOCATABLE:: TriAreas(:)
      TYPE(DistributedForce),ALLOCATABLE:: Forces(:)

    END TYPE
  
  CONTAINS

  !WELL MODEL.
  !Reads the base file (.geo).
  SUBROUTINE ReadBaseFile(Path,BaseFile,Chanel_base,Problem,WellData,nthreads)

    USE IO
    USE MeshAuxiliares
    USE OMP_LIB  !bibliotecas OpenMP(R)

    !VARIABLES.
    IMPLICIT NONE
    CHARACTER(LEN=320),INTENT(IN):: BaseFile, &   !file with basic information about the case.
                                    Path          !work directory where the case file is located and where temporary files will
                                                  !be created.
    INTEGER,INTENT(IN)::  Chanel_base
    TYPE (PD), INTENT(IN OUT):: Problem
    TYPE (WellboreData), INTENT(IN OUT):: WellData
    INTEGER, INTENT(IN OUT):: nthreads

    INTEGER:: layerCount,i,j,k1,k2,n,nModels,node,counter
    REAL(8):: length,aperture
    REAL(8):: rate
    LOGICAL:: success,found
    CHARACTER(LEN=320):: text
    REAL(8):: z
      
    REAL(8):: Cwater
    INTEGER:: nlayers
    
    CHARACTER(LEN=50),ALLOCATABLE:: table(:,:)

    !***************************************************************************
    WRITE(text,*) 'Reading the problem data.'
    CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )

    !READING BASE FILE.
    WRITE(text,*)'path: ',TRIM(Path)
    CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
    WRITE(text,*)'BaseFile: ',TRIM(BaseFile)
    CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )

    Problem%BaseFile = BaseFile

    !Work path.
    !CALL SearchInFile( Chanel_base, '#WORKPATH', success )
    !IF(success) THEN
    !    READ( Chanel_base, * ) Problem%Path
    !ELSE
    !    CALL Error( "command #WORKPATH not found." )
    !END IF

    !N�CLEOS DE PROCESSAMENTO.
    CALL SearchInFile( Chanel_base, '#NTHREADS', success )
    IF(success) THEN
        READ( Chanel_base, * ) nthreads
    ELSE
        nthreads = 2
    END IF
    Problem%nthreads = nthreads
    CALL OMP_SET_NUM_THREADS( nthreads )

    !***************************************************************************
    ! PRESSURE INITIALIZATION.
    !***************************************************************************
    CALL SearchInFile( Chanel_base, '#FLUIDINITIALIZATION', success )
    IF(success)THEN
      READ(Chanel_base,*) Problem%OWdepth, Problem%OWpressure
      READ(Chanel_base,*) Problem%nwGrad,Problem%wGrad
    END IF

    !***************************************************************************
    ! WELLBORE DATA.
    !***************************************************************************
    CALL SearchInFile( Chanel_base, '#NWB_STRUCTUREDMESH', success )
    IF(.NOT.success) CALL SearchInFile( Chanel_base, '#STRUCTUREDMESH', success )
    IF(success) THEN

          READ(Chanel_base,*) text  ! >
          IF(TRIM(text).NE.">")THEN
              CALL Error( "> not found after #NWB_STRUCTUREDMESH." )
          END IF

          !Geological layers.
          READ(Chanel_base,*) WellData%WellName
          Problem%ModelName = WellData%WellName
          READ(Chanel_base,*) WellData%MeshType
          !READ(Chanel_base,*) WellData%NTHETA,WellData%NRAD,WellData%NZ,WellData%Ratio
          READ(Chanel_base,*) WellData%NRAD,WellData%NTHETA,WellData%NZ,WellData%Ratio
          READ(Chanel_base,*) WellData%Radius(1), WellData%Radius(2)
          READ(Chanel_base,*) text  !#TABLEofProperties.
          IF(TRIM(text).NE."#TABLEofProperties")THEN
              CALL Error( "#TABLEofProperties not found." )
          END IF
          Problem%MeshType = WellData%MeshType
            
          WRITE(text,*) '...number of layers:',WellData%NZ
          CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )

          ALLOCATE( WellData%LayersNames( WellData%NZ ),&
              WellData%LayersTops( WellData%NZ ),&
              WellData%Porosity( WellData%NZ ),&
              WellData%EffectivePerms( 2, WellData%NZ ),&
              WellData%YoungModulus( 2, WellData%NZ ),&
              WellData%Poisson( 2, WellData%NZ ),&
              WellData%Kg( WellData%NZ ),&
              WellData%VerticalStress( WellData%NZ ),&
              WellData%MaxHorizontalStress( WellData%NZ ),&
              WellData%MinHorizontalStress( WellData%NZ ),&
              WellData%Porepressure( WellData%NZ ))!,&
              !WellData%MCparameters(3,WellData%NZ) )
          DO layerCount = 1, WellData%NZ
              READ(Chanel_base,*) WellData%LayersNames(layerCount),&
                  WellData%LayersTops(layerCount),&
                  WellData%Porosity(layerCount),&
                  WellData%EffectivePerms(1,layerCount),&
                  WellData%EffectivePerms(2,layerCount),&
                  WellData%YoungModulus(1,layerCount),&
                  WellData%YoungModulus(2,layerCount),&
                  WellData%Poisson(1,layerCount),&
                  WellData%Poisson(2,layerCount),&
                  WellData%Kg(layerCount),&
                  WellData%VerticalStress(layerCount),&
                  WellData%MaxHorizontalStress(layerCount),&
                  WellData%MinHorizontalStress(layerCount)!,&
                  !WellData%MCparameters(1,layerCount),&
                  !WellData%MCparameters(2,layerCount),&
                  !WellData%MCparameters(3,layerCount)
          END DO  !loop nos layers.
            
          !Lendo dados do peso do fluido de perfuracao.
          CALL SearchInFile( Chanel_base, '#MudHistoryPPG', success )
          IF(success) THEN
            READ(Chanel_base,*) WellData%MudWeight                    !ppg
            WellData%MudWeight = WellData%MudWeight * (9810.0/8.34)   !Pa/m
          ELSE
            WellData%MudWeight = 1.028 * 9810.0                       !Pa/m
          END IF

          !Calculo da poropressao inicial nos layers do modelo.
          DO layerCount = 1, WellData%NZ
              IF(layerCount<WellData%NZ)THEN
                  z = 0.5* (WellData%LayersTops(layerCount) + WellData%LayersTops(layerCount+1) )
              ELSE
                  z = WellData%LayersTops(layerCount)
              END IF
              IF(z>Problem%OWdepth)THEN
                  !acima do contato OW.
                  WellData%Porepressure(layerCount) = Problem%OWpressure - &
                      (z-Problem%OWdepth) * Problem%nwGrad
              ELSE
                  !abaixo do contato OW.
                  WellData%Porepressure(layerCount) = Problem%OWpressure - &
                      (z-Problem%OWdepth) * Problem%wGrad
              END IF
          END DO

    END IF      !if do comando #NWB_STRUCTUREDMESH.
    
    !GENERAL WELL GRID-------------------------------------------------------
    CALL SearchInFile( Chanel_base, '#WELLBOREGRID', success )
    IF(success)THEN
        
        READ(Chanel_base,*) text  ! >
        IF(TRIM(text).NE.">")THEN
            CALL Error( "> not found after #WELLBOREGRID." )
        END IF

        !Geological layers.
        READ(Chanel_base,*) WellData%WellName
        Problem%ModelName = WellData%WellName
        READ(Chanel_base,*) WellData%NRAD,WellData%NTHETA,WellData%NZ,WellData%Ratio
        READ(Chanel_base,*) WellData%Radius(1), WellData%Radius(2)
        READ(Chanel_base,*) text  !#TABLEofProperties.
        IF(TRIM(text).NE."#TABLEofProperties")THEN
            CALL Error( "#TABLEofProperties not found." )
        END IF
        WellData%MeshType = 'RAD'
        Problem%MeshType = WellData%MeshType
            
        WRITE(text,*) '...number of layers:',WellData%NZ
        CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )

        ALLOCATE( WellData%LayersNames( WellData%NZ ),&
            WellData%LayersTops( WellData%NZ ),&
            WellData%Porosity( WellData%NZ ),&
            WellData%EffectivePerms( 2, WellData%NZ ),&
            WellData%YoungModulus( 2, WellData%NZ ),&
            WellData%Poisson( 2, WellData%NZ ),&
            WellData%Kg( WellData%NZ ),&
            WellData%StressTensor(6,WellData%NZ ),&
            WellData%Porepressure( WellData%NZ ) )
        DO layerCount = 1, WellData%NZ
            READ(Chanel_base,*) WellData%LayersNames(layerCount),&
                WellData%LayersTops(layerCount),&
                WellData%Porosity(layerCount),&
                WellData%EffectivePerms(1,layerCount),&
                WellData%EffectivePerms(2,layerCount),&
                WellData%YoungModulus(1,layerCount),&
                WellData%Poisson(1,layerCount),&
                WellData%Kg(layerCount),&
                WellData%StressTensor(1,layerCount),&
                WellData%StressTensor(2,layerCount),&
                WellData%StressTensor(3,layerCount),&
                WellData%StressTensor(4,layerCount),&
                WellData%StressTensor(5,layerCount),&
                WellData%StressTensor(6,layerCount),&
                WellData%Porepressure( layerCount )
        END DO  !loop nos layers.
        WellData%YoungModulus(2,:) = WellData%YoungModulus(1,:)
        WellData%Poisson(2,:) = WellData%Poisson(1,:)
            
        !Lendo dados do peso do fluido de perfuracao.
        CALL SearchInFile( Chanel_base, '#MudHistoryPPG', success )
        IF(success) THEN
          READ(Chanel_base,*) WellData%MudWeight                    !ppg
          WellData%MudWeight = WellData%MudWeight * (9810.0/8.34)   !Pa/m
        ELSE
          WellData%MudWeight = 1.028 * 9810.0                       !Pa/m
        END IF

        !Calculo da poropressao inicial nos layers do modelo.
        CALL SearchInFile( Chanel_base, '#FLUIDINITIALIZATION', success )
        IF( success )THEN
          DO layerCount = 1, WellData%NZ
            IF(layerCount<WellData%NZ)THEN
                z = 0.5* (WellData%LayersTops(layerCount) + WellData%LayersTops(layerCount+1) )
            ELSE
                z = WellData%LayersTops(layerCount)
            END IF
            IF(z>Problem%OWdepth)THEN
                !acima do contato OW.
                WellData%Porepressure(layerCount) = Problem%OWpressure - &
                    (z-Problem%OWdepth) * Problem%nwGrad
            ELSE
                !abaixo do contato OW.
                WellData%Porepressure(layerCount) = Problem%OWpressure - &
                    (z-Problem%OWdepth) * Problem%wGrad
            END IF
          END DO
        END IF

    END IF      !if do comando #NWB_STRUCTUREDMESH.

    !***************************************************************************
    ! DADOS DE PROPRIEDADES DOS FLUIDOS.
    !***************************************************************************
    CALL SearchInFile( Chanel_base, '#WATERPHASEPROPS', success )
    IF(success) THEN
        READ(Chanel_base,*) Problem%water_bwI
        READ(Chanel_base,*) Problem%water_pref                !kPa
        Problem%water_pref = Problem%water_pref * 1000.0      !kPa->Pa
        READ(Chanel_base,*) Problem%density_pref
        Problem%density_pref = Problem%density_pref * 1000.0  !kg/m³
        READ(Chanel_base,*) Problem%water_cw                  !1/kPa.
        Problem%water_cw  = Problem%water_cw  / 1000.0
        READ(Chanel_base,*) Problem%water_viscI               !cP.
        Problem%water_viscI = Problem%water_viscI * 0.001     !Pa.s
    ELSE
    END IF

    CALL SearchInFile( Chanel_base, '#OILPHASEPROPS', success )
    IF(success) THEN
        READ(Chanel_base,*) Problem%oil_bwI
        READ(Chanel_base,*) Problem%oil_viscI             !cP.
        Problem%oil_viscI = Problem%oil_viscI * 0.001     !Pa.s
    ELSE
    END IF

    !***************************************************************************
    ! DADOS DE PROPRIEDADES DOS FLUIDOS.
    !***************************************************************************
    CALL SearchInFile( Chanel_base, '#PRESSUREBC', success )

    IF(success)THEN
        ! #PRESSUREBC
        ! 9   15  11.667  <-internal: {K1,K2,Grad}
        ! 9   15  10.000  <-external: {K1,K2,Grad}
        READ(Chanel_base,*) WellData%Kint1,WellData%Kint2,WellData%gradInt
        READ(Chanel_base,*) WellData%Kext1,WellData%Kext2,WellData%gradExt
    ENDIF

    !***************************************************************************
    !PRE-FRACTURE DATA.
    !***************************************************************************
    !LOGICAL:: PreFractureExists
    !REAL(8),ALLOCATABLE:: PreFractureLengh(:) !<-length for each k_layer.
    !REAL(8),ALLOCATABLE:: PreFractureAperture(:) !<-aperture for each k_layer.
    !
    !#PREFRACTURE
    !{number of lines to be read}
    !{k1(most top),k2,length,aperture} <-line 1
    !...
    !
    CALL SearchInFile( Chanel_base, '#PREFRACTURE', success )
    IF(success)THEN

        WellData%PreFractureExists = .TRUE.
        ALLOCATE(WellData%PreFractureLengh(WellData%NZ),WellData%PreFractureAperture(WellData%NZ))
        WellData%PreFractureLengh=0.0
        WellData%PreFractureAperture = 0.0

        READ(Chanel_base,*)layerCount
        DO i=1,layerCount
          READ(Chanel_base,*)k1,k2,length,aperture
          WellData%PreFractureLengh(k1:k2) = length
          WellData%PreFractureAperture(k1:k2) = aperture
        END DO

    END IF

    !***************************************************************************
    !FRACTURE MODEL.
    !***************************************************************************

    !INTEGER:: K1,K2
    !REAL(8):: ST
    !REAL(8):: EyyMax
    
    CALL SearchInFile(Chanel_base,'#CDM_FRACTURE',success)
    IF(success) THEN
        READ(Chanel_base,*)nModels
        
        ALLOCATE( Problem%CDM(nModels) )
        
        DO i=1,nModels
          READ(Chanel_base,*)Problem%CDM(i)%K1,Problem%CDM(i)%K2,&
                             Problem%CDM(i)%ST,&
                             Problem%CDM(i)%EyyMax,&
                             Problem%CDM(i)%Gc
          Problem%CDM(i)%ST = Problem%CDM(i)%ST * 1000.0                  !kPa->Pa.
        END DO
        
    END IF
    
    CALL SearchInFile( Chanel_base, '#FRAC_CLOUSURE', success )
    IF(success) Problem%FracClosure = .TRUE.
    
    CALL SearchInFile( Chanel_base, '#FRAC_UNBALVOLUME', success )
    IF(success) Problem%FracUnbalVolume = .TRUE.

    !***************************************************************************
    !CONSTITUTIVE MODEL.
    !***************************************************************************
    !INTEGER,ALLOCATABLE:: ConsModel(:,:)              !(1,:) - model number.
    !                                                  !(2,:) - K1
    !                                                  !(3,:) - K2
    !REAL(8),ALLOCATABLE:: ConsModelParameters(:,:)    !(:,i) - Cons. model parameters.
    !
    !model types.
    !2 - Mohr-Coulomb
    !3 - CAM-CLAY (modified)
    CALL SearchInFile(Chanel_base,'#ELASTOPLASTIC_MODEL',success )
    IF(success) THEN
        READ(Chanel_base,*)nModels
          
        ALLOCATE(WellData%ConsModel(3,nModels),WellData%ConsModelParameters(10,nModels))    !10-arbitrário.
        WellData%ConsModel = 0
        WellData%ConsModelParameters = 0.0
          
        DO i=1,nModels
          READ(Chanel_base,*)text
          READ(Chanel_base,*)WellData%ConsModel(1,i)
          READ(Chanel_base,*)WellData%ConsModel(2,i),WellData%ConsModel(3,i)
          IF(WellData%ConsModel(1,i)==1) THEN                      !Damage (EXPLICITY SOLVER).
            WellData%ConsModelParameters(1,i) = 0.0
          END IF
          IF(WellData%ConsModel(1,i)==2) THEN
            READ(Chanel_base,*)WellData%ConsModelParameters(1:3,i) !Mohr-Coulomb.
            WellData%ConsModelParameters(1,i) = 1.0E3*WellData%ConsModelParameters(1,i)
          END IF
          IF(WellData%ConsModel(1,i)==3) THEN
            READ(Chanel_base,*)WellData%ConsModelParameters(1:6,i) !Modified Cam-clay.
            WellData%ConsModelParameters(1,i) = 1.0E3*WellData%ConsModelParameters(1,i)
            WellData%ConsModelParameters(6,i) = 1.0E3*WellData%ConsModelParameters(6,i)
          END IF
          IF(WellData%ConsModel(1,i)==4) THEN
            READ(Chanel_base,*)WellData%ConsModelParameters(1:5,i) !Damage-shear.
            WellData%ConsModelParameters(1,i) = 1.0E3*WellData%ConsModelParameters(1,i)
            WellData%ConsModelParameters(4,i) = 1.0E3*WellData%ConsModelParameters(4,i)
          END IF
          IF(WellData%ConsModel(1,i)==5) THEN
            READ(Chanel_base,*)WellData%ConsModelParameters(1:5,i) !Hyperbolic.
            WellData%ConsModelParameters(1,i) = 1.0E3*WellData%ConsModelParameters(1,i)
            WellData%ConsModelParameters(3,i) = 1.0E3*WellData%ConsModelParameters(3,i)
          END IF
        END DO
          
    END IF
      
    !SPECIAL PORPOUSE MODELS (undrained, liquefaction).
    !%ListOfUndrainedLayers(0:)
    !%UndrainedLayerParameters(1,:)
    CALL SearchInFile(Chanel_base,'#UNDRAINEDLAYERS',success)
    IF(success) THEN
        READ(Chanel_base,*)Cwater
        READ(Chanel_base,*)nlayers
        ALLOCATE(WellData%ListOfUndrainedLayers(0:WellData%NZ),&
                  WellData%UndrainedLayerParameters(1,WellData%NZ))
        WellData%ListOfUndrainedLayers = 0
        WellData%UndrainedLayerParameters(1,:) = Cwater/1000.0
        counter = 0
        DO i=1,nlayers
          READ(Chanel_base,*)k1,k2
          DO j=k1,k2
            counter = counter + 1
            WellData%ListOfUndrainedLayers(counter) = j
          END DO
        END DO
        WellData%ListOfUndrainedLayers(0) = counter
    END IF
      
    !TRANSVERSE ISOTROPIC MODEL
    !      LOGICAL:: TrasverseIsotropicModel
    !      REAL(8),ALLOCATABLE:: Dip(:),Strike(:)
    CALL SearchInFile(Chanel_base,'#TRASVERSEISOTROPIC',success)
    IF(success) THEN
          !Z1   Z2   DIP(DEG)   STRIKE(DEG)
          CALL ReadTable( Chanel_base, 3, table )
          ALLOCATE( WellData%InjPressureTable(3,SIZE(table(1,:))) )
          Problem%nSteps = SIZE(table(1,:))
          DO i=1,Problem%nSteps
              READ(table(1,i),*) WellData%InjPressureTable(1,i)
          END DO
    END IF
    
    !CREEP MODEL
    CALL SearchInFile(Chanel_base,'#SALTSSCREEPMODEL',success)
    IF(success) THEN
          !I1  I2  J1  J2   K1  K2  A  B  Ess  mu0   Qc  N2
          CALL ReadTable( Chanel_base, 9, table )
          n = SIZE(table(1,:))
          ALLOCATE( Problem%CreepTable(9,n) )
          DO i=1,n
              READ(table(1,i),*) Problem%CreepTable(1,i)
              READ(table(2,i),*) Problem%CreepTable(2,i)
              READ(table(3,i),*) Problem%CreepTable(3,i)
              READ(table(4,i),*) Problem%CreepTable(4,i)
              READ(table(5,i),*) Problem%CreepTable(5,i)
              READ(table(6,i),*) Problem%CreepTable(6,i)
              READ(table(7,i),*) Problem%CreepTable(7,i)
              READ(table(8,i),*) Problem%CreepTable(8,i)
              READ(table(9,i),*) Problem%CreepTable(9,i)
          END DO
    END IF
    
    CALL SearchInFile(Chanel_base,'#TIMESTEPS',success)
    IF(success) THEN
          !Time steps in minutes
          CALL ReadTable( Chanel_base, 1, table )
          n = SIZE(table(1,:))
          ALLOCATE( Problem%TimeSteps(n) )
          DO i=1,n
              READ(table(1,i),*) Problem%TimeSteps(i)
          END DO
    END IF

    CALL SearchInFile(Chanel_base,'#CRITICALTIMESTEPFACTOR',success)
    IF(success) THEN
        READ(Chanel_base,*)Problem%CriticalTimeStepFactor
    ELSE
        Problem%CriticalTimeStepFactor = 1.0
    END IF
  
    !***************************************************************************
    !TRANSIENT FLOW DATA.
    !***************************************************************************
      
    CALL SearchInFile(Chanel_base,'#WELLPRESSURECONSTRAINT',success)
    IF(success) THEN
        WellData%InjectionPressure = .TRUE.
        READ(Chanel_base,*)WellData%InjZ1, WellData%InjZ2
        !Steps of pressure.
        !DEALLOCATE( table )
        CALL ReadTable( Chanel_base, 3, table )
        ALLOCATE( WellData%InjPressureTable(3,SIZE(table(1,:))) )
        Problem%nSteps = SIZE(table(1,:))
        DO i=1,Problem%nSteps
            READ(table(1,i),*) WellData%InjPressureTable(1,i)
            READ(table(2,i),*) WellData%InjPressureTable(2,i)
            READ(table(3,i),*) WellData%InjPressureTable(3,i)
        END DO
        Problem%TotalTime = WellData%InjPressureTable(1,Problem%nSteps)
        WellData%InjPressureTable(1,:) = WellData%InjPressureTable(1,:) * 60
        WellData%InjPressureTable(2,:) = WellData%InjPressureTable(2,:) * 1000.0  !kPa -> Pa
        
        !WellData%InjX
        CALL SearchInFile(Chanel_base,'#PERFORATIONDEPTH',found)
        IF(found) THEN
            READ(Chanel_base,*)WellData%InjX
        ELSE
            IF( WellData%MeshType=='RAD' )THEN
                WellData%InjX = WellData%Radius(1) * (1.001)
            ELSE
                WellData%InjX = 1.0E-4
            END IF
        END IF
        
    END IF
    
    CALL SearchInFile(Chanel_base,'#WELLRATECONSTRAINT',success)
    IF(success) THEN
        WellData%InjectionRate = .TRUE.
        READ(Chanel_base,*)WellData%InjZ1, WellData%InjZ2
        !Steps of rates.
        !DEALLOCATE( table )
        CALL ReadTable( Chanel_base, 3, table )
        ALLOCATE( WellData%InjRateTable(3,SIZE(table(1,:))) )
        Problem%nSteps = SIZE(table(1,:))
        DO i=1,Problem%nSteps
            READ(table(1,i),*) WellData%InjRateTable(1,i)
            READ(table(2,i),*) WellData%InjRateTable(2,i)
            READ(table(3,i),*) WellData%InjRateTable(3,i)
        END DO
        Problem%TotalTime = WellData%InjRateTable(1,Problem%nSteps)
        WellData%InjRateTable(1,:) = WellData%InjRateTable(1,:) * 60
        WellData%InjRateTable(2,:) = WellData%InjRateTable(2,:) / (24*60*60)  !m³/d -> m³/sec
        
        !WellData%InjX
        CALL SearchInFile(Chanel_base,'#PERFORATIONDEPTH',found)
        IF(found) THEN
            READ(Chanel_base,*)WellData%InjX
        ELSE
            IF( WellData%MeshType=='RAD' )THEN
                WellData%InjX = WellData%Radius(1) * (1.001)
            ELSE
                WellData%InjX = 1.0E-4
            END IF
        END IF
        
        Problem%FracClosure = .TRUE.
        
    END IF
    
    CALL SearchInFile(Chanel_base,'#TRANSIENTFLOW',success )
    IF(success) THEN
        READ(Chanel_base,*)Problem%nSteps
        READ(Chanel_base,*)Problem%DeltaTime    !ratio
        READ(Chanel_base,*)Problem%TotalTime
        Problem%TotalTime = Problem%TotalTime*60.0
    END IF
      
    CALL SearchInFile(Chanel_base,'#WATERINJECTIONRATE',success)
    IF(success) THEN
        WellData%InjectionRate = .TRUE.
        READ(Chanel_base,*)WellData%InjX
        READ(Chanel_base,*)WellData%InjZ1, WellData%InjZ2
        READ(Chanel_base,*)WellData%Rate
        WellData%Rate = WellData%Rate / (24.0*60.0*60) !m³/dia to m³/sec.
    END IF

    CALL SearchInFile(Chanel_base,'#WATERINJECTIONPRESSURE',success)
    IF(success) THEN
        WellData%InjectionPressure = .TRUE.
        READ(Chanel_base,*)WellData%InjX
        READ(Chanel_base,*)WellData%InjZ1, WellData%InjZ2
        READ(Chanel_base,*)WellData%Pinj,WellData%Zinj
        WellData%Pinj = WellData%Pinj * 1000.0  !kPa -> Pa
    END IF
    
    ALLOCATE( Problem%MaxRate(Problem%nSteps) )
    Problem%MaxRate = 1.0E30
    CALL SearchInFile(Chanel_base,'#MAXRATELIMIT',success)
    IF(success) THEN
        !Steps of pressure.
        !DEALLOCATE( table )
        CALL ReadTable( Chanel_base, 2, table )
        DO i=1,SIZE(table(1,:))
            READ(table(:,i),*) j, rate
            Problem%MaxRate(j) = rate/24/60/60
        END DO
        
    END IF
    
    !STATIC PRESSURE
    ! FarFieldFlow(:)
    !                 INTEGER:: I1,I2,J1,J2,K1,K2
    !                 REAL(8):: J !(m3/d / Pa)
    !                 REAL(8):: Pest,Zref_Pest        !static pressure at a given Z depth.
    CALL SearchInFile(Chanel_base,'#OUTFLOW',success)
    IF(success) THEN
          CALL ReadTable( Chanel_base, 10, table )
          ALLOCATE( Problem%FarFieldFlow(SIZE(table(1,:))) )
          DO i=1,SIZE(table(1,:))
            READ(table(:,i),*)  Problem%FarFieldFlow(i)%I1,Problem%FarFieldFlow(i)%I2,&
                                Problem%FarFieldFlow(i)%J1,Problem%FarFieldFlow(i)%J2,&
                                Problem%FarFieldFlow(i)%K1,Problem%FarFieldFlow(i)%K2,&
                                Problem%FarFieldFlow(i)%DIR,&
                                Problem%FarFieldFlow(i)%J,&
                                Problem%FarFieldFlow(i)%Pest,Problem%FarFieldFlow(i)%Zref_Pest
            END DO
    END IF
    
    !CALL SearchInFile(Chanel_base,'#WATERINJECTIONSTAGES',WellData%InjectionPressure)
    !IF(WellData%InjectionPressure) THEN
    !  WellData%InjectionPressure = .TRUE.
    !  READ(Chanel_base,*)WellData%InjX
    !  READ(Chanel_base,*)WellData%InjZ1, WellData%InjZ2
    !  !Steps of pressure.
    !  !DEALLOCATE( table )
    !  CALL ReadTable( Chanel_base, 3, table )
    !  ALLOCATE( WellData%InjPressureTable(3,SIZE(table(1,:))) )
    !  DO i=1,SIZE(table(1,:))
    !      READ(table(1,i),*) WellData%InjPressureTable(1,i)
    !      READ(table(2,i),*) WellData%InjPressureTable(2,i)
    !      READ(table(3,i),*) WellData%InjPressureTable(3,i)
    !  END DO
    !  WellData%InjPressureTable(1,:) = WellData%InjPressureTable(1,:) * 60
    !  WellData%InjPressureTable(2,:) = WellData%InjPressureTable(2,:) * 1000.0  !kPa -> Pa
    !END IF
          
    !INITIAL STATE VARIATIONS
    CALL SearchInFile(Chanel_base,'#InitialStateVariation',success)
    IF(success) THEN
        CALL ReadTable( Chanel_base, 8, table )
        n = SIZE(table(1,:))
        ALLOCATE( Problem%InitialVariations(n) )
        DO i=1,n
            READ(table(1,i),*) Problem%InitialVariations(i)%I1
            READ(table(2,i),*) Problem%InitialVariations(i)%I2
            READ(table(3,i),*) Problem%InitialVariations(i)%J1
            READ(table(4,i),*) Problem%InitialVariations(i)%J2
            READ(table(5,i),*) Problem%InitialVariations(i)%K1
            READ(table(6,i),*) Problem%InitialVariations(i)%K2
            READ(table(7,i),*) Problem%InitialVariations(i)%DeltaP
            Problem%InitialVariations(i)%DeltaP = Problem%InitialVariations(i)%DeltaP * 1000.0
            READ(table(8,i),*) Problem%InitialVariations(i)%DeltaT
        END DO
    END IF
    
    !***************************************************************************
    !PHYSICS FLAG.
    !***************************************************************************        
      
    CALL SearchInFile(Chanel_base,'#FLOWFLAG',success )
    IF(success) Problem%Flow_flag = .TRUE.
    CALL SearchInFile(Chanel_base,'#STRUCTURALFLAG',success )
    IF(success) Problem%Structural_flag = .TRUE.
    CALL SearchInFile(Chanel_base,'#THERMALFLAG',success )
    IF(success) Problem%Thermal_flag = .TRUE.
      
    !***************************************************************************
    !OUTPUT CONTROL.
    !***************************************************************************        
    !INTEGER:: StructuralSolverPrintFreq
    !INTEGER:: ElementMonitor               !IX,IY,IZ
    !INTEGER:: NodeMonitor
      
    CALL SearchInFile(Chanel_base,'#OUTPUT_MONITOR',success )
    IF(success) THEN
          
        !Output Frequency (number of cycles).
        !Element Monitor.
        !Node Monitor.
          
        READ(Chanel_base,*)Problem%StructuralSolverPrintFreq
        READ(Chanel_base,*)Problem%ElementMonitor
        READ(Chanel_base,*)Problem%NodeMonitor
          
    ELSE
          
        Problem%StructuralSolverPrintFreq = 5000
        Problem%ElementMonitor = 1
        Problem%NodeMonitor = 1
          
    END IF
      
    !***************************************************************************
    !WELL BORE DATA. ( #DirectionalWell )
    !***************************************************************************
    !REAL(8):: WellZ,WellRadius
      
    CALL SearchInFile(Chanel_base,'#HorizontallWell',success )
    IF(success) THEN
        READ(Chanel_base,*)WellData%WellZ,WellData%WellRadius,WellData%WellLength,&
                           WellData%WellPermMult
        WellData%HorWellFlag = .TRUE.
    END IF
    
    CALL SearchInFile(Chanel_base,'#ProductionCasingDepth',success )
    IF(success) THEN
        READ(Chanel_base,*)WellData%Z_phase1
    END IF
    
    CALL SearchInFile(Chanel_base,'#VerticaDisplacementRestriction',success )
    IF(success)THEN
        WellData%VerticaDisplacementRestriction = .TRUE.
    ELSE
        WellData%VerticaDisplacementRestriction = .FALSE.
    END IF
    
    ! #ProductionCompletion
    !CHARACTER(LEN=10):: ProductionCompletionType			!<- OPEN , CASED or CASED-SLIDE
		!REAL(8):: ProductionCompletionPerfLength         !<- PERFORATION LENGTH (m)
		!INTEGER:: ProductionCompletionDistribution			  !<- 1-ONE, 2-TWO or 0-DISTRIBUTED
    CALL SearchInFile(Chanel_base,'#ProductionCompletion',success )
    IF(success) THEN
        READ(Chanel_base,*)WellData%ProductionCompletionType
        READ(Chanel_base,*)WellData%ProductionCompletionPerfLength
        READ(Chanel_base,*)text
        WellData%ProductionCompletionDistribution = 0
        IF(TRIM(text)=='ONE') WellData%ProductionCompletionDistribution = 1
        IF(TRIM(text)=='TWO') WellData%ProductionCompletionDistribution = 2
    ELSE
        WellData%ProductionCompletionType = 'OPEN'
        WellData%ProductionCompletionPerfLength = 0.01
        WellData%ProductionCompletionDistribution = 0
    END IF
    
    !DRILLING DATA for drilling workflow----------------------------------------
    !REAL(8):: OpenIntervalZ(2)
    !REAL(8),ALLOCATABLE:: MudWeightHistory(:)
    
    CALL SearchInFile(Chanel_base,'#DrillingInterval',success )
    IF(success) READ(Chanel_base,*) WellData%OpenIntervalZ(1), WellData%OpenIntervalZ(2)
    
    CALL SearchInFile(Chanel_base,'#MUDWEIGHTppg',success )
    IF(success) THEN
        READ(Chanel_base,*) WellData%DrillingSteps
        ALLOCATE( WellData%MudWeightHistory(WellData%DrillingSteps) )
        DO i=1,WellData%DrillingSteps
          READ(Chanel_base,*) WellData%MudWeightHistory(i)
        END DO
        WellData%MudWeightHistory = WellData%MudWeightHistory * (9810.0/8.34)   !Pa/m
    END IF
    
    !***************************************************************************
    !FLUID BOUNDARY CONDITION.
    !***************************************************************************
    Problem%InfiniteFluidBC = .TRUE. !Constant pressure B.C.
      
    CALL SearchInFile(Chanel_base,'#NoFlowBC',success)
    IF(success) Problem%InfiniteFluidBC = .FALSE.
      
    !***************************************************************************
    !PROPERTIES EDITOR.
    !***************************************************************************
    !REAL(8),ALLOCATABLE:: EditProps_values(:)
    !REAL(8),ALLOCATABLE:: EditProps_IJKinterval(:,:)
    !CHARACTER(LEN=40),ALLOCATABLE:: EditProps_PropsNames(:)
    CALL SearchInFile(Chanel_base,'#EDITprops',success)
    IF(success) THEN
        READ(Chanel_base,*)n
        ALLOCATE( Problem%EditProps_values(n), Problem%EditProps_IJKinterval(6,n), &
                  Problem%EditProps_PropsNames(n) )
        DO i=1,n
          READ(Chanel_base,*)Problem%EditProps_PropsNames(i),Problem%EditProps_values(i),&
                             Problem%EditProps_IJKinterval(1,i),Problem%EditProps_IJKinterval(2,i),&
                             Problem%EditProps_IJKinterval(5,i),Problem%EditProps_IJKinterval(6,i),&
                             Problem%EditProps_IJKinterval(3,i),Problem%EditProps_IJKinterval(4,i)
        END DO
    END IF
    
    !***************************************************************************
    !DAMAGE.
    !***************************************************************************
    !TYPE (DamageShear):: ShearStrainDamage
    CALL SearchInFile(Chanel_base,'#DamageShearStrain',success)
    IF(success) THEN
        READ(Chanel_base,*)n
        ALLOCATE( Problem%ShearStrainDamage(n) )
        DO i=1,n
          READ(Chanel_base,*) Problem%ShearStrainDamage(i)%K1,Problem%ShearStrainDamage(i)%K2,&
                              Problem%ShearStrainDamage(i)%LimitPlasticDeviatoricStrain,&
                              Problem%ShearStrainDamage(i)%ResidualDamage
        END DO
    END IF
    
    !TYPE(DamageStiffness):: Liquefaction
    CALL SearchInFile(Chanel_base,'#DamageLiquefaction',success)
    IF(success) THEN
        READ(Chanel_base,*)n
        ALLOCATE( Problem%Liquefaction(n) )
        DO i=1,n
          READ(Chanel_base,*) Problem%Liquefaction(i)%K1,Problem%Liquefaction(i)%K2,&
                              Problem%Liquefaction(i)%MaxMeanEffectiveStress,&
                              Problem%Liquefaction(i)%ShearDamage,&
                              Problem%Liquefaction(i)%ElasticDamage
          Problem%Liquefaction(i)%MaxMeanEffectiveStress = 1000.0*Problem%Liquefaction(i)%MaxMeanEffectiveStress !kPa->Pa
        END DO
    END IF
        
    CALL SearchInFile(Chanel_base,'#PERMENHANCEMENT',success)
    IF(success) THEN
        READ(Chanel_base,*)nModels
        
        ALLOCATE( Problem%Frac(nModels) )
        
        DO i=1,nModels
          READ(Chanel_base,*)Problem%Frac(i)%K1,Problem%Frac(i)%K2,&
                             Problem%Frac(i)%PermIncrementByPlasticDilation,&
                             Problem%Frac(i)%PermIncrementByPlasticShear
          Problem%Frac(i)%PermIncrementByPlasticDilation = Problem%Frac(i)%PermIncrementByPlasticDilation * 9.869D-16/0.0001
          Problem%Frac(i)%PermIncrementByPlasticShear = Problem%Frac(i)%PermIncrementByPlasticShear * 9.869D-16/0.0001
        END DO
    END IF
    
    CALL SearchInFile(Chanel_base,'#FORMATIONFLUIDDAMAGE',success)
    IF(success) THEN
          CALL ReadTable( Chanel_base, 4, Problem%FluidPermDamage )
    END IF
    
    !***************************************************************************
    !THERMAL MODEL.
    !***************************************************************************
    
    CALL SearchInFile(Chanel_base,'#ThermalModel',success)
    IF(success) THEN
        READ(Chanel_base,*)n
        ALLOCATE( Problem%ThermalProps(5,n),Problem%ThermalPropsKinterval(2,n)  )
        Problem%ThermalProps = 0.0
        DO i=1,n
          READ(Chanel_base,*) Problem%ThermalPropsKinterval(1,i),Problem%ThermalPropsKinterval(2,i),&
                              Problem%ThermalProps(1,i),Problem%ThermalProps(2,i),&
                              Problem%ThermalProps(3,i),Problem%ThermalProps(4,i),&
                              Problem%ThermalProps(5,i)
        END DO
    END IF
    
    CALL SearchInFile(Chanel_base,'#PerkinsGonzales',Problem%PerkinsGonzales)
    
    CALL SearchInFile(Chanel_base,'#InitialTemperature',success)
    IF(success) THEN
        READ(Chanel_base,*)Problem%TopTemperatureZ, Problem%TopTemperature, Problem%TemperatureGrad
    END IF

    CALL SearchInFile(Chanel_base,'#WaterInjectionTemperature',success)
    IF(success) THEN
        READ(Chanel_base,*)Problem%TinjZ, Problem%Tinj, Problem%TinjGrad
    ELSE
        Problem%Tinj = 0.0
        Problem%TinjZ = 0.0
        Problem%TinjGrad = 0.0
    END IF
    
    !EXPLICIT TEMPERATURE VARIATION
    !TYPE(TemperatureVariation),ALLOCATABLE:: ExplicitTemperatureChange(:)    !INTEGER:: Step
                                                                              !INTEGER:: I1,I2,J1,J2,K1,K2
                                                                              !REAL(8):: DeltaT
    CALL SearchInFile(Chanel_base,'#TemperatureVariation',success)
    IF(success) THEN
        CALL ReadTable( Chanel_base, 8, table )
        n = SIZE(table(1,:))
        ALLOCATE( Problem%ExplicitTemperatureChange(n) )
        DO i=1,n
            READ(table(1,i),*) Problem%ExplicitTemperatureChange(i)%Step
            READ(table(2,i),*) Problem%ExplicitTemperatureChange(i)%DeltaT
            READ(table(3,i),*) Problem%ExplicitTemperatureChange(i)%I1
            READ(table(4,i),*) Problem%ExplicitTemperatureChange(i)%I2
            READ(table(5,i),*) Problem%ExplicitTemperatureChange(i)%J1
            READ(table(6,i),*) Problem%ExplicitTemperatureChange(i)%J2
            READ(table(7,i),*) Problem%ExplicitTemperatureChange(i)%K1
            READ(table(8,i),*) Problem%ExplicitTemperatureChange(i)%K2
        END DO
    END IF
    
    !***************************************************************************
    !THERMAL MODEL.
    !***************************************************************************
    !LOGICAL,INTENT(IN OUT):: LongTermFlag
    !INTEGER,INTENT(IN OUT):: LongTermNumStep
    !REAL(8),INTENT(IN OUT):: LongTermDt
    
    CALL SearchInFile(Chanel_base,'#LONGTERMFLOW',success)
    IF(success) THEN
        Problem%LongTermFlag = .TRUE.
        READ(Chanel_base,*)Problem%LongTermNumStep
        READ(Chanel_base,*)Problem%LongTermDt       !days
        Problem%LongTermDt = Problem%LongTermDt * 24*60*60          !sec
    END IF
    
    !***************************************************************************
    
    
    !***************************************************************************
    !PIPE DATA
    !***************************************************************************
    !TYPE WellPipe
    !    REAL(8):: Diameter
    !    INTEGER,ALLOCATABLE:: WellPath(:,:)         !(:,# pipe node) - (4) DIR,i,j,k
                                                     !DIR: 1:X, 2:Y, 3:Z
    !END IF
    
    CALL SearchInFile(Chanel_base,'#VERTICALWELLPIPE',success)
    IF(success) THEN
          Problem%WellPipeModel = .TRUE.
          READ(Chanel_base,*)Problem%Pipe%Diameter
          READ(Chanel_base,*)k1,k2
          ALLOCATE( Problem%Pipe%WellPath(4,k2-k1+1) )
          Problem%Pipe%WellPath(1,:) = 3
          Problem%Pipe%WellPath(2,:) = 1
          Problem%Pipe%WellPath(3,:) = 1
          DO i=1,k2-k1+1
              Problem%Pipe%WellPath(4,i) = k1+(i-1)
          END DO
    END IF
        
    !***************************************************************************
    !TRANSIENT FLOW DATA.
    !***************************************************************************
    !Data for multiphase coupled solver

    !DISTRIBUTED FORCES (OVER AREA)
    ! TYPE(DistributedForce),ALLOCATABLE:: Forces(:)
    CALL SearchInFile(Chanel_base,'#FORCE',success)
    IF(success) THEN
        !force name, time series' name, surface's name, magnitude (kN/m²)
        CALL ReadTable( Chanel_base, 4, table )
        n = SIZE(table(1,:))
        ALLOCATE( Problem%Forces(n) )
        DO i=1,n
          READ(table(:,i),*) Problem%Forces(i)%Name,&
                              Problem%Forces(i)%TimeSerieName,&
                              Problem%Forces(i)%TriangularSurfaceName,&
                              Problem%Forces(i)%Amplitude
        END DO
    END IF

    !GEOMETRICAL (TRIANGULAR) SURFACES
    ! TYPE(TriangularSurface),ALLOCATABLE:: TriAreas(:)
    CALL SearchInFile(Chanel_base,'#TRISURFACES',success)
    IF(success) THEN
        !name, (X0,Y0,Z0), (X1,Y1,Z1), (X2,Y2,Z2)
        CALL ReadTable( Chanel_base, 10, table )
        n = SIZE(table(1,:))
        ALLOCATE( Problem%TriAreas(n) )
        DO i=1,n
          READ(table(:,i),*) Problem%TriAreas(i)%Name,&
                              Problem%TriAreas(i)%X1,&
                              Problem%TriAreas(i)%Y1,&
                              Problem%TriAreas(i)%Z1,&
                              Problem%TriAreas(i)%X2,&
                              Problem%TriAreas(i)%Y2,&
                              Problem%TriAreas(i)%Z2,&
                              Problem%TriAreas(i)%X3,&
                              Problem%TriAreas(i)%Y3,&
                              Problem%TriAreas(i)%Z3
        END DO
    END IF

    !TIME SERIES
    ! TYPE(TimeSerie),ALLOCATABLE:: TimeSeries(:)
    CALL SearchInFile(Chanel_base,'#TIMESERIES',success)
    IF(success) THEN
        n=0
        DO WHILE (TRIM(text) .NE. "<<")
          READ(Chanel_base,*)text
          IF(TRIM(text)=="<")n = n + 1
        END DO
        ALLOCATE( Problem%TimeSeries(n) )
        CALL SearchInFile(Chanel_base,'#TIMESERIES',success)
        READ(Chanel_base,*)text  !">>"
        DO i=1,n
          READ(Chanel_base,*)Problem%TimeSeries(i)%Name
          !time(sec), magnitude (relative)
          CALL ReadTable( Chanel_base, 2, table )
          ALLOCATE( Problem%TimeSeries(i)%Table( 2,SIZE(table(1,:)) ) )
          DO j=1,SIZE(table(1,:))
            READ(table(:,j),*) Problem%TimeSeries(i)%Table(1,j),&
                                Problem%TimeSeries(i)%Table(2,j)
          END DO
        END DO
    END IF

    !***************************************************************************
    
    Problem%WellData = WellData
          
    !***************************************************************************
    !OTHER OPTIONS.
    !***************************************************************************
    
    !#UPDATEGEOMETRY
    CALL SearchInFile(Chanel_base,'#UPDATEGEOMETRY',Problem%UpdateGeometry)
    
    !#REPORT
    CALL SearchInFile(Chanel_base,'#REPORT',Problem%Report)
    
    RETURN
  END SUBROUTINE ReadBaseFile

  SUBROUTINE ReadModelType(Chanel_base,ModelType)

      USE IO
    
      IMPLICIT NONE
      INTEGER,INTENT(IN)::  Chanel_base
      CHARACTER(LEN=*),INTENT(IN OUT):: ModelType
      
      LOGICAL:: wellModel,wellModel1,geotechnicalModel
    
      !WELL MODEL--------------------------------------------------------------
      CALL SearchInFile(Chanel_base,'#NWB_STRUCTUREDMESH',wellModel)
      CALL SearchInFile(Chanel_base,'#WELLBOREGRID',wellModel1)
      
      !GEOTECHNICAL------------------------------------------------------------
      CALL SearchInFile(Chanel_base,'#GEOTECHMODEL',geotechnicalModel)
      
      !MODEL TYPE--------------------------------------------------------------
      IF(wellModel)THEN
        ModelType = 'WellModel0'
      ELSEIF(wellModel1)THEN
        ModelType = 'WellModel1'
      ELSEIF(geotechnicalModel)THEN
        ModelType = 'Geotechnic'
      ELSE
        !Exception.
        
      END IF
      
  RETURN
  END SUBROUTINE
      
  !Imprime mensagem de erro e finaliza o codigo.
  SUBROUTINE Error( Message )

      USE IO

      CHARACTER(LEN=*),INTENT(IN):: Message

      CALL PrintMessage( .TRUE., .TRUE. , "ERROR MESSAGE:" )
      CALL PrintMessage( .FALSE., .FALSE., Message )
      CALL PrintMessage( .FALSE., .FALSE. , "FINILIZING" )
      STOP

  END SUBROUTINE


END MODULE ProblemDefinition
