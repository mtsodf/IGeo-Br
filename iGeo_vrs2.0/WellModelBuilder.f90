  
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE WellModelBuilder
  
  CONTAINS
  
  SUBROUTINE BuildWellModel(Chanel_base,Problem,GeomechanicalModel)

    USE IO
    USE MeshAuxiliares
    USE OMP_LIB  !bibliotecas OpenMP(R)
    USE DiscretizedModel
    USE ProblemDefinition

    !VARIABLES.
    IMPLICIT NONE
    INTEGER,INTENT(IN)::  Chanel_base
    TYPE(PD),INTENT(IN OUT):: Problem
    TYPE(Model),INTENT(IN OUT):: GeomechanicalModel

    INTEGER:: layerCount,i,k1,k2,n,nModels,node,counter,j
    REAL(8):: length,aperture
    LOGICAL:: success
    CHARACTER(LEN=320):: text
    REAL(8):: z
        
    REAL(8):: Cwater
    INTEGER:: nlayers
    
    CHARACTER(LEN=50),ALLOCATABLE:: table(:,:)

    !***************************************************************************
    ! WELLBORE DATA.
    !***************************************************************************  
    !GENERAL WELL GRID
    
    CALL SearchInFile( Chanel_base, '#WELLBOREGRID', success )
    IF(success)THEN
        
      READ(Chanel_base,*) text  ! >
      IF(TRIM(text).NE.">")THEN
          CALL Error( "> not found after #WELLBOREGRID." )
      END IF

      !Geological layers.
      READ(Chanel_base,*) Problem%ModelName
      READ(Chanel_base,*) Problem%WellData%NRAD,Problem%WellData%NTHETA,Problem%WellData%NZ,Problem%WellData%Ratio
      READ(Chanel_base,*) Problem%WellData%Radius(1), Problem%WellData%Radius(2)
      READ(Chanel_base,*) text  !#TABLEofProperties.
      IF(TRIM(text).NE."#TABLEofProperties")THEN
          CALL Error( "#TABLEofProperties not found." )
      END IF
      Problem%WellData%MeshType = 'RAD'
      Problem%MeshType = Problem%WellData%MeshType
            
      WRITE(text,*) '...number of layers:',Problem%WellData%NZ
      CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )

      ALLOCATE( Problem%WellData%LayersNames( Problem%WellData%NZ ),&
          Problem%WellData%LayersTops( Problem%WellData%NZ ),&
          Problem%WellData%Porosity( Problem%WellData%NZ ),&
          Problem%WellData%EffectivePerms( 2, Problem%WellData%NZ ),&
          Problem%WellData%YoungModulus( 2, Problem%WellData%NZ ),&
          Problem%WellData%Poisson( 2, Problem%WellData%NZ ),&
          Problem%WellData%Kg( Problem%WellData%NZ ),&
          Problem%WellData%StressTensor(6,Problem%WellData%NZ ),&
          Problem%WellData%Porepressure( Problem%WellData%NZ ) )
      DO layerCount = 1, Problem%WellData%NZ
          READ(Chanel_base,*) Problem%WellData%LayersNames(layerCount),&
              Problem%WellData%LayersTops(layerCount),&
              Problem%WellData%Porosity(layerCount),&
              Problem%WellData%EffectivePerms(1,layerCount),&
              Problem%WellData%EffectivePerms(2,layerCount),&
              Problem%WellData%YoungModulus(1,layerCount),&
              Problem%WellData%Poisson(1,layerCount),&
              Problem%WellData%Kg(layerCount),&
              Problem%WellData%StressTensor(1,layerCount),&
              Problem%WellData%StressTensor(2,layerCount),&
              Problem%WellData%StressTensor(3,layerCount),&
              Problem%WellData%StressTensor(4,layerCount),&
              Problem%WellData%StressTensor(5,layerCount),&
              Problem%WellData%StressTensor(6,layerCount),&
              Problem%WellData%Porepressure( layerCount )
      END DO  !loop nos layers.
      Problem%WellData%YoungModulus(2,:) = Problem%WellData%YoungModulus(1,:)
      Problem%WellData%Poisson(2,:) = Problem%WellData%Poisson(1,:)
            
      !Lendo dados do peso do fluido de perfuracao.
      CALL SearchInFile( Chanel_base, '#MudHistoryPPG', success )
      IF(success) THEN
        READ(Chanel_base,*) Problem%WellData%MudWeight                    !ppg
        Problem%WellData%MudWeight = Problem%WellData%MudWeight * (9810.0/8.34)   !Pa/m
      ELSE
        Problem%WellData%MudWeight = 1.028 * 9810.0                       !Pa/m
      END IF

      !Calculo da poropressao inicial nos layers do modelo.
      CALL SearchInFile( Chanel_base, '#FLUIDINITIALIZATION', success )
      IF( success )THEN
        DO layerCount = 1, Problem%WellData%NZ
          IF(layerCount<Problem%WellData%NZ)THEN
              z = 0.5* (Problem%WellData%LayersTops(layerCount) + Problem%WellData%LayersTops(layerCount+1) )
          ELSE
              z = Problem%WellData%LayersTops(layerCount)
          END IF
          IF(z>Problem%OWdepth)THEN
              !acima do contato OW.
              Problem%WellData%Porepressure(layerCount) = Problem%OWpressure - &
                  (z-Problem%OWdepth) * Problem%nwGrad
          ELSE
              !abaixo do contato OW.
              Problem%WellData%Porepressure(layerCount) = Problem%OWpressure - &
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
      READ(Chanel_base,*) Problem%WellData%Kint1,Problem%WellData%Kint2,Problem%WellData%gradInt
      READ(Chanel_base,*) Problem%WellData%Kext1,Problem%WellData%Kext2,Problem%WellData%gradExt
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

      Problem%WellData%PreFractureExists = .TRUE.
      ALLOCATE(Problem%WellData%PreFractureLengh(Problem%WellData%NZ),Problem%WellData%PreFractureAperture(Problem%WellData%NZ))
      Problem%WellData%PreFractureLengh=0.0
      Problem%WellData%PreFractureAperture = 0.0

      READ(Chanel_base,*)layerCount
      DO i=1,layerCount
        READ(Chanel_base,*)k1,k2,length,aperture
        Problem%WellData%PreFractureLengh(k1:k2) = length
        Problem%WellData%PreFractureAperture(k1:k2) = aperture
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
          
      ALLOCATE(Problem%WellData%ConsModel(3,nModels),Problem%WellData%ConsModelParameters(10,nModels))    !10-arbitrário.
      Problem%WellData%ConsModel = 0
      Problem%WellData%ConsModelParameters = 0.0
          
      DO i=1,nModels
        READ(Chanel_base,*)text
        READ(Chanel_base,*)Problem%WellData%ConsModel(1,i)
        READ(Chanel_base,*)Problem%WellData%ConsModel(2,i),Problem%WellData%ConsModel(3,i)
        IF(Problem%WellData%ConsModel(1,i)==1) THEN                      !Damage (EXPLICITY SOLVER).
          Problem%WellData%ConsModelParameters(1,i) = 0.0
        END IF
        IF(Problem%WellData%ConsModel(1,i)==2) THEN
          READ(Chanel_base,*)Problem%WellData%ConsModelParameters(1:3,i) !Mohr-Coulomb.
          Problem%WellData%ConsModelParameters(1,i) = 1.0E3*Problem%WellData%ConsModelParameters(1,i)
        END IF
        IF(Problem%WellData%ConsModel(1,i)==3) THEN
          READ(Chanel_base,*)Problem%WellData%ConsModelParameters(1:6,i) !Modified Cam-clay.
          Problem%WellData%ConsModelParameters(1,i) = 1.0E3*Problem%WellData%ConsModelParameters(1,i)
          Problem%WellData%ConsModelParameters(6,i) = 1.0E3*Problem%WellData%ConsModelParameters(6,i)
        END IF
        IF(Problem%WellData%ConsModel(1,i)==4) THEN
          READ(Chanel_base,*)Problem%WellData%ConsModelParameters(1:5,i) !Damage-shear.
          Problem%WellData%ConsModelParameters(1,i) = 1.0E3*Problem%WellData%ConsModelParameters(1,i)
          Problem%WellData%ConsModelParameters(4,i) = 1.0E3*Problem%WellData%ConsModelParameters(4,i)
        END IF
        IF(Problem%WellData%ConsModel(1,i)==5) THEN
          READ(Chanel_base,*)Problem%WellData%ConsModelParameters(1:5,i) !Hyperbolic.
          Problem%WellData%ConsModelParameters(1,i) = 1.0E3*Problem%WellData%ConsModelParameters(1,i)
          Problem%WellData%ConsModelParameters(3,i) = 1.0E3*Problem%WellData%ConsModelParameters(3,i)
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
      ALLOCATE(Problem%WellData%ListOfUndrainedLayers(0:Problem%WellData%NZ),&
                Problem%WellData%UndrainedLayerParameters(1,Problem%WellData%NZ))
      Problem%WellData%ListOfUndrainedLayers = 0
      Problem%WellData%UndrainedLayerParameters(1,:) = Cwater/1000.0
      counter = 0
      DO i=1,nlayers
        READ(Chanel_base,*)k1,k2
        DO j=k1,k2
          counter = counter + 1
          Problem%WellData%ListOfUndrainedLayers(counter) = j
        END DO
      END DO
      Problem%WellData%ListOfUndrainedLayers(0) = counter
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
        
    !***************************************************************************
    !TRANSIENT FLOW DATA.
    !***************************************************************************
        
    CALL SearchInFile(Chanel_base,'#TRANSIENTFLOW',success )
    IF(success) THEN
      READ(Chanel_base,*)Problem%nSteps
      READ(Chanel_base,*)Problem%DeltaTime    !ratio
      READ(Chanel_base,*)Problem%TotalTime
      Problem%TotalTime = Problem%TotalTime*60.0
    END IF
        
    CALL SearchInFile(Chanel_base,'#WATERINJECTIONRATE',Problem%WellData%InjectionRate)
    IF(Problem%WellData%InjectionRate) THEN
      Problem%WellData%InjectionRate = .TRUE.
      READ(Chanel_base,*)Problem%WellData%InjX
      READ(Chanel_base,*)Problem%WellData%InjZ1, Problem%WellData%InjZ2
      READ(Chanel_base,*)Problem%WellData%Rate
      Problem%WellData%Rate = Problem%WellData%Rate / (24.0*60.0*60) !m³/dia to m³/sec.
    END IF

    CALL SearchInFile(Chanel_base,'#WATERINJECTIONPRESSURE',Problem%WellData%InjectionPressure)
    IF(Problem%WellData%InjectionPressure) THEN
      Problem%WellData%InjectionPressure = .TRUE.
      READ(Chanel_base,*)Problem%WellData%InjX
      READ(Chanel_base,*)Problem%WellData%InjZ1, Problem%WellData%InjZ2
      READ(Chanel_base,*)Problem%WellData%Pinj,Problem%WellData%Zinj
      Problem%WellData%Pinj = Problem%WellData%Pinj * 1000.0  !kPa -> Pa
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
      READ(Chanel_base,*)Problem%WellData%WellZ,Problem%WellData%WellRadius,Problem%WellData%WellLength,&
                          Problem%WellData%WellPermMult
      Problem%WellData%HorWellFlag = .TRUE.
    END IF
      
    CALL SearchInFile(Chanel_base,'#ProductionCasingDepth',success )
    IF(success) THEN
      READ(Chanel_base,*)Problem%WellData%Z_phase1
    END IF
      
    CALL SearchInFile(Chanel_base,'#VerticaDisplacementRestriction',success )
    IF(success)THEN
      Problem%WellData%VerticaDisplacementRestriction = .TRUE.
    ELSE
      Problem%WellData%VerticaDisplacementRestriction = .FALSE.
    END IF
           
    ! #ProductionCompletion
    !CHARACTER(LEN=10):: ProductionCompletionType			!<- OPEN , CASED or CASED-SLIDE
		!REAL(8):: ProductionCompletionPerfLength         !<- PERFORATION LENGTH (m)
		!INTEGER:: ProductionCompletionDistribution			  !<- 1-ONE, 2-TWO or 0-DISTRIBUTED
    CALL SearchInFile(Chanel_base,'#ProductionCompletion',success )
    IF(success) THEN
      READ(Chanel_base,*)Problem%WellData%ProductionCompletionType
      READ(Chanel_base,*)Problem%WellData%ProductionCompletionPerfLength
      READ(Chanel_base,*)text
      Problem%WellData%ProductionCompletionDistribution = 0
      IF(TRIM(text)=='ONE') Problem%WellData%ProductionCompletionDistribution = 1
      IF(TRIM(text)=='TWO') Problem%WellData%ProductionCompletionDistribution = 2
    ELSE
      Problem%WellData%ProductionCompletionType = 'OPEN'
      Problem%WellData%ProductionCompletionPerfLength = 0.01
      Problem%WellData%ProductionCompletionDistribution = 0
    END IF
      
    !DRILLING DATA for drilling workflow----------------------------------------
    !REAL(8):: OpenIntervalZ(2)
    !REAL(8),ALLOCATABLE:: MudWeightHistory(:)
      
    CALL SearchInFile(Chanel_base,'#DrillingInterval',success )
    IF(success) READ(Chanel_base,*) Problem%WellData%OpenIntervalZ(1), Problem%WellData%OpenIntervalZ(2)
      
    CALL SearchInFile(Chanel_base,'#MUDWEIGHTppg',success )
    IF(success) THEN
      READ(Chanel_base,*) Problem%WellData%DrillingSteps
      ALLOCATE( Problem%WellData%MudWeightHistory(Problem%WellData%DrillingSteps) )
      DO i=1,Problem%WellData%DrillingSteps
        READ(Chanel_base,*) Problem%WellData%MudWeightHistory(i)
      END DO
      Problem%WellData%MudWeightHistory = Problem%WellData%MudWeightHistory * (9810.0/8.34)   !Pa/m
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
                            Problem%EditProps_IJKinterval(1:6,i)
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
      
    Problem%WellData = Problem%WellData
      
            
    !***************************************************************************
    !OTHER OPTIONS.
    !***************************************************************************
      
    !#UPDATEGEOMETRY
    CALL SearchInFile(Chanel_base,'#UPDATEGEOMETRY',Problem%UpdateGeometry)
      
    !#REPORT
    CALL SearchInFile(Chanel_base,'#REPORT',Problem%Report)
      
    RETURN
END SUBROUTINE

  END MODULE
  