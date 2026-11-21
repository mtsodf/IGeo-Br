!*****************************************************************************
!COMPUTATIONAL GEOMECHANICS
!Eng. Luis Carlos de Sousa Junior
!Computational Geomechanics
!*****************************************************************************

MODULE WF_MEDSolver

  CONTAINS

  ! #MEDSolver
  ! Runs a multiphysics dynamic solver
  SUBROUTINE RunPoromechanicsProblem( GeomechanicalModel,Matrices,&
                                      SolverParameters,&
                                      Problem,BaseFile,Path,VTKfileName,Chanel_base )
    USE omp_lib
    USE IO
    USE DiscretizedModel
    USE ProblemDefinition
    USE StiffnessMatrix
    USE NumericalUtilities
    USE VTK
    USE Report
    USE DynamicSolver

    IMPLICIT NONE

    TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
    TYPE (SolverMatricesContainer):: Matrices
    TYPE(DynamicSolverParameters):: SolverParameters
    TYPE(PD),INTENT(IN):: Problem
    CHARACTER(LEN=320),INTENT(IN):: BaseFile, &   !file with basic information about the case.
                                    Path          !work directory where the case file is located and where temporary files will
                                                  !be created.
    CHARACTER(LEN=*),INTENT(IN):: VTKfileName
    INTEGER,INTENT(IN):: Chanel_base
    
    TYPE(VariablesContainer):: solution
    TYPE(PoreFlowModel):: resFlowModel
    
    INTEGER:: vtk_fileChanel=30
    
    TYPE(DisplacementBC):: DispBC
    TYPE(StepReport):: reportData

    INTEGER:: nSteps                          !total number of steps
    INTEGER:: timeStepCounter                 !time step counter
    REAL(8):: time                            !cumulated time
    
    CHARACTER(LEN=320):: text
    INTEGER:: Chanel_log = 9
    
    INTEGER:: i,j    
    LOGICAL:: found

    
    INTEGER:: nst = 6                         !number of stress tensor.
    INTEGER:: nip                             !number of integration points.
    INTEGER:: numElements, numEquations, nunOpenFracElements, numFluidComponents
    
    LOGICAL:: Xl,Xp,Yl,Yp,Zl,Zp

    REAL(8),ALLOCATABLE:: fractureTIPsReports(:,:)                              !(1,:) - Top
                                                                                !(2,:) - Base
                                                                                !(3,:) - X
                                                                                !(4.:) - Fracture volume (m3)
    
    !***************************************************************************
    
    CALL PrintMessage( .FALSE., .FALSE., 'WORKFLOW: MULTIPHYSICS POROMECHANICS' )
    
    CLOSE(vtk_fileChanel)  !each time step will go on a separeted file.
    
    reportData%time = 0.0
    reportData%globalStep = 0
    reportData%cycles = 0
    reportData%stepResult = ''
    
    nip = GeomechanicalModel%NumIntegrationPoints
    numElements = GeomechanicalModel%NumElements
    numEquations = GeomechanicalModel%NumEquations
    
    ! Elastic stiffness
    IF(Problem%Structural_flag.AND..NOT.Matrices%Built) THEN
      CALL PrintMessage(.FALSE.,.FALSE.,'Building elastic stiffness matrix')
      CALL BuildElasticStiffness(GeomechanicalModel,Matrices,GeomechanicalModel%ActiveElements)
    END IF
    
    ! Memory Allocation
    CALL BuildVariablesContainer(GeomechanicalModel,Problem,solution)

    !**************************************************************************
    !BUILD POROUS FLOW MODEL                                                  *
    !**************************************************************************
    CALL PrintMessage(.FALSE.,.FALSE.,'Building porous flow model')

    ! Porous flow model
    ! Builds resFlow
    numFluidComponents = 1
    CALL BuildResFlowModel(GeomechanicalModel,numFluidComponents,ResFlowModel)

 !   
 !   !***************************************************************************
 !   ! INITIALIZE FLOW MODEL                                                    *
 !   !***************************************************************************
 !   CALL PrintMessage(.FALSE.,.FALSE.,'INITIALIZING MODELS')
 !   
 !   IF(Problem%Flow_flag) CALL InitializeFlowModel(GeomechanicalModel,Problem,FlowModel)
 !   
 !   IF(Problem%Thermal_flag) CALL BuildThermalModel(GeomechanicalModel,Problem,FlowModel,TempModel)
 !   
 !   !***************************************************************************
 !   ! FLOW BOUNDARY CONDITIONS                                                 *
 !   !***************************************************************************
 !   CALL PrintMessage(.FALSE.,.FALSE.,'PRESSURE B.C.')
 !   
 !   IF(Problem%Flow_flag)THEN
 !     
 !     !For Radial:
 !     !   X -> RAD
 !     !   Y -> THETA
 !     !   .TRUE. -> INFINITI B.C.    
 !     Xl = .FALSE.
 !     IF(FlowModel%RAD)THEN
 !       Xp = .TRUE.
 !     ELSE
 !       Xp = .FALSE.
 !       !Xp = .TRUE.
 !     END IF
 !     Yl = .FALSE.
 !     IF(FlowModel%RAD)THEN
 !       Yp = .FALSE.
 !     ELSE
 !       Yp = .TRUE.
 !     END IF
 !     Zl = .FALSE.
 !     Zp = .FALSE.
 !     
 !     IF( .NOT.ALLOCATED(Problem%FarFieldFlow) ) CALL ExternalPressureBC(GeomechanicalModel,Xl,Xp,Yl,Yp,Zl,Zp,FlowModel)
 !     
 !     IF(WellData%InjectionPressure) THEN
 !       IF(FlowModel%RAD.AND.WellData%ProductionCompletionType=='OPEN')THEN
 !         !Source at well wall.
 !         CALL InternalPressureBC(FlowModel,WellData,GeomechanicalModel,0)
 !         !Mechanical BC.
 !         !CALL FreeMechanicalBC(GeomechanicalModel,Problem,StructuralBC)
 ! !stop !not implemented.
 !       ELSEIF(FlowModel%RAD.AND.WellData%ProductionCompletionType=='CASED')THEN 
 !         !Perforation.
 !         CALL ProductionCompletionAndInjection(GeomechanicalModel,Problem,FlowModel)
 !       ELSE
 !         !Source at well wall.
 !         CALL InternalPressureBC(FlowModel,WellData,GeomechanicalModel,0)
 !       END IF
 !     END IF
 !   
 !   END IF
 !   
 !   IF(Problem%Thermal_flag) THEN
 !     !CALL PrintMessage(.FALSE.,.FALSE.,'TEMPERATURE B.C.')
 !     CALL TemperatureBC(TempModel,Problem,FlowModel,GeomechanicalModel)
 !     ALLOCATE( CumFlowVolume(GeomechanicalModel%GridDivision(2)+1) )
 !     CumFlowVolume = 0.0
 !   END IF
 !   
 !   !***************************************************************************
 !   !INITIAL STATE VARIATIONS                                                  *
 !   !***************************************************************************
 !   
 !   IF( Problem%Structural_flag .AND. ALLOCATED(Problem%InitialVariations) )THEN
 !       CALL PrintMessage(.FALSE.,.FALSE.,'INITIAL STATE VARIATIONS')
 !       
 !       CALL ApplyInitialStateVariations(GeomechanicalModel,Problem,&
 !           Solution%Porepressure,Solution%Temperature)
 !       TempModel%RockTemperature = Solution%Temperature
 !       
 !       Load = 0.0
 !       deltaP = Solution%Porepressure - Solution%Porepressure0
 !       CALL PoroelasticLoad(GeomechanicalModel,deltaP,Load)                              ! -> Load(deltaP)

 !       deltaT = Solution%Temperature - Solution%Temperature0
 !       CALL ThermalStrainLoad(GeomechanicalModel,deltaT,Load(0:GeomechanicalModel%NumEquations))
 !       
 !       StructuralBC_modif = StructuralBC
 !       !Unlok model top
 !       CALL UnlockModelTop(GeomechanicalModel,StructuralBC_modif)
 !       
 !       IF(TRIM(GeomechanicalModel%SolverType)=='DynamicFEM')THEN
 !           !
 !           !CALL DynamicSolver( GeomechanicalModel,Matrices,&                               ! -> u,eps,sigma
 !           !                    SolverParameters, StructuralBC_modif%ListOfRestrictedDOFs,&
 !           !                    DispBC,Load(1:),displacement,&
 !           !                    Problem%StructuralSolverPrintFreq,&
 !           !                    Problem%ElementMonitor,&
 !           !                    Problem%NodeMonitor,&
 !           !                    Chanel_log,GeomechanicalModel%ActiveElements)
 !           !!converged:   0 - success
 !           !!             1 - not converged
 !           !!             2 - failed
 !           !IF((SolverParameters%SolverConvergence==2))THEN                                 ! -> failed to converge -> break time step
 !           !  converged = .FALSE.
 !           !END IF          !=0 or =1
 !           !converged = .TRUE.
 !        ELSE  !GeomechanicalModel%SolverType=='Viscoplastic'
 !           cg_tol = SolverParameters%linTol
 !           cg_limit = SolverParameters%maxLinIters
 !           tol = SolverParameters%viscoTol
 !           limit = SolverParameters%maxViscoIters
 !           CALL ViscoSolver(GeomechanicalModel,reportData,Matrices%Stiff%a, Matrices%Stiff%ia, &
 !                            Matrices%Stiff%ja,&
 !                            StructuralBC_modif%ListOfRestrictedDOFs,Load(1:),DeltaT,&
 !                            cg_tol,cg_limit,tol,limit,&
 !                            displacement)
 !        END IF
 !       
 !        CALL UpdateModels(Solution,GeomechanicalModel,FlowModel,HFmodel,WellData)

 !   END IF
 !   
 !   !***************************************************************************
 !   ! BUILD HYDRAULIC FRACTURE FLOW MODEL.                                     *
 !   !***************************************************************************
 !   
 !   CALL PrintMessage( .FALSE., .FALSE., 'BUILDING HYDRAULIC FRACTURE MODEL.' )
 !   
 !   CALL BuildHydraulicFractureModel( HFmodel,Problem,GeomechanicalModel,&
 !                                     StructuralBC%ListOfRestrictedDOFs)
 !   
 !   !***************************************************************************
 !   !WATER INJECTION.                                                          *
 !   !***************************************************************************
 !   
 !   CALL PrintMessage(.FALSE.,.FALSE.,'INJECTION RATE')
 !   IF(WellData%InjectionRate) CALL InjectionRate(FlowModel,WellData,GeomechanicalModel)
 !   
 !   !***************************************************************************
 !   ! OUTPUT INITIAL STATE.                                                    *
 !   !***************************************************************************
 !   
 !   OutStep=OutStep+1
 !   CALL WriteStepResults2VTK(vtk_fileChanel,VTKfileName,Problem,GeomechanicalModel,&
 !                                TRIM(Problem%ModelName),OutStep)
 !   
 !   IF( WellData%MeshType=='RAD' .AND. Problem%Report ) THEN
 !       CALL WriteReport(vtk_fileChanel,TRIM(Path)//'/'//TRIM(WellData%WellName),&
 !                                   Problem,GeomechanicalModel,FlowModel,OutStep)
 !   END IF
 !   
 !   !***************************************************************************
 !   ! COUPLED MULTIPHASE SIMULATION.                                           *
 !   !***************************************************************************
 !   
 !   CALL PrintMessage( .FALSE., .FALSE., 'TRANSIENT SOLVER.' )
 !   
 !   du = 0.0
 !   tensor = 0.
 !   
 !   nunOpenFracElements = 0
 !   oldNum = 0
 !   timeStepCounter = 0
 !   time = 0.0
 !   
 !   !Print Time Steps.
 !   totalTime = Problem%TotalTime
 !   nSteps = Problem%nSteps
 !   ALLOCATE(printTimes(0:nSteps))
 !   printTimes = 0.
 !   IF( SIZE(WellData%InjPressureTable(1,:))>0 )THEN
 !       DO i=1,nSteps
 !           printTimes(i) = WellData%InjPressureTable(1,i)
 !       END DO
 !   ELSEIF( SIZE(WellData%InjRateTable(1,:))>0 )THEN
 !       DO i=1,nSteps
 !           printTimes(i) = WellData%InjRateTable(1,i)
 !       END DO
 !   ELSE   
 !       sn = totalTime
 !       n = nSteps - 1
 !       q = Problem%DeltaTime
 !       a1 = sn * (q-1.0) / (q**n-1)
 !       DO i=1,nSteps
 !           printTimes(i) = a1 * q**(i-1)
 !       END DO
 !   END IF
 !   
 !   ALLOCATE(fractureTIPsReports(4,nSteps))
 !   fractureTIPsReports = 0.0
 !   
 !   ALLOCATE(stepRate(nSteps))
 !   stepRate = 0.0

 !   ! GLOBAL TIME STEP ========================================================
 !   CALL SaveVariables(GeomechanicalModel,variablesInitial)
 !   DO
 !     CALL PrintStepReportHEAD()
 !     
 !     timeStepCounter = timeStepCounter + 1

 !     reportData%globalStep = timeStepCounter
 !     reportData%flowStep = 0
 !     reportData%structuralStep = 0
 !     reportData%fractureStep = 0
 !     reportData%linIters = 0
 !     reportData%nonLinIters = 0
 !     reportData%stepResult = ''
 !   
 !     !TIME STEP---------------------------------------------------------------
 !     couplingTimeStep = printTimes(timeStepCounter)-printTimes(timeStepCounter-1)
 !     
 !     ! UPDATE BOUNDARY CONDITIONS---------------------------------------------
 !     IF(WellData%InjectionPressure)THEN
 !         CALL InternalPressureBC(FlowModel,WellData,GeomechanicalModel,&
 !                                   timeStepCounter)
 !     ELSE
 !         CALL InternalRateBC(FlowModel,WellData,GeomechanicalModel,&
 !                               timeStepCounter)
 !     END IF
 !     Solution%Porepressure = GeomechanicalModel%Porepressure
 !     
 !     preLoad = 0.0
 !     !EXPLICIT TEMPERATURE CHANGE
 !     !CALL ExplicitTemperatureChange(GeomechanicalModel, Problem, timeStepCounter,&
 !     !                                  one, Solution%Temperature)
 !     !deltaT = Solution%Temperature - Solution%Temperature0
 !     !CALL ThermalStrainLoad(GeomechanicalModel,deltaT,preLoad(0:GeomechanicalModel%NumEquations))
 !     !Solution%Temperature0 = Solution%Temperature

 !     ! TRANSIENT COUPLED SOLUTION---------------------------------------------
 !     coupledTime = 0.

 !     iters = 0
 !     DO WHILE( coupledTime < couplingTimeStep ) !COUPLING LOOPING------
 !       iters = iters + 1

 !       !CRITICAL TIME STEP####################################################
 !       !Flow
 !       dt = FlowCriticalTimeStep(FlowModel,HFmodel)
 !       dt = dt * Problem%CriticalTimeStepFactor

 !       !Thermal
 !       !...

 !       IF( dt >= couplingTimeStep) dt = couplingTimeStep/2.0

 !       !Update time step
 !       time = time + dt
 !       coupledTime = coupledTime + dt

 !       IF( coupledTime > couplingTimeStep )THEN
 !         time = time - dt
 !         coupledTime = coupledTime - dt 
 !         
 !         dt = couplingTimeStep - coupledTime
 !         
 !         time = time + dt
 !         coupledTime = coupledTime + dt
 !       END IF

 !       !FLUID FLOW TIME STEP##################################################
 !       IF(Problem%Flow_flag)THEN

 !         reportData%flowStep = reportData%flowStep + 1          
 !         CALL UpdateTransmissibility(Problem,dt,FlowModel,GeomechanicalModel)

 !         IF(HFmodel%ContainsFracModel) CALL UpdateHFtransmissibilityMatrix(HFmodel,GeomechanicalModel)
 !         CALL FlowStepSolution(Problem%nthreads,FlowModel,HFmodel,GeomechanicalModel,WellData,&             ! ->Porepressure, WaterSaturation
 !                               dt,Problem%MaxRate(timeStepCounter),reportData,Solution%Porepressure,&
 !                               Solution%Sw)
 !         stepRate(timeStepCounter) = stepRate(timeStepCounter) + FlowModel%AverageTotalRate * dt
 !         
 !         ! Updating hydraulic fracture fluid pressure - catch from nodes--------
 !         IF(HFmodel%ContainsFracModel) CALL UpdateHFfluidPressure(GeomechanicalModel,HFmodel,Solution%Porepressure)        ! P -> Hydraulic Fracture Model
 !       
 !       END IF

 !       !THERMAL TIME STEP#####################################################
 !       IF(Problem%Thermal_flag)THEN
 !           CALL TemperatureBC(TempModel,Problem,FlowModel,GeomechanicalModel)
 !           !CALL ComputeFlowVelocities(GeomechanicalModel,FlowModel,&
 !           !                            Solution%Porepressure,Solution%Porepressure0)
 !           IF( TempModel%Numerical )THEN
 !             !CALL UpdateTemperature( TempModel,GeomechanicalModel,FlowModel,&
 !             !                        couplingTimeStep )
 !             CALL ThermalStepSolution( Problem%nthreads,GeomechanicalModel,FlowModel,TempModel,&
 !                                       couplingTimeStep,TempModel%FluidTemperature,&
 !                                       reportData )
 !             Solution%Temperature = TempModel%RockTemperature
 !           ELSE
 !             CALL PerkinsGonzalesThermalZone( GeomechanicalModel, TempModel, FlowModel, &
 !                                               couplingTimeStep, cumFlowVolume )
 !             Solution%Temperature = TempModel%RockTemperature
 !           END IF
 !       ELSE
 !           CALL ExplicitTemperatureChange(GeomechanicalModel, Problem, timeStepCounter,&
 !                                            dt/couplingTimeStep, Solution%Temperature)
 !           TempModel%RockTemperature = Solution%Temperature
 !       END IF
 !       
 !       !STRUCTURAL SIMULATION#################################################

 !       converged = .TRUE.
 !       IF(Problem%Structural_flag)THEN
 !     
 !         Load = preLoad
 !           
 !         !UPDATE POREPRESSURE LOAD--------------------------------------------
 !         !Load = 0.0
 !         deltaP = Solution%Porepressure - Solution%Porepressure0
 !         CALL PoroelasticLoad(GeomechanicalModel,deltaP,Load)                              ! -> Load(deltaP)

 !         !TEMPERATURE--------------------------------------------------------
 !         deltaT = Solution%Temperature - Solution%Temperature0
 !         CALL ThermalStrainLoad(GeomechanicalModel,deltaT,Load(0:GeomechanicalModel%NumEquations))

 !         !HYDRAULIC FRACTURE GEOMETRY AND FRACTURE LOAD-----------------------
 !         IF(HFmodel%ContainsFracModel)THEN
 !           FractureLoad(:,2) = 0.0
 !           CALL UpdateHydraulicFracture(HFmodel,GeomechanicalModel,FractureLoad(:,2))                  ! -> + Load(frac pressure)
 !           Load = Load + FractureLoad(:,2)
 !           FractureLoad(:,1) = FractureLoad(:,2)
 !         END IF

 !         !WELL WALL LOAD------------------------------------------------------
 !         !Completion: OPEN
 !         !IF(WellData%ProductionCompletionType=='OPEN') CALL WellwallLoad(GeomechanicalModel,Load)
 !         !STRUCTURAL SOLVER---------------------------------------------------
 !         j = 0                            
 !         StructuralBC_modif%ListOfRestrictedDOFs = StructuralBC%ListOfRestrictedDOFs   !StructuralBC_modif  used only in the next loop
 !         StructuralBC_modif%ListOfSupportedDOFs = StructuralBC%ListOfSupportedDOFs
 !         fractureIteration = .TRUE.
 !         DO WHILE( fractureIteration )                                  !To avoid material transpassing for fracture wall
 !           j = j + 1         
 !           reportData%fractureStep = reportData%fractureStep + 1
 !           
 !           IF(HFmodel%ContainsFracModel) nunOpenFracElements = COUNT(GeomechanicalModel%IsFractured)
 !           CALL SaveVariables(GeomechanicalModel,variables)                                                                      
 !           
 !           IF(MAXVAL(ABS(Load(1:)))<0.01)EXIT
 !           
 !           !ELATIC RUN -> adjust fracture boundary conditions--------------
 !           DO WHILE(fractureIteration)
 !               reportData%structuralStep =  reportData%structuralStep + 1
 !               
 !               IF(TRIM(GeomechanicalModel%SolverType)=='DynamicFEM')THEN
 !                 !CALL DynamicSolver( ...
 !               ELSE  !GeomechanicalModel%SolverType=='Viscoplastic'
 !                   cg_tol = SolverParameters%linTol
 !                   cg_limit = SolverParameters%maxLinIters
 !                   tol = SolverParameters%viscoTol
 !                   limit = 1    !SolverParameters%maxViscoIters
 !                   CALL ViscoSolver(GeomechanicalModel,reportData,Matrices%Stiff%a, Matrices%Stiff%ia, &
 !                                    Matrices%Stiff%ja,&
 !                                    StructuralBC_modif%ListOfRestrictedDOFs,Load(1:),DeltaT,&
 !                                    cg_tol,cg_limit,tol,limit,&
 !                                    displacement)
 !               END IF
 !               !CHECK FRACTURE BOUNDARY CONDITION
 !               CALL CheckFractureDisplacement(GeomechanicalModel,HFmodel,StructuralBC_modif,fractureIteration)
 !               IF(fractureIteration)THEN
 !                   CALL RecoverVariables(GeomechanicalModel,variables)
 !                   j = StructuralBC_modif%ListOfRestrictedDOFs(0)
 !                   displacement(StructuralBC_modif%ListOfRestrictedDOFs(1:j)) = 0.0
 !               END IF
 !           END DO
 !           CALL RecoverVariables(GeomechanicalModel,variables)
 !           GeomechanicalModel%Displacement( StructuralBC_modif%ListOfRestrictedDOFs(1:StructuralBC_modif%ListOfRestrictedDOFs(0)) ) = 0.0      !(08.08.2021) force zero displacements for new b.c.
 !           !VISCOPLASTIC RUN-----------------------------------------------
 !           fractureIteration = .TRUE.
 !           DO WHILE(fractureIteration)
 !               IF(TRIM(GeomechanicalModel%SolverType)=='DynamicFEM')THEN
 !                 !CALL DynamicSolver( ...
 !               ELSE  !GeomechanicalModel%SolverType=='Viscoplastic'
 !                   cg_tol = SolverParameters%linTol
 !                   cg_limit = SolverParameters%maxLinIters
 !                   tol = SolverParameters%viscoTol
 !                   limit = SolverParameters%maxViscoIters
 !                   CALL ViscoSolver(GeomechanicalModel,reportData,Matrices%Stiff%a, Matrices%Stiff%ia, &
 !                                    Matrices%Stiff%ja,&
 !                                    StructuralBC_modif%ListOfRestrictedDOFs,Load(1:),DeltaT,&
 !                                    cg_tol,cg_limit,tol,limit,&
 !                                    displacement)
 !               END IF
 !               !CHECK FRACTURE BOUNDARY CONDITION
 !               IF(HFmodel%ContainsFracModel)THEN
 !                   CALL CheckFractureDisplacement(GeomechanicalModel,HFmodel,StructuralBC_modif,fractureIteration)
 !                   IF(fractureIteration)THEN
 !                       CALL RecoverVariables(GeomechanicalModel,variables)
 !                       j = StructuralBC_modif%ListOfRestrictedDOFs(0)
 !                       displacement(StructuralBC_modif%ListOfRestrictedDOFs(1:j)) = 0.
 !                   END IF
 !               ELSE
 !                   fractureIteration = .FALSE.
 !               END IF
 !           END DO
 !           ! FRACTURE PROPAGATION---------------------------------------------
 !           fractureIteration = .FALSE.
 !           IF(HFmodel%ContainsFracModel)THEN
 !             !  CALL CalculateInternalForceVector(GeomechanicalModel)
 !             du(1:) = GeomechanicalModel%Displacement - variables%Displacement
 !             tensor = GeomechanicalModel%EffectiveStressTensor - variables%EffectiveStressTensor                   !effective stress increment since previous step
 !             !tensor = GeomechanicalModel%EffectiveStressTensor
 !             CALL UpdateStrainWork(GeomechanicalModel,tensor,du)
 !             tensor = GeomechanicalModel%EffectiveStressTensor - variablesInitial%EffectiveStressTensor            !effective stress increment since initial
 !             CALL PropagateFracture_TractionDamage(HFmodel,GeomechanicalModel,&
 !                                                   StructuralBC%ListOfRestrictedDOFs,tensor)!,&                    !Obs.: não usar StructuralBC_modif para não restringir os nós. O StructuralBC_modif pára o nó
 !                                                   !Load(0:))                                                      !      quando o deslocamento fica nulo, mas não impõe a restrição, liberando a fratura para dilatar.
 !             !du = displacement
 !           END IF
 !           IF( COUNT(GeomechanicalModel%IsFractured) > nunOpenFracElements ) THEN
 !               fractureIteration=.TRUE.
 !               CALL RecoverVariables(GeomechanicalModel,variables)
 !           END IF
 !                 
 !           preLoad = 0.0   !due explicity temperature change.
 !           
 !         END DO!fracture propagation loop-------------------------------------

 !         !UPDATE PORE COMPRESSIBILITY-----------------------------------------
 !         deltaVolumetricStrain = deltaVolumetricStrain - volumetricStrain

 !       END IF!STRUCTURAL SIMULATION###########################################

 !       maxDisp = MAXVAL(ABS(GeomechanicalModel%Displacement))
 !       MaxDT = MAXVAL(ABS(Solution%Temperature - Solution%Temperature0))
 !       maxDP = MAXVAL(ABS(Solution%Porepressure - Solution%Porepressure0))
 !       
 !       CALL UpdateVariablesContainers(Solution)

 !       !######################################################################
 !       IF(converged)THEN
 !         CALL UpdateModels(Solution,GeomechanicalModel,FlowModel,HFmodel,WellData)
 !       ELSE
 !         !Structural model did not converged -> Break time step.
 !         CALL PrintMessage(.FALSE.,.FALSE.,"...structural model failed to converge.")
 !         CALL PrintMessage( .FALSE., .FALSE., 'SIMULATION ENDED.' )
 !         CALL SLEEP(3)

 !       END IF
 !       
 !     END DO!COUPLING LOOPING--------------------------------------------------

 !     ! OUTPUT-----------------------------------------------------------------
 !     IF(WellData%InjectionPressure)THEN
 !         stepRate(timeStepCounter) = stepRate(timeStepCounter) / coupledTime
 !     ELSE !RATE CONTROLED
 !         IF( WellData%InjRateTable(2,timeStepCounter) .NE. 0.0 )THEN
 !             stepRate(timeStepCounter) = WellData%InjRateTable(2,timeStepCounter)
 !         ELSE
 !             stepRate(timeStepCounter) = 0.0
 !         END IF
 !     END IF
 !     
 !     OutStep=OutStep+1
 !     CALL WriteStepResults2VTK(vtk_fileChanel,VTKfileName,Problem,GeomechanicalModel,&
 !                                   TRIM(Problem%ModelName),OutStep)
 !     
 !     !General FEM output file
 !     WRITE(100)GeomechanicalModel%Porepressure
 !     WRITE(100)GeomechanicalModel%Displacement
 !     WRITE(100)GeomechanicalModel%EffectiveStressTensor
 !     WRITE(100)GeomechanicalModel%StrainTensor
 !     WRITE(100)GeomechanicalModel%PlasticStrainTensor

 !     CALL ReportFractureTips(GeomechanicalModel,HFmodel,fractureTIPsReports(:,timeStepCounter))
 !     CALL SimulationReport(MaxDT,maxDP,MaxDisp,StepRate(timeStepCounter))
 !     
 !     ! STOP CRITERIA----------------------------------------------------------
 !     ! Global time step
 !     IF( timeStepCounter>=nSteps )THEN
 !       CALL PrintMessage(.FALSE.,.FALSE.,'...Simulation ended.')
 !       EXIT
 !     END IF
 !     
 !     !------------------------------------------------------------------------
 !     
 !   END DO  !global time step
 !   
 !   !==========================================================================
 !   !Fracture TIP report
 !   CALL ReportFractureTipsAllSteps(fractureTIPsReports)
 !   IF( WellData%MeshType=='RAD' .AND. Problem%Report ) THEN
 !     CALL WriteReport(vtk_fileChanel,TRIM(Path)//'/'//TRIM(WellData%WellName),&
 !                                 Problem,GeomechanicalModel,FlowModel,OutStep)
 !   END IF
 !   
 !   !==========================================================================
    
    CALL PrintMessage( .FALSE., .FALSE., 'WORKFLOW FINISHED.' )
    
  RETURN
  END SUBROUTINE
       
  SUBROUTINE BuildVariablesContainer(GeomechanicalModel,Problem,Solution)
  
    USE DiscretizedModel
    USE ProblemDefinition
  
    IMPLICIT NONE

    TYPE(Model),INTENT(IN):: GeomechanicalModel
    TYPE(PD),INTENT(IN):: Problem
    TYPE(VariablesContainer),INTENT(IN OUT):: Solution
  
    INTEGER:: numNodes,numElements,numFluidComponents
    INTEGER:: num(8),i
    
    numNodes = GeomechanicalModel%NumNodes
    numElements = GeomechanicalModel%numElements
    numFluidComponents = 1
    
    ALLOCATE(&
            Solution%Porepressure0(numElements),&
            Solution%Porepressure1(numElements),&
            Solution%Temperature0(numElements),&
            Solution%Temperature1(numElements),&
            Solution%Displacement0(3,numNodes),&
            Solution%Displacement1(3,numNodes)&
            )

    ! Initial porepressure and temperature initialized from nodes
    DO i=1,numElements
      num = GeomechanicalModel%ElementNodes(:,i)
      Solution%Porepressure0(i) = 1./8.*SUM(GeomechanicalModel%Porepressure(num))
      Solution%Temperature0(i) = 1./8.*SUM(GeomechanicalModel%Temperature(num))
    END DO
  
    Solution%Porepressure1 = Solution%Porepressure0
    Solution%Temperature1 = Solution%Temperature0
    Solution%Displacement0 = 0.0
    Solution%Displacement1 = 0.0

    !     REAL(8),ALLOCATABLE:: Sw0(:,:),Sw(:,:)                      !cells - for each componnent
    !     REAL(8),ALLOCATABLE:: FractureNodes(:,:)                    !4 nodes each
    !     REAL(8),ALLOCATABLE:: FracturePressure(:,:)                 !fracture cell  
    !     REAL(8),ALLOCATABLE:: FractureSat(:,:)                      !fracture cell - for each componnent
    !     CHARACTER(LEN=150),ALLOCATABLE:: ComponentName(:)

  RETURN
  END SUBROUTINE
  
  SUBROUTINE BuildResFlowModel(GeomechanicalModel,NumFluidComponents,ResFlowModel)
    
    USE DiscretizedModel
    USE DynamicSolver
    USE NumericalUtilities
    USE FEMlib
  
    IMPLICIT NONE
    TYPE(Model),INTENT(IN):: GeomechanicalModel
    INTEGER,INTENT(IN):: NumFluidComponents
    TYPE(PoreFlowModel),INTENT(IN OUT):: ResFlowModel

    INTEGER:: elementNum,numElements,numNodes,numInterfaces,i
    INTEGER:: ndim = 3                      !number of dimensions: 3D -> 3.
    INTEGER:: nod = 8                       !number of nodes per element.
    INTEGER:: nip = 8                       !number of integration points.
    CHARACTER(LEN=320)::element             !'hexahedron' ou 'tetrahedron'.
    INTEGER,ALLOCATABLE:: num(:)            !8 nodes per hex element
    REAL(8):: det
    REAL(8),ALLOCATABLE:: g_coord(:,:),&    !nodal coordenates.
                          coord(:,:)        !element nodal coordenates.
    REAL(8),ALLOCATABLE:: points(:,:),&     !coordenadas dos pontos de integração.
                          weights(:)        !peso dos pontos de integração.
    REAL(8),ALLOCATABLE:: jac(:,:),der(:,:),deriv(:,:)

    REAL(8):: elemVOL

    nod = GeomechanicalModel%NumNodesPerElement
    ndim = GeomechanicalModel%NumDimensions

    numNodes = GeomechanicalModel%NumNodes
    numElements = GeomechanicalModel%NumElements
    numInterfaces = 3*numElements

    ! Dimensioning
    ALLOCATE( ResFlowModel%FlowElementPoreVolume(numElements),&  !current pore volume
              ResFlowModel%FlowElementTransmissibilities(3,numElements),&  !I,J,K transmissibilities
              ResFlowModel%FlowElementKrel(NumFluidComponents,numElements),&  !phases relative perms
              ResFlowModel%FlowElementIndex(numElements),&  !index for each element in the grid
              ResFlowModel%FlowInterfaceNodes(4,numInterfaces),&  !four nodes interface
              ResFlowModel%FlowInterfaceElements(2,numInterfaces),&  !adjacent flow elements
              ResFlowModel%FlowInterfaceArea(numInterfaces),&  !area normal to flow
              ResFlowModel%FlowInterfaceNormal(3,numInterfaces)&  !normal vector to area
              )
    ResFlowModel%FlowElementPoreVolume=0.0
    ResFlowModel%FlowElementTransmissibilities=0.0
    ResFlowModel%FlowElementKrel=0.0
    ResFlowModel%FlowElementIndex=0
    ResFlowModel%FlowInterfaceNodes=0
    ResFlowModel%FlowInterfaceElements=0
    ResFlowModel%FlowInterfaceArea=0.0
    ResFlowModel%FlowInterfaceNormal=0.0

    ! Mesh processing preparation
    ALLOCATE(g_coord(3,numNodes))
    g_coord = 0.0
    g_coord(1,:) = GeomechanicalModel%X
    g_coord(2,:) = GeomechanicalModel%Y
    g_coord(3,:) = GeomechanicalModel%Z
    ALLOCATE(points(nip,ndim),coord(nod,ndim),  &
              jac(ndim,ndim),weights(nip),num(nod),der(ndim,nod) )
    element = GeomechanicalModel%element
    CALL sample(element,points,weights)

    ! Looping on HEXAHEDRON MESH
    ! Interfaces transmissibility
    ! Conection index
    DO elementNum=1,numElements

      IF(GeomechanicalModel%ActiveElements(elementNum))THEN
        
        num = GeomechanicalModel%ElementNodes(:,elementNum)

        !ELEMENT VOLUME
        elemVOL = 0.0
        DO i=1,nip
          CALL shape_der(der,points,i)
          jac = MATMUL(der,coord)
          CALL determinante(jac,det)
          elemVOL = elemVOL + det * weights(i)
        END DO
        ResFlowModel%FlowElementPoreVolume(elementNum) = elemVOL * GeomechanicalModel%POR(elementNum)

        !FACE I+
        

        !FACE J+

        !FACE K+



      END IF

    END DO


  RETURN
  END SUBROUTINE

  SUBROUTINE FreeMechanicalBC(GeomechanicalModel,Problem,StructuralBC)
  
    USE DiscretizedModel
    USE ProblemDefinition
    USE InterfaceMechanics
    
    IMPLICIT NONE

    TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
    TYPE(InterfaceStiffnessVector),INTENT(IN OUT):: StructuralBC
    TYPE(PD),INTENT(IN):: Problem
    
    INTEGER:: NX,NY,NZ,ix,iy,iz,i,j,node,dofs(3)
    
    NX = GeomechanicalModel%GridDivision(1)
    NZ = GeomechanicalModel%GridDivision(2)
    NY = GeomechanicalModel%GridDivision(3)
    
    node=0
    DO iy=1,NY+1
      DO iz=1,NZ+1
        DO ix=1,NX+1
          node=node+1
          
          IF( (GeomechanicalModel%Z(node) < Problem%WellData%Z_phase1) ) THEN
            
            dofs = GeomechanicalModel%NodeDOFs(:,node)
            
            DO i=1,StructuralBC%ListOfRestrictedDOFs(0)
              IF( (dofs(1) == StructuralBC%ListOfRestrictedDOFs(i)).AND.(iy<NY+1) ) StructuralBC%ListOfRestrictedDOFs(i) = 0
              IF( (dofs(2) == StructuralBC%ListOfRestrictedDOFs(i)).AND.(iy>1) )    StructuralBC%ListOfRestrictedDOFs(i) = 0
              IF(  dofs(3) == StructuralBC%ListOfRestrictedDOFs(i))                 StructuralBC%ListOfRestrictedDOFs(i) = 0
            END DO
            
          END IF
        END DO
      END DO
    END DO
    
    !Updating restricted dofs.
    j = 0
    DO i=1,StructuralBC%ListOfRestrictedDOFs(0)
      IF( StructuralBC%ListOfRestrictedDOFs(i).NE.0 )THEN
        j=j+1
        StructuralBC%ListOfRestrictedDOFs(j) = StructuralBC%ListOfRestrictedDOFs(i)
        IF(j.NE.i)StructuralBC%ListOfRestrictedDOFs(i) = 0
      END IF
    END DO
    StructuralBC%ListOfRestrictedDOFs(0) = j
  
  RETURN
  END SUBROUTINE
  
  SUBROUTINE WellWallLoad(GeomechanicalModel,Problem)
  
    USE DiscretizedModel
    USE ProblemDefinition
    
    IMPLICIT NONE

    TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
    TYPE(PD),INTENT(IN):: Problem
    
    INTEGER:: NX,NY,NZ,ix,iy,iz,i,j,node,dofs(3)
    
    NX = GeomechanicalModel%GridDivision(1)
    NZ = GeomechanicalModel%GridDivision(2)
    NY = GeomechanicalModel%GridDivision(3)
    
    node=0
    DO iy=1,NY+1
      DO iz=1,NZ+1
        DO ix=1,NX+1
          node=node+1
          
          IF( (GeomechanicalModel%Z(node) < Problem%WellData%Z_phase1) ) THEN
            
            dofs = GeomechanicalModel%NodeDOFs(:,node)
            
            !Well load.
            
            
          END IF
        END DO
      END DO
    END DO
  
  RETURN
  END SUBROUTINE

  FUNCTION FlowCriticalTimeStep(FlowModel,HFmodel)
          
      USE FluidFlowSolver
      USE HydraulicFractureModel
      
      IMPLICIT NONE
      REAL(8)::FlowCriticalTimeStep
      TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
      TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel

      REAL(8):: dt,dt_min,cw,ct
      REAL(8),ALLOCATABLE:: sumT(:),sumVC(:)
      INTEGER:: i,iel,n,dof1,dof2

      n = SIZE(FlowModel%NodeDepth)
      ALLOCATE( sumT(n), sumVC(n) )
      sumT = 0.
      sumVC = 0.

      DO i=1,FlowModel%ActiveToFlowConections(0)
        iel = FlowModel%ActiveToFlowConections(i) 
        dof1 = FlowModel%TransPointer(1,iel)
        dof2 = FlowModel%TransPointer(2,iel)
        sumT(dof1) = sumT(dof1) + ABS(FlowModel%TransArray(iel))
        sumT(dof2) = sumT(dof2) + ABS(FlowModel%TransArray(iel))
      END DO

      !n = SIZE(HFmodel%Trans)
      !DO i=1,n
      !  dof1 = HFmodel%TransDOFs(1,i)
      !  dof2 = HFmodel%TransDOFs(2,i)
      !  sumT(dof1) = sumT(dof1) + ABS(HFmodel%Trans(i))
      !  sumT(dof2) = sumT(dof2) + ABS(HFmodel%Trans(i))
      !END DO

      n = SIZE(FlowModel%NodeDepth)
      cw = FlowModel%Cw
      DO i=1,n
        ct = cw + FlowModel%NodalPoreCompressbility(i)
        sumVC(i) = sumVC(i) + FlowModel%NodeVolume(i) * ct
      END DO

      dt_min = 1.0E30
      DO i=1,n
        IF(sumVC(i)>0 .AND. sumT(i)>0)THEN
          dt = sumVC(i) / sumT(i)
          IF( dt<dt_min ) dt_min = dt
        END IF
      END DO
      FlowCriticalTimeStep = dt_min

    RETURN
  END FUNCTION
  
  SUBROUTINE UnlockModelTop(GeomechanicalModel,BC)
    
    USE DiscretizedModel
    USE ProblemDefinition
    USE InterfaceMechanics
    USE IO
    
    IMPLICIT NONE
    TYPE(Model),INTENT(IN):: GeomechanicalModel
    TYPE(InterfaceStiffnessVector):: BC
    
    INTEGER:: i,j,iel,k
    LOGICAL,ALLOCATABLE:: freeDofs(:)
    INTEGER:: num(8),g(4)
    
    ALLOCATE( freeDofs(GeomechanicalModel%NumEquations) )
    freeDofs = .FALSE.
    
    DO iel=1,GeomechanicalModel%NumElements
        k = GeomechanicalModel%ElementIJK(2,iel)
        IF( k==1 )THEN
            num = GeomechanicalModel%ElementNodes(:,iel)
            g(1) = GeomechanicalModel%NodeDOFs(3,num(2))
            g(2) = GeomechanicalModel%NodeDOFs(3,num(3))
            g(3) = GeomechanicalModel%NodeDOFs(3,num(6))
            g(4) = GeomechanicalModel%NodeDOFs(3,num(7))
            freeDofs(g) = .TRUE.
        END IF
    END DO
    
    DO i=1,BC%ListOfRestrictedDOFs(0)
        g(1) = BC%ListOfRestrictedDOFs(i)
        IF( freeDofs(g(1)) ) BC%ListOfRestrictedDOFs(i) = 0
    END DO
    
    !ORDENATE ListOfSupportedDOFs
    j = 0
    DO i=1,BC%ListOfRestrictedDOFs(0)
      IF( BC%ListOfRestrictedDOFs(i).NE.0 )THEN
        j=j+1
        BC%ListOfRestrictedDOFs(j) = BC%ListOfRestrictedDOFs(i)
        IF(j.NE.i) BC%ListOfRestrictedDOFs(i) = 0
      END IF
    END DO
    BC%ListOfRestrictedDOFs(0) = j
    
    RETURN
  END SUBROUTINE

  SUBROUTINE ReportFractureTips(GeomechanicalModel,HFmodel,fractureTIPsReports)
  
    USE DiscretizedModel
    USE IO
    USE NumericalLib
    USE HydraulicFractureModel
    
    IMPLICIT NONE

    TYPE(Model),INTENT(IN):: GeomechanicalModel
    TYPE(HydraulicFracture),INTENT(IN):: HFmodel
    REAL(8),INTENT(IN OUT):: fractureTIPsReports(:)
    
    REAL(8):: Xmax,Zmax,Zmin
    REAL(8):: z1,z2,x
    INTEGER:: iel,i
    
    REAL(8):: vol,area,disp(4)
    INTEGER:: num(8), g(4)
    
    CHARACTER(LEN=320):: text
    
    Xmax = 0
    Zmax = -1.E+30
    Zmin = +1.E+30

    fractureTIPsReports = 0.0
    
    !FRACTURE TIPS AND VOLUMES.
    fractureTIPsReports(4) = 0.0
    DO iel=1,GeomechanicalModel%NumElements
        
        i = HFmodel%SolidToFracElement(iel)
        IF( GeomechanicalModel%IsFractured(iel) .OR. GeomechanicalModel%FracAperture(iel)>0.0 )THEN
            !IF(HFmodel%Aperture(i)>0.0)THEN
                
                !fracture volume
                num = GeomechanicalModel%ElementNodes(:,iel) 
                g(1) = GeomechanicalModel%NodeDOFs(2,num(1))
                g(2) = GeomechanicalModel%NodeDOFs(2,num(2))
                g(3) = GeomechanicalModel%NodeDOFs(2,num(3))
                g(4) = GeomechanicalModel%NodeDOFs(2,num(4))
                area = AreaQuad( GeomechanicalModel, num(1:4) )
                disp = GeomechanicalModel%Displacement(g)
                vol = area * (SUM(disp)/4.0)
                fractureTIPsReports(4) = fractureTIPsReports(4) + vol
                
                !fracture tip
                z1 = MAXVAL( GeomechanicalModel%Z(  GeomechanicalModel%ElementNodes(:,iel) ) )
                z2 = MINVAL( GeomechanicalModel%Z(  GeomechanicalModel%ElementNodes(:,iel) ) )
                x = MAXVAL( GeomechanicalModel%X(  GeomechanicalModel%ElementNodes(:,iel) ) )
                IF( z1 > Zmax .AND. vol>0.0 ) Zmax = z1
                IF( z2 < Zmin .AND. vol>0.0 ) Zmin = z2
                IF( x  > Xmax .AND. vol>0.0 ) Xmax = x
            
            !END IF
        END IF
    END DO
    
    !IF( Zmax == -1.E+30 ) Zmax = 0.0
    !IF( Zmin == +1.E+30 ) Zmin = 0.0
    
    CALL PrintMessage(.FALSE.,.FALSE.,'FRACTURE')
    WRITE(text,'(A,F12.2,A)')"...                             Length:  ", Xmax," m"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    WRITE(text,'(A,F12.2,A)')"...                             Z top :  ", Zmax," m"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    WRITE(text,'(A,F12.2,A)')"...                             Z down:  ", Zmin," m"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    WRITE(text,'(A,E12.2,A)')"...                             Volume:  ", fractureTIPsReports(4)," m3"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
  
    fractureTIPsReports(1) = Zmax
    fractureTIPsReports(2) = Zmin
    fractureTIPsReports(3) = Xmax
    
  RETURN
  END SUBROUTINE
  
  SUBROUTINE ReportFractureTipsAllSteps(fractureTIPsReports)
  
    USE DiscretizedModel
    USE IO
    
    IMPLICIT NONE

    REAL(8),INTENT(IN OUT):: fractureTIPsReports(:,:)
    
    INTEGER:: i,nSteps
    REAL(8):: a1,a2
    CHARACTER(LEN=320):: text
    
    nSteps = SIZE( fractureTIPsReports(1,:) )

    CALL PrintMessage(.TRUE.,.TRUE.,'FRACTURE TIPS TRACKING')
    CALL PrintMessage(.FALSE.,.FALSE.,'  STEP   ZTOP(m)    ZDOWN(m)   LENGTH(m)   Frac.VOLUME(m3)')
    
    DO i=1,nSteps
        IF( fractureTIPsReports(1,i)>-1.0e30 ) a1 = 0.0
        IF( fractureTIPsReports(2,i)<1.0E30 ) a2 = 0.0
        WRITE(text,'(I4,3F12.2,E17.4)')i,a1,a2,&
                                       fractureTIPsReports(3,i),fractureTIPsReports(4,i)
        CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    END DO
    
  RETURN
  END SUBROUTINE
            
  SUBROUTINE SimulationReport(MaxDT,maxDP,MaxDisp,StepsRate)
  
    USE IO
    
    IMPLICIT NONE
    REAL(8),INTENT(IN):: MaxDT,maxDP,MaxDisp,StepsRate

    CHARACTER(LEN=320):: text
    
    !Solution
    !   Porepressure0(:),Porepressure(:)
    !   Temperature0(:),Temperature(:)
    !   Sw0(:),Sw(:)
        
    CALL PrintMessage(.FALSE.,.FALSE.,'SIMULATION REPORT')

    WRITE(text,'(A,F12.2,A)')"...                        Injection Rate:  ", StepsRate*60*60*24," m3/d"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    WRITE(text,'(A,F12.2,A)')"...   Max. Pressure Variation (increment):  ", maxDP/1000.0," kPa"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    WRITE(text,'(A,F12.2,A)')"...Max. Temperature Variation (increment):  ", maxDT," deg."
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    WRITE(text,'(A,F12.2,A)')"...          Max. Displacement (absolute):  ", MaxDisp," meter"
    CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
    
  END SUBROUTINE
  
  FUNCTION AreaQuad( GeomechanicalModel, num )
    
    USE DiscretizedModel
    USE NumericalLib
    
    IMPLICIT NONE
    REAL(8):: AreaQuad
    TYPE(Model),INTENT(IN):: GeomechanicalModel
    INTEGER,INTENT(IN):: num(4)
    REAL(8):: coord(3,3)          !3 node coordenates (node,coordenate)
    REAL(8):: area1,area2,&       !triagle area
              normal(3)           !normal vect
    
    INTEGER:: nodes1(3),nodes2(3)
    
    nodes1(1) = num(1)
    nodes1(2) = num(2)
    nodes1(3) = num(3)
    coord(:,1) = GeomechanicalModel%X(nodes1)
    coord(:,2) = GeomechanicalModel%Y(nodes1)
    coord(:,3) = GeomechanicalModel%Z(nodes1)
    CALL Triangle( coord, area1, normal )
    
    nodes2(1) = num(1)
    nodes2(2) = num(3)
    nodes2(3) = num(4)
    coord(:,1) = GeomechanicalModel%X(nodes2)
    coord(:,2) = GeomechanicalModel%Y(nodes2)
    coord(:,3) = GeomechanicalModel%Z(nodes2)
    CALL Triangle( coord, area2, normal )
    
    AreaQuad = area1 + area2
    
  END FUNCTION

END MODULE