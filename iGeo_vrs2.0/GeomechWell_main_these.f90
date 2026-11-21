  
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Geotechnical Engineering
  !University of Brasilia
  !*****************************************************************************
  
  MODULE GeomechWell_main_these
  
  CONTAINS
  
  !THESES MAIN PROGRAM
  SUBROUTINE GeomecW(Path,BaseFile)
  
    USE StructuralModel         !modelo estrutural em Elementos Finitos.
    USE IO                      !métodos auxiliares para IO.
    USE ProblemDefinition       !objeto que guarda dos dados de entrada do problema.
    USE DiscretizedModel        !modelo geomecanico discreto.
    USE InterfaceMechanics
    USE DynamicSolver
    USE FluidFlowSolver
      
    USE WF_Excavation
    !USE WF_UniaxialModeling
    !USE WF_TriaxialModeling
    !USE WF_DisplacementLoad
      
    !USE WF_Construction
    !USE WF_LateralLoad
    !USE WF_BeamTractionTest
    !USE WF_TherzaguiColumn
      
    USE VTK
    USE VTKgeneral
    USE Report
    USE OMP_LIB                 !bibliotecas OpenMP(R)
    USE IFPORT
      
    USE NearWellboreModel       !objeto para um modelo "near wellbore".
    USE HydraulicFractureModel
    USE WF_WaterInjection
    !USE WF_WaterInjectionRAD
    !USE WF_HydraulicFractureConstPressure

    USE ConstitutiveGeoModels
    USE ViscoplasticExplicitySolver
    USE StiffnessMatrix
    USE WF_MEDSolver

    !***************************************************************************
    !VARIABLES
    !***************************************************************************
    IMPLICIT NONE
    
    CHARACTER(LEN=320), INTENT(IN)::  BaseFile, &   !file with basic information about the case.
                                      Path          !work directory where the case file is located
                                                    !and where temporary files will be created.
                                                    
    TYPE (PD):: Problem
    TYPE (Model):: GeomechanicalModel
    TYPE (PorousFlowModel):: RESflowModel
    TYPE (StructuralSolver):: StrucSolver
    TYPE (HydraulicFracture):: HFmodel
    TYPE (InterfaceStiffnessVector):: BC
    TYPE (ViscoSolverParameters):: SolverParameters
    TYPE (DynamicSolverParameters):: DynamicParameters
    TYPE (SolverMatricesContainer):: Matrices
      
    !Well.
    TYPE (WellboreData):: WellData
    TYPE (NWB_Model):: WellModel
    
    !Variables.
  
    INTEGER:: Chanel_log=9, Chanel_base=10, Chanel_vtk = 30, FreeChanel = 99

    CHARACTER(LEN=320)::buffer  !variável para armazenamento de argumentos passados ao programa.
    CHARACTER(LEN=320)::VTKfileName
      
    CHARACTER(LEN=25):: ModelType       !Well
                                        !Soil

    INTEGER:: nthreads,&
              numARG,&
              step

    LOGICAL:: found,&
              test=.TRUE.

    REAL(8):: temp
    
    !***************************************************************************
    !BUILD MODELS
    !***************************************************************************
      
    CALL ReadModelType(Chanel_base,ModelType)
      
    !READING BASE FILE----------------------------------------------------------

    !Read base file and building the Problem.
    CALL PrintMessage( .FALSE., .FALSE., 'READING BASE FILE.' )
    CALL ReadBaseFile( Path, BaseFile, Chanel_base, Problem, WellData, nthreads )

    !Build well models.
    CALL PrintMessage( .FALSE., .FALSE., 'BUILDING WELL MODEL.' )
    CALL BuildNearWellBoreModel( WellData, Chanel_base, WellModel )
    !BuildNearWellBoreModel completa WellData (ex.: WellBoreWallFlag); mantém a
    !cópia em Problem%WellData atualizada, como no fim de ReadBaseFile.
    Problem%WellData = WellData

    !Geomechanical Model.
    CALL PrintMessage( .FALSE., .FALSE., 'BUILDING FE GEOMECHANICAL MODEL.' )
    CALL BuildWellGeomechanicalModel( GeomechanicalModel, Path, BaseFile, WellModel,WellData, Problem )
    CALL EditProperties( GeomechanicalModel, Problem )
  
    !Structural Boundary Conditions.
    BC%NumEquations = GeomechanicalModel%NumEquations
    !BC%Values = 0.0
    CALL SupportBoundaryConditions(WellModel%RestrictedDOFs,WellModel%SupportedNodeDOFs,&
                                    GeomechanicalModel,&
                                    BC%ListOfRestrictedDOFs,BC%ListOfSupportedDOFs,BC%Kvalues) 
      
    !PREPARING STRUCTURAL SOLVER------------------------------------------------
    GeomechanicalModel%NumPartitions = nthreads
    CALL ReadSolverParameters( Path, BaseFile, Chanel_base, StrucSolver )
    GeomechanicalModel%SolverType = TRIM(StrucSolver%StructuralSolver)
    IF(StrucSolver%StructuralSolver=='Viscoplastic')THEN
      CALL PrepareViscoplasticSolver(Chanel_base,GeomechanicalModel,GeomechanicalModel%ActiveElements,&
                                            SolverParameters,Matrices%Stiff)
    ELSEIF(StrucSolver%StructuralSolver=='DynamicFEM')THEN
      CALL PrintMessage( .FALSE., .FALSE., 'Structutal explicity solver.' )
      CALL PrepareDynamicExplicitSolver(GeomechanicalModel,Matrices,DynamicParameters,&
                                        BC%ListOfRestrictedDOFs,Chanel_base)
    ELSE
    END IF
      
    !BUILD GENERAL FEM OUTPUT FILE----------------------------------------------
    OPEN(100,FILE=TRIM(TRIM(Path)//'/'//TRIM(BaseFile)//'.femOutput'),FORM = 'UNFORMATTED')
    WRITE(100)GeomechanicalModel%GridDivision               !NX,NZ,NY for structured meshes
    WRITE(100)GeomechanicalModel%NumNodes
    WRITE(100)GeomechanicalModel%NumElements
    WRITE(100)GeomechanicalModel%NumIntegrationPoints
    WRITE(100)GeomechanicalModel%NumEquations
    WRITE(100)GeomechanicalModel%X
    WRITE(100)GeomechanicalModel%Y
    WRITE(100)GeomechanicalModel%Z
    WRITE(100)GeomechanicalModel%ElementNodes
    WRITE(100)GeomechanicalModel%NodeDofs
    WRITE(100)GeomechanicalModel%BiotCoef
    !Transient data.
    WRITE(100)GeomechanicalModel%Porepressure
    WRITE(100)GeomechanicalModel%Displacement
    WRITE(100)GeomechanicalModel%EffectiveStressTensor
    WRITE(100)GeomechanicalModel%StrainTensor
    WRITE(100)GeomechanicalModel%PlasticStrainTensor
      
    !***************************************************************************
    !INITIALIZING POST-PROCESSING FILES.
    !***************************************************************************
  
    !Post-processing.
    CALL PrintMessage( .FALSE., .FALSE., 'INITIALIZING VTK FILE.' )

    VTKfileName = TRIM( TRIM(Path)//'/'//TRIM(WellData%WellName)//'.vtk' )
    OPEN(Chanel_vtk,FILE = TRIM(VTKfileName) )
    CALL HexLin_MeshToVTKfile( Problem%Path, WellData%WellName, Chanel_vtk,&
                                WellModel%Coords, WellModel%ElementIncidence)
    CALL WriteStaticDataTo2VTK(Chanel_vtk,Problem,WellData,WellModel,GeomechanicalModel)
      
    CLOSE(Chanel_vtk)
      
    !***************************************************************************
    !WORKFLOWS.
    !***************************************************************************
    
    step = 0
      
    !!DRILLING-------------------------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_Drilling',found)
    !IF(found)THEN
    !  
    !  CALL PrintMessage( .FALSE., .FALSE., 'DRILLING.' )
    !    
    !  !Write results prior to excavation.
    !  CALL PrintMessage( .FALSE., .FALSE., '...Report' )
    !  CALL WriteStepResults2VTK(Chanel_vtk,VTKfileName,Problem,GeomechanicalModel,&
    !                            TRIM(Problem%ModelName),step)
    !  CALL WriteDrillingReport(FreeChanel,TRIM(Path)//'/'//TRIM(WellData%WellName),&
    !                            Problem,GeomechanicalModel,step)
    !    
    !  CALL Excavation( WellData,WellModel,GeomechanicalModel,Matrices,&
    !                      BC,SolverParameters,Problem,.FALSE.,ModelType )
    !    
    !  !Report
    !  CALL PrintMessage( .FALSE., .FALSE., '...Report' )
    !  step = step + 1
    !  CALL WriteStepResults2VTK(Chanel_vtk,VTKfileName,Problem,GeomechanicalModel,&
    !                            TRIM(Problem%ModelName),step)
    !  CALL WriteDrillingReport(FreeChanel,TRIM(Path)//'/'//TRIM(WellData%WellName),&
    !                            Problem,GeomechanicalModel,step)
    !  !General FEM output file
    !  WRITE(100)GeomechanicalModel%Porepressure
    !  WRITE(100)GeomechanicalModel%Displacement
    !  WRITE(100)GeomechanicalModel%EffectiveStressTensor
    !  WRITE(100)GeomechanicalModel%StrainTensor
    !  WRITE(100)GeomechanicalModel%PlasticStrainTensor
    !    
    !  CALL SLEEP(3)
    !  STOP
    !    
    !END IF
    
    !!UNIAXIAL STRESS MODELLING--------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_UniaxialLoad',found)
    !IF(found)THEN
    !  CALL RunUniaxialLoad(step,GeomechanicalModel,Matrices,&
    !                        BC,SolverParameters,&
    !                        Problem,BaseFile,Path,Chanel_vtk,Chanel_base)
    !  STOP
    !END IF
    !
    !!TRIAXIAL STRESS MODELLING--------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_TriaxialLoad',found)
    !IF(found)THEN
    !  CALL RunTriaxialLoad(step,GeomechanicalModel,Matrices,&
    !                        BC,SolverParameters,&
    !                        Problem,BaseFile,Path,Chanel_vtk,Chanel_base,VTKfileName)
    !  STOP
    !END IF

    !TRANSIENT POROELASTIC FRACTURE PROPAGATION-------------------------------------------
    CALL SearchInFile(Chanel_base,'#WF_MEDSolver',found)
    IF(found)THEN
      step = 0
      temp = 0.0

      CALL RunPoromechanicsProblem( GeomechanicalModel,Matrices,&
                                    DynamicParameters,&
                                    Problem,BaseFile,Path,VTKfileName,Chanel_base )

      CALL SLEEP(3)
      STOP
        
    END IF
      
    !TRANSIENT POROELASTIC FRACTURE PROPAGATION-------------------------------------------
    CALL SearchInFile(Chanel_base,'#WF_WaterInj',found)
    IF(found)THEN
      step = 0
      temp = 0.0
      
      IF( Problem%Structural_flag .AND. WellData%MeshType=='RAD' .AND. Problem%Structural_flag ) &
                  CALL Excavation( WellData,WellModel,GeomechanicalModel,Matrices,&
                                    BC,SolverParameters,Problem,.FALSE.,ModelType,step,temp )
      
      CALL RunFracWaterInjection(step,WellData,GeomechanicalModel,Matrices,&
                                  BC,SolverParameters,&
                                  Problem,BaseFile,Path,VTKfileName,Chanel_base)
      CALL SLEEP(3)
      STOP
        
    END IF
      
    !!DRILLING-------------------------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_Excavation',found)
    !IF(found)THEN
    !    
    !  step = 0
    !  IF( Problem%Structural_flag .AND. WellData%MeshType=='RAD' .AND. Problem%Structural_flag ) &
    !              CALL Excavation( WellData,WellModel,GeomechanicalModel,Matrices,&
    !                                BC,SolverParameters,Problem,.FALSE.,ModelType )
    !  CALL SLEEP(3)
    !  STOP
    !    
    !END IF
    !
    !!!CONSTANT PRESSURE FRACTURE PROPAGATION (NO FILTRATION)--------------------
    !!CALL SearchInFile(Chanel_base,'#WF_HydraulicFractureConstPressure',found)
    !!IF(found)THEN
    !!  CALL RunConstPressFracturePropagation(step,GeomechanicalModel,Matrices,&
    !!                                        BC,SolverParameters,&
    !!                                        Problem,BaseFile,Path,Chanel_vtk,Chanel_base,VTKfileName)
    !!  STOP
    !!END IF   
    !  
    !!UNIAXIAL DISPLACEMENT MODELLING--------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_PrescribedDisplacement',found)
    !IF(found)THEN
    !  CALL RunUniaxialDisplacement(step,GeomechanicalModel,Matrices,&
    !                                BC,SolverParameters,&
    !                                Problem,BaseFile,Path,Chanel_vtk,Chanel_base)
    !  STOP
    !END IF
    !  
    !!LATERAL LOAD---------------------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_LateralLoad',found)
    !IF(found)THEN
    !  CALL RunLateralLoad(step,GeomechanicalModel,Matrices,&
    !                      BC,SolverParameters,&
    !                      Problem,BaseFile,Path,Chanel_vtk,Chanel_base,VTKfileName)
    !  STOP
    !END IF
    !  
    !!!BEAM TRACTION TEST--------------------------------------------------------
    !!CALL SearchInFile(Chanel_base,'#WF_BeamTractionTest',found)
    !!IF(found)THEN
    !!  CALL RunBeamTractionLoad( step,GeomechanicalModel,Matrices,&
    !!                            BC,SolverParameters,&
    !!                            Problem,BaseFile,Path,VTKfileName,Chanel_base)
    !!  STOP
    !!END IF
    !  
    !!THERZAGUI COLUMN-----------------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_TherzaguiColumn',found)
    !IF(found)THEN
    !  CALL RunTherzaguiColumn( step,WellData,GeomechanicalModel,Matrices,&
    !                            BC,SolverParameters,&
    !                            Problem,BaseFile,Path,VTKfileName,Chanel_base )
    !  STOP
    !END IF
      
    !!SOIL GEOMECHANICS---------------------------------------------------------
    !CALL SearchInFile(Chanel_base,'#WF_SoilGeomechanics',found)
    !IF(found)THEN
    !  CALL RunSoilMechanics( step,GeomechanicalModel,Matrices,&
    !                          BC,SolverParameters,&
    !                          Problem,BaseFile,Path,Chanel_vtk,Chanel_base )
    !END IF
      
    !***************************************************************************
    !FINALIZING.
    !***************************************************************************
    CALL PrintMessage( .FALSE., .FALSE., 'SIMULATION ENDED.' )
    CALL SLEEP(3)

  END SUBROUTINE
  
  
  END MODULE