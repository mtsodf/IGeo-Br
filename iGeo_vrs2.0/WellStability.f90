  
  !*****************************************************************************
  !Geotechnical Modelling
  !Eng. Luis Carlos de Sousa Junior
  !Geotechnical Engineering
  !Civil and Environmental Engineering
  !University of Brasilia - BRA
  !*****************************************************************************
   
  MODULE WellStability
  
  CONTAINS
  
  !THESES MAIN PROGRAM
  SUBROUTINE WellStabilityModule(Path,BaseFile)
  
    USE StructuralModel         !fe structural model
    USE IO
    USE ProblemDefinition       !saves the entry data
    USE WellModelBuilder
    USE DiscretizedModel        !discrete model
    USE InterfaceMechanics
    !USE DynamicExplicitSolver
    USE NearWellboreModel
    USE ConstitutiveGeoModels
    USE WellModelBuilder
    USE WF_Excavation
    USE NumericalUtilities
    USE StiffnessMatrix
    USE ViscoplasticExplicitySolver
    
    USE VTK
    USE VTKgeneral
    USE Report
    USE OMP_LIB
    USE IFPORT

    !***************************************************************************
    !VARIABLES
    !***************************************************************************
    IMPLICIT NONE
    
    CHARACTER(LEN=320), INTENT(IN)::  BaseFile, &   !file with basic information about the case.
                                      Path          !work directory where the case file is located
                                                    !and where temporary files will be created.

    TYPE (PD):: Problem
    TYPE (Model):: GeomechanicalModel
    TYPE (StructuralSolver):: StrucSolver
    TYPE (InterfaceStiffnessVector):: BC
    TYPE (ViscoSolverParameters):: SolverParameters
    TYPE (SolverMatricesContainer):: Matrices

    TYPE (NWB_Model):: WellModel

    INTEGER:: Chanel_log=9, Chanel_base=10, Chanel_vtk = 30, FreeChanel = 99
    CHARACTER(LEN=320)::VTKfileName

    INTEGER:: step,i

    LOGICAL:: found
    CHARACTER(LEN=25):: ModelType = 'WellModel1'
    
    REAL(8):: dt
    
    Problem%Path = TRIM(Path)
    Problem%WellData%WellName = TRIM(BaseFile)
    
    CALL SearchInFile( Chanel_base, '#NTHREADS', found )
    IF(found) THEN
        READ( Chanel_base, * ) GeomechanicalModel%NumPartitions
    ELSE
        GeomechanicalModel%NumPartitions = 2
    END IF
    CALL OMP_SET_NUM_THREADS( GeomechanicalModel%NumPartitions )
    
    
    !***************************************************************************
    !BUILD MODELS
    !***************************************************************************
    !General wellbore model.

    !Build well data model.
    CALL PrintMessage( .FALSE., .FALSE., 'BUILDING WELL MODEL.' )
    CALL BuildWellModel(Chanel_base,Problem, GeomechanicalModel)
    
    !Build well models.
    CALL PrintMessage( .FALSE., .FALSE., 'BUILDING WELL MODEL.' )
    !CALL BuildNearWellBoreModel( Problem%WellData, Chanel_base, WellModel )
    CALL BuildNearWellBoreModelFull( Problem%WellData, Chanel_base, WellModel )
    
    !Geomechanical Model.
    CALL PrintMessage( .FALSE., .FALSE., 'BUILDING FE GEOMECHANICAL MODEL.' )
    CALL BuildWellGeomechanicalModel( GeomechanicalModel, Path, BaseFile,&
                                      WellModel,Problem%WellData, Problem )
    CALL EditProperties( GeomechanicalModel, Problem )
  
    !Structural Boundary Conditions.
    BC%NumEquations = GeomechanicalModel%NumEquations
    CALL SupportBoundaryConditions(WellModel%RestrictedDOFs,&
                                    WellModel%SupportedNodeDOFs,&
                                    GeomechanicalModel,&
                                    BC%ListOfRestrictedDOFs,BC%ListOfSupportedDOFs,&
                                    BC%Kvalues)
    
    !PREPARING STRUCTURAL SOLVER------------------------------------------------
    CALL ReadSolverParameters( Path, BaseFile, Chanel_base, StrucSolver )
    GeomechanicalModel%SolverType = TRIM(StrucSolver%StructuralSolver)
    IF(StrucSolver%StructuralSolver=='Viscoplastic')THEN
      CALL PrepareViscoplasticSolver(Chanel_base,GeomechanicalModel,GeomechanicalModel%ActiveElements,&
                                            SolverParameters,Matrices%Stiff)
    ELSEIF(StrucSolver%StructuralSolver=='DynamicFEM')THEN
      !CALL PrintMessage( .FALSE., .FALSE., 'Structutal explicity solver.' )
      !CALL PrepareDynamicExplicitSolver(GeomechanicalModel,Matrices,SolverParameters,&
      !                                  BC%ListOfRestrictedDOFs,Chanel_base)
    ELSE
    END IF
      
    !BUILD GENERAL FEM OUTPUT FILE----------------------------------------------
    OPEN(100,FILE=TRIM(TRIM(Path)//'/'//TRIM(BaseFile)//'.femOutput'),FORM = 'UNFORMATTED')
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
    step = 0
    
    VTKfileName = TRIM( TRIM(Path)//'/'//TRIM(Problem%ModelName)//'.vtk' )
    OPEN(Chanel_vtk,FILE = TRIM(VTKfileName) )
    CALL HexLin_MeshToVTKfile( Problem%Path, Problem%ModelName, Chanel_vtk,&
                               WellModel%Coords, WellModel%ElementIncidence)
    CALL WriteStaticDataTo2VTK(Chanel_vtk,Problem,Problem%WellData,WellModel,GeomechanicalModel)
    CALL WriteStepResults2VTK(Chanel_vtk,VTKfileName,Problem,GeomechanicalModel,&
                                  TRIM(Problem%ModelName),step)
        
    !OPEN(Chanel_vtk,FILE = TRIM(VTKfileName) )
      
    !***************************************************************************
    !WORKFLOWS.
    !***************************************************************************

    !EXCATAVION-----------------------------------------------------------------
    DO i=1,SIZE(Problem%TimeSteps)
        step = step + 1
          
        READ(Problem%TimeSteps(i),*)dt
        dt = 60 * dt
    
        CALL SearchInFile(Chanel_base,'#WF_Excavation',found)
        IF(found)THEN
          CALL SearchInFile( Chanel_base, '#MudHistoryPPG', found )
          IF(found) THEN
            READ(Chanel_base,*) Problem%WellData%MudWeight                              !ppg
            Problem%WellData%MudWeight = Problem%WellData%MudWeight * (9810.0/8.34)     !Pa/m
          ELSE
            Problem%WellData%MudWeight = 1.028 * 9810.0                                 !Pa/m
          END IF
          IF( Problem%Structural_flag .AND. Problem%WellData%MeshType=='RAD' .AND. &
                                                  Problem%Structural_flag ) &
            CALL Excavation( Problem%WellData,WellModel,GeomechanicalModel,Matrices,&
                                        BC,SolverParameters,Problem,.FALSE.,ModelType,step,dt )
            CALL WriteStepResults2VTK(Chanel_vtk,VTKfileName,Problem,GeomechanicalModel,&
                                      TRIM(Problem%ModelName),step)
            CALL SectionReport(Chanel_base,GeomechanicalModel)
        END IF
        
        !Transient effects: temperature and porepressure----------------------------
        IF(ANY(GeomechanicalModel%KH>0.))THEN
        !  CALL RunFracWaterInjection( step,Problem%WellData,GeomechanicalModel,Matrices,&
        !                              BC,SolverParameters,&
        !                              Problem,BaseFile,Path,VTKfileName,Chanel_base)
        END IF
        
    END DO

    RETURN
  END SUBROUTINE
  
  SUBROUTINE SectionReport(Chanel_base,GeomechanicalModel)
    
    USE DiscretizedModel
    USE SolidMechanicsLib
    USE NumericalLib
    USE IO
    
    IMPLICIT NONE
    TYPE (Model),INTENT(IN):: GeomechanicalModel
    INTEGER,INTENT(IN):: Chanel_base
    
    REAL(8):: Srad,St,Sv
    REAL(8):: radDisplacement, verticalDisplacement,ux,uy
    REAL(8):: coordsTRI(3,3),area,normal(3),stress_rot(3)
    INTEGER:: KK,i,j,k,elementNumber,num(8)
    LOGICAL:: found
    CHARACTER(LEN=320):: text
    
    CALL SearchInFile(Chanel_base,'#SECTIONREPORT',found)
    IF(found)THEN
        
        READ(Chanel_base,*)KK
        
        CALL PrintMessage( .TRUE., .TRUE., 'SECTION DATA.' )
                
        CALL PrintMessage( .FALSE., .FALSE., 'Cilindrical Effective Stress (kgf/cm2) and Displacement (mm)' )
        CALL PrintMessage( .FALSE., .FALSE., '   I   J   K      Srad      Stan     Svert       Urad      Uvert     ' )
        
        elementNumber = 0
        DO j=1, GeomechanicalModel%GridDivision(3)
            DO k=1, GeomechanicalModel%GridDivision(2)
                DO i=1, GeomechanicalModel%GridDivision(1)
                    elementNumber = elementNumber + 1
                    IF(k==KK.AND.i==1)THEN
                        num = GeomechanicalModel%ElementNodes(:,elementNumber)
                        
                        !RADIAL DISPLACEMENT
                        ux = 1.0/8.0 * SUM( GeomechanicalModel%Displacement( GeomechanicalModel%NodeDofs(1,num) ))
                        uy = 1.0/8.0 * SUM( GeomechanicalModel%Displacement( GeomechanicalModel%NodeDofs(2,num) ))
                        radDisplacement = SQRT( ux**2.0 + uy**2.0 )
                        verticalDisplacement = 1.0/8.0 * SUM(GeomechanicalModel%Displacement( GeomechanicalModel%NodeDofs(3,num) ))
                        
                        !RADIAL STRESS
                        coordsTRI(1,1)=GeomechanicalModel%X(num(5))
                        coordsTRI(1,2)=GeomechanicalModel%Y(num(5))
                        coordsTRI(1,3)=GeomechanicalModel%Z(num(5))
                        coordsTRI(2,1)=GeomechanicalModel%X(num(1))
                        coordsTRI(2,2)=GeomechanicalModel%Y(num(1))
                        coordsTRI(2,3)=GeomechanicalModel%Z(num(1))
                        coordsTRI(3,1)=GeomechanicalModel%X(num(2))
                        coordsTRI(3,2)=GeomechanicalModel%Y(num(2))
                        coordsTRI(3,3)=GeomechanicalModel%Z(num(2))
                        CALL Triangle( coordsTRI, area, normal )
                        CALL rot_tensao( GeomechanicalModel%EffectiveStressTensor(:,1,elementNumber),&
                                            coordsTRI(1,:), coordsTRI(2,:), coordsTRI(3,:),stress_rot )    !(Sn,Sxy,Szx)
                        Srad = stress_rot(1)
                        
                        !TANGENTIAL STRESS
                        coordsTRI(1,1)=0.5*( GeomechanicalModel%X(num(1)) + GeomechanicalModel%X(num(5)) )
                        coordsTRI(1,2)=0.5*( GeomechanicalModel%Y(num(1)) + GeomechanicalModel%Y(num(5)) )
                        coordsTRI(1,3)=0.5*( GeomechanicalModel%Z(num(1)) + GeomechanicalModel%Z(num(5)) )
                        coordsTRI(2,1)=0.5*( GeomechanicalModel%X(num(3)) + GeomechanicalModel%X(num(7)) )
                        coordsTRI(2,2)=0.5*( GeomechanicalModel%Y(num(3)) + GeomechanicalModel%Y(num(7)) )
                        coordsTRI(2,3)=0.5*( GeomechanicalModel%Z(num(3)) + GeomechanicalModel%Z(num(7)) )
                        coordsTRI(3,1)=0.5*( GeomechanicalModel%X(num(2)) + GeomechanicalModel%X(num(6)) )
                        coordsTRI(3,2)=0.5*( GeomechanicalModel%Y(num(2)) + GeomechanicalModel%Y(num(6)) )
                        coordsTRI(3,3)=0.5*( GeomechanicalModel%Z(num(2)) + GeomechanicalModel%Z(num(6)) )
                        CALL Triangle( coordsTRI, area, normal )
                        CALL rot_tensao( GeomechanicalModel%EffectiveStressTensor(:,1,elementNumber),&
                                            coordsTRI(1,:), coordsTRI(2,:), coordsTRI(3,:),stress_rot )    !(Sn,Sxy,Szx)
                        St = stress_rot(1)
                        
                        !VERTICAL STRESS
                        coordsTRI(1,1)=GeomechanicalModel%X(num(1))
                        coordsTRI(1,2)=GeomechanicalModel%Y(num(1))
                        coordsTRI(1,3)=GeomechanicalModel%Z(num(1))
                        coordsTRI(2,1)=GeomechanicalModel%X(num(4))
                        coordsTRI(2,2)=GeomechanicalModel%Y(num(4))
                        coordsTRI(2,3)=GeomechanicalModel%Z(num(4))
                        coordsTRI(3,1)=GeomechanicalModel%X(num(5))
                        coordsTRI(3,2)=GeomechanicalModel%Y(num(5))
                        coordsTRI(3,3)=GeomechanicalModel%Z(num(5))
                        CALL Triangle( coordsTRI, area, normal )
                        CALL rot_tensao( GeomechanicalModel%EffectiveStressTensor(:,1,elementNumber),&
                                            coordsTRI(1,:), coordsTRI(2,:), coordsTRI(3,:),stress_rot )    !(Sn,Sxy,Szx)
                        Sv = stress_rot(1)
                        
                        WRITE(text,'(I4,I4,I4,3(F10.2),2(E11.3))')i,j,k, Srad*1.01972E-5, St*1.01972E-5, Sv*1.01972E-5, radDisplacement*1000, verticalDisplacement*1000
                        CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
                        
                    END IF
                END DO
            END DO
        END DO

    END IF
    
  RETURN
  END SUBROUTINE
  
  END MODULE