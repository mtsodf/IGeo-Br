
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia
  !*****************************************************************************
  
  MODULE WF_Excavation

  CONTAINS
  
    SUBROUTINE Excavation( WellData,WellModel,GeomechanicalModel,Matrices,&
                           StructuralBC,SolverParameters,Problem,UpdateBC,&
                           ModelType,Step,TimeIncrement )
      USE omp_lib
      USE IO
      USE DiscretizedModel
      USE InterfaceMechanics
      USE ProblemDefinition
      USE NearWellboreModel
      !USE DynamicExplicitSolver
      USE PorepressureLoad
      USE GravityLoad
      USE HydraulicFractureModel
      USE FluidFlowSolver
      USE ViscoplasticExplicitySolver
      USE StiffnessMatrix
      USE NumericalUtilities
      USE CreepMod
      USE VTK

      IMPLICIT NONE

      TYPE(WellboreData),INTENT(IN):: WellData
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      TYPE (SolverMatricesContainer):: Matrices
      TYPE(InterfaceStiffnessVector),INTENT(IN OUT):: StructuralBC
      TYPE(ViscoSolverParameters):: SolverParameters
      TYPE(PD),INTENT(IN):: Problem
      LOGICAL,INTENT(IN):: UpdateBC
      CHARACTER(LEN=25):: ModelType
      INTEGER,INTENT(IN):: Step
      REAL(8),INTENT(IN):: TimeIncrement
    
      TYPE(InterfaceStiffnessVector):: StructuralBC_modif
    
      TYPE(DisplacementBC):: DispBC
      TYPE(HydraulicFracture):: HFmodel
      TYPE(StepReport):: reportData
    
      !Viscoplastic solver parameters.
      INTEGER:: i,j,k,iel
      REAL(8)::cg_tol,&                        !tolerância do solver linear de equações - razão de máximos deslocamentos
               tol                             !tolerância do solver viscoplástico - razão de máximas forças desbalanceadas
      INTEGER::limit,&                         !limite de iterações para o solver viscoplástico
               cg_limit                        !limete de iterações para o solver linear de equações
    
      REAL(8):: z

      REAL(8),ALLOCATABLE:: Load(:)
      REAL(8),ALLOCATABLE:: deltaP(:),deltaT(:)
    
      INTEGER,ALLOCATABLE:: ListOfElements(:)

      INTEGER:: nst = 6                       !number of stress tensor.
      INTEGER:: nip = 8                       !number of integration points.
      INTEGER:: numElements, numEquations, numFluidEquations
    
      REAL(8),ALLOCATABLE:: displacement0(:),displacement(:)
    
      LOGICAL:: converged
    
      INTEGER:: Chanel_log = 9
    
      LOGICAL:: cementElement
      
      !***************************************************************************
    
      CALL PrintMessage( .FALSE., .FALSE., 'WORKFLOW: EXCAVATION.' )

      reportData%time = 0.0
      reportData%globalStep = Step
      reportData%flowStep = 0
      reportData%structuralStep = 0
      reportData%fractureStep = 0
      reportData%linIters = 0
      reportData%nonLinIters = 0
      reportData%stepResult = ''
      CALL PrintStepReportHEAD()
      
      ALLOCATE( Load(0:GeomechanicalModel%NumEquations),&
                ListOfElements(GeomechanicalModel%NumElements) )
      Load = 0.0
      ListOfElements = 0
      DO i=1,GeomechanicalModel%NumElements
        ListOfElements(i) = i
      END DO
    
      numElements = GeomechanicalModel%NumElements
      numEquations = GeomechanicalModel%NumEquations
    
      numFluidEquations = GeomechanicalModel%NumNodes
      ALLOCATE( deltaP(numFluidEquations),deltaT(numFluidEquations))
      deltaP = 0.0
      deltaT = 0.0
        
      !Y displacemente restriction for negative displacement.
      ALLOCATE(DispBC%ListYnodalRestriction(0:GeomechanicalModel%NumNodes))
      DispBC%ListYnodalRestriction = 0
      
      !***************************************************************************
      ! STRUCTURAL LOAD AND B.C.                                                 *
      !***************************************************************************
    
      StructuralBC_modif = StructuralBC
      
      !Displacement Boundary Conditions
      CALL OpenWellBC(GeomechanicalModel,StructuralBC_modif%ListOfRestrictedDOFs)
      
      !Excavation equivalent nodal load
      Load = 0.
      !CALL NumericalExcavationLoad(WellData,WellModel,GeomechanicalModel,Load)
      CALL GeneralExcavationLoad( WellData,WellModel,GeomechanicalModel,Load )
      
      !***************************************************************************
      ! ELASTIC STIFFNESS                                                        *
      !***************************************************************************

      !Elastic stiffness.
      IF(.NOT.Matrices%Built) THEN
        Matrices%Stiff%a = 0.0
        IF(TRIM(GeomechanicalModel%SolverType)=='DynamicFEM')THEN
          !CALL BuildElasticStiffness(GeomechanicalModel,Matrices,GeomechanicalModel%ActiveElements)
        ELSE
          CALL BuildStiffness(GeomechanicalModel,Matrices%Stiff,GeomechanicalModel%ActiveElements)
        END IF
      END IF
    
      !***************************************************************************
      ! STRUCTURAL SOLVER                                                        *
      !***************************************************************************

      ALLOCATE( displacement0(numEquations), displacement(0:numEquations) )
      displacement0 = 0.0
      displacement = 0.0                                                                ! <- initial displacement
      
      IF(GeomechanicalModel%Creep)THEN
          CALL CreepSolver( GeomechanicalModel,Matrices,&
                            StructuralBC_modif,SolverParameters,ReportData,DeltaT,Load,TimeIncrement)
          
      ELSEIF(TRIM(GeomechanicalModel%SolverType)=='DynamicFEM')THEN
        !CALL DynamicSolver( GeomechanicalModel,Matrices,&                               ! -> u,eps,sigma
        !                    SolverParameters,StructuralBC_modif%ListOfRestrictedDOFs,&
        !                    DispBC,Load(1:),displacement,&
        !                    Problem%StructuralSolverPrintFreq,&
        !                    Problem%ElementMonitor,&
        !                    Problem%NodeMonitor,&
        !                    Chanel_log,GeomechanicalModel%ActiveElements)
        !!converged:   0 - success
        !!             1 - not converged
        !!             2 - failed
        !IF((SolverParameters%SolverConvergence==2))THEN                                 ! -> failed to converge -> break time step
        !  converged = .FALSE.
        !END IF          !=0 or =1
        !converged = .TRUE.
          
      ELSEIF(GeomechanicalModel%SolverType=='Viscoplastic')THEN
        cg_tol = 1.0E-6
        cg_limit = 1.0E6
        tol = 1.0E-2
        limit = 20
        CALL ViscoSolver(GeomechanicalModel,reportData,Matrices%Stiff%a, Matrices%Stiff%ia, &
                          Matrices%Stiff%ja,&
                          StructuralBC_modif%ListOfRestrictedDOFs,Load(1:),DeltaT,&
                          cg_tol,cg_limit,tol,limit,&
                          displacement)
        
      ELSE
      END IF
      
      CALL UpdateStrainWork(GeomechanicalModel,&
                            GeomechanicalModel%EffectiveStressTensor,&
                            displacement)
      
      !***************************************************************************
      !UPDATING BOUNDARY CONDITION AT WELL WALL - UpdateBC
      IF(UpdateBC) StructuralBC = StructuralBC_modif
      !***************************************************************************
      
      !General FEM output file.
      WRITE(100)GeomechanicalModel%Porepressure
      WRITE(100)GeomechanicalModel%Displacement
      WRITE(100)GeomechanicalModel%EffectiveStressTensor
      WRITE(100)GeomechanicalModel%StrainTensor
      WRITE(100)GeomechanicalModel%PlasticStrainTensor

    RETURN
    END SUBROUTINE

    !Numerial routine for wellbore excavation LOAD.
    !Build Tri3D FE inside the well.
    SUBROUTINE NumericalExcavationLoad(WellData,WellModel,GeomechanicalModel,Load)

        USE MainConstants
        USE FEMlib
        USE NumericalLib
        USE ProblemDefinition
        USE IO
        USE MeshAuxiliares
        USE MeshIO
        USE DiscretizedModel
        USE NearWellboreModel
        USE StressLoad

        IMPLICIT NONE
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(NWB_Model),INTENT(IN OUT):: WellModel
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        REAL(8),INTENT(IN OUT):: Load(0:)

        INTEGER:: numIntegrationPoints
        PARAMETER(numIntegrationPoints = 1 )

        INTEGER:: i,k,j,m,elementNumber

        INTEGER,ALLOCATABLE:: num(:)

        REAL(8):: Zavg

        REAL(8):: theta
        REAL(8):: force(3)
        REAL(8):: normal, shear, sxx, sxy, syy, stress_rot(3)

        REAL(8):: stress(6,numIntegrationPoints)   !pressao hidrostatica de fluido.
        REAL(8),ALLOCATABLE:: HydrostaticStress(:)
        
        REAL(8):: z
        
        REAL(8),ALLOCATABLE:: points(:,:),&  !coordenadas dos pontos de integração
                              weights(:),&  !peso dos pontos de integração
                              coord(:,:),&  !coordenada dos nós de um elemento
                              elemNodalForce(:)
        INTEGER:: globalDOFs(12),localDOFs(12)
        
        ndim = GeomechanicalModel%NumDimensions
        nst = GeomechanicalModel%NumTensorComponents
        nod = GeomechanicalModel%NumNodesPerElement
        ndof = GeomechanicalModel%NumDOFUperElement
        element = GeomechanicalModel%element
        
        ALLOCATE( num(nod), points(numIntegrationPoints,ndim),&
            weights(numIntegrationPoints),coord(nod,ndim),elemNodalForce(ndof),&
            HydrostaticStress(nst) )
        
        !INTEGRATION POINTS SETUP
        CALL sample(element,points,weights)
        
        localDOFs(1) = 7 !3
        localDOFs(2) = 8
        localDOFs(3) = 9
        localDOFs(4) = 10 !4
        localDOFs(5) = 11
        localDOFs(6) = 12
        localDOFs(7) = 19 !7
        localDOFs(8) = 20
        localDOFs(9) = 21
        localDOFs(10) = 22 !8
        localDOFs(11) = 23
        localDOFs(12) = 24

        !*************************************************************************
        !STRUCTURAL LOAD
        !*************************************************************************
    
        !Looping on elements
        elementNumber = 0
        DO j=1, WellModel%NTHETA
          DO k=1, WellModel%NZ
            DO i=1, WellModel%NRAD
              elementNumber = elementNumber + 1
              
              IF( i==1 ) THEN
                  
                  num = GeomechanicalModel%ElementNodes(:,elementNumber)
                  coord(:,1) = GeomechanicalModel%X( num ) 
                  coord(:,2) = GeomechanicalModel%Y( num )
                  coord(:,3) = GeomechanicalModel%Z( num )
                  
                  !tri-hex element inside the well (>>>VERTICAL)
                  coord(3,:) = coord(2,:)
                  coord(4,:) = coord(1,:)
                  coord(7,:) = coord(6,:)
                  coord(8,:) = coord(5,:)
                  coord(1,1:2) = 0.
                  coord(2,1:2) = 0.
                  coord(5,1:2) = 0.
                  coord(6,1:2) = 0.
                  
                  globalDOFs(1:3) = GeomechanicalModel%NodeDOFs(:,num(1))
                  globalDOFs(4:6) = GeomechanicalModel%NodeDOFs(:,num(2))
                  globalDOFs(7:9) = GeomechanicalModel%NodeDOFs(:,num(5))
                  globalDOFs(10:12) = GeomechanicalModel%NodeDOFs(:,num(6))
                                
                  !TOTAL STRESS AT FREE SURFACE--------------------------------------------------------
                  stress(:,1) = GeomechanicalModel%EffectiveStressTensor(1:6,1,elementNumber)
                  stress(1:3,1) = stress(1:3,1) - GeomechanicalModel%BiotCoef(elementNumber)*&
                                                  GeomechanicalModel%Porepressure(elementNumber)
                  
                  !FLUID HYDROSTATIC STRESS------------------------------------------------------------
                  HydrostaticStress = 0.
                  Zavg = SUM( GeomechanicalModel%Z( num ) ) / 8
                  HydrostaticStress(1:3) = WellData%MudWeight * Zavg  !(-)
                  stress(:,1) = stress(:,1) - HydrostaticStress

                  !ELEMENT FORCE-----------------------------------------------------------------------
                  CALL ElementStressForce(numIntegrationPoints,points,weights,coord,stress,elemNodalForce)
                  
                  !GLOBAL FORCE------------------------------------------------------------------------
                  Load( globalDOFs ) = Load( globalDOFs ) + elemNodalForce( localDOFs )
              
              END IF
              
              !Hydrostatic stress-------------------------------------------------
              IF( .NOT.GeomechanicalModel%ActiveElements(elementNumber)) THEN
                GeomechanicalModel%EffectiveStressTensor(1:3,:,elementNumber) = WellData%MudWeight * Zavg
              END IF
              
            END DO
          END DO
        END DO
      
    RETURN
    END SUBROUTINE
      
    !Routine for wellbore excavation LOAD.
    !Calculates explicitally the load at well wall.
    SUBROUTINE GeneralExcavationLoad(WellData,WellModel,GeomechanicalModel,Load)

      USE FEMlib
      USE NumericalLib
      USE SolidMechanicsLib
      USE ProblemDefinition
      USE IO
      USE MeshAuxiliares
      USE MeshIO
      USE DiscretizedModel
      USE NearWellboreModel

      IMPLICIT NONE
      TYPE(WellboreData),INTENT(IN):: WellData
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),INTENT(IN OUT):: Load(0:)

      INTEGER:: numIntegrationPoints
      PARAMETER(numIntegrationPoints = 1 )

      INTEGER:: i,k,j,m,elementNumber

      INTEGER:: numEquations
      INTEGER:: nodes(4)
      INTEGER,ALLOCATABLE:: num(:),g(:)

      REAL(8):: area,n(3),coordsTRI(3,3),area1,area2,normal1(3),normal2(3)
      REAL(8):: Zavg

      REAL(8):: theta
      REAL(8):: force(3),unbalForce(3),hydrostaticForce(3)
      REAL(8):: normal, shear, sxx, sxy, syy, stress_rot(3)

      REAL(8):: stress(6,numIntegrationPoints), HydrostaticStress   !pressao hidrostatica de fluido.

      LOGICAL:: computeExcavationLoad
      REAL(8):: z

      numEquations = GeomechanicalModel%NumEquations
      ALLOCATE( num(8), g(24) )
      Load = 0.0

      !*************************************************************************
      !STRUCTURAL LOAD.
      !*************************************************************************

      Load = 0.0
    
      !Looping na interface da parede do poço.
      !   2*******3
      !   *       *
      !   *       *
      !   *       *
      !   1*******4
      elementNumber = 0
      DO j=1, WellModel%NTHETA
        DO k=1, WellModel%NZ
          DO i=1, WellModel%NRAD
            elementNumber = elementNumber + 1
                        
            num = GeomechanicalModel%ElementNodes(:,elementNumber)
            Zavg = SUM( GeomechanicalModel%Z( num ) ) / 8

            computeExcavationLoad = .FALSE.
            IF( i==1 ) computeExcavationLoad= .TRUE.
            
            IF(computeExcavationLoad)THEN
              
              g = GeomechanicalModel%ElementDOFs(:,elementNumber)

              !FREE WALL SURFACE AREA AND NORMAL VECTOR--------------------------------------------
              coordsTRI(1,1)=GeomechanicalModel%X(num(2))
              coordsTRI(1,2)=GeomechanicalModel%Y(num(2))
              coordsTRI(1,3)=GeomechanicalModel%Z(num(2))
              coordsTRI(2,1)=GeomechanicalModel%X(num(6))
              coordsTRI(2,2)=GeomechanicalModel%Y(num(6))
              coordsTRI(2,3)=GeomechanicalModel%Z(num(6))
              coordsTRI(3,1)=GeomechanicalModel%X(num(5))
              coordsTRI(3,2)=GeomechanicalModel%Y(num(5))
              coordsTRI(3,3)=GeomechanicalModel%Z(num(5))
              CALL Triangle( coordsTRI, area1, normal1 )
              coordsTRI(1,1)=GeomechanicalModel%X(num(5))
              coordsTRI(1,2)=GeomechanicalModel%Y(num(5))
              coordsTRI(1,3)=GeomechanicalModel%Z(num(5))
              coordsTRI(2,1)=GeomechanicalModel%X(num(1))
              coordsTRI(2,2)=GeomechanicalModel%Y(num(1))
              coordsTRI(2,3)=GeomechanicalModel%Z(num(1))
              coordsTRI(3,1)=GeomechanicalModel%X(num(2))
              coordsTRI(3,2)=GeomechanicalModel%Y(num(2))
              coordsTRI(3,3)=GeomechanicalModel%Z(num(2))
              CALL Triangle( coordsTRI, area2, normal2 )
              area = area1+area2
              n = (normal1 + normal2) / 2.0
              
              !STRESS AT FREE SURFACE--------------------------------------------------------------
              stress(:,1) = GeomechanicalModel%EffectiveStressTensor(1:6,1,elementNumber)
              stress(1:3,1) = stress(1:3,1) - GeomechanicalModel%BiotCoef(elementNumber)*&
                                              GeomechanicalModel%Porepressure(elementNumber)
              CALL rot_tensao( stress(:,1), coordsTRI(1,:), coordsTRI(2,:), coordsTRI(3,:),&
                               stress_rot )    !(Sn,Sxy,Szx)
              
              !VECTORS NORMAL AND TANGENTIAL TO FREE SURFACE---------------------------------------
              !"horizontal tangent"
              coordsTRI(1,1)=0.5*( GeomechanicalModel%X(num(1)) + GeomechanicalModel%X(num(5)) )
              coordsTRI(1,2)=0.5*( GeomechanicalModel%Y(num(1)) + GeomechanicalModel%Y(num(5)) )
              coordsTRI(1,3)=0.5*( GeomechanicalModel%Z(num(1)) + GeomechanicalModel%Z(num(5)) )
              coordsTRI(2,1)=0.5*( GeomechanicalModel%X(num(3)) + GeomechanicalModel%X(num(7)) )
              coordsTRI(2,2)=0.5*( GeomechanicalModel%Y(num(3)) + GeomechanicalModel%Y(num(7)) )
              coordsTRI(2,3)=0.5*( GeomechanicalModel%Z(num(3)) + GeomechanicalModel%Z(num(7)) )
              coordsTRI(3,1)=0.5*( GeomechanicalModel%X(num(2)) + GeomechanicalModel%X(num(6)) )
              coordsTRI(3,2)=0.5*( GeomechanicalModel%Y(num(2)) + GeomechanicalModel%Y(num(6)) )
              coordsTRI(3,3)=0.5*( GeomechanicalModel%Z(num(2)) + GeomechanicalModel%Z(num(6)) )
              CALL Triangle( coordsTRI, area1, normal1 )
              
              !"vertical" tangent
              coordsTRI(1,1)=GeomechanicalModel%X(num(1))
              coordsTRI(1,2)=GeomechanicalModel%Y(num(1))
              coordsTRI(1,3)=GeomechanicalModel%Z(num(1))
              coordsTRI(2,1)=GeomechanicalModel%X(num(4))
              coordsTRI(2,2)=GeomechanicalModel%Y(num(4))
              coordsTRI(2,3)=GeomechanicalModel%Z(num(4))
              coordsTRI(3,1)=GeomechanicalModel%X(num(5))
              coordsTRI(3,2)=GeomechanicalModel%Y(num(5))
              coordsTRI(3,3)=GeomechanicalModel%Z(num(5))
              CALL Triangle( coordsTRI, area2, normal2 )  
              normal2 = -normal2

              !ELEMENT FORCE-----------------------------------------------------------------------
              unbalForce = stress_rot(1) * area * n
              unbalForce = unbalForce + stress_rot(2) * area * normal1
              unbalForce = unbalForce + stress_rot(3) * area * normal2
              
              !FLUID HYDROSTATIC STRESS------------------------------------------------------------
              HydrostaticStress = WellData%MudWeight * Zavg
              stress(:,1) = 0
              stress(1:3,1) = HydrostaticStress
              coordsTRI(1,1)=GeomechanicalModel%X(num(5))
              coordsTRI(1,2)=GeomechanicalModel%Y(num(5))
              coordsTRI(1,3)=GeomechanicalModel%Z(num(5))
              coordsTRI(2,1)=GeomechanicalModel%X(num(1))
              coordsTRI(2,2)=GeomechanicalModel%Y(num(1))
              coordsTRI(2,3)=GeomechanicalModel%Z(num(1))
              coordsTRI(3,1)=GeomechanicalModel%X(num(2))
              coordsTRI(3,2)=GeomechanicalModel%Y(num(2))
              coordsTRI(3,3)=GeomechanicalModel%Z(num(2))
              CALL rot_tensao( stress(:,1), coordsTRI(1,:), coordsTRI(2,:), coordsTRI(3,:),&
                               stress_rot )    !(Sn,Sxy,Szx)
              hydrostaticForce = stress_rot(1) * area * n
              hydrostaticForce = hydrostaticForce + stress_rot(2) * area * normal1
              hydrostaticForce = hydrostaticForce + stress_rot(3) * area * normal2
                            
              !GLOBAL FORCE------------------------------------------------------------------------
              force = hydrostaticForce - unbalForce
              nodes(1) = num(1)
              nodes(2) = num(2)
              nodes(3) = num(6)
              nodes(4) = num(5)
              Load( GeomechanicalModel%NodeDOFs(1,nodes) ) = &
                  Load( GeomechanicalModel%NodeDOFs(1,nodes) ) + 0.25 * force(1)
              Load( GeomechanicalModel%NodeDOFs(2,nodes) ) = &
                  Load( GeomechanicalModel%NodeDOFs(2,nodes) ) + 0.25 * force(2)
              Load( GeomechanicalModel%NodeDOFs(3,nodes) ) = &
                  Load( GeomechanicalModel%NodeDOFs(3,nodes) ) + 0.25 * force(3)
            
            END IF
            
            !Hydrostatic stress-------------------------------------------------
            IF( .NOT.GeomechanicalModel%ActiveElements(elementNumber)) THEN
              GeomechanicalModel%EffectiveStressTensor(1:3,:,elementNumber) = WellData%MudWeight * Zavg
            END IF
            
          END DO
        END DO
      END DO
      
    RETURN
    END SUBROUTINE
                              
    !2D formularion.
    !D.Sc.
    SUBROUTINE ExcavationLoad(WellData,WellModel,GeomechanicalModel,Load)

      USE FEMlib
      USE NumericalLib
      USE ProblemDefinition
      USE IO
      USE MeshAuxiliares
      USE MeshIO
      USE DiscretizedModel
      USE NearWellboreModel

      IMPLICIT NONE
      TYPE(WellboreData),INTENT(IN):: WellData
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),INTENT(IN OUT):: Load(0:)

      INTEGER:: numIntegrationPoints
      PARAMETER(numIntegrationPoints = 1 )

      INTEGER:: i,k,j,m,elementNumber

      INTEGER:: numEquations
      INTEGER:: nodes(4)
      INTEGER,ALLOCATABLE:: num(:),g(:)
      REAL(8),ALLOCATABLE :: coord(:,:),g_coord(:,:)  !coordenadas dos nós da malha.

      REAL(8):: area,n(3),coordsTRI(3,3),area1,area2,normal1(3),normal2(3)
      REAL(8):: Zavg

      REAL(8):: theta
      !REAL(8):: area, h
      REAL(8):: force(2) !,Fx,Fy
      REAL(8):: normal, shear, sxx, sxy, syy

      REAL(8):: stress(6,numIntegrationPoints), HydrostaticStress   !pressao hidrostatica de fluido.
      
      LOGICAL:: computeExcavationLoad
      REAL(8):: z

      numEquations = GeomechanicalModel%NumEquations
      ALLOCATE( coord(8,3), num(8), g(24), g_coord(3,GeomechanicalModel%NumNodes) )
      Load = 0.0
      
      g_coord=0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z

      !*************************************************************************
      !STRUCTURAL LOAD.
      !*************************************************************************

      Load = 0.0
    
      !Looping na interface da parede do poço.
      !   2*******3
      !   *       *
      !   *       *
      !   *       *
      !   1*******4
      elementNumber = 0
      DO j=1, WellModel%NTHETA
        DO k=1, WellModel%NZ
          DO i=1, WellModel%NRAD
            elementNumber = elementNumber + 1
            
            z = SUM( GeomechanicalModel%Z( GeomechanicalModel%ElementNodes(:,elementNumber) ) ) / 8

            computeExcavationLoad = .FALSE.
            IF( i==1 ) computeExcavationLoad= .TRUE.
            
            num = GeomechanicalModel%ElementNodes(:,elementNumber)
            coord = TRANSPOSE(g_coord(:,num))
            
            IF(computeExcavationLoad)THEN
              
              g = GeomechanicalModel%ElementDOFs(:,elementNumber)

              !Surface area and normal vector.
              coordsTRI(1,:)=coord(1,:)
              coordsTRI(2,:)=coord(2,:)
              coordsTRI(3,:)=coord(5,:)
              CALL Triangle( coordsTRI, area1, normal1 )
              coordsTRI(1,:)=coord(2,:)
              coordsTRI(2,:)=coord(6,:)
              coordsTRI(3,:)=coord(5,:)
              CALL Triangle( coordsTRI, area2, normal2 )
              area = area1+area2
              n = (normal1 + normal2) / 2.0

              !Ângulo theta em relação ao eixo X das duas bordas (esquerda e direita).
              theta = ATAN(n(2)/n(1))
              !                IF(coord(1,1)>0.001)THEN
              !                    theta1 = ATAN( coord(1,2)/coord(1,1) )
              !                ELSE
              !                    theta1 = 3.1415926535897932384626433832795 / 2.0
              !                END IF
              !                IF(coord(3,1)>0.001)THEN
              !                    theta2 = ATAN( coord(3,2)/coord(3,1) )
              !                ELSE
              !                    theta2 = pi / 2.0
              !                END IF
              !                theta = 0.5* (theta2 + theta1)

              !Leitura do tensor de tensões no sólido.
              stress(:,1) = GeomechanicalModel%EffectiveStressTensor(1:6,1,elementNumber)
              sxx = stress(1,1)
              syy = stress(2,1)
              sxy = stress(4,1)
              normal = (sxx+syy)/2.0 + (sxx-syy)/2.0*COS(2.0*theta) + sxy*SIN(2.0*theta)
              shear = - (sxx-syy)/2.0*SIN(2.0*theta) + sxy*COS(2.0*theta)
              normal = normal - GeomechanicalModel%BiotCoef(elementNumber)*&
                                GeomechanicalModel%Porepressure(elementNumber)

              !Calculando a pressão do fluido interno do poço.
              Zavg = SUM( coord(:,3) ) / 8
              HydrostaticStress = WellData%MudWeight * Zavg

              !Force.
              force(1) = (normal-HydrostaticStress) * area * COS(theta)             !componente x da força.
              force(2) = (normal-HydrostaticStress) * area * SIN(theta)             !componente y da força.
              force(1) = force(1) - shear * area * SIN(theta)                       !componente x da força.
              force(2) = force(2) + shear * area * COS(theta)                       !componente y da força.
              nodes(1) = num(1)
              nodes(2) = num(2)
              nodes(3) = num(6)
              nodes(4) = num(5)
              Load( GeomechanicalModel%NodeDOFs(1,nodes) ) = &
                  Load( GeomechanicalModel%NodeDOFs(1,nodes) ) + 0.25 * force(1)
              Load( GeomechanicalModel%NodeDOFs(2,nodes) ) = &
                  Load( GeomechanicalModel%NodeDOFs(2,nodes) ) + 0.25 * force(2)
            
            END IF
            
            !Hydrostatic stress-------------------------------------------------
            IF( .NOT.GeomechanicalModel%ActiveElements(elementNumber)) THEN
              Zavg = SUM( coord(:,3) ) / 8.0
              GeomechanicalModel%EffectiveStressTensor(1:3,:,elementNumber) = WellData%MudWeight * Zavg
            END IF
            
          END DO
        END DO
      END DO
      
    RETURN
    END SUBROUTINE
    
    
    SUBROUTINE OpenWellBC(GeomechanicalModel,ListOfRestrictedDOFs)
    
      USE DiscretizedModel

      IMPLICIT NONE

      TYPE(Model),INTENT(IN):: GeomechanicalModel
      INTEGER,INTENT(IN OUT):: ListOfRestrictedDOFs(0:)
    
      LOGICAL,ALLOCATABLE:: freeDOF(:)
        
      INTEGER::i,j,k,node,dof,numRestrictedDOFs
      
      ALLOCATE( freeDOF(0:GeomechanicalModel%NumEquations) )
      freeDOF = .TRUE.
    
      !DISPLACEMENT BOUNDARY CONDITIONS.
      node = 0
      DO j=1, GeomechanicalModel%GridDivision(3)+1
        DO k=1, GeomechanicalModel%GridDivision(2)+1
          DO i=1, GeomechanicalModel%GridDivision(1)+1
            node = node + 1
            !X face
            IF(j==1)THEN
                dof = GeomechanicalModel%NodeDOFs(2,node)
                freeDOF(dof) = .FALSE.
            END IF
            !Y face
            IF(j==(GeomechanicalModel%GridDivision(3)+1))THEN
                dof = GeomechanicalModel%NodeDOFs(1,node)
                freeDOF(dof) = .FALSE.
            END IF
            !External face
            IF(i==(GeomechanicalModel%GridDivision(1)+1))THEN
                dof = GeomechanicalModel%NodeDOFs(1,node)
                freeDOF(dof) = .FALSE.
                dof = GeomechanicalModel%NodeDOFs(2,node)
                freeDOF(dof) = .FALSE.
            END IF
            !Top face
            IF(k==1)THEN
                dof = GeomechanicalModel%NodeDOFs(3,node)
                freeDOF(dof) = .FALSE.
            END IF
            !Base face
            IF(k==(GeomechanicalModel%GridDivision(2)+1))THEN
                dof = GeomechanicalModel%NodeDOFs(3,node)
                freeDOF(dof) = .FALSE.
            END IF
          END DO
        END DO
      END DO
    
      !CORRECTING DOFs IN THE CURRENT LIST
      ListOfRestrictedDOFs = 0
      numRestrictedDOFs = 0
      DO i=1,GeomechanicalModel%NumEquations
          IF(.NOT.freeDOF(i))THEN
              numRestrictedDOFs = numRestrictedDOFs + 1
              ListOfRestrictedDOFs(numRestrictedDOFs) = i
          END IF
      END DO
      ListOfRestrictedDOFs(0) = numRestrictedDOFs

    RETURN
    END SUBROUTINE
    

  END MODULE
