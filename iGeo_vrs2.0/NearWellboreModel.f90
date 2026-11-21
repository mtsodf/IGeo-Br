
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************


  MODULE NearWellboreModel

    TYPE NWB_Model

      !***************************************************************************
      !MALHA DE ELEMENTOS FINITOS.
      !***************************************************************************
      !Malha de Elementos Finitos.
      INTEGER:: NumElements,&
                NumNodes

      !Divisão da malha estrutural.
      INTEGER:: NTHETA, NRAD, NZ

      INTEGER:: ElementType=3   !1 - tetraedro de 4 nós e 3 pontos de Gauss.
                                !2 - tetraedro de 8 nós e 6 pontos de Gauss.
                                !3 - hexaédro de 8 nós e 8 pontos de Gauss.
                                !4 - hexaédro de 14 nós e 27 pontos de Gauss.
                                !5 - hexaédro de 20 nós e 27 pontos de Gauss.

      REAL(8),ALLOCATABLE:: Coords(:,:)             !coordenadas nodais.
      INTEGER,ALLOCATABLE:: ElementIncidence(:,:)   !incidência dos elementos.

      !***************************************************************************
      !STRUCTURAL MODEL.
      !***************************************************************************
      INTEGER,ALLOCATABLE:: ElementDOFs(:,:)          !(24,element)
      INTEGER,ALLOCATABLE:: NodeDOFs(:,:)             !(3,node)
      LOGICAL,ALLOCATABLE:: RestrictedDOFs(:,:)       !(3,node)
      LOGICAL,ALLOCATABLE:: SupportedNodeDOFs(:,:)    !(3,node) Spring Support.
          
      !Estado de tensoes.
      REAL(8),ALLOCATABLE:: Porepressure(:)
      REAL(8),ALLOCATABLE:: Sxx(:),Syy(:),Szz(:),Sxy(:),Syz(:),Szx(:)
      REAL(8),ALLOCATABLE:: RadialStress(:,:),TangStress(:,:),VerticalStress(:,:)

      !Energia.
      REAL(8),ALLOCATABLE:: PlasticVolumetricEnergy(:)
      REAL(8),ALLOCATABLE:: PlasticShearEnergy(:)

      !Carga nodal equivalente a pressao de fluido do poco.
      REAL(8),ALLOCATABLE:: WellBoreNodalLoad(:)

      !Arquivos.
      CHARACTER(LEN=320)::  StressFile      !Salva o estado de tensoes nos pontos de Gauss.
          
      !Casing phases
      REAL(8):: Z_phase1

    END TYPE NWB_Model

    CONTAINS

    !Builds 'near wellbore model' object.
    SUBROUTINE BuildNearWellBoreModel(Well,FreeChanelNum,WellModel)

      USE FEMlib
      USE ProblemDefinition
      USE IO
      USE MeshAuxiliares
      USE NumericalUtilities

      IMPLICIT NONE
      TYPE(WellboreData),INTENT(IN OUT):: Well
      INTEGER,INTENT(IN)::  FreeChanelNum
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel

      INTEGER:: numIntegrationPoints
      PARAMETER( numIntegrationPoints = 1 )

      INTEGER:: i,j,k, nodeNumber, elementNumber
      INTEGER:: numElements, numNodes
      INTEGER,ALLOCATABLE:: num(:),g(:)
      INTEGER:: xconstrain, yconstrain, zconstrain

      REAL(8):: x,y,z,theta,xlength,d
      REAL(8),ALLOCATABLE :: xPillars(:), yPillars(:), zPillars(:), &
          coord(:,:)

      CHARACTER(LEN=15)::element='hexahedron'

      !Variáveis da progressão geométrica para formar os raios da malha.
      INTEGER:: n
      REAL(8):: sn,a1,q,radius
      REAL(8),ALLOCATABLE:: r(:)
      REAL(8):: pi=3.1415926535897932384626433832795

      CHARACTER(LEN=320):: text
      
      REAL(8),ALLOCATABLE:: delta(:,:)

      !***************************************************************************
      !MESH.
      !***************************************************************************
      !ALLOCATE( xPillars(Well%NTHETA+1), yPillars(Well%NRAD+1), &
      !    zPillars(Well%NZ+1) )
      ALLOCATE( xPillars(Well%NRAD+1), yPillars(Well%NTHETA+1), zPillars(Well%NZ+1) )
      xPillars=0.0
      yPillars=0.0
      zPillars=0.0

      !Profundidade das camadas.
      DO k=1 , (Well%NZ)
        zPillars(k) = Well%LayersTops(k)
      END DO
      IF(Well%NZ>1)THEN
        zPillars( Well%NZ+1 ) =  zPillars( Well%NZ ) + zPillars( Well%NZ ) - zPillars( Well%NZ-1 )
      ELSE
        zPillars( Well%NZ+1 ) =  zPillars( Well%NZ ) - 1.0
      END IF
          
      !Criando falsos pilares em x,y para posterior correção de coordenadas.
      DO i=1 , (Well%NRAD+1)
        xPillars(i) = i
      END DO
      DO j=1 , (Well%NTHETA+1)
        yPillars(j) = j
      END DO
      !Gerando malha de elementos finitos.
      ALLOCATE( coord(8,3), num(8), g(24) )

      CALL mesh_size( element,8,numElements,numNodes,Well%NRAD,Well%NTHETA,Well%NZ )
        
      WRITE(text,'(A40,I15)') 'number of FE:', numElements
      CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
      WRITE(text,'(A40,I15)') 'number of nodes:', numNodes
      CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
        
      WellModel%NumNodes = numNodes
      WellModel%NumElements = numElements

      ALLOCATE( WellModel%ElementIncidence(8,numElements),&
          WellModel%Coords(3,numNodes),&
          WellModel%ElementDOFs(24,numElements),&
          WellModel%NodeDOFs(3,numNodes), WellModel%RestrictedDOFs(3,numNodes), &
          WellModel%SupportedNodeDOFs(3,numNodes),&
          WellModel%PlasticVolumetricEnergy(numElements), &
          WellModel%PlasticShearEnergy(numElements) )
      WellModel%ElementIncidence = 0
      WellModel%Coords = 0.0
      WellModel%ElementDOFs = 0
      WellModel%PlasticVolumetricEnergy = 0.0
      WellModel%PlasticShearEnergy = 0.0
      WellModel%RestrictedDOFs = .FALSE.
      WellModel%SupportedNodeDOFs = .FALSE.

      !coordenadas nodais e incidência.
      DO elementNumber = 1, numElements
        CALL hexahedron_xz(elementNumber,xPillars,yPillars,zPillars,coord,num)
        WellModel%Coords(:,num) = TRANSPOSE(coord)
        WellModel%ElementIncidence(:,elementNumber) = num
      END DO
        
      !graus de liberdade.
      WellModel%Z_phase1 = Well%Z_phase1
      WellModel%NodeDOFs = 1
      DO nodeNumber = 1, numNodes
        x = WellModel%Coords(1,nodeNumber)
        y = WellModel%Coords(2,nodeNumber)
        z = WellModel%Coords(3,nodeNumber)

        xconstrain = 1
        yconstrain = 1
        zconstrain = 1
              
        IF(Well%MeshType=='RAD')THEN
            
          IF(y<1.01) yconstrain = 0                                     !fracture border
          
          IF(y>(Well%NTHETA+1)-1.0e-4) xconstrain = 0                   !normal to fracture border
          
          IF(x<1.01)THEN                                                !well wall.
            xconstrain = 0
            yconstrain = 0
            IF(Well%VerticaDisplacementRestriction) zconstrain = 0
          END IF
          
          IF( x>(Well%NRAD+1)-1.0e-4 ) THEN                             !external border.
            xconstrain = 0
            yconstrain = 0
          END IF
          
          IF( z < (Well%LayersTops(Well%NZ)-1.0e-4) ) THEN              !base.
            xconstrain = 0
            yconstrain = 0  
            zconstrain = 0
          END IF
          
          !>>> Model TOP.
          
          !!Casing and cemment restriction at cap rock.
          !!   Casing phases
          !!   REAL(8):: k_phase1
          !IF( (x<2.01).AND.(z>=(WellModel%Z_phase1) ) )THEN
          !  xconstrain = 0
          !  yconstrain = 0
          !  zconstrain = 0
          !END IF
          
        ELSEIF(Well%MeshType=='BOX')THEN
            
          IF(x<1.01) xconstrain = 0                                     !fronteira normal à tensão horizontal máxima.
          IF(x>(Well%NRAD+1)-1.0e-4) xconstrain = 0           
          IF(y<1.01) yconstrain = 0                                     !face normal à tensão horizontal mínima (face da propagação da fratura).
          IF(y>(Well%NTHETA+1)-1.0e-4) yconstrain = 0                   !face oposta à fratura hidráulica
          IF(z<(Well%LayersTops(Well%NZ)-1.0e-4) ) zconstrain = 0       !base
          IF(z>zPillars(1)-1.0e-4) zconstrain = 0                       !top
          
        ELSE
          CALL PrintMessage( .FALSE., .FALSE., 'ERROR: WRONG MESHTYPE' )
          CALL PrintMessage( .FALSE., .FALSE., Well%MeshType )
          CALL PrintMessage( .FALSE., .FALSE., 'FINALIZING' )
          STOP
        END IF
              
        IF(xconstrain==0) WellModel%RestrictedDOFs(1,nodeNumber) = .TRUE.
        IF(yconstrain==0) WellModel%RestrictedDOFs(2,nodeNumber) = .TRUE.
        IF(zconstrain==0) WellModel%RestrictedDOFs(3,nodeNumber) = .TRUE.
        IF(x>(Well%NRAD+1)-1.0e-4) WellModel%SupportedNodeDOFs(2,nodeNumber) = .TRUE.
              
        !all dofs will be initialy free.
        !they will be restricted with springs (SupportBoundaryConditions)
        WellModel%NodeDOFs(1,nodeNumber) = 1
        WellModel%NodeDOFs(2,nodeNumber) = 1
        WellModel%NodeDOFs(3,nodeNumber) = 1
          
      END DO
          
      CALL formnf(WellModel%NodeDOFs)
          
      DO elementNumber = 1, numElements
        CALL hexahedron_xz(elementNumber,xPillars,yPillars,zPillars,coord,num)
        CALL num_to_g(num,WellModel%NodeDOFs,g)
        WellModel%ElementDOFs(:,elementNumber) = g
      END DO

      !adjusting x,y coordinates.
      IF(Well%MeshType=='RAD')THEN
          
        ALLOCATE( r( Well%NRAD+1 ) )
        IF( Well%Radius(1)<0.001 )THEN
          r(1) = 193.68E-3/2.0                 !7 5/8"            Well%Radius(1)
                                               !r(2) = 244.48E-3/2.0          !9 5/8"
                                               !r(3) = 311.15E-3/2.0          !12 1/4"
        ELSE
          r(1) = Well%Radius(1)
        END IF
        
        sn = Well%Radius(2) - r(1)           !Well%Radius(2)-Well%Radius(1)   !311.15E-3/2.0    !9 5/8"
        n = Well%NRAD - 1                    !- 3
        q = Well%Ratio                       !1.15 !1.08 !1.18
        a1 = sn * (q-1.0) / (q**n-1)
        
        DO j=1,Well%NRAD
          r(j+1) = r(j) + a1 * q**(j-1)
        END DO
              
        DO nodeNumber = 1 , numNodes
          x = WellModel%Coords(1,nodeNumber)
          y = WellModel%Coords(2,nodeNumber)
          z = WellModel%Coords(3,nodeNumber)
          
          radius = r( NINT(x) )
          IF( (z>WellModel%Z_phase1).AND.(NINT(x).NE.(Well%NRAD+1)) .AND. ( Well%Radius(1)<0.001) )THEN
            radius = radius + (311.15E-3-193.68E-3)/2.0
          END IF
        
          theta = (y-1.0)/Well%NTHETA * pi/2.0
          WellModel%Coords(1,nodeNumber) = radius * COS(theta)
          WellModel%Coords(2,nodeNumber) = radius * SIN(theta)
        END DO
          
      ELSEIF(Well%MeshType=='BOX')THEN
        
        ALLOCATE( r( Well%NRAD+1 ) )
        IF( Well%NRAD>1 )THEN
          sn = Well%Radius(1)
          n = Well%NRAD
          !q = Well%Ratio  !1.15 !1.08 !1.18
          a1 = sn/n  !sn * (q-1.0) / (q**n-1)
          r(1) = 0.0
          DO j=1,Well%NRAD
            r(j+1) = r(j) + a1 !a1 * q**(j-1)
          END DO
        ELSE
          r(1) = 0.0
          r(2) = Well%Radius(1)
        END IF
        DO nodeNumber = 1 , numNodes
          x = WellModel%Coords(1,nodeNumber)
          WellModel%Coords(1,nodeNumber) =  r( NINT(x) )
        END DO
        
        DEALLOCATE(r)
        ALLOCATE( r( Well%NTHETA+1 ) )
        IF(Well%NTHETA>1)THEN
          sn = Well%Radius(2)
          n = Well%NTHETA
          q = Well%Ratio  !1.15 !1.08 !1.18
          a1 = sn * (q-1.0) / (q**n-1)
          r(1) = 0.0
          DO j=1,Well%NTHETA
            r(j+1) = r(j) + a1 * q**(j-1)
          END DO
        ELSE
          r(1) = 0.0
          r(2) = Well%Radius(2)
        END IF
        DO nodeNumber = 1 , numNodes
          y = WellModel%Coords(2,nodeNumber)
          WellModel%Coords(2,nodeNumber) = r( NINT(y) )
        END DO

      ELSE
            
      END IF
        
      WellModel%NTHETA = Well%NTHETA
      WellModel%NRAD = Well%NRAD
      WellModel%NZ = Well%NZ

      !WELL BORE DATA. ( #HorizontallWell )
      !REAL(8):: WellZ,WellRadius
      ALLOCATE(Well%WellBoreWallFlag(numNodes))
      Well%WellBoreWallFlag = .FALSE.
      IF(Well%HorWellFlag)THEN
        ALLOCATE(delta((Well%NRAD+1),(Well%NZ+1)))
        delta = 0.0
        Well%WellBoreWallFlag = .FALSE.
        xlength = MAXVAL(WellModel%Coords(1,:))
        nodeNumber = 0
        DO k=1,(Well%NZ+1)
          DO i=1,(Well%NRAD+1)
            nodeNumber = nodeNumber + 1
            x = WellModel%Coords(1,nodeNumber)
            z = WellModel%Coords(3,nodeNumber)
            IF( (ABS(z-Well%WellZ)<Well%WellRadius).AND.&
                (x>=(xlength-Well%WellLength)/2.0).AND.&
                (x<=xlength-(xlength-Well%WellLength)/2.0))THEN
        !      theta = ACOS(ABS(z-Well%WellZ)/Well%WellRadius)
        !      WellModel%Coords(2,nodeNumber) = Well%WellRadius * SIN(theta)
              delta(i,k) = WellModel%Coords(2,nodeNumber)
              Well%WellBoreWallFlag(nodeNumber) = .TRUE.
            END IF
          END DO
        END DO
        
        !KEEP nodeNumber
        !Make a circle for 5xWell%WellRadius 
        !DO j=2 , (Well%NTHETA+1)
        !  DO k=1,(Well%NZ+1)
        !    DO i=1,(Well%NRAD+1)
        !      nodeNumber = nodeNumber + 1
        !      y = WellModel%Coords(2,nodeNumber)
        !      z = WellModel%Coords(3,nodeNumber)
        !      d = SQRT( y**2.0 + (z-Well%WellZ)**2.0 )
        !      IF( (d < 5*Well%WellRadius ) )THEN
        !        WellModel%Coords(2,nodeNumber) = WellModel%Coords(2,nodeNumber) + delta(i,k)
        !      END IF
        !    END DO
        !  END DO
        !END DO
        
      END IF
      
      RETURN
    END SUBROUTINE BuildNearWellBoreModel

    !Builds 'near wellbore model' object.
    SUBROUTINE BuildNearWellBoreModelFull(Well,FreeChanelNum,WellModel)

      USE FEMlib
      USE ProblemDefinition
      USE IO
      USE MeshAuxiliares
      USE NumericalUtilities

      IMPLICIT NONE
      TYPE(WellboreData),INTENT(IN OUT):: Well
      INTEGER,INTENT(IN)::  FreeChanelNum
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel

      INTEGER:: numIntegrationPoints
      PARAMETER( numIntegrationPoints = 1 )

      INTEGER:: i,j,k, nodeNumber, elementNumber
      INTEGER:: numElements, numNodes
      INTEGER,ALLOCATABLE:: num(:),g(:)
      INTEGER:: xconstrain, yconstrain, zconstrain

      REAL(8):: x,y,z,theta,xlength,d
      REAL(8),ALLOCATABLE :: xPillars(:), yPillars(:), zPillars(:), &
          coord(:,:)

      CHARACTER(LEN=15)::element='hexahedron'

      INTEGER:: n
      REAL(8):: sn,a1,q,radius
      REAL(8),ALLOCATABLE:: r(:)
      REAL(8):: pi=3.1415926535897932384626433832795
      
      CHARACTER(LEN=320):: text
      
      REAL(8),ALLOCATABLE:: delta(:,:)

      !***************************************************************************
      !MESH.
      !***************************************************************************
      !ALLOCATE( xPillars(Well%NTHETA+1), yPillars(Well%NRAD+1), &
      !    zPillars(Well%NZ+1) )
      ALLOCATE( xPillars(Well%NRAD+1), yPillars(Well%NTHETA+1), zPillars(Well%NZ+1) )
      xPillars=0.0
      yPillars=0.0
      zPillars=0.0

      !Profundidade das camadas.
      DO k=1 , (Well%NZ)
        zPillars(k) = Well%LayersTops(k)
      END DO
      IF(Well%NZ>1)THEN
        zPillars( Well%NZ+1 ) =  zPillars( Well%NZ ) + zPillars( Well%NZ ) - zPillars( Well%NZ-1 )
      ELSE
        zPillars( Well%NZ+1 ) =  zPillars( Well%NZ ) - 1.0
      END IF
          
      !Criando falsos pilares em x,y para posterior correção de coordenadas.
      DO i=1 , (Well%NRAD+1)
        xPillars(i) = i
      END DO
      DO j=1 , (Well%NTHETA+1)
        yPillars(j) = j
      END DO
      !Gerando malha de elementos finitos.
      ALLOCATE( coord(8,3), num(8), g(24) )

      CALL mesh_size( element,8,numElements,numNodes,Well%NRAD,Well%NTHETA,Well%NZ )
        
      WRITE(text,'(A40,I15)') 'number of FE:', numElements
      CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
      WRITE(text,'(A40,I15)') 'number of nodes:', numNodes
      CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
        
      WellModel%NumNodes = numNodes
      WellModel%NumElements = numElements

      ALLOCATE( WellModel%ElementIncidence(8,numElements),&
          WellModel%Coords(3,numNodes),&
          WellModel%ElementDOFs(24,numElements),&
          WellModel%NodeDOFs(3,numNodes), WellModel%RestrictedDOFs(3,numNodes), &
          WellModel%SupportedNodeDOFs(3,numNodes),&
          WellModel%PlasticVolumetricEnergy(numElements), &
          WellModel%PlasticShearEnergy(numElements) )
      WellModel%ElementIncidence = 0
      WellModel%Coords = 0.0
      WellModel%ElementDOFs = 0
      WellModel%PlasticVolumetricEnergy = 0.0
      WellModel%PlasticShearEnergy = 0.0
      WellModel%RestrictedDOFs = .FALSE.
      WellModel%SupportedNodeDOFs = .FALSE.

      !coordenadas nodais e incidência.
      DO elementNumber = 1, numElements
        CALL hexahedron_xz(elementNumber,xPillars,yPillars,zPillars,coord,num)
        WellModel%Coords(:,num) = TRANSPOSE(coord)
        WellModel%ElementIncidence(:,elementNumber) = num
      END DO
      !Closing the model - radial conection.
      DO elementNumber = 1, numElements
        num = WellModel%ElementIncidence(:,elementNumber)
        y = SUM(WellModel%Coords(2,num)) / 8
        IF( y>Well%NTHETA )THEN
          num(5:8) = num(5:8) - Well%NTHETA * (Well%NRAD+1)*(Well%NZ+1)
          WellModel%ElementIncidence(:,elementNumber) = num
        END IF
      END DO
        
      !graus de liberdade.
      WellModel%Z_phase1 = Well%Z_phase1
      WellModel%NodeDOFs = 1
      DO nodeNumber = 1, numNodes
        x = WellModel%Coords(1,nodeNumber)
        y = WellModel%Coords(2,nodeNumber)
        z = WellModel%Coords(3,nodeNumber)

        xconstrain = 1
        yconstrain = 1
        zconstrain = 1

        IF( x<1.01 )THEN                                              !well wall.
          xconstrain = 0
          yconstrain = 0
          IF(Well%VerticaDisplacementRestriction) zconstrain = 0
        END IF
          
        IF( x>(Well%NRAD+1)-1.0e-4 ) THEN                             !external border.
          xconstrain = 0
          yconstrain = 0
        END IF
          
        IF( z < (Well%LayersTops(Well%NZ)-1.0e-4) ) THEN              !base.
          xconstrain = 0
          yconstrain = 0  
          zconstrain = 0
        END IF
        
        IF(y>(Well%NTHETA+1)-1.0e-4)THEN                              !non-used nodes
          xconstrain = 0
          yconstrain = 0  
          zconstrain = 0
        END IF
        
        IF( (z>=Well%OpenIntervalZ(1)).OR.(z<=Well%OpenIntervalZ(2)) ) THEN              !model top.
          xconstrain = 0
          yconstrain = 0 
          zconstrain = 0
        END IF
              
        IF(xconstrain==0) WellModel%RestrictedDOFs(1,nodeNumber) = .TRUE.
        IF(yconstrain==0) WellModel%RestrictedDOFs(2,nodeNumber) = .TRUE.
        IF(zconstrain==0) WellModel%RestrictedDOFs(3,nodeNumber) = .TRUE.
        IF(x>(Well%NRAD+1)-1.0e-4) WellModel%SupportedNodeDOFs(2,nodeNumber) = .TRUE.
          
      END DO
          
      !all dofs will be initialy free.
      !they will be restricted with springs (SupportBoundaryConditions)
      WellModel%NodeDOFs = 1
      CALL formnf(WellModel%NodeDOFs)
      DO nodeNumber = 1, numNodes
        y = WellModel%Coords(2,nodeNumber)
        IF(y>(Well%NTHETA+1)-1.0e-4)THEN                              !non-used nodes
          WellModel%NodeDOFs(:,nodeNumber) = 0
        END IF
      END DO
      
      !!Element nodes
      !DO elementNumber = 1, numElements
      !  CALL hexahedron_xz(elementNumber,xPillars,yPillars,zPillars,coord,num)
      !  CALL num_to_g(num,WellModel%NodeDOFs,g)
      !  WellModel%ElementDOFs(:,elementNumber) = g
      !END DO
      
      !Attatching nodes (360o).
      !DO nodeNumber = 1, numNodes
      !  y = WellModel%Coords(2,nodeNumber)
      !  IF( y>(Well%NTHETA+1)-1.0e-4 )THEN
      !    i = nodeNumber - Well%NTHETA * (Well%NRAD+1)*(Well%NZ+1)
      !    WellModel%NodeDOFs(:,nodeNumber) = WellModel%NodeDOFs(:,i)
      !  END IF
      !END DO
      DO elementNumber = 1, numElements
        num = WellModel%ElementIncidence(:,elementNumber)
        DO i=1,8
          g(i*3-2) = WellModel%NodeDOFs(1,num(i))
          g(i*3-1) = WellModel%NodeDOFs(2,num(i))
          g(i*3) = WellModel%NodeDOFs(3,num(i))
        END DO
        WellModel%ElementDOFs(:,elementNumber) = g
      END DO
      
      !adjusting x,y coordinates.  
      ALLOCATE( r( Well%NRAD+1 ) )
      IF( Well%Radius(1)<0.001 )THEN
        r(1) = 193.68E-3/2.0               !7 5/8"            Well%Radius(1)
                                           !r(2) = 244.48E-3/2.0          !9 5/8"
                                           !r(3) = 311.15E-3/2.0          !12 1/4"
      ELSE
        r(1) = Well%Radius(1)
      END IF
        
      sn = Well%Radius(2) - r(1)           !Well%Radius(2)-Well%Radius(1)   !311.15E-3/2.0    !9 5/8"
      n = Well%NRAD - 1                    !- 3
      q = Well%Ratio                       !1.15 !1.08 !1.18
      a1 = sn * (q-1.0) / (q**n-1)
        
      DO j=1,Well%NRAD
        r(j+1) = r(j) + a1 * q**(j-1)
      END DO
              
      DO nodeNumber = 1 , numNodes
        x = WellModel%Coords(1,nodeNumber)
        y = WellModel%Coords(2,nodeNumber)
        z = WellModel%Coords(3,nodeNumber)
          
        radius = r( NINT(x) )
        IF( (z>WellModel%Z_phase1).AND.(NINT(x).NE.(Well%NRAD+1)) .AND. ( Well%Radius(1)<0.001) )THEN
          radius = radius + (311.15E-3-193.68E-3)/2.0
        END IF
        
        theta = (y-1.0)/Well%NTHETA * 2.0*pi
        WellModel%Coords(1,nodeNumber) = radius * COS(theta)
        WellModel%Coords(2,nodeNumber) = radius * SIN(theta)
      END DO
      
      WellModel%NTHETA = Well%NTHETA
      WellModel%NRAD = Well%NRAD
      WellModel%NZ = Well%NZ
      
      ALLOCATE(Well%WellBoreWallFlag(numNodes))
      Well%WellBoreWallFlag = .FALSE.

      RETURN
    END SUBROUTINE BuildNearWellBoreModelFull

    !Escreve tensões cilíndricas na parede do poço no arquivo .log.
    SUBROUTINE WriteWellboreAnalysis( FileChanel, Path, StageName, WellModel, StressTensor )

        USE IO

        IMPLICIT NONE

        INTEGER,INTENT(IN)::  FileChanel   !Canal livre para leitura dos dados de tensões.
        TYPE (NWB_Model):: WellModel
        CHARACTER(LEN=*):: Path,StageName
      
        REAL(8),ALLOCATABLE,INTENT(IN):: StressTensor(:,:,:)

        !INTEGER:: ntheta, nk, elementNumber
        !INTEGER:: itheta,k
        !REAL(8):: theta,r,z,x,y
        !
        !INTEGER:: numIntegrationPoints
        !PARAMETER(numIntegrationPoints = 1)
        !
        !REAL(8):: sxx,syy,sxy,szz,sr,st!,str   !tensões na parede do poço.
        !
        !INTEGER:: num(8)
        !
        !REAL(8):: pi = 3.1415926535897932384626433832795
        !
        !CALL PrintMessage( .FALSE., .FALSE., '...wellbore stress analysis.' )
        !
        !!Output file.
        !OPEN(FileChanel, FILE = TRIM( TRIM(Path)//'/'//TRIM(StageName)//'.wellbore.out' ))
        !WRITE(FileChanel,*)'Wellbore Stress Analysis'
        !WRITE(FileChanel,*)''
        !
        !ntheta = WellModel%NTHETA
        !nk = WellModel%NZ
        !
        !elementNumber = 0
        !DO k=1, nk
        !    WRITE(FileChanel,*)'    k        theta              szz              sr              st'
        !    WRITE(FileChanel,*)'              deg.              kPa              kPa             kPa'
        !    DO itheta=1, ntheta
        !
        !        elementNumber = elementNumber + 1
        !        num = WellModel%ElementIncidence(:,elementNumber)
        !
        !        x = SUM( WellModel%Coords(1,num) ) / 8
        !        y = SUM( WellModel%Coords(2,num) ) / 8
        !        z = SUM( WellModel%Coords(3,num) ) / 8
        !        r = sqrt( x*x + y*y )
        !
        !        theta = ATAN(y/x)  !ângulo com o eixo x (horizontal).
        !
        !        sxx =  StressTensor(1,1,elementNumber)
        !        syy =  StressTensor(2,1,elementNumber)
        !        szz =  StressTensor(3,1,elementNumber)
        !        sxy =  StressTensor(4,1,elementNumber)
        !
        !        !Tensões cilíndricas.
        !        sr = (sxx+syy)/2.0 + (sxx-syy)/2.0*COS(2.0*theta) + sxy*SIN(2.0*theta)
        !        theta = theta + pi/2.0
        !        st = (sxx+syy)/2.0 + (sxx-syy)/2.0*COS(2.0*theta) + sxy*SIN(2.0*theta)
        !        theta = theta - pi/2.0
        !        !str = - (sxx-syy)/2.0*SIN(2.0*theta) + sxy*COS(2.0*theta)
        !
        !        !Imprimindo no Arquivo.
        !        WRITE(FileChanel,'(I6,F12.1,F20.1,F16.1,F16.1)') k,(theta/pi*180.0),szz/1.0E3,sr/1.0E3,st/1.0E3
        !
        !    END DO
        !    WRITE(FileChanel,*)'---------------------------------------------------------------------'
        !END DO
        !WRITE(FileChanel,'(A,F20.1,A)')'Radius:', r,' meter'
        !
        !CLOSE(FileChanel)

        RETURN
    END SUBROUTINE WriteWellboreAnalysis


  END MODULE NearWellboreModel

