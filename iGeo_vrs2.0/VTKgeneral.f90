
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************


  MODULE VTKgeneral


  CONTAINS

    SUBROUTINE Results2VTK(FreeChanelNum,VTKfileName,Problem,GeomechanicalModel,&
                                        Name,step)

      USE IO
      USE ProblemDefinition
      USE NearWellboreModel
      USE FEMlib
      USE SolidMechanicsLib
      USE DiscretizedModel        !modelo geomecanico discreto.
      USE VTK

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

      INTEGER,INTENT(IN):: FreeChanelNum
      CHARACTER(LEN=*):: VTKfileName
      TYPE(PD),INTENT(IN):: Problem
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      CHARACTER(*),INTENT(IN):: Name
      INTEGER,INTENT(IN):: step

      INTEGER:: vtkFileChanel
      INTEGER:: neq
      INTEGER:: i,ip,nip,dof1,dof2,dof3,nod,numElements
      PARAMETER(nod=8)

      REAL(8),ALLOCATABLE:: points(:,:), weights(:)
      REAL(8):: coord(8,3), det
      REAL(8):: der(3,8),fun(8),jac(3,3)
      INTEGER:: num(8)

      REAL(8),ALLOCATABLE::tensor_s(:,:),tensor_ep(:,:)  !tensors.
      REAL(iwp):: Sxx,Syy,Szz,Sxy,Syz,Szx !stress components.
      REAL(iwp),ALLOCATABLE:: defVol(:)
      REAL(8):: sigm, dsbar, theta, esh, deformation,S3
      REAL(8):: damage     
      
      REAL(8),ALLOCATABLE:: g_coord(:,:)
      
      REAL(8):: pp

      REAL(8):: biot,temp,pi = 3.1415926535897932384626433832795
      
      REAL(8):: f
      
      INTEGER,ALLOCATABLE:: activeNodes(:)

      CHARACTER(len=320):: stepChar,text,file

      vtkFileChanel = FreeChanelNum
      
      nip = GeomechanicalModel%NumIntegrationPoints
      neq = 3*GeomechanicalModel%NumNodes

      ALLOCATE( defVol(GeomechanicalModel%NumElements),points(nip,3), weights(nip),&
                tensor_s(6,nip),tensor_ep(6,nip))
      defVol=0.0
      
      !ACTIVE NODES++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ALLOCATE( activeNodes(GeomechanicalModel%NumNodes) )
      CALL GetActiveNodes( GeomechanicalModel, activeNodes )
      
      WRITE(stepChar,*)step
      
      !OPENING AND INITIALIZING VTK FILE+++++++++++++++++++++++++++++++++++++++++++++++++
      i = LEN_TRIM(VTKfileName)
      WRITE(text,*)VTKfileName(1:(i-4)),'_',TRIM(ADJUSTL(stepChar)),'.vtk'
      OPEN(vtkFileChanel,FILE=TRIM(ADJUSTL(text)))

      ALLOCATE( g_coord(3,GeomechanicalModel%NumNodes) )
      g_coord = 0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
      CALL HexLin_MeshToVTKfile_general( Problem%Path, Name, vtkFileChanel,&
                                          g_coord, GeomechanicalModel%ElementNodes,&
                                          GeomechanicalModel, activeNodes )
      
      !NODAL DATA++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      WRITE(vtkFileChanel,*)'POINT_DATA', MAXVAL(activeNodes)
      
      !DISPLACEMENT----------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        WRITE(text,*)'displacement' !,TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'VECTORS '//TRIM(text)//' float'

        DO i=1, GeomechanicalModel%NumNodes
          IF(activeNodes(i).NE.0)THEN
            dof1 = (i-1)*3+1  !GeomechanicalModel%NodeDOFs(1,i)
            dof2 = dof1+1     !GeomechanicalModel%NodeDOFs(2,i)
            dof3 = dof2+1     !GeomechanicalModel%NodeDOFs(3,i)
            WRITE(vtkFileChanel,'(F20.10,F20.10,F20.10)')&
                  GeomechanicalModel%Displacement(dof1),&
                  GeomechanicalModel%Displacement(dof2),&
                  GeomechanicalModel%Displacement(dof3)
          END IF
        END DO

      END IF
      
      !PRESSURE--------------------------------------------------------------------------
      WRITE(text,*)'Pp_kgf_cm2' !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumNodes
        IF(activeNodes(i).NE.0)&
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Porepressure(i) * 1.01972E-5
      END DO

      !TEMPERATURE-----------------------------------------------------------------------
      !WRITE(text,*)'RockTemperature' !,TRIM(ADJUSTL(stepChar))
      !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !DO i=1,GeomechanicalModel%NumNodes
      !  WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Temperature(i)
      !END DO
      
      !WRITE(text,*)'FluidTemperature' !,TRIM(ADJUSTL(stepChar))
      !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !DO i=1,GeomechanicalModel%NumNodes
      !  WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Temperature(i)
      !END DO

      !CELL DATA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      numElements= 0
      DO i=1,GeomechanicalModel%NumElements
        IF(GeomechanicalModel%ActiveElements(i)) numElements = numElements + 1
      END DO
      WRITE(vtkFileChanel,*)'CELL_DATA', numElements
        
      !ELEMENT POREPRESSURE--------------------------------------------------------------
      WRITE(text,*)'PpElement_kgf_cm2'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
        IF(GeomechanicalModel%ActiveElements(i))&
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%ElementAVGporepressure(i) * 1.01972E-5
      END DO
      
      !CONSTITUTIVE GEOMODEL-------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        WRITE(text,*)'ConstitutiveGeoModel' !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(I20)') GeomechanicalModel%ConstitutiveGeoModel(i)
        END DO
      END IF
      
      !HYDRAULIC FRACTURE----------------------------------------------------------------
      !DAMAGE
      IF(Problem%Structural_flag)THEN
        WRITE(text,*)'DamageShear'  !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%DamageShear(i)
        END DO
      END IF
      !IF(Problem%Structural_flag)THEN
      !  WRITE(text,*)'DamageTraction'  !_',TRIM(ADJUSTL(stepChar))
      !  WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !  WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !  DO i=1,GeomechanicalModel%NumElements
      !    IF(GeomechanicalModel%ActiveElements(i))&
      !      WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%DamageTraction(i)
      !  END DO
      !END IF
      
      !FRAC PRESSURE
      !WRITE(text,*)'FracPressure'  !_',TRIM(ADJUSTL(stepChar))
      !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !DO i=1,GeomechanicalModel%NumElements
      !  IF(GeomechanicalModel%ActiveElements(i))&
      !    WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%Fracpressure(i) * 1.01972E-5
      !END DO
      
      !STRESS----------------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        !SXX STRESS
        WRITE(text,*)'SxxeffecAVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(1,:,i)/nip) * 1.01972E-5
        END DO      
      
        !SYY STRESS
        WRITE(text,*)'SyyeffecAVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)/nip) * 1.01972E-5
        END DO 
        WRITE(text,*)'SyyAVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          biot = GeomechanicalModel%BiotCoef(i)
          num = GeomechanicalModel%ElementNodes(:,i)
          pp = GeomechanicalModel%ElementAVGporepressure(i)
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(F20.5)')&
              ( SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)/nip) - biot * pp ) * 1.01972E-5
        END DO
        !WRITE(text,*)'SyyeffecConcAVG_kgf_cm2'
        !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !DO i=1,GeomechanicalModel%NumElements
        !  damage = GeomechanicalModel%Damage(i)
        !  IF(damage>0.0)&
        !    WRITE(vtkFileChanel,'(F20.5)')&
        !                      SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)) / nip * 1.01972E-5 / (1.0-damage)
        !END DO
      
        !SZZ STRESS
        WRITE(text,*)'SzzeffecAVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(3,:,i)) / nip * 1.01972E-5
        END DO 
      
        !DEVIATORIC STRESS.
        WRITE(text,*)'J2AVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN
            temp = 0.0
            tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
            DO ip = 1, nip
              CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
              temp = temp + dsbar
            END DO
            temp = temp / nip
            WRITE(vtkFileChanel,'(F20.5)') temp * 1.01972E-5
          END IF
        END DO 
        
        !MEAN EFFECTIVE STRESS.
        WRITE(text,*)'Peffective_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN
            tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
            temp = 0.0
            DO ip = 1, nip
              CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
              temp = temp + sigm
            END DO
            temp = temp / nip
            WRITE(vtkFileChanel,'(F20.5)') temp * 1.01972E-5
          END IF
        END DO
      
        !MINIMUM EFFECTIVE STRESS.
        WRITE(text,*)'S3effec_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN            
            S3 = -1.0E30
            tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
            DO ip = 1, nip
              CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
              temp = sigm + 2.0/3.0 * dsbar * SIN( theta + 2.0 * pi / 3.0)
              IF( temp > S3 ) S3 = temp
            END DO
            WRITE(vtkFileChanel,'(F20.5)') S3 * 1.01972E-5
          END IF
        END DO
      
        !MINIMUM STRESS.
        WRITE(text,*)'S3_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN
            S3 = -1.0E30
            tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
            DO ip = 1, nip
              CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
              temp = sigm + 2.0/3.0 * dsbar * SIN( theta + 2.0 * pi / 3.0)
              IF( temp > S3 ) S3 = temp
            END DO
            biot = GeomechanicalModel%BiotCoef(i)
            pp = GeomechanicalModel%ElementAVGporepressure(i)
            WRITE(vtkFileChanel,'(F20.5)') (S3 - biot*pp) * 1.01972E-5
          END IF
        END DO
        
        !YELD FUNCTION.
        WRITE(text,*)'Yield_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN
            f =  MAXVAL(GeomechanicalModel%YieldFunc(:,i))
            WRITE(vtkFileChanel,'(F20.5)') f * 1.01972E-5
          END IF
        END DO
        
        !YOUNG MODULUS.
        !WRITE(text,*)'Young_kPa'
        !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !DO i=1,GeomechanicalModel%NumElements
        !  IF(GeomechanicalModel%ActiveElements(i))THEN
        !    WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%YoungMod(i) * 1.01972E-3
        !  END IF
        !END DO

      END IF
      
      !STRAIN----------------------------------------------------------------------------
      
      IF(Problem%Structural_flag)THEN
        
        !STRAIN YY.
        WRITE(text,*)'EyyAVG_kgf_cm2'  !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(E20.5)') SUM(GeomechanicalModel%StrainTensor(2,:,i)) / nip
        END DO
      
        !VOLUMETRIC.
        WRITE(text,*)'evol_plastic'  !,TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN
            deformation = 0.0
            tensor_s(:,:) = GeomechanicalModel%PlasticStrainTensor(:,:,i)
            DO ip = 1, nip
              CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
              deformation = deformation + sigm
            END DO
            deformation = deformation / nip
            WRITE(vtkFileChanel,'(E20.5)') deformation
          END IF
        END DO
      
        !DEVIATORIC.
        WRITE(text,*)'edev_plastic'  !,TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))THEN
            deformation = 0.0
            tensor_s(:,:) = GeomechanicalModel%PlasticStrainTensor(:,:,i)
            DO ip = 1, nip
              CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
              deformation = deformation + dsbar
            END DO
            deformation = deformation / nip
            WRITE(vtkFileChanel,'(E20.5)') deformation
          END IF
        END DO
      
      END IF

      !ENERGY----------------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        !PLASTIC STRAIN WORK.
        WRITE(text,*)'PlasticStrainWork_J'  !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          IF(GeomechanicalModel%ActiveElements(i))&
            WRITE(vtkFileChanel,'(E20.5)') GeomechanicalModel%PlasticStrainWork(i)
        END DO
      
      END IF

      CLOSE(vtkFileChanel)

      RETURN
    END SUBROUTINE
    
    SUBROUTINE HexLin_MeshToVTKfile_general( Path, FileName, FileNum, Coords, &
                                             ElementIncidence, GeomechanicalModel, &
                                             ActiveNodes )

      USE IO
      USE MeshAuxiliares
      USE DiscretizedModel
      USE VTK

      IMPLICIT NONE
      CHARACTER(LEN=*):: Path, FileName
      INTEGER,INTENT(IN):: FileNum

      REAL(8),INTENT(IN):: Coords(:,:)
      INTEGER,INTENT(IN):: ElementIncidence(:,:)
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      INTEGER,INTENT(IN):: ActiveNodes(:)

      INTEGER:: numNodes, numElements
      INTEGER:: i,j
      INTEGER:: VTKincidence(8)

      !OPEN(FileNum,FILE = TRIM( TRIM(Path)//'/'//TRIM(FileName)//'.vtk' ) )

      numNodes = MAXVAL(ActiveNodes) !SIZE( Coords(1,:) )
      numElements = SIZE( ElementIncidence(1,:) )

      WRITE(FileNum,'(A)')'# vtk DataFile Version 1.0'
      WRITE(FileNum,'(A)')'Unstructured Grid of hexahedrons'
      WRITE(FileNum,'(A)')'ASCII'
      WRITE(FileNum,'(A)')' '
      WRITE(FileNum,'(A)')'DATASET UNSTRUCTURED_GRID'
      WRITE(FileNum,'(A,I10,A)')'POINTS ', numNodes, ' float'

      !NODAL COODENATES--------------------------------------------------------
      DO i=1, GeomechanicalModel%NumNodes
        IF(ActiveNodes(i).NE.0) &
            WRITE(FileNum,'(F20.5,F20.5,F20.5)')Coords(1,i),Coords(2,i),Coords(3,i)
      END DO

      numElements= 0
      DO i=1,GeomechanicalModel%NumElements
        IF(GeomechanicalModel%ActiveElements(i)) numElements = numElements + 1
      END DO
      
      !INCIDENCE---------------------------------------------------------------
      WRITE(FileNum,'(A,I10,I10)')'CELLS ', numElements, 9*numElements
      DO i=1, GeomechanicalModel%NumElements
        IF(GeomechanicalModel%ActiveElements(i))THEN            
          VTKincidence = ElementIncidence(:,i)
          DO j=1,8
            VTKincidence(j) = ActiveNodes( VTKincidence(j) )
          END DO
          CALL Hex8_IncidenceCorrection( VTKincidence )
          VTKincidence = VTKincidence - 1
          WRITE(FileNum,'(I10,I10,I10,I10,I10,I10,I10,I10,I10)')&
                          8,&
                          VTKincidence(1),VTKincidence(2),VTKincidence(3),&
                          VTKincidence(4),VTKincidence(5),VTKincidence(6),&
                          VTKincidence(7),VTKincidence(8)
        END IF
      END DO

      !Tipo de elemento finito.
      WRITE(FileNum,'(A,I10)')'CELL_TYPES ',numElements
      DO i=1, GeomechanicalModel%NumElements
        IF(GeomechanicalModel%ActiveElements(i))&
          WRITE(FileNum,'(A)')'12'
      END DO
        
      RETURN
    END SUBROUTINE

    SUBROUTINE GetActiveNodes( GeomechanicalModel, NodeToValidNode )

      USE IO
      USE MeshAuxiliares
      USE DiscretizedModel
      USE VTK

      IMPLICIT NONE
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      INTEGER,INTENT(IN OUT):: NodeToValidNode(:)

      INTEGER:: numNodes, numElements
      INTEGER:: i,counter,num(8)
      
      LOGICAL,ALLOCATABLE:: activeNodes(:)

      numNodes = GeomechanicalModel%NumNodes
      numElements = GeomechanicalModel%NumElements
      
      NodeToValidNode = 0

      ALLOCATE( activeNodes(numNodes) )
      activeNodes = .FALSE.
      
      !Activate nodes following active elements.
      DO i=1,numElements
        IF(GeomechanicalModel%ActiveElements(i))THEN
          num = GeomechanicalModel%ElementNodes(:,i)
          activeNodes(num) = .TRUE.
        END IF
      END DO
      
      counter = 0
      DO i=1,numNodes
        IF( activeNodes(i) )THEN
          counter = counter + 1
          NodeToValidNode(i) = counter
        END IF
      END DO
      
      RETURN
    END SUBROUTINE
                                         
  END MODULE
