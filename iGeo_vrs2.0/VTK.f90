
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************


  MODULE VTK


  CONTAINS

    SUBROUTINE WriteStepResults2VTK(FreeChanelNum,VTKfileName,Problem,GeomechanicalModel,&
                                        Name,step)

      USE IO
      USE ProblemDefinition
      USE NearWellboreModel
      USE FEMlib
      USE SolidMechanicsLib
      USE DiscretizedModel        !modelo geomecanico discreto.

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
      INTEGER:: i,ix,iy,iz,ip,nip,dof1,dof2,dof3,nod,numElements
      PARAMETER(nod=8)

      REAL(8),ALLOCATABLE:: points(:,:), weights(:)
      REAL(8):: coord(8,3), det
      REAL(8):: der(3,8),fun(8),jac(3,3)
      INTEGER:: num(8),dofs(4)

      REAL(8),ALLOCATABLE::tensor_s(:,:),tensor_ep(:,:)  !tensors.
      REAL(iwp):: Sxx,Syy,Szz,Sxy,Syz,Szx !stress components.
      REAL(iwp),ALLOCATABLE:: defVol(:)
      REAL(8):: sigm, dsbar, theta, esh, deformation,S3
      REAL(8):: damage   
      REAL(8):: ux,uy,uz
      
      REAL(8),ALLOCATABLE:: g_coord(:,:)

      REAL(8),ALLOCATABLE:: fracAperture(:)
      
      REAL(8):: pp

      REAL(8):: biot,temp,pi = 3.1415926535897932384626433832795
      
      REAL(8):: f

      CHARACTER(len=320):: stepChar,text,file

      vtkFileChanel = FreeChanelNum
      
      nip = GeomechanicalModel%NumIntegrationPoints
      neq = 3*GeomechanicalModel%NumNodes

      ALLOCATE( defVol(GeomechanicalModel%NumElements),points(nip,3), weights(nip),&
                tensor_s(6,nip),tensor_ep(6,nip), fracAperture(GeomechanicalModel%NumElements))
      defVol=0.0
      fracAperture = 0.0
      
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
      CALL HexLin_MeshToVTKfile( Problem%Path, Name, vtkFileChanel,&
                                  g_coord, GeomechanicalModel%ElementNodes)
      
      !NODAL DATA++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      WRITE(vtkFileChanel,*)'POINT_DATA', GeomechanicalModel%numNodes
      
      !DISPLACEMENT----------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        WRITE(text,*)'displacement' !,TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'VECTORS '//TRIM(text)//' float'

        DO i=1, GeomechanicalModel%NumNodes
          dof1 = GeomechanicalModel%NodeDOFs(1,i) !(i-1)*3+1  
          dof2 = GeomechanicalModel%NodeDOFs(2,i) !dof1+1     
          dof3 = GeomechanicalModel%NodeDOFs(3,i) !dof2+1  
          ux=0
          uy=0
          uz=0
          IF(dof1.NE.0) ux=GeomechanicalModel%Displacement(dof1)
          IF(dof2.NE.0) uy=GeomechanicalModel%Displacement(dof2)
          IF(dof3.NE.0) uz=GeomechanicalModel%Displacement(dof3)
          WRITE(vtkFileChanel,'(F20.10,F20.10,F20.10)')ux,uy,uz
        END DO

      END IF
      
      !PRESSURE--------------------------------------------------------------------------
      WRITE(text,*)'Pp_kgf_cm2' !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumNodes
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Porepressure(i) * 1.01972E-5
      END DO

      !TEMPERATURE-----------------------------------------------------------------------
      WRITE(text,*)'RockTemperature' !,TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumNodes
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Temperature(i)
      END DO
      
      WRITE(text,*)'FluidTemperature' !,TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumNodes
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Temperature(i)
      END DO
      
      !WRITE(text,*)'FluidRates_m3_d' !,TRIM(ADJUSTL(stepChar))
      !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !DO i=1,GeomechanicalModel%NumNodes
      !    WRITE(vtkFileChanel,'(F20.3)') -GeomechanicalModel%Flow(i)*60*24       !days
      !END DO
      !
      !WRITE(text,*)'InjectionRate_m3_d' !,TRIM(ADJUSTL(stepChar))
      !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !DO i=1,GeomechanicalModel%NumNodes
      !    WRITE(vtkFileChanel,'(F20.3)') -GeomechanicalModel%InjRate(i)*60*24    !days
      !END DO

      !CELL DATA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      numElements= 0
      DO i=1,GeomechanicalModel%NumElements
        numElements = numElements + 1
      END DO
      WRITE(vtkFileChanel,*)'CELL_DATA', numElements
        
      !ELEMENT POREPRESSURE--------------------------------------------------------------
      WRITE(text,*)'PpElement_kgf_cm2'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%ElementAVGporepressure(i) * 1.01972E-5
      END DO
      
      !CONSTITUTIVE GEOMODEL-------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        WRITE(text,*)'ConstitutiveGeoModel' !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
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
          WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%DamageShear(i)
        END DO
      END IF
      !IF(Problem%Structural_flag)THEN
      !  WRITE(text,*)'DamageTraction'  !_',TRIM(ADJUSTL(stepChar))
      !  WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !  WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !  DO i=1,GeomechanicalModel%NumElements
      !    WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%DamageTraction(i)
      !  END DO
      !END IF
      
      !FRAC APERTURE
      fracAperture = 0.0
      WRITE(text,*)'FractureAperture_mm'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
          !IF( GeomechanicalModel%IsFractured(i) )THEN
          !    num = GeomechanicalModel%ElementNodes(:,i)
          !    dofs = GeomechanicalModel%NodeDofs(2,num(1:4))
          !    fracAperture(i) = SUM(GeomechanicalModel%Displacement(dofs))/4.0
          !END IF
          fracAperture(i) = GeomechanicalModel%FracAperture(i)
          WRITE(vtkFileChanel,'(F20.5)') fracAperture(i) * 1000.0
      END DO
      
      !FRAC PRESSURE
      WRITE(text,*)'FracPressure'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
          f = 0
          IF(fracAperture(i)>0.0)  f = GeomechanicalModel%Fracpressure(i)
          WRITE(vtkFileChanel,'(F20.5)') f * 1.01972E-5
      END DO
      
      !CUMULATIVE J-INTEGRAL
      !Z-UP
      WRITE(text,*)'JintegralCum_Zup'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
          iy = GeomechanicalModel%ElementIJK(3,i)
          IF(iy==1)THEN
              ix = GeomechanicalModel%ElementIJK(1,i)
              iz = GeomechanicalModel%ElementIJK(2,i)
              f = GeomechanicalModel%JintegralCum(1,i)
          ELSE
              f = 0.
          END IF
          WRITE(vtkFileChanel,'(F20.5)') f
      END DO
      !Z-DOWN
      WRITE(text,*)'JintegralCum_Zdown'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
          iy = GeomechanicalModel%ElementIJK(3,i)
          IF(iy==1)THEN
              ix = GeomechanicalModel%ElementIJK(1,i)
              iz = GeomechanicalModel%ElementIJK(2,i)
              f = GeomechanicalModel%JintegralCum(2,i)
          ELSE
              f = 0.
          END IF
          WRITE(vtkFileChanel,'(F20.5)') f
      END DO
      !X
      WRITE(text,*)'JintegralCum_X'  !_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
          iy = GeomechanicalModel%ElementIJK(3,i)
          IF(iy==1)THEN
              ix = GeomechanicalModel%ElementIJK(1,i)
              iz = GeomechanicalModel%ElementIJK(2,i)
              f = GeomechanicalModel%JintegralCum(3,i)
          ELSE
              f = 0.
          END IF
          WRITE(vtkFileChanel,'(F20.5)') f
      END DO
      
      !STRESS----------------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        !SXX STRESS
        WRITE(text,*)'SxxeffecAVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
                WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(1,:,i)/nip) * 1.01972E-5
        END DO      
      
        !SYY STRESS
        WRITE(text,*)'SyyeffecAVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)/nip) * 1.01972E-5
        END DO 
        !WRITE(text,*)'SyyAVG_kgf_cm2'
        !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !DO i=1,GeomechanicalModel%NumElements
        !  biot = GeomechanicalModel%BiotCoef(i)
        !  num = GeomechanicalModel%ElementNodes(:,i)
        !  pp = GeomechanicalModel%ElementAVGporepressure(i)
        !  WRITE(vtkFileChanel,'(F20.5)')&
        !      ( SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)/nip) - biot * pp ) * 1.01972E-5
        !END DO
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
              WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(3,:,i)) / nip * 1.01972E-5
        END DO 
      
        !DEVIATORIC STRESS.
        WRITE(text,*)'J2AVG_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          temp = 0.0
          tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            temp = temp + dsbar
          END DO
          temp = temp / nip
          WRITE(vtkFileChanel,'(F20.5)') temp * 1.01972E-5
        END DO 
        
        !MEAN EFFECTIVE STRESS.
        WRITE(text,*)'Peffective_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
          temp = 0.0
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            temp = temp + sigm
          END DO
          temp = temp / nip
          WRITE(vtkFileChanel,'(F20.5)') temp * 1.01972E-5
        END DO
      
        !MINIMUM EFFECTIVE STRESS.
        WRITE(text,*)'S3effec_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          S3 = -1.0E30
          tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            temp = sigm + 2.0/3.0 * dsbar * SIN( theta + 2.0 * pi / 3.0)
            IF( temp > S3 ) S3 = temp
          END DO
          WRITE(vtkFileChanel,'(F20.5)') S3 * 1.01972E-5
        END DO
      
        !MINIMUM STRESS.
        WRITE(text,*)'S3_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
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
        END DO
        
        !YELD FUNCTION.
        WRITE(text,*)'Yield_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          f =  MAXVAL(GeomechanicalModel%YieldFunc(:,i))
          WRITE(vtkFileChanel,'(F20.5)') f * 1.01972E-5
        END DO
        
        !IF(Problem%MeshType=='RAD')THEN
        !
        !    !CILINDRICAL STRESS.
        !    WRITE(text,*)'RadEffec_kgf_cm2'
        !    WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !    WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !    DO i=1,GeomechanicalModel%NumElements
        !      biot = GeomechanicalModel%BiotCoef(i)
        !      num = GeomechanicalModel%ElementNodes(:,i)
        !      pp = GeomechanicalModel%ElementAVGporepressure(i)
        !      WRITE(vtkFileChanel,'(F20.5)') &
        !              -(GeomechanicalModel%CilindricalStress(1,1,i) + biot*pp) * 1.01972E-5
        !    END DO
        !    WRITE(text,*)'CircEffec_kgf_cm2'
        !    WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !    WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !    DO i=1,GeomechanicalModel%NumElements
        !      biot = GeomechanicalModel%BiotCoef(i)
        !      num = GeomechanicalModel%ElementNodes(:,i)
        !      pp = GeomechanicalModel%ElementAVGporepressure(i)
        !      WRITE(vtkFileChanel,'(F20.5)')&
        !              -(GeomechanicalModel%CilindricalStress(2,1,i) + biot*pp) * 1.01972E-5
        !    END DO
        !    WRITE(text,*)'CircShear_kgf_cm2'
        !    WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !    WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !    DO i=1,GeomechanicalModel%NumElements
        !      biot = GeomechanicalModel%BiotCoef(i)
        !      num = GeomechanicalModel%ElementNodes(:,i)
        !      pp = GeomechanicalModel%ElementAVGporepressure(i)
        !      WRITE(vtkFileChanel,'(F20.5)') &
        !          -(GeomechanicalModel%CilindricalStress(3,1,i) + biot*pp) * 1.01972E-5
        !    END DO
        !
        !END IF

      END IF
      
      !STRAIN----------------------------------------------------------------------------
      
      IF(Problem%Structural_flag)THEN
        
        !STRAIN YY.
        WRITE(text,*)'EyyAVG_kgf_cm2'  !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
            WRITE(vtkFileChanel,'(E20.5)') SUM(GeomechanicalModel%StrainTensor(2,:,i)) / nip
        END DO
      
        !VOLUMETRIC.
        WRITE(text,*)'evol_plastic'  !,TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          deformation = 0.0
          tensor_s(:,:) = GeomechanicalModel%PlasticStrainTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            deformation = deformation + sigm
          END DO
          deformation = deformation / nip
          WRITE(vtkFileChanel,'(E20.5)') deformation
        END DO
      
        !DEVIATORIC.
        WRITE(text,*)'edev_plastic'  !,TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          deformation = 0.0
          tensor_s(:,:) = GeomechanicalModel%PlasticStrainTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            deformation = deformation + dsbar
          END DO
          deformation = deformation / nip
          WRITE(vtkFileChanel,'(E20.5)') deformation
        END DO
      
      END IF

      !ENERGY----------------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        !PLASTIC STRAIN WORK.
        WRITE(text,*)'PlasticStrainWork_J_m3'  !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(E20.5)') MAXVAL(GeomechanicalModel%StrainWorkDensity(:,i))   !GeomechanicalModel%PlasticStrainWork(i)
        END DO
        
        !J-INTEGRAL.
        WRITE(text,*)'J_integral'  !_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(E20.5)') GeomechanicalModel%Jintegral(i)
        END DO
      
      END IF
      
      !PROPS-----------------------------------------------------------------------------
      WRITE(text,*)'PermH_mD'
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
        WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%KH(i) / 9.869D-16
      END DO
      
      WRITE(text,*)'PermV_mD'
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
        WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%KV(i) / 9.869D-16
      END DO
      
      WRITE(text,*)'DamagePerm'
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
        WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%DamagePerm(i)
      END DO
      
      CLOSE(vtkFileChanel)

      RETURN
    END SUBROUTINE
                                  
    !ESCREVE RESULTADOS DO step 0 PARA O VTK.
    SUBROUTINE WriteStaticDataTo2VTK(FreeChanelNum,Problem,WellData,WellModel,GeomechanicalModel)

      USE ProblemDefinition
      USE IO
      USE NearWellboreModel
      USE DiscretizedModel
      USE SolidMechanicsLib
      USE NumericalLib
      USE FEMlib

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

      INTEGER,INTENT(IN):: FreeChanelNum
      TYPE(PD),INTENT(IN):: Problem
      TYPE(WellboreData),INTENT(IN):: WellData
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel

      INTEGER:: vtkFileChanel
      INTEGER:: i,ip

      INTEGER:: nip = 1
      REAL(8):: biot, pi = 3.1415926535897932384626433832795
      REAL(8):: tensor_s(6,1),tensor_ep(6,1)   !tensor de tensões e deformações plásticas.
      REAL(8):: Sxx,Syy,Szz,Sxy,Syz,Szx     !componentes de tensão.
      REAL(8),ALLOCATABLE:: defVol(:),S3(:),S30(:)
      REAL(8):: sigm, dsbar, theta, esh

      REAL(8):: points(1,3), weights(1), coord(8,3), det
      REAL(8):: der(3,8),fun(8),jac(3,3),temp
      INTEGER:: num(8)

      CHARACTER(LEN=320):: file, text
  
      ALLOCATE( defVol(GeomechanicalModel%NumElements),S3(GeomechanicalModel%NumElements),S30(GeomechanicalModel%NumElements))
      defVol=0.0
      S3=-1.0E30
      S30 = 0.0
      !ABRINDO O ARQUIVO .vtk DO MDOELO.
      vtkFileChanel = FreeChanelNum
        
      WRITE(vtkFileChanel,*)'CELL_DATA', WellModel%numElements
        
      !POROSITY.
      WRITE(vtkFileChanel,*)'SCALARS Porosity float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%POR(i)
      END DO

      !YOUNG.
      WRITE(vtkFileChanel,*)'SCALARS Young1_GPa float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          !WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%YoungMod(i)/(1.0E9)
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%TransverIsotropicProps(i)%E1/(1.0E9)
      END DO
      WRITE(vtkFileChanel,*)'SCALARS Young2_GPa float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          !WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%YoungMod(i)/(1.0E9)
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%TransverIsotropicProps(i)%E2/(1.0E9)
      END DO

      !POISSON.
      WRITE(vtkFileChanel,*)'SCALARS Poisson1 float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          !WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Poisson(i)
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%TransverIsotropicProps(i)%n1
      END DO
      WRITE(vtkFileChanel,*)'SCALARS Poisson2 float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          !WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Poisson(i)
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%TransverIsotropicProps(i)%n2
      END DO

      !BIOT.
      WRITE(vtkFileChanel,*)'SCALARS Biot float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          biot = GeomechanicalModel%BiotCoef(i)    !Coef. de Biot.
          WRITE(vtkFileChanel,'(F20.3)') biot
      END DO
      
      !PERM. H.
      WRITE(vtkFileChanel,*)'SCALARS KH_md float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%KH(i)/9.869D-16
      END DO
      
      !PERM V.
      WRITE(vtkFileChanel,*)'SCALARS KV_md float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,WellModel%NumElements
          WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%KV(i)/9.869D-16
      END DO

      RETURN
    END SUBROUTINE
    
    !Constrói arquivo VTK para uma malha de elementos finitos hexaédrixos quadráticos.
    SUBROUTINE HexLin_MeshToVTKfile( Path, FileName, FileNum, Coords, ElementIncidence )

        USE IO
        USE MeshAuxiliares

        IMPLICIT NONE
        CHARACTER(LEN=*):: Path, FileName
        INTEGER,INTENT(IN):: FileNum

        REAL(8),INTENT(IN):: Coords(:,:)
        INTEGER,INTENT(IN):: ElementIncidence(:,:)

        INTEGER:: numNodes, numElements
        INTEGER:: i
        INTEGER:: VTKincidence(8)

        !OPEN(FileNum,FILE = TRIM( TRIM(Path)//'/'//TRIM(FileName)//'.vtk' ) )

        numNodes = SIZE( Coords(1,:) )
        numElements = SIZE( ElementIncidence(1,:) )

        WRITE(FileNum,'(A)')'# vtk DataFile Version 1.0'
        WRITE(FileNum,'(A)')'Unstructured Grid of hexahedrons'
        WRITE(FileNum,'(A)')'ASCII'
        WRITE(FileNum,'(A)')' '
        WRITE(FileNum,'(A)')'DATASET UNSTRUCTURED_GRID'
        WRITE(FileNum,'(A,I10,A)')'POINTS ', numNodes, ' float'

        !Coordenadas nodais.
        DO i=1, numNodes
          WRITE(FileNum,'(F20.5,F20.5,F20.5)')Coords(1,i),Coords(2,i),Coords(3,i)
        END DO

        !Incidência.
        numElements= SIZE(ElementIncidence(1,:))
        
        WRITE(FileNum,'(A,I10,I10)')'CELLS ', numElements, 9*numElements
        DO i=1, numElements
          VTKincidence = ElementIncidence(:,i)
          CALL Hex8_IncidenceCorrection( VTKincidence )
          VTKincidence = VTKincidence - 1
          WRITE(FileNum,'(I10,I10,I10,I10,I10,I10,I10,I10,I10)')&
                          8,&
                          VTKincidence(1),VTKincidence(2),VTKincidence(3),&
                          VTKincidence(4),VTKincidence(5),VTKincidence(6),&
                          VTKincidence(7),VTKincidence(8)
        END DO

        !Tipo de elemento finito.
        WRITE(FileNum,'(A,I10)')'CELL_TYPES ',numElements
        DO i=1, numElements
          WRITE(FileNum,'(A)')'12'
        END DO
        
        RETURN
    END SUBROUTINE HexLin_MeshToVTKfile

    !Converte a incidência nodal para formato VTK para hexaédro de 20 nós.
    SUBROUTINE Hex20_IncidenceCorrection( Incidence )

        IMPLICIT NONE
        INTEGER,INTENT(IN OUT):: Incidence(20)
        INTEGER:: oldIncidence(20)

        !Passing current values.
        oldIncidence(1) = Incidence(1)
        oldIncidence(2) = Incidence(2)
        oldIncidence(3) = Incidence(3)
        oldIncidence(4) = Incidence(4)
        oldIncidence(5) = Incidence(5)
        oldIncidence(6) = Incidence(6)
        oldIncidence(7) = Incidence(7)
        oldIncidence(8) = Incidence(8)
        oldIncidence(9) = Incidence(9)
        oldIncidence(10) = Incidence(10)
        oldIncidence(11) = Incidence(11)
        oldIncidence(12) = Incidence(12)
        oldIncidence(13) = Incidence(13)
        oldIncidence(14) = Incidence(14)
        oldIncidence(15) = Incidence(15)
        oldIncidence(16) = Incidence(16)
        oldIncidence(17) = Incidence(17)
        oldIncidence(18) = Incidence(18)
        oldIncidence(19) = Incidence(19)
        oldIncidence(20) = Incidence(20)

        !Passing back the nodes in the right order.
        Incidence(1) = oldIncidence(1)
        Incidence(2) = oldIncidence(7)
        Incidence(3) = oldIncidence(19)
        Incidence(4) = oldIncidence(13)
        Incidence(5) = oldIncidence(3)
        Incidence(6) = oldIncidence(5)
        Incidence(7) = oldIncidence(17)
        Incidence(8) = oldIncidence(15)
        Incidence(9) = oldIncidence(8)
        Incidence(10) = oldIncidence(12)
        Incidence(11) = oldIncidence(20)
        Incidence(12) = oldIncidence(9)
        Incidence(13) = oldIncidence(4)
        Incidence(14) = oldIncidence(11)
        Incidence(15) = oldIncidence(16)
        Incidence(16) = oldIncidence(10)
        Incidence(17) = oldIncidence(2)
        Incidence(18) = oldIncidence(6)
        Incidence(19) = oldIncidence(18)
        Incidence(20) = oldIncidence(14)

    END SUBROUTINE Hex20_IncidenceCorrection

    !Corrige a incidência do formato VTK para o elemento Hexaédrico de 8 nodes.
    SUBROUTINE Hex8_IncidenceCorrection( Incidence )

        IMPLICIT NONE
        INTEGER,INTENT(IN OUT):: Incidence(8)
        INTEGER:: oldIncidence(8)

        !Passing current values.
        oldIncidence(1) = Incidence(1)
        oldIncidence(2) = Incidence(2)
        oldIncidence(3) = Incidence(3)
        oldIncidence(4) = Incidence(4)
        oldIncidence(5) = Incidence(5)
        oldIncidence(6) = Incidence(6)
        oldIncidence(7) = Incidence(7)
        oldIncidence(8) = Incidence(8)

        !Passing back the nodes in the right order.
        Incidence(1) = oldIncidence(5)
        Incidence(2) = oldIncidence(6)
        Incidence(3) = oldIncidence(2)
        Incidence(4) = oldIncidence(1)
        Incidence(5) = oldIncidence(8)
        Incidence(6) = oldIncidence(7)
        Incidence(7) = oldIncidence(3)
        Incidence(8) = oldIncidence(4)

    END SUBROUTINE Hex8_IncidenceCorrection

    !Imprime resultados de célula no arquivo VTK.
    SUBROUTINE WriteCellResults2VTK(FreeChanelNum,Problem,WellData,VariableName,Vector)

        USE ProblemDefinition
        USE IO
        !USE msflib

        IMPLICIT NONE

        INTEGER,INTENT(IN):: FreeChanelNum
        TYPE(PD),INTENT(IN):: Problem
        TYPE(WellboreData),INTENT(IN):: WellData
        CHARACTER(LEN=32),INTENT(IN):: VariableName
        REAL(8),INTENT(IN):: Vector(:)

        INTEGER:: i, num
        CHARACTER(len=320):: file

        !ABRINDO ARQUIVO VTK.
        !file =TRIM(Problem%Path)//'\'//TRIM(WellData%WellName)//'.vtk'
        !OPEN( FreeChanelNum, FILE = TRIM(file), POSITION="APPEND" )

        num = SIZE(Vector)

        !ESCREVENDO RESULTADOS.
        !    WRITE(FreeChanelNum,*)'CELL_DATA', num
        WRITE(FreeChanelNum,*)'SCALARS ', TRIM(VariableName) ,' float 1'
        WRITE(FreeChanelNum,*)'LOOKUP_TABLE default'
        DO i=1, num
            WRITE(FreeChanelNum,'(E20.6)') Vector(i)
        END DO

        !CLOSE( FreeChanelNum )

    END SUBROUTINE WriteCellResults2VTK

    
    SUBROUTINE WriteStepResults2VTK_mod(FreeChanelNum,Problem,GeomechanicalModel,&
                                        Name,step)

      USE IO
      USE ProblemDefinition
      USE NearWellboreModel
      USE FEMlib
      USE SolidMechanicsLib
      USE DiscretizedModel        !modelo geomecanico discreto.

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

      INTEGER,INTENT(IN):: FreeChanelNum
      TYPE(PD),INTENT(IN):: Problem
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      CHARACTER(*),INTENT(IN):: Name
      INTEGER,INTENT(IN):: step

      INTEGER:: vtkFileChanel
      INTEGER:: neq
      INTEGER:: i,ip,nip,dof1,dof2,dof3,nod
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
      
      REAL(8):: pp

      REAL(8):: biot,temp,pi = 3.1415926535897932384626433832795      
      
      REAL(8):: f

      CHARACTER(len=320):: stepChar,text,file

      vtkFileChanel = FreeChanelNum
      
      nip = GeomechanicalModel%NumIntegrationPoints
      neq = 3*GeomechanicalModel%NumNodes

      ALLOCATE( defVol(GeomechanicalModel%NumElements),points(nip,3), weights(nip),&
                tensor_s(6,nip),tensor_ep(6,nip))
      defVol=0.0

      !NODAL DATA++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      WRITE(stepChar,*)step
      WRITE(vtkFileChanel,*)'POINT_DATA', GeomechanicalModel%numNodes
      
      !DISPLACEMENT----------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        WRITE(text,*)'displacement',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'VECTORS '//TRIM(text)//' float'

        DO i=1, GeomechanicalModel%NumNodes
          dof1 = GeomechanicalModel%NodeDOFs(1,i)   !(i-1)*3+1  !WellModel%NodeDOFs(1,i)
          dof2 = GeomechanicalModel%NodeDOFs(2,i)   !dof1+1     !WellModel%NodeDOFs(2,i)
          dof3 = GeomechanicalModel%NodeDOFs(3,i)   !dof2+1     !WellModel%NodeDOFs(3,i)
          WRITE(vtkFileChanel,'(F20.10,F20.10,F20.10)')&
                GeomechanicalModel%Displacement(dof1),&
                GeomechanicalModel%Displacement(dof2),&
                GeomechanicalModel%Displacement(dof3)
        END DO

      END IF
      
      !PRESSURE--------------------------------------------------------------------------
      WRITE(text,*)'Pp_kgf_cm2_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumNodes
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Porepressure(i) * 1.01972E-5
      END DO

      !TEMPERATURE-----------------------------------------------------------------------
      WRITE(text,*)'T',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumNodes
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%Temperature(i)
      END DO

      !CELL DATA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(vtkFileChanel,*)'CELL_DATA', GeomechanicalModel%numElements
        
      !ELEMENT POREPRESSURE--------------------------------------------------------------
      WRITE(text,*)'PpElement_kgf_cm2_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
        WRITE(vtkFileChanel,'(F20.3)') GeomechanicalModel%ElementAVGporepressure(i) * 1.01972E-5
      END DO
      
      !CONSTITUTIVE GEOMODEL-------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        WRITE(text,*)'ConstitutiveGeoModel_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
            WRITE(vtkFileChanel,'(I20)') GeomechanicalModel%ConstitutiveGeoModel(i)
        END DO
      END IF
      
      !HYDRAULIC FRACTURE----------------------------------------------------------------
      !DAMAGE
      !IF(Problem%Structural_flag)THEN
      !  WRITE(text,*)'Damage_',TRIM(ADJUSTL(stepChar))
      !  WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      !  WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      !  DO i=1,GeomechanicalModel%NumElements
      !    WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%DamageTraction(i)
      !  END DO
      !END IF
      
      !FRAC PRESSURE
      WRITE(text,*)'FracPressure_',TRIM(ADJUSTL(stepChar))
      WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
      WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
      DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(F20.5)') GeomechanicalModel%Fracpressure(i) * 1.01972E-5
      END DO
      
      !STRESS----------------------------------------------------------------------------
      IF(Problem%Structural_flag)THEN
        
        !SXX STRESS
        WRITE(text,*)'SxxeffecAVG_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(1,:,i)/nip) * 1.01972E-5
        END DO      
      
        !SYY STRESS
        WRITE(text,*)'SyyeffecAVG_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)/nip) * 1.01972E-5
        END DO 
        WRITE(text,*)'SyyAVG_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          biot = GeomechanicalModel%BiotCoef(i)
          num = GeomechanicalModel%ElementNodes(:,i)
          pp = GeomechanicalModel%ElementAVGporepressure(i)
          WRITE(vtkFileChanel,'(F20.5)')&
              ( SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)/nip) - biot * pp ) * 1.01972E-5
        END DO
        !WRITE(text,*)'SyyeffecConcAVG_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        !WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        !WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        !DO i=1,GeomechanicalModel%NumElements
        !  damage = GeomechanicalModel%Damage(i)
        !  IF(damage>0.0)&
        !   WRITE(vtkFileChanel,'(F20.5)')&
        !                   SUM(GeomechanicalModel%EffectiveStressTensor(2,:,i)) / nip * 1.01972E-5 / (1.0-damage)
        !END DO
      
        !SZZ STRESS
        WRITE(text,*)'SzzeffecAVG_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(F20.5)') SUM(GeomechanicalModel%EffectiveStressTensor(3,:,i)) / nip * 1.01972E-5
        END DO 
      
        !DEVIATORIC STRESS.
        WRITE(text,*)'J2AVG_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          temp = 0.0
          tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            temp = temp + dsbar
          END DO
          temp = temp / nip
          WRITE(vtkFileChanel,'(F20.5)') temp * 1.01972E-5
        END DO      
      
        !MINIMUM EFFECTIVE STRESS.
        WRITE(text,*)'S3effec_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          S3 = -1.0E30
          tensor_s(:,:)=GeomechanicalModel%EffectiveStressTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            temp = sigm + 2.0/3.0 * dsbar * SIN( theta + 2.0 * pi / 3.0)
            IF( temp > S3 ) S3 = temp
          END DO
          WRITE(vtkFileChanel,'(F20.5)') S3 * 1.01972E-5
        END DO
      
        !MINIMUM STRESS.
        WRITE(text,*)'S3_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
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
        END DO
        
        !YELD FUNCTION.
        WRITE(text,*)'Yield_kgf_cm2'
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          f =  MINVAL(GeomechanicalModel%YieldFunc(:,i))
          WRITE(vtkFileChanel,'(F20.5)') f * 1.01972E-5
        END DO
        
        IF(Problem%MeshType=='RAD')THEN
        
        !CILINDRICAL STRESS.
        WRITE(text,*)'RadEffec_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          biot = GeomechanicalModel%BiotCoef(i)
          num = GeomechanicalModel%ElementNodes(:,i)
          pp = GeomechanicalModel%ElementAVGporepressure(i)
          WRITE(vtkFileChanel,'(F20.5)') &
                  -(GeomechanicalModel%CilindricalStress(1,1,i) + biot*pp) * 1.01972E-5
        END DO
        WRITE(text,*)'CircEffec_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          biot = GeomechanicalModel%BiotCoef(i)
          num = GeomechanicalModel%ElementNodes(:,i)
          pp = GeomechanicalModel%ElementAVGporepressure(i)
          WRITE(vtkFileChanel,'(F20.5)')&
                  -(GeomechanicalModel%CilindricalStress(2,1,i) + biot*pp) * 1.01972E-5
        END DO
        WRITE(text,*)'CircShear_kgf_cm2_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          biot = GeomechanicalModel%BiotCoef(i)
          num = GeomechanicalModel%ElementNodes(:,i)
          pp = GeomechanicalModel%ElementAVGporepressure(i)
          WRITE(vtkFileChanel,'(F20.5)') &
              -(GeomechanicalModel%CilindricalStress(3,1,i) + biot*pp) * 1.01972E-5
        END DO
      
        END IF

      END IF
      !STRAIN----------------------------------------------------------------------------
      
      IF(Problem%Structural_flag)THEN
        
        !STRAIN YY.
        WRITE(text,*)'Eyy_vp_AVG_adm_',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          WRITE(vtkFileChanel,'(E20.5)') SUM(GeomechanicalModel%PlasticStrainTensor(2,:,i)) / nip
        END DO
      
        !VOLUMETRIC.
        WRITE(text,*)'evol_plastic',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          deformation = SUM( GeomechanicalModel%PlasticStrainTensor(1:3,:,i) ) / nip
          WRITE(vtkFileChanel,'(E20.5)') deformation
        END DO
      
        !DEVIATORIC.
        WRITE(text,*)'edev_plastic',TRIM(ADJUSTL(stepChar))
        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
        DO i=1,GeomechanicalModel%NumElements
          deformation = 0.0
          tensor_s(:,:) = GeomechanicalModel%PlasticStrainTensor(:,:,i)
          DO ip = 1, nip
            CALL invar(tensor_s(:,ip),sigm,dsbar,theta)
            deformation = deformation + dsbar
          END DO
          deformation = deformation / nip
          WRITE(vtkFileChanel,'(E20.5)') deformation
        END DO
      
      END IF

!        !PLASTIC VOLUMETRIC WORK.
!        WRITE(text,*)'Wp_vol',TRIM(ADJUSTL(stepChar))
!        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
!        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
!        DO i=1,GeomechanicalModel%NumElements
!            WRITE(vtkFileChanel,'(F20.5)') WellModel%PlasticVolumetricEnergy(i)
!        END DO
!
!        !PLASTIC SHEAR WORK.
!        WRITE(text,*)'Wp_shear',TRIM(ADJUSTL(stepChar))
!        WRITE(vtkFileChanel,*)'SCALARS ',TRIM(text),' float 1'
!        WRITE(vtkFileChanel,*)'LOOKUP_TABLE default'
!        DO i=1,GeomechanicalModel%NumElements
!            WRITE(vtkFileChanel,'(F20.5)') WellModel%PlasticShearEnergy(i)
!        END DO

      !CLOSE(vtkFileChanel)

      RETURN
    END SUBROUTINE
    
    
END MODULE VTK
