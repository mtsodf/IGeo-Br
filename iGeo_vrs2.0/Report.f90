
    !*****************************************************************************
    !COMPUTATIONAL GEOMECHANICS
    !Eng. Luis Carlos de Sousa Junior
    !Programa de Pos-graduação em Geotecnia.
    !University of Brasilia

  
    !*****************************************************************************
  
  MODULE Report
    
    
  CONTAINS
  
  SUBROUTINE WriteReport(FreeChanelNum,fileName,Problem,GeomechanicalModel,FlowModel,step)
  
    USE DiscretizedModel
    USE ProblemDefinition
    USE NumericalLib
    USE FluidFlowSolver
    
    IMPLICIT NONE
    INTEGER,INTENT(IN):: FreeChanelNum
    CHARACTER(LEN=*):: FileName
    TYPE(PD),INTENT(IN):: Problem
    TYPE(Model),INTENT(IN):: GeomechanicalModel
    TYPE(PorousFlowModel),INTENT(IN):: FlowModel
    INTEGER,INTENT(IN):: step
    
    INTEGER:: i,j,irad,ik,itheta,Nrad,Nk,Ntheta,numEls
    INTEGER:: nels,nip,num(8),g(8)
    
    CHARACTER(len=320):: text
    
    REAL(8):: Zavg,Pp,Q,Trad,Umax,FracRad,ShearRad,Jint,rad
    
    REAL(8):: coords1(3),coords2(3),dy
    
    nip = GeomechanicalModel%NumIntegrationPoints
    numEls = GeomechanicalModel%NumElements
    
    Nk = GeomechanicalModel%GridDivision(2)
    
    i = LEN_TRIM(FileName)
    WRITE(text,*) FileName(1:i),'.REPORT.txt'
    IF(step==0)THEN
      OPEN(FreeChanelNum,FILE=TRIM(text))
      WRITE(FreeChanelNum,*)'OUTPUT REPORT'
    ELSE
      OPEN(FreeChanelNum,FILE=TRIM(text),ACCESS='APPEND')
    END IF
    
    WRITE(FreeChanelNum,*)' '
    WRITE(FreeChanelNum,*)'Step:', step
    WRITE(FreeChanelNum,*)' '
    
    !LAYER MECHANICAL DATA.
    !--------------------------------------------------------------------------------------------------------
    !Zavg  Max.PP  Temp.Front  Max.Disp  Min.s3  Min.s3effec  TracFracLength  ShearFracLength  Max.J-integral
    !--------------------------------------------------------------------------------------------------------
    WRITE(FreeChanelNum,*)'-------------------------------------------------------------------------------'
    WRITE(FreeChanelNum,*)'   k    Zavg    Max.PP      Rate   T.Front Max.Zdisp   F.Front   S.Front  Max.J'
    WRITE(FreeChanelNum,*)'           m       kPa      m3/d         m        cm         m         m      J'
    WRITE(FreeChanelNum,*)'-------------------------------------------------------------------------------'
    DO ik=1,Nk
      
      Zavg = 0.0
      Pp = 0.0
      Q = 0.0
      Trad = 0.0
      Umax = 0.0
      FracRad = 0.0
      ShearRad = 0.0
      Jint = 0.0
      
      j = 0
      DO i=1,numEls
        IF(GeomechanicalModel%ElementIJK(2,i)==ik)THEN 
          
          j = j + 1
          num = GeomechanicalModel%ElementNodes(:,i)
          g = GeomechanicalModel%NodeDofs(3,num)
          
          rad = ( (SUM(GeomechanicalModel%X(num))/8.0)**2.0 + (SUM(GeomechanicalModel%Y(num))/8.0)**2.0 )
          
          Zavg = Zavg + SUM(GeomechanicalModel%Z(num))/8.0
          IF(GeomechanicalModel%ElementAVGporepressure(i)>Pp) Pp = GeomechanicalModel%ElementAVGporepressure(i)
          
          IF(Problem%Flow_flag)THEN
            IF(ANY(FlowModel%InjNode(num)))  Q = SUM( FlowModel%NodeFlows(num) )
          END IF
          
          Trad = 0.0
          
          IF( ABS(SUM(GeomechanicalModel%Displacement(g))/8.0)>Umax) Umax = ABS(SUM(GeomechanicalModel%Displacement(g))/8.0)
          
          IF( GeomechanicalModel%DamageTraction(i)>0.0 .AND. rad>FracRad ) FracRad = rad
          IF( ANY(ABS(GeomechanicalModel%PlasticStrainTensor(:,:,i))>0.0) .AND. rad>ShearRad ) ShearRad = rad
          
          coords1(1) = GeomechanicalModel%X( num(1) )
          coords1(2) = GeomechanicalModel%Y( num(1) )
          coords1(3) = GeomechanicalModel%Z( num(1) )
          coords2(1) = GeomechanicalModel%X( num(2) )
          coords2(2) = GeomechanicalModel%Y( num(2) )
          coords2(3) = GeomechanicalModel%Z( num(2) )
          dy = Distance(coords1,coords2)
          IF( GeomechanicalModel%YYstrainDesityEnergy(i)*dy>Jint ) Jint = GeomechanicalModel%YYstrainDesityEnergy(i)*dy
          
        END IF
      END DO
      
      Zavg = Zavg / j
      
      WRITE(FreeChanelNum,'(I5,F8.1,F10.0,F10.1,F10.1,F10.1,F10.1,F10.1,F7.0)')&
                          ik,Zavg,Pp/1E3,Q/(3600*24),Trad,Umax*1000.0,FracRad,ShearRad,Jint
      
    END DO
    
    WRITE(FreeChanelNum,*)'-------------------------------------------------------------------------------'
    
    CLOSE(FreeChanelNum)

    RETURN
  END SUBROUTINE
  
  SUBROUTINE WriteDrillingReport(FreeChanelNum,fileName,Problem,GeomechanicalModel,step)
  
    USE DiscretizedModel
    USE ProblemDefinition
    USE NumericalLib
    USE FluidFlowSolver
    USE SolidMechanicsLib
    
    IMPLICIT NONE
    INTEGER,INTENT(IN):: FreeChanelNum
    CHARACTER(LEN=*):: FileName
    TYPE(PD),INTENT(IN):: Problem
    TYPE(Model),INTENT(IN):: GeomechanicalModel
    INTEGER,INTENT(IN):: step
    
    INTEGER:: iel,j,ip,irad,ik,itheta,Nrad,Nk,Ntheta,numEls
    INTEGER:: nels,nip,num(8),g(8)
    
    CHARACTER(len=320):: text
    
    REAL(8):: Zavg,Pp,Umax,Ws,Jint,rad,FracRad
    REAL(8):: S3,dsbar,sigm,theta
    REAL(8):: temp, pi=3.1415926535897932384626433832795
    
    REAL(8):: coords1(3),coords2(3),dy
    
    nip = GeomechanicalModel%NumIntegrationPoints
    numEls = GeomechanicalModel%NumElements
    
    Nk = GeomechanicalModel%GridDivision(2)
    
    iel = LEN_TRIM(FileName)
    WRITE(text,*) FileName(1:iel),'.REPORT.txt'
    IF(step==0)THEN
      OPEN(FreeChanelNum,FILE=TRIM(text))
      WRITE(FreeChanelNum,*)'OUTPUT REPORT'
    ELSE
      OPEN(FreeChanelNum,FILE=TRIM(text),ACCESS='APPEND')
    END IF
    
    WRITE(FreeChanelNum,*)' '
    WRITE(FreeChanelNum,*)'Step:', step
    WRITE(FreeChanelNum,*)' '
    
    !LAYER MECHANICAL DATA.
    WRITE(FreeChanelNum,*)'-------------------------------------------------------------------------------'
    WRITE(FreeChanelNum,*)'   k    Zavg      Max.PP   Max.Disp    S3    Max.DESV    Ws      J-int   Frac.L'
    WRITE(FreeChanelNum,*)'          m         kPa       mm     kgf/cm²  kgf/cm²    J         J         cm'
    WRITE(FreeChanelNum,*)'-------------------------------------------------------------------------------'
    DO ik=1,Nk
      
      Zavg = 0.0
      Pp = 0.0
      Umax = 0.0
      Ws = 0.0
      S3 = -1.0E30
      dsbar = 0.0
      Jint = 0.0
      FracRad = 0.0
      
      j = 0
      DO iel=1,numEls
        IF(GeomechanicalModel%ElementIJK(2,iel)==ik)THEN     !loop o layer.
          
          j = j + 1
          num = GeomechanicalModel%ElementNodes(:,iel)
          g = GeomechanicalModel%NodeDofs(3,num)
          rad = ( (SUM(GeomechanicalModel%X(num))/8.0)**2.0 + (SUM(GeomechanicalModel%Y(num))/8.0)**2.0 )
          
          !Maximum displacemenet.
          Zavg = Zavg + SUM(GeomechanicalModel%Z(num))/8.0
          IF(GeomechanicalModel%ElementAVGporepressure(iel)>Pp) Pp = GeomechanicalModel%ElementAVGporepressure(iel)
          
          !Maximum displacement.
          IF( ABS(SUM(GeomechanicalModel%Displacement(g))/8.0)>Umax) Umax = ABS(SUM(GeomechanicalModel%Displacement(g))/8.0)
          
          !Fracture or plastic radius.
          IF( GeomechanicalModel%DamageTraction(iel)>0.0 .AND. rad>FracRad ) FracRad = rad
          IF( ANY(ABS(GeomechanicalModel%PlasticStrainTensor(:,:,iel))>0.0) .AND. rad>FracRad ) FracRad = rad
          
          !Stress invariants.
          DO ip = 1, nip
            CALL invar(GeomechanicalModel%EffectiveStressTensor(:,ip,iel),sigm,temp,theta)
            IF(temp>dsbar) dsbar = temp
            temp = sigm + 2.0/3.0 * temp * SIN( theta + 2.0 * pi / 3.0)
            IF(temp>S3) S3 = temp
          END DO
          
          !Equivalent J-Integral
          coords1(1) = GeomechanicalModel%X( num(1) )
          coords1(2) = GeomechanicalModel%Y( num(1) )
          coords1(3) = GeomechanicalModel%Z( num(1) )
          coords2(1) = GeomechanicalModel%X( num(2) )
          coords2(2) = GeomechanicalModel%Y( num(2) )
          coords2(3) = GeomechanicalModel%Z( num(2) )
          dy = Distance(coords1,coords2)
          Jint = Jint + GeomechanicalModel%YYstrainDesityEnergy(iel) * dy
          
          !Strain Energy density.
          Ws = GeomechanicalModel%PlasticStrainWork(iel)
          
        END IF
      END DO
      
      Zavg = Zavg / j
      
      WRITE(FreeChanelNum,'(I5,F10.2,F10.0,F8.0,F10.0,F10.0,E10.2,E10.2,F7.1)')&
                            ik,Zavg,Pp/1E3,Umax*1000.0,&
                            -S3*1.01972E-5,dsbar*1.01972E-5,&
                            Ws,Jint,Fracrad*100.0
      
    END DO
    
    WRITE(FreeChanelNum,*)'-------------------------------------------------------------------------------'
    
    CLOSE(FreeChanelNum)

    RETURN
  END SUBROUTINE
  
  END MODULE