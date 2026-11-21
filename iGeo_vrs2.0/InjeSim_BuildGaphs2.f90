!  InjeSim_BuildGaphs.f90 
!
!  FUNCTIONS:
!  InjeSim_BuildGaphs - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: InjeSim_BuildGaphs
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

  PROGRAM InjeSim_BuildGaphs

    IMPLICIT NONE
    
    INTEGER:: iel,ix,iy,iz,GridDivision(3),step,gp,num(8)
    INTEGER:: direction,index1_user,index2_user,nPoints,index(3)
    REAL(8),ALLOCATABLE:: arrayCells(:,:)
    INTEGER:: NumNodes
    INTEGER:: NumElements
    INTEGER:: NumIntegrationPoints
    INTEGER:: NumEquations
    REAL(8),ALLOCATABLE:: X(:)
    REAL(8),ALLOCATABLE:: Y(:)
    REAL(8),ALLOCATABLE:: Z(:)
    INTEGER,ALLOCATABLE:: ElementNodes(:,:)
    INTEGER,ALLOCATABLE:: NodeDofs(:,:)
    REAL(8)::BiotCoef
    !Transient data.
    REAL(8),ALLOCATABLE:: Porepressure(:)
    REAL(8),ALLOCATABLE:: Displacement(:)
    REAL(8),ALLOCATABLE:: EffectiveStressTensor(:,:,:)
    REAL(8),ALLOCATABLE:: StrainTensor(:,:,:)
    REAL(8),ALLOCATABLE:: PlasticStrainTensor(:,:,:)
    
    REAL(8):: sigm,dsbar,theta,e_dsbar
    
    CHARACTER(LEN=320)InputFile
    
    WRITE(*,*)'INPUT FILE:'
    READ(*,*)InputFile
    
    OPEN(100,FILE=TRIM(InputFile),FORM='UNFORMATTED')
    OPEN(10,FILE='OUTPUT.txt')
    
    READ(100)GridDivision               !NX,NZ,NY for structured meshes
    READ(100)NumNodes
    READ(100)NumElements
    READ(100)NumIntegrationPoints
    READ(100)NumEquations

    ALLOCATE( X(NumNodes) )
    ALLOCATE( Y(NumNodes) )
    ALLOCATE( Z(NumNodes) )
    ALLOCATE( ElementNodes(8,NumElements) )
    ALLOCATE( NodeDofs(3,NumNodes) )
    ALLOCATE( Porepressure(NumNodes) )
    ALLOCATE( Displacement(3*NumNodes) )
    ALLOCATE( EffectiveStressTensor(6,NumIntegrationPoints,NumElements) )
    ALLOCATE( StrainTensor(6,NumIntegrationPoints,NumElements) )
    ALLOCATE( PlasticStrainTensor(6,NumIntegrationPoints,NumElements) )

    READ(100)X
    READ(100)Y
    READ(100)Z
    READ(100)ElementNodes
    READ(100)NodeDofs
    READ(100)BiotCoef

    WRITE(*,*)'NX,NZ,NY: ',GridDivision  
    WRITE(*,*)'Number of Gauss Points: ',NumIntegrationPoints

    !Transient data.
    !WRITE(*,*)'GAUSS POINT:'
    !READ(*,*)gp
    gp=1
    WRITE(*,*)'Direction'
    WRITE(*,*)'1-x, 2-z, 3-y'
    READ(*,*)direction
    WRITE(*,*)'PAIR of IX or IZ or IY'
    READ(*,*)index1_user,index2_user

    IF(direction==1)THEN
        nPoints=GridDivision(1)
        index(1)=0
        index(2)=index1_user
        index(3)=index2_user
    ELSE IF(direction==2)THEN
        nPoints=GridDivision(2)
        index(1)=index1_user
        index(2)=0
        index(3)=index2_user
    ELSE IF(direction==3)THEN
        nPoints=GridDivision(3)
        index(1)=index1_user
        index(2)=index2_user
        index(3)=0
    ELSE
    END IF
    ALLOCATE(arrayCells(3,nPoints))
    
    WRITE(10,*)'      X(m)      Y(m)    Z(m)     p(kPa)    p_ef(kPa)  q(kPa)    Disp(m)'
    WRITE(10,*)'-----------------------------------------------------------------------'
    step = 0
    DO
      step = step + 1
      WRITE(*,*)step

      READ(100,END=200)Porepressure
      READ(100)Displacement
      READ(100)EffectiveStressTensor
      READ(100)StrainTensor
      READ(100)PlasticStrainTensor
    
      WRITE(10,*)'step:',step

      IF(direction==1)THEN

        iel= 0
        DO iy=1,GridDivision(3)
          DO iz=1,GridDivision(2)
            DO ix=1,GridDivision(1)
              iel = iel + 1

              IF( (iz==index1_user ).AND.&
                  (iy==index2_user ) )THEN
                num = ElementNodes(:,iel)
                CALL invar(PlasticStrainTensor(:,gp,iel),sigm,e_dsbar,theta)
                CALL invar(EffectiveStressTensor(:,gp,iel),sigm,dsbar,theta)
                WRITE(10,'(3F10.3,3F13.1,F13.3)')SUM(X(num))/8.,SUM(Y(num))/8.,SUM(Z(num))/8.,&
                                    SUM(Porepressure(num))/8/1e3,&
                                    sigm/1e3,&
                                    dsbar/1e3,&
                                    SUM(Displacement(num*3))/8/1e3
              END IF

            END DO
          END DO
        END DO
      
      ELSE IF(direction==2)THEN

        iel= 0
        DO iy=1,GridDivision(3)
          DO iz=1,GridDivision(2)
            DO ix=1,GridDivision(1)
              iel = iel + 1

              IF( (ix==index1_user ).AND.&
                  (iy==index2_user ) )THEN
                num = ElementNodes(:,iel)
                CALL invar(PlasticStrainTensor(:,gp,iel),sigm,e_dsbar,theta)
                CALL invar(EffectiveStressTensor(:,gp,iel),sigm,dsbar,theta)
                WRITE(10,'(3F10.3,3F10.1,F10.3)')SUM(X(num))/8.,SUM(Y(num))/8.,SUM(Z(num))/8.,&
                                    SUM(Porepressure(num))/8/1e3,&
                                    sigm/1e3,&
                                    dsbar/1e3,&
                                    SUM(Displacement(num*3))/8/1e3
              END IF

            END DO
          END DO
        END DO
      
      ELSE IF(direction==3)THEN

        iel= 0
        DO iy=1,GridDivision(3)
          DO iz=1,GridDivision(2)
            DO ix=1,GridDivision(1)
              iel = iel + 1

              IF( (ix==index1_user ).AND.&
                  (iz==index2_user ) )THEN
                num = ElementNodes(:,iel)
                CALL invar(PlasticStrainTensor(:,gp,iel),sigm,e_dsbar,theta)
                CALL invar(EffectiveStressTensor(:,gp,iel),sigm,dsbar,theta)
                WRITE(10,'(3F10.3,3F10.1,F10.3)')SUM(X(num))/8.,SUM(Y(num))/8.,SUM(Z(num))/8.,&
                                    SUM(Porepressure(num))/8/1e3,&
                                    sigm/1e3,&
                                    dsbar/1e3,&
                                    SUM(Displacement(num*3))/8/1e3
              END IF

            END DO
          END DO
        END DO

      ELSE
      END IF
      
    END DO
    
200 WRITE(*,*) 'END'

    
  CONTAINS
  
    SUBROUTINE invar(stress,sigm,dsbar,theta)
      ! Rotina de Smith e Griffiths para cálculo dos inveriantes de tensão.
      ! This subroutine forms the stress invariants in 2- or 3-d.
      !
      IMPLICIT NONE
      REAL(8),INTENT(IN)::stress(:)
      REAL(8),INTENT(OUT),OPTIONAL::sigm,dsbar,theta
      REAL(8)::sx,sy,sz,txy,dx,dy,dz,xj3,sine,s1,s2,s3,s4,s5,s6,ds1,ds2,ds3, &
          d2,d3,sq3,zero=0.0,small=1.e-10,one=1.0,two=2.0,       &
          three=3.0,six=6.0,thpt5=13.5
      INTEGER::nst
      nst=UBOUND(stress,1)
      SELECT CASE(nst)
          CASE(4)
              sx=stress(1)
              sy=stress(2)
              txy=stress(3)
              sz=stress(4)
              sigm=(sx+sy+sz)/three
              dsbar=SQRT((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+six*txy**2)/SQRT(two)
              IF(dsbar<small)THEN
                  theta=zero
              ELSE
                  dx=(two*sx-sy-sz)/three
                  dy=(two*sy-sz-sx)/three
                  dz=(two*sz-sx-sy)/three
                  xj3=dx*dy*dz-dz*txy**2
                  sine=-thpt5*xj3/dsbar**3
                  IF(sine>=one)sine=one
                  IF(sine<-one)sine=-one
                  theta=ASIN(sine)/three
              END IF
          CASE(6)
              sq3=SQRT(three)
              s1=stress(1)
              s2=stress(2)
              s3=stress(3)
              s4=stress(4)
              s5=stress(5)
              s6=stress(6)
              sigm=(s1+s2+s3)/three
              d2=((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/six+s4*s4+s5*s5+s6*s6
              ds1=s1-sigm
              ds2=s2-sigm
              ds3=s3-sigm
              d3=ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+two*s4*s5*s6
              dsbar=sq3*SQRT(d2)
              IF(dsbar<small)THEN
                  theta=zero
              ELSE
                  sine=-three*sq3*d3/(two*SQRT(d2)**3)
                  IF(sine>=one)sine=one
                  IF(sine<-one)sine=-one
                  theta=ASIN(sine)/three
              END IF
          CASE DEFAULT
              WRITE(*,*)"wrong size for nst in invar"
      END SELECT
      RETURN
    END SUBROUTINE invar

  END PROGRAM InjeSim_BuildGaphs

