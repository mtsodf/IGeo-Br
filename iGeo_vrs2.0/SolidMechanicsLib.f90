
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE SolidMechanicsLib


  CONTAINS
  
  SUBROUTINE TensorDeformacaoTermica(CoefThermalExpansion,DeltaT,eterm)
    !Rotina para cálculo do tensor de deformação térmica.
    !Os arquivos de propriedades dinâmicas (14) e (28) precisam estar posicionado no passo de tempo correto e no elemento correto.
    
    !VARIÁVEIS.
    IMPLICIT NONE
    REAL(8),INTENT(IN)::CoefThermalExpansion
    REAL(8),INTENT(OUT)::eterm(6)  !deformação de origem térmica
    REAL(8):: DeltaT
    

    REAL(8):: alfa,&                            !coeficiente de dilatação térmica
              zero=0.0,one=1.0
    
    !Calculando a deformação térmica
  
    alfa = CoefThermalExpansion                 !coeficiente de dilatação linear
    eterm = zero
    eterm = (/one,one,one,zero,zero,zero/)
    eterm = eterm*alfa*DeltaT
  
  RETURN
  END SUBROUTINE TensorDeformacaoTermica

    !Calcula Estado de Tensões na parede do poço.
    FUNCTION WellboreStress(angle,sxx,syy,sxy)

        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp):: WellboreStress(2)  !1 - tensão radial, 2 - tensão tangencial.
        REAL(iwp),INTENT(IN):: angle                  !ângulo em graus.
        REAL(iwp),INTENT(IN):: sxx,syy,sxy            !tensões cartesianas na parede do poço.

        REAL(iwp):: pi = 3.1415926535897932384626433832795
        REAL(iwp):: theta

        !Rotacionando tensões.
        !"angle" é o ângulo com o eixo x.
        theta = angle
        WellboreStress(1) = (sxx+syy)/2.0 + (sxx-syy)/2.0*COS(2.0*theta) + sxy*SIN(2.0*theta)
        theta = angle + pi/2.0
        WellboreStress(2) = (sxx+syy)/2.0 + (sxx-syy)/2.0*COS(2.0*theta) + sxy*SIN(2.0*theta)

    END FUNCTION

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

    SUBROUTINE TensorEig(tensorVector,eigenValues,eigenVectors)
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: tensorVector(6)
      REAL(8),INTENT(OUT):: eigenValues(3),eigenVectors(3,3)
      
      REAL(8),ALLOCATABLE:: tensor(:,:)

      INTEGER::        N=3
      INTEGER::        LDA,LDVL,LDVR
      INTEGER::        LWMAX = 1000
      INTEGER::        INFO, LWORK
      REAL(8),ALLOCATABLE:: VL(:,:),VR(:,:),WR(:),WI(:),WORK(:)
      
      LDA = N
      LDVL = N
      LDVR = N
       
      ALLOCATE( tensor( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),WR( N ), WI( N ), WORK( LWMAX ))
      
      tensor(1,1) = tensorVector(1)
      tensor(1,2) = tensorVector(4)
      tensor(1,3) = tensorVector(6)
      tensor(2,1) = tensorVector(4)
      tensor(2,2) = tensorVector(2)
      tensor(2,3) = tensorVector(5)
      tensor(3,1) = tensorVector(6)
      tensor(3,2) = tensorVector(5)
      tensor(3,3) = tensorVector(3)

      !Query the optimal workspace.
      LWORK = -1
      CALL DGEEV( 'Vectors', 'Vectors', N, tensor, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      
      !Solve eigenproblem.
      CALL DGEEV( 'Vectors', 'Vectors', N, tensor, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
      
      eigenValues = WR
      eigenVectors(:,1) = VR(:,1)
      eigenVectors(:,2) = VR(:,2)
      eigenVectors(:,3) = VR(:,3)
     
      RETURN
    END SUBROUTINE

    SUBROUTINE TensorEigenvalues(tensorVector,eigenValues)
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: tensorVector(6)
      REAL(8),INTENT(OUT):: eigenValues(3)
      
      REAL(8),ALLOCATABLE:: tensor(:,:)

      INTEGER::        N=3
      INTEGER::        LDA,LDVL,LDVR
      INTEGER::        LWMAX = 1000
      INTEGER::        INFO, LWORK
      REAL(8),ALLOCATABLE:: VL(:,:),VR(:,:),WR(:),WI(:),WORK(:)
      
      LDA = N
      LDVL = N
      LDVR = N
       
      ALLOCATE( tensor( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),WR( N ), WI( N ), WORK( LWMAX ))
      
      tensor(1,1) = tensorVector(1)
      tensor(1,2) = tensorVector(4)
      tensor(1,3) = tensorVector(6)
      tensor(2,1) = tensorVector(4)
      tensor(2,2) = tensorVector(2)
      tensor(2,3) = tensorVector(5)
      tensor(3,1) = tensorVector(6)
      tensor(3,2) = tensorVector(5)
      tensor(3,3) = tensorVector(3)

      !Query the optimal workspace.
      LWORK = -1
      CALL DGEEV( 'Vectors', 'Vectors', N, tensor, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      
      !Solve eigenproblem.
      CALL DGEEV( 'Vectors', 'Vectors', N, tensor, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
      
      eigenValues = WR
     
      RETURN
    END SUBROUTINE
    
    SUBROUTINE RotateStress(tensorVector,normal,outTensorVector)
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: tensorVector(6)
      REAL(8),INTENT(IN):: normal(3,3)
      REAL(8),INTENT(IN OUT):: outTensorVector(6)
      
      REAL(8):: tensor(3,3)
      
      ! Matriz de rotação
      REAL(8), dimension(3,3) :: R

      tensor(1,1) = tensorVector(1)
      tensor(1,2) = tensorVector(4)
      tensor(1,3) = tensorVector(6)
      tensor(2,1) = tensorVector(4)
      tensor(2,2) = tensorVector(2)
      tensor(2,3) = tensorVector(5)
      tensor(3,1) = tensorVector(6)
      tensor(3,2) = tensorVector(5)
      tensor(3,3) = tensorVector(3)
      
      !Matriz de rotação 
      R(:,1) = normal(:,1)
      R(:,2) = normal(:,2)
      R(:,3) = normal(:,3)

      !Matriz constitutiva elástica global 
      tensor =  MATMUL(TRANSPOSE(R),MATMUL(tensor,R))
      
      outTensorVector(1) = tensor(1,1)
      outTensorVector(2) = tensor(2,2)
      outTensorVector(3) = tensor(3,3)
      outTensorVector(4) = tensor(1,2)
      outTensorVector(5) = tensor(2,3)
      outTensorVector(6) = tensor(1,3)
      
    RETURN
    END SUBROUTINE
        
    !***********************************************************************
    !																	   !
    ! Tecnologia em Computação Gráfica - TECGRAF						   !
    ! Pontifícia Universidade Católica do Rio de Janeiro - PUC-Rio         !
    !                                                                      !
    ! Grupo de Geomecânica Computacional                                   !
    ! Coordenadora: Deane Roehl                                            !
    ! Contato: geomec@tecgraf.puc-rio.br                                   !
    !                                                                      !
    !***********************************************************************
    !																	   !
    ! Rotina:  rot_tensao											       !
    ! Arquivo: rot_tensao.f90 											   !
    ! Data:    05/08/2011												   !
    ! Autor:   Marianna A. de O. Coelho									   !
    !																	   !
    ! Modificações													       !
    !   Data:															   !
    !																	   !
    !***********************************************************************
    ! 																	   !
    ! Descrição: Rotaciona tensões para um plano de de falha			   !
    ! 																	   !
    ! Argumentos:														   !
    !   TENSAO_G   - Vetor de tensões no sistema global de coordenadas com !
    !                 6 componentes (Sxx, Syy, Szz, Sxy, Syz, Szx)         !
    !   P1, P2, P3 - Vetores de coordenadas (x,y) que definem um triângulo !
    !                 no de falha, no sentido anti-horário                 !
    !   TENSAO_L   - Vetor de tensões no sistema local de coordenadas com  !
    !                 3 componentes (Sn,Sxy,Szx), normal e cisalhantes     !
    !                 respectivamente                                      !
    !																	   !
    !***********************************************************************
    subroutine rot_tensao( TENSAO_G, P1, P2, P3, TENSAO_L )

      implicit none

      ! Variáveis
      ! Tensões no sistema global e coordenadas do triângulo
      real(8), dimension(:), intent(in) :: TENSAO_G(6), P1(3), P2(3), P3(3)
      ! Tensões no sistema local da falha
      real(8), dimension(:), intent(out) :: TENSAO_L(3)
      ! Vetores unitários e vetores auxiliares
      real(8), dimension(3) :: X, Y, Z, XL, YL, ZL, XL_AUX, YL_AUX, ZL_AUX,& 
						       COORD_AUX, COORD_AUX2
      ! Tensores
      real(8), dimension(3,3) :: R, TENSAO_S, TENSAO_F


      ! Vetores unitários do sistema global
      X = (/1.00_8, 0.00_8, 0.00_8/)   !Vetor base x
      Y = (/0.00_8, 1.00_8, 0.00_8/)   !Vetor base y
      Z = (/0.00_8, 0.00_8, 1.00_8/)   !Vetor base z

      ! Vetores normalizados no plano do triângulo
      ! Vetor local x normalizado no plano da falha
      XL_AUX = P2 - P1  
      XL = XL_AUX / DSQRT(XL_AUX(1)**2 + XL_AUX(2)**2 + XL_AUX(3)**2)  
      ! Vetor local auxiliar normalizado no plano do triângulo
      COORD_AUX2 = P3 - P1
      COORD_AUX  = COORD_AUX2 / DSQRT(COORD_AUX2(1)**2 + COORD_AUX2(2)**2 +&
                    COORD_AUX2(3)**2)  

      ! Produto vetorial dos vetores normalizados no plano (XL x COORD_AUX)
      ZL_AUX(1) = XL(2)*COORD_AUX(3) - XL(3)*COORD_AUX(2)
      ZL_AUX(2) = XL(3)*COORD_AUX(1) - XL(1)*COORD_AUX(3)
      ZL_AUX(3) = XL(1)*COORD_AUX(2) - XL(2)*COORD_AUX(1)
      ! Vetor local Z normalizado no plano da falha
      ZL = ZL_AUX / DSQRT(ZL_AUX(1)**2 + ZL_AUX(2)**2 + ZL_AUX(3)**2)  
    
      ! Produto vetorial dos vetores locais (ZL x XL)
      YL_AUX(1) = ZL(2)*XL(3) - ZL(3)*XL(2)
      YL_AUX(2) = ZL(3)*XL(1) - ZL(1)*XL(3)
      YL_AUX(3) = ZL(1)*XL(2) - ZL(2)*XL(1)

      ! Vetor local y normalizado no plano da falha
      YL = YL_AUX / DSQRT(YL_AUX(1)**2 + YL_AUX(2)**2 + YL_AUX(3)**2)  

      ! Matriz de rotação entre o sistema local da falha e o sistema global
      R(1,1) = DOT_PRODUCT(XL,X)
      R(1,2) = DOT_PRODUCT(XL,Y)
      R(1,3) = DOT_PRODUCT(XL,Z)
      R(2,1) = DOT_PRODUCT(YL,X)
      R(2,2) = DOT_PRODUCT(YL,Y)
      R(2,3) = DOT_PRODUCT(YL,Z)
      R(3,1) = DOT_PRODUCT(ZL,X)
      R(3,2) = DOT_PRODUCT(ZL,Y)
      R(3,3) = DOT_PRODUCT(ZL,Z)

      ! Tensões no sistema global (TENSAO_G) em forma matricial
      TENSAO_S(1,1) = TENSAO_G(1)   !xx
      TENSAO_S(1,2) = TENSAO_G(4)   !xy
      TENSAO_S(1,3) = TENSAO_G(6)   !xz
      TENSAO_S(2,1) = TENSAO_G(4)   !yx
      TENSAO_S(2,2) = TENSAO_G(2)   !yy
      TENSAO_S(2,3) = TENSAO_G(5)   !yz
      TENSAO_S(3,1) = TENSAO_G(6)   !zx
      TENSAO_S(3,2) = TENSAO_G(5)   !zy
      TENSAO_S(3,3) = TENSAO_G(3)   !zz

      ! Rotação do tensor de tensões para o plano do triângulo
      TENSAO_F = MATMUL(R,MATMUL(TENSAO_S,TRANSPOSE(R)))

      ! Componentes de tensão normal e cisalhante no plano do triângulo
      TENSAO_L(1) = TENSAO_F(3,3)
      !TENSAO_L(2) = DSQRT(TENSAO_F(1,3)**2 + TENSAO_F(2,3)**2)
      TENSAO_L(2) = TENSAO_F(1,3)
      TENSAO_L(3) = TENSAO_F(2,3)

    end subroutine

    !SUBROUTINE S3direction(tensorVector,S3,S3_normal)
    !
    !  IMPLICIT NONE
    !  REAL(8),INTENT(IN):: tensorVector(6),S3
    !  REAL(8),INTENT(OUT):: S3_normal(3)
    !  
    !  REAL(8):: a(6), rhs(3)=0.0
    !  INTEGER:: ia(4), ja(6)
    !  
    !  a = 0.0
    !  
    !  a(1) = tensorVector(1) - S3
    !  ia(1) = 1
    !  ja(1) = 1
    !  a(2) = tensorVector(4)
    !  ja(1) = 2
    !  a(3) = tensorVector(6)
    !  ja(3) = 3
    !  
    !  a(4) = tensorVector(2) - S3
    !  ia(2) = 4
    !  ja(1) = 2
    !  a(5) = tensorVector(5)
    !  ja(5) = 3
    !  
    !  a(6) = tensorVector(3) - S3
    !  ia(3) = 4
    !  ja(6) = 3
    !  
    !  ia(4) = 4
    !  
    !  CALL LinearSolverPAR(3,6,a,ja,ia,rhs)
    !  S3_normal = rhs
    !  
    !  RETURN
    !END SUBROUTINE
    !
    SUBROUTINE D_HomogeneousElasticMatrix(dee,e,v)
        !Subrotina para cálculo da matriz elástica [D] para materiais elásticos homogêneos.
        !
        !Variáveis
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::e,v
        REAL(iwp),INTENT(IN OUT)::dee(:,:)
        REAL(iwp)::v1,v2,c,vv,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp
        INTEGER::i
        !
        !Qualidade
        IF((v.GT.0.5_iwp).OR.(v.LT.zero))THEN
            WRITE(*,*)"ERRO NA MONTAGEM DA MATRIZ ELASTICA [D]"
            WRITE(*,*)"POISSON FORA DO LIMITE ACEITAVEL {0.0 , 0.5}"
            WRITE(*,*)"subrotina: D_HomogeneousElasticMatrix (D_HomogeneousElasticMatrix.dll)"
            WRITE(*,*)"FINALIZANDO..."
            STOP
        END IF
        IF(e.LT.zero)THEN
            WRITE(*,*)"ERRO NA MONTAGEM DA MATRIZ ELASTICA [D]"
            WRITE(*,*)"MODULO DE YOUNG NEGATIVO"
            WRITE(*,*)"subrotina: D_HomogeneousElasticMatrix (D_HomogeneousElasticMatrix.dll)"
            WRITE(*,*)"FINALIZANDO..."
            STOP
        END IF
        !
        !Construção da matriz D
        dee=zero
        v1=one-v
        c=e/((one+v)*(one-two*v))
        v2=v/(one-v)
        vv=(one-two*v)/(one-v)*pt5
        DO i=1,3
            dee(i,i)=one
        END DO
        DO i=4,6
            dee(i,i)=vv
        END DO
        dee(1,2)=v2
        dee(2,1)=v2
        dee(1,3)=v2
        dee(3,1)=v2
        dee(2,3)=v2
        dee(3,2)=v2
        dee=dee*e/(two*(one+v)*vv)
        !
        RETURN
    END SUBROUTINE D_HomogeneousElasticMatrix

    SUBROUTINE ThermalStrain(alfa,DeltaT,eterm)

        IMPLICIT NONE
        REAL(8),INTENT(IN):: DeltaT
        REAL(8),INTENT(IN)::alfa
        REAL(8),INTENT(OUT)::eterm(6)  !deformação de origem térmica

        REAL(8):: zero = 0.0, one = 1.0

        !calculando acréscimo de temperatura e deformação térmica
        eterm=zero
        eterm=(/one,one,one,zero,zero,zero/)
        eterm=eterm*alfa*DeltaT

        RETURN
    END SUBROUTINE ThermalStrain

    SUBROUTINE LinearSolverPAR(neq,a,ja,ia,loads)
      
      !---------------------------------------------------------------------------------------------------------------
      !Solver linear direto PARDISO(R): parallel direct sparse solver. Interl Fortran.
      !Matrized não-simétricas no formato CSR.
      !---------------------------------------------------------------------------------------------------------------
      
      USE OMP_LIB  !bibliotecas OpenMP(R)
      IMPLICIT NONE
      
      !variáveis estáticas
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER::iparm(64), pt(64),&  !parâmetros do PARDISO
                idum,maxfct, mnum, mtype, phase, n, nrhs, error, msglvl  !parâmetros do PARDISO
      INTEGER,INTENT(IN)::neq  !número de equações
      REAL(iwp)::ddum,&  !parâmetros do PARDISO
                  zero=0.0_iwp
      !variáveis dinâmicas
      INTEGER,INTENT(IN)::ia(:),&  !identificação da posição de início das linhas no vetor a(:)
                          ja(:)  !identificação das colunas dos elementos não-zero da matriz
      REAL(iwp),ALLOCATABLE::x(:),&  !vetor solução
                              b(:)  !vetor de cargas e resultado do sistema A*x=b
      REAL(iwp),INTENT(IN)::a(:)  !vetor de elementos não zero da matriz A do sistema A*x=b
      REAL(iwp),INTENT(IN OUT)::loads(:)  !vetor de carregamento estrutural (IN) e resultados (OUT)

      !alocação de memória
      ALLOCATE(b(neq),x(neq))
      
      !Ajustando parâmetros de controle do PARDISO
      n=neq  !number of equations
      nrhs=1  !number of right-hand side to be solved
      maxfct=1  !maxumum number of factors in the memory at the same time
      mnum=1  !actual matriz for solution phase ( 1 ≤ mnum ≤ maxfct )
      
      iparm = 0 !inicializando
      iparm(1) = 1  !=0 for solver default
      iparm(2) = 3  !fill-in reducing ordering for the input matrix
                    !=0 for minimum degree algorithm (VERY SLOW)
                    !=2 nested dissection algorithm from the METIS package
                    !=3 parallel (OpenMP) version of the nested dissection algorithm is used
      iparm(3) = 0  !currently is not in use
      iparm(4) = 31  !Preconditioned CGS. It´s possible use CG instead of LU as pre-conditionator.
                      !=0 no iterative-direct algorithm. 
                      !=31 LU-preconditioned CGS iteration with a stopping criterion of 1.0E-3 for unsymmetric matrices.
                      !=61 LU-preconditioned CGS iteration with a stopping criterion of 1.0E-6 for unsymmetric matrices.
                      !=62 LU-preconditioned CGS iteration with a stopping criterion of 1.0E-6 for symmetric matrices 
      iparm(5) = 0  !no user fill-in reducing permutation
      iparm(6) = 0  !=0 solution on the first n compoments of x
                    !=1 a solução é escrita em b ( Ax=b )
      iparm(7) = 0  !not in use
      iparm(8) = 0  !iterative refinement step.
      iparm(9) = 0  !not in use
      iparm(10) = 13  !pivoting perturbation.
                      !The default value of iparm(10) is 13 and therefore eps = 1.0E-13 for unsymmetric matrices 
      iparm(11) = 1  !scaling vectors.
                      !use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy <<<<<<<<<<?
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1  !Output: number of nonzeros in the factor LU
      iparm(19) = -1  !Output: Mflops for LU factorization
      iparm(20) = 0  !Output: Numbers of CG Iterations
      iparm(60) = 0  !=0 in-core
                      !=2 for out-of-core OOC algorithm (=2 fica muito lento)
      
      error = 0 ! initialize error flag
      msglvl = 0  !Message level information
                  !=0 for NO output information
                  !=1 for print statistical information
      mtype = 11  !Kind of the matrix
                  ! = 2 for symmetric, positive definite matrix
                  ! =11 for real and unsymmetric matrix
      
      !Initiliaze the internal solver memory pointer (for memory access). This is only
      !necessary for the FIRST call of the PARDISO solver.
      pt = 0 
      
      !Reordering and Symbolic Factorization. This step also allocates
      !all memory that is necessary for the factorization

      !WRITE(*,*)'  ...reordering and symbolic factorization'
      phase = 11 ! only reordering and symbolic factorization <<<<<CONTROLS CALL pardiso
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,&
                    idum, nrhs, iparm, msglvl, ddum, ddum, error)
      IF (error .NE. 0) THEN
        WRITE(*,*)'  the following ERROR was detected: ', error
        STOP
      END IF
      
      !Factorization.
      WRITE(*,*)'  ...factorization'
      phase = 22 ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,&
                    idum, nrhs, iparm, msglvl, ddum, ddum, error)
      IF (error .NE. 0) THEN
        WRITE(*,*) '  the following ERROR was detected: ', error
        STOP
      ENDIF

      !Back substitution and iterative refinement
      WRITE(*,*)'  ...back substitution and iterative refinement'
      iparm(8) = 2 ! max numbers of iterative refinement steps
      phase = 33 ! solving the system
      b = loads(2:neq+1)
      x = zero
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,&
                    idum, nrhs, iparm, msglvl, b, x, error)
      loads(2:neq+1)=x
      WRITE(*,*)'  ...CGS diagnostics = ',iparm(20)

      !Termination and release of memory
      phase = -1 ! release internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
                    idum, nrhs, iparm, msglvl, ddum, ddum, error)
      
      !WRITE(*,*) '  PARDISO(R) solve completed'
      
    RETURN
    END SUBROUTINE LinearSolverPAR


END MODULE SolidMechanicsLib
