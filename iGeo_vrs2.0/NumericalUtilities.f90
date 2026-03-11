  
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE NumericalUtilities

  CONTAINS
  


  SUBROUTINE MatrizDerivadaFuncForma(i,ndim,nod,coord,points,det,deriv,bee)
  !Rotina para c�lculo da matriz de [B] dee derivadas das
  !fun��es de forma em um dado ponto de Gauss.
  !{deforma��o}=[B]*{deslocamentos nodais}
  !
  !USE NumericalUtilities
  !vari�veis
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  INTEGER,INTENT(IN)::i,&  !n�mero do ponto de Gauss
  ndim,&  !n�mero de dimens�es
  nod  !n�mero de n�s de deslocamento por elemento
  REAL(iwp),INTENT(IN OUT)::det  !determinante de matriz
  REAL(iwp),ALLOCATABLE::der(:,:),&  !derivadas das fun��es de forma com rela��o �s coordenadas locais
  jac(:,:)    !matriz jacobiana das fun��es de forma
  REAL(iwp),ALLOCATABLE,INTENT(IN)::coord(:,:),&  !coordenada dos n�s de um elemento
  points(:,:)  !coordenadas dos pontos de integra��o
  REAL(iwp),ALLOCATABLE,INTENT(IN OUT)::bee(:,:),&  !matriz de derivadas das fun��es de forma
  deriv(:,:)  !derivadas das fun��es de forma com rela��o �s coordenadas globais
  ALLOCATE(der(ndim,nod),jac(ndim,ndim))
  !
  CALL shape_der(der,points,i)  !fun��es de forma [N]
  jac=MATMUL(der,coord)         !jacobiano da transforma��o de sistema de coordenadas [J]
  CALL determinante(jac,det)    !determinante do Jacobiano det[J]
  CALL invert(jac)              !inverso do Jacobiano: 1/[J]
  deriv=MATMUL(jac,der)         !matriz de derivadas
  CALL beemat(bee,deriv)        !matriz [B]
  !
  RETURN
  END SUBROUTINE MatrizDerivadaFuncForma


  !rotina para formar matriz de graus de liberdade nf(nodof,nn)
  ! Smith & Griffiths 2004
  SUBROUTINE formnf(nf)
  !
  ! This subroutine forms the nf matrix.
  !
  !IMPLICIT NONE
  INTEGER,INTENT(IN OUT)::nf(:,:)
  INTEGER::i,j,m
  m=0
  DO j=1,UBOUND(nf,2)
    DO i=1,UBOUND(nf,1)
      IF(nf(i,j)/=0)THEN
        m=m+1
        nf(i,j)=m
      END IF
    END DO
  END DO
  RETURN
  END SUBROUTINE formnf

  SUBROUTINE num_to_g(num,nf,g)
  !
  ! This subroutine finds the g vector from num and nf.
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN)::num(:),nf(:,:)
  INTEGER,INTENT(OUT)::g(:)
  INTEGER::i,k,nod,nodof
  nod=UBOUND(num,1)
  nodof=UBOUND(nf,1)
  DO i=1,nod
    k=i*nodof
    g(k-nodof+1:k)=nf(:,num(i))
  END DO
  RETURN
  END SUBROUTINE num_to_g

  SUBROUTINE determinante(jac,det)
  !
  ! This function returns the determinant of a 1x1, 2x2 or 3x3
  ! Jacobian matrix.
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp),INTENT(IN)::jac(:,:)
  REAL(iwp),INTENT(IN OUT)::det
  INTEGER::it
  it=UBOUND(jac,1)
  det=0.0
  SELECT CASE(it)
  CASE(1)
    det=1.0_iwp
  CASE(2)
    det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
  CASE(3)
    det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
    det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
    det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
    CASE DEFAULT
    WRITE(*,*)' wrong dimension for Jacobian matrix'
  END SELECT
  RETURN
  END SUBROUTINE determinante
  !
  SUBROUTINE invert(matrix)
  !
  ! This subroutine inverts a small square matrix onto itself.
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp),INTENT(IN OUT)::matrix(:,:)
  REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
  INTEGER::ndim,i,k
  ndim=UBOUND(matrix,1)
  IF(ndim==2)THEN
    det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
    j11=matrix(1,1)
    matrix(1,1)=matrix(2,2)
    matrix(2,2)=j11
    matrix(1,2)=-matrix(1,2)
    matrix(2,1)=-matrix(2,1)
    matrix=matrix/det
  ELSE IF(ndim==3)THEN
    det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
    det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
    det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
    j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
    j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
    j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
    j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
    j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
    j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
    j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
    j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
    j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
    matrix(1,1)=j11
    matrix(1,2)=j12
    matrix(1,3)=j13
    matrix(2,1)=j21
    matrix(2,2)=j22
    matrix(2,3)=j23
    matrix(3,1)=j31
    matrix(3,2)=j32
    matrix(3,3)=j33
    matrix=matrix/det
  ELSE
    DO k=1,ndim
      con=matrix(k,k)
      matrix(k,k)=1.0_iwp
      matrix(k,:)=matrix(k,:)/con
      DO i=1,ndim
        IF(i/=k)THEN
          con=matrix(i,k)
          matrix(i,k)=0.0_iwp
          matrix(i,:)=matrix(i,:)-matrix(k,:)*con
        END IF
      END DO
    END DO
  END IF
  RETURN
  END SUBROUTINE invert
  !
  SUBROUTINE beemat(bee,deriv)
  !
  ! This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp),INTENT(IN)::deriv(:,:)
  REAL(iwp),INTENT(OUT)::bee(:,:)
  INTEGER::k,l,m,n,ih,nod
  REAL(iwp)::x,y,z
  bee=0.0_iwp
  ih=UBOUND(bee,1)
  nod=UBOUND(deriv,2)
  SELECT CASE (ih)
  CASE(3,4)
    DO m=1,nod
      k=2*m
      l=k-1
      x=deriv(1,m)
      y=deriv(2,m)
      bee(1,l)=x
      bee(3,k)=x
      bee(2,k)=y
      bee(3,l)=y
    END DO
  CASE(6)
    DO m=1,nod
      n=3*m
      k=n-1
      l=k-1
      x=deriv(1,m)
      y=deriv(2,m)
      z=deriv(3,m)
      bee(1,l)=x
      bee(4,k)=x
      bee(6,n)=x
      bee(2,k)=y
      bee(4,l)=y
      bee(5,n)=y
      bee(3,n)=z
      bee(5,k)=z
      bee(6,l)=z
    END DO
    CASE DEFAULT
    WRITE(*,*)'wrong dimension for nst in bee matrix'
  END SELECT
  RETURN
  END SUBROUTINE beemat
  !
  SUBROUTINE cross_product(b,c,a)
  !
  ! This subroutine forms the cross product of two vectors, a = b x c
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp),INTENT(IN)::b(:),c(:)
  REAL(iwp),INTENT(OUT)::a(:,:)
  INTEGER::ib,ic,i,j
  ib=SIZE(b)
  ic=SIZE(c)
  DO i=1,ib
    DO j=1,ic
      a(i,j)=b(i)*c(j)
    END DO
  END DO
  RETURN
  END SUBROUTINE cross_product
  !
  SUBROUTINE shape_fun(fun,points,i)
  !
  !   This subroutine computes the values of the shape functions.
  !   to local coordinates
  !
  !   points(:,:) - coordenadas naturais dos pontos de integra��o
  !   i - n�mero do ponto de integra��o
  !   fun(:) - retorna as fun��es de forma no ponto de integra��o i
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  INTEGER,INTENT(in)::i
  REAL(iwp),INTENT(IN)::points(:,:)
  REAL(iwp),INTENT(OUT)::fun(:)
  REAL(iwp)::eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3
  REAL(iwp)::L1,L2,L3,L4 !coordenadas naturais para o tetraedro
  REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9
  REAL(iwp)::zeta,xi0,eta0,zeta0
  INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
  REAL,PARAMETER::pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,pt75=0.75_iwp, &
                  one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d8=8.0_iwp,d9=9.0_iwp,   &
                  d16=16.0_iwp,d27=27.0_iwp,d32=32.0_iwp,d64=64.0_iwp,d128=128.0_iwp
  ndim=UBOUND(points,2)
  nod=UBOUND(fun,1)
  SELECT CASE(ndim)
  CASE(1) ! one dimensional case
    xi=points(i,1)
    SELECT CASE(nod)
    CASE(2)
      t1=-one-xi
      t2= one-xi
      fun(1)=t2/two
      fun(2)=-t1/two
    CASE(3)
      t1=-one-xi
      t2=-xi
      t3=one-xi
      fun(1)=t2*t3/two
      fun(2)=-t1*t3
      fun(3)=t1*t2/two
    CASE(4)
      t1=-one-xi
      t2=-one/d3-xi
      t3=one/d3-xi
      t4=one-xi
      fun(1)=t2*t3*t4*d9/d16
      fun(2)=-t1*t3*t4*d27/d16
      fun(3)=t1*t2*t4*d27/d16
      fun(4)=-t1*t2*t3*d9/d16
    CASE(5)
      t1=-one -xi
      t2=-pt5-xi
      t3=-xi
      t4=pt5-xi
      t5=one-xi
      fun(1)=t2*t3*t4*t5*two/d3
      fun(2)=-t1*t3*t4*t5*d8/d3
      fun(3)=t1*t2*t4*t5*d4
      fun(4)=-t1*t2*t3*t5*d8/d3
      fun(5)=t1*t2*t3*t4*two/d3
      CASE DEFAULT
      WRITE(*,*)"wrong number of nodes in shape_fun"
    END SELECT
  CASE(2) ! two dimensional case
    c1=points(i,1)
    c2=points(i,2)
    c3=one-c1-c2
    xi=points(i,1)
    eta=points(i,2)
    etam=pt25*(one-eta)
    etap=pt25*(one+eta)
    xim=pt25*(one-xi)
    xip=pt25*(one+xi)
    SELECT CASE(nod)
    CASE(3)
      fun = (/c1,c3,c2/)
    CASE(6)
      fun(1)=(two*c1-one)*c1
      fun(2)=d4*c3*c1
      fun(3)=(two*c3-one)*c3
      fun(4)=d4*c2*c3
      fun(5)=(two*c2-one)*c2
      fun(6)=d4*c1*c2
    CASE(10)
      fun(1)= ((d3*c1-one)*(d3*c1-two)*c1)/two
      fun(2)= -(d9*(d3*c1-one)*(c1+c2-one)*c1)/two
      fun(3)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c1)/two
      fun(4)=-((d3*c1+d3*c2-one)*(d3*c1+d3*c2-two)*(c1+c2-one))/two
      fun(5)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c2)/two
      fun(6)= -(d9*(c1+c2-one)*(d3*c2-one)*c2)/two
      fun(7)= ((d3*c2-one)*(d3*c2-two)*c2)/two
      fun(8)=  (d9*(d3*c2-one)*c1*c2)/two
      fun(9)=  (d9*(d3*c1-one)*c1*c2)/two
      fun(10)=-d27*((c2-one)+c1)*c1*c2
    CASE(15)
      t1=c1-pt25
      t2=c1-pt5
      t3=c1-pt75
      t4=c2-pt25
      t5=c2-pt5
      t6=c2-pt75
      t7=c3-pt25
      t8=c3-pt5
      t9=c3-pt75
      fun(1)=d32/d3*c1*t1*t2*t3
      fun(2)=d128/d3*c3*c1*t1*t2
      fun(3)=d64*c3*c1*t1*t7
      fun(4)=d128/d3*c3*c1*t7*t8
      fun(5)=d32/d3*c3*t7*t8*t9
      fun(6)=d128/d3*c2*c3*t7*t8
      fun(7)=d64*c2*c3*t4*t7
      fun(8)=d128/d3*c2*c3*t4*t5
      fun(9)=d32/d3*c2*t4*t5*t6
      fun(10)=d128/d3*c1*c2*t4*t5
      fun(11)=d64*c1*c2*t1*t4
      fun(12)=d128/d3*c1*c2*t1*t2
      fun(13)=d128*c1*c2*t1*c3
      fun(15)=d128*c1*c2*c3*t4
      fun(14)=d128*c1*c2*c3*t7
    CASE(4)
      fun=(/d4*xim*etam,d4*xim*etap,d4*xip*etap,d4*xip*etam/)
    CASE(5)
      fun=(/d4*xim*etam-pt25*(one-xi**2)*(one-eta**2),             &
      d4*xim*etap-pt25*(one-xi**2)*(one-eta**2),             &
      d4*xip*etap-pt25*(one-xi**2)*(one-eta**2),             &
      d4*xip*etam-pt25*(one-xi**2)*(one-eta**2),             &
      (one-xi**2)*(one-eta**2)/)
    CASE(8)
      fun=(/d4*etam*xim*(-xi-eta-one),d32*etam*xim*etap,                   &
      d4*etap*xim*(-xi+eta-one),d32*xim*xip*etap,                    &
      d4*etap*xip*(xi+eta-one), d32*etap*xip*etam,                   &
      d4*xip*etam*(xi-eta-one), d32*xim*xip*etam/)
    CASE(9)
      etam=eta-one
      etap=eta+one
      xim=xi-one
      xip=xi+one
      fun=(/pt25*xi*xim*eta*etam,-pt5*xi*xim*etap*etam,                    &
      pt25*xi*xim*eta*etap,-pt5*xip*xim*eta*etap,                    &
      pt25*xi*xip*eta*etap,-pt5*xi*xip*etap*etam,                    &
      pt25*xi*xip*eta*etam,-pt5*xip*xim*eta*etam,                    &
      xip*xim*etap*etam/)
      CASE DEFAULT
      WRITE(*,*)"wrong number of nodes in shape_fun"
    END SELECT
  CASE(3) ! d3 dimensional case
    xi=points(i,1)
    eta=points(i,2)
    zeta=points(i,3)
    etam=one-eta
    xim=one-xi
    zetam=one-zeta
    etap=eta+one
    xip=xi+one
    zetap=zeta+one
    SELECT CASE(nod)
    CASE(4)
      fun(1)=xi
      fun(2)=eta
      fun(3)=zeta
      fun(4)=one-fun(1)-fun(2)-fun(3)
    CASE(8)
      fun=(/pt125*xim*etam*zetam,pt125*xim*etam*zetap,                     &
      pt125*xip*etam*zetap,pt125*xip*etam*zetam,                     &
      pt125*xim*etap*zetam,pt125*xim*etap*zetap,                     &
      pt125*xip*etap*zetap,pt125*xip*etap*zetam/)
    CASE(10)!ELEMENTO TETRA�DRICO QUADR�TICO--------------------------------
      L1=points(i,1);L2=points(i,2);L3=points(i,3);L4=one-L1-L2-L3
      fun(1) = (2*L1-1)*L1
      fun(2) = (2*L2-1)*L2
      fun(3) = (2*L3-1)*L3
      fun(4) = (2*L4-1)*L4
      fun(5) =  4*L1*L2
      fun(6) =  4*L2*L3
      fun(7) =  4*L1*L3
      fun(8) =  4*L1*L4
      fun(9) =  4*L2*L4
      fun(10) =  4*L3*L4
      !Zienkiewicz e Taylor vol.1 pp 188
      !-----------------------------------------------------------------------
    CASE(14) ! type 6 element
      fun(1) = (xi*eta+xi*zeta+two*xi+eta*zeta+two*eta+two*zeta+two)*      &
      (xi-one)*(eta-one)*(zeta-one)/d8
      fun(2) =-(xi*eta-xi*zeta+two*xi-eta*zeta+two*eta-two*zeta+two)*      &
      (xi-one)*(eta-one)*(zeta+one)/d8
      fun(3) =-(xi*eta-xi*zeta+two*xi+eta*zeta-two*eta+two*zeta-two)*      &
      (xi+one)*(eta-one)*(zeta+one)/d8
      fun(4) = (xi*eta+xi*zeta+two*xi-eta*zeta-two*eta-two*zeta-two)*      &
      (xi+one)*(eta-one)*(zeta-one)/d8
      fun(5) =-(xi+one)*(xi-one)*(eta-one)*(zeta+one)*(zeta-one)/two
      fun(6) =-(xi-one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
      fun(7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta+one)/two
      fun(8) = (xi+one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
      fun(9) =-(xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta-one)/two
      fun(10)= (xi*eta-xi*zeta-two*xi+eta*zeta+two*eta-two*zeta-two)*      &
      (xi-one)*(eta+one)*(zeta-one)/d8
      fun(11)=-(xi*eta+xi*zeta-two*xi-eta*zeta+two*eta+two*zeta-two)*      &
      (xi-one)*(eta+one)*(zeta+one)/d8
      fun(12)=-(xi*eta+xi*zeta-two*xi+eta*zeta-two*eta-two*zeta+two)*      &
      (xi+one)*(eta+one)*(zeta+one)/d8
      fun(13)= (xi*eta-xi*zeta-two*xi-eta*zeta-two*eta+two*zeta+two)*      &
      (xi+one)*(eta+one)*(zeta-one)/d8
      fun(14)= (xi+one)*(xi-one)*(eta+one)*(zeta+one)*(zeta-one)/two
    CASE(20)
      xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
      etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
      zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
      DO l=1,20
        xi0=xi*xii(l)
        eta0=eta*etai(l)
        zeta0=zeta*zetai(l)
        IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
          fun(l)=pt25*(one-xi*xi)*(one+eta0)*(one+zeta0)
        ELSE IF(l>=9.AND.l<=12)THEN
          fun(l)=pt25*(one+xi0)*(one-eta*eta)*(one+zeta0)
        ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18)THEN
          fun(l)=pt25*(one+xi0)*(one+eta0)*(one-zeta*zeta)
        ELSE
          fun(l)=pt125*(one+xi0)*(one+eta0)*(one+zeta0)*(xi0+eta0+zeta0-2)
        END IF
      END DO
      CASE DEFAULT
      WRITE(*,*)"wrong number of nodes in shape_fun"
    END SELECT
    CASE DEFAULT
    WRITE(*,*)"wrong number of dimensions in shape_fun"
  END SELECT
  RETURN
  END SUBROUTINE shape_fun
  !
  SUBROUTINE shape_der(der,points,i)
  !
  !   This subroutine produces derivatives of shape functions withe respect
  !   to local coordinates.
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  INTEGER,INTENT(IN)::i
  REAL(iwp),INTENT(IN)::points(:,:)
  REAL(iwp),INTENT(OUT)::der(:,:)
  REAL(iwp)::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3
  REAL(iwp)::L1,L2,L3,L4 !coordenadas naturais para o tetraedro
  REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9,x2p1,x2m1,e2p1,e2m1,zetam,zetap
  REAL,PARAMETER::zero=0.0_iwp,pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,  &
  pt75=0.75_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,&
  d6=6.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d10=10.0_iwp,d11=11.0_iwp,            &
  d12=12.0_iwp,d16=16.0_iwp,d18=18.0_iwp,d27=27.0_iwp,d32=32.0_iwp,      &
  d36=36.0_iwp,d54=54.0_iwp,d64=64.0_iwp,d128=128.0_iwp
  INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
  ndim=UBOUND(der,1)
  nod= UBOUND(der,2)
  SELECT CASE(ndim)
  CASE(1)   ! one dimensional elements
    xi=points(i,1)
    SELECT CASE(nod)
    CASE(2)
      der(1,1)=-pt5
      der(1,2)= pt5
    CASE(3)
      t1=-one-xi
      t2=-xi
      t3=one-xi
      der(1,1)=-(t3+t2)/two
      der(1,2)=(t3+t1)
      der(1,3)=-(t2+t1)/two
    CASE(4)
      t1=-one-xi
      t2=-one/d3-xi
      t3=one/d3-xi
      t4=one-xi
      der(1,1)=-(t3*t4+t2*t4+t2*t3)*d9/d16
      der(1,2)=(t3*t4+t1*t4+t1*t3)*d27/d16
      der(1,3)=-(t2*t4+t1*t4+t1*t2)*d27/d16
      der(1,4)=(t2*t3+t1*t3+t1*t2)*d9/d16
    CASE(5)
      t1=-one-xi
      t2=-pt5-xi
      t3=-xi
      t4=pt5-xi
      t5=one-xi
      der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*two/d3
      der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*d8/d3
      der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*d4
      der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*d8/d3
      der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*two/d3
      CASE DEFAULT
      WRITE(*,*)"wrong number of nodes in shape_der"
    END SELECT
  CASE(2)      ! two dimensional elements
    xi=points(i,1)
    eta=points(i,2)
    c1=xi
    c2=eta
    c3=one-c1-c2
    etam=pt25*(one-eta)
    etap=pt25*(one+eta)
    xim= pt25*(one-xi)
    xip= pt25*(one+xi)
    x2p1=two*xi+one
    x2m1=two*xi-one
    e2p1=two*eta+one
    e2m1=two*eta-one
    SELECT CASE(nod)
    CASE(3)
      der(1,1)=one
      der(1,3)=zero
      der(1,2)=-one
      der(2,1)=zero
      der(2,3)=one
      der(2,2)=-one
    CASE(6)
      der(1,1)=d4*c1-one
      der(1,6)=d4*c2
      der(1,5)=zero
      der(1,4)=-d4*c2
      der(1,3)=-(d4*c3-one)
      der(1,2)=d4*(c3-c1)
      der(2,1)=zero
      der(2,6)=d4*c1
      der(2,5)=d4*c2-one
      der(2,4)=d4*(c3-c2)
      der(2,3)=-(d4*c3-one)
      der(2,2)=-d4*c1
    CASE(10)
      der(1,1)=(d27*c1**2-d18*c1+two)/two
      der(1,9)=(d9*(d6*c1-one)*c2)/two
      der(1,8)=(d9*(d3*c2-one)*c2)/two
      der(1,7)=zero
      der(1,6)=-(d9*(d3*c2-one)*c2)/two
      der(1,5)= (d9*(d6*c1+d6*c2-d5)*c2)/two
      der(1,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
      der(1,3)= (d9*(d9*c1**2+d12*c1*c2-d10*c1+d3*c2**2-d5*c2+two))/two
      der(1,2)=-(d9*(d9*c1**2+d6*c1*c2-d8*c1-c2+one))/two
      der(1,10)=-d27*(((c2-one)+c1)+c1)*c2
      der(2,1)=zero
      der(2,9)= (d9*(d3*c1-one)*c1)/two
      der(2,8)= (d9*(d6*c2-one)*c1)/two
      der(2,7)=(d27*c2**2-d18*c2+two)/two
      der(2,6)=-(d9*((c1+c2-one)*(d6*c2-one)+(d3*c2-one)*c2))/two
      der(2,5)= (d9*(d3*c1**2+d12*c1*c2-d5*c1+d9*c2**2-d10*c2+two))/two
      der(2,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
      der(2,3)= (d9*(d6*c1+d6*c2-d5)*c1)/two
      der(2,2)=-(d9*(d3*c1-one)*c1)/two
      der(2,10)=-d27*(((c2-one)+c1)+c2)*c1
    CASE(15)
      t1=c1-pt25
      t2=c1-pt5
      t3=c1-pt75
      t4=c2-pt25
      t5=c2-pt5
      t6=c2-pt75
      t7=c3-pt25
      t8=c3-pt5
      t9=c3-pt75
      der(1,1)=d32/d3*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
      der(1,12)=d128/d3*c2*(t2*(t1+c1)+c1*t1)
      der(1,11)=d64*c2*t4*(t1+c1)
      der(1,10)=d128/d3*c2*t4*t5
      der(1,9)=zero
      der(1,8)=-d128/d3*c2*t4*t5
      der(1,7)=-d64*c2*t4*(t7+c3)
      der(1,6)=-d128/d3*c2*(t8*(t7+c3)+c3*t7)
      der(1,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
      der(1,4)=d128/d3*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
      der(1,3)=d64*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
      der(1,2)=d128/d3*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
      der(1,13)=d128*c2*(c3*(t1+c1)-c1*t1)
      der(1,15)=d128*c2*t4*(c3-c1)
      der(1,14)=d128*c2*(c3*t7-c1*(t7+c3))
      der(2,1)=zero
      der(2,12)=d128/d3*c1*t1*t2
      der(2,11)=d64*c1*t1*(t4+c2)
      der(2,10)=d128/d3*c1*(t5*(t4+c2)+c2*t4)
      der(2,9)=d32/d3*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
      der(2,8)=d128/d3*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
      der(2,7)=d64*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
      der(2,6)=d128/d3*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
      der(2,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
      der(2,4)=-d128/d3*c1*(t8*(t7+c3)+c3*t7)
      der(2,3)=-d64*c1*t1*(t7+c3)
      der(2,2)=-d128/d3*c1*t1*t2
      der(2,13)=d128*c1*t1*(c3-c2)
      der(2,15)=d128*c1*(c3*(t4+c2)-c2*t4)
      der(2,14)=d128*c1*(c3*t7-c2*(c3+t7))
    CASE (4)
      der(1,1)=-etam
      der(1,2)=-etap
      der(1,3)=etap
      der(1,4)=etam
      der(2,1)=-xim
      der(2,2)=xim
      der(2,3)=xip
      der(2,4)=-xip
    CASE(5)
      der(1,1)=-etam+pt5*xi*(one-eta**2)
      der(1,2)=-etap+pt5*xi*(one-eta**2)
      der(1,3)=etap+pt5*xi*(one-eta**2)
      der(1,4)=etam+pt5*xi*(one-eta**2)
      der(1,5)=-two*xi*(one-eta**2)
      der(2,1)=-xim+pt5*eta*(one-xi**2)
      der(2,2)=xim+pt5*eta*(one-xi**2)
      der(2,3)=xip+pt5*eta*(one-xi**2)
      der(2,4)=-xip+pt5*eta*(one-xi**2)
      der(2,5)=-two*eta*(one-xi**2)
    CASE(8)
      der(1,1)=etam*(two*xi+eta)
      der(1,2)=-d8*etam*etap
      der(1,3)=etap*(two*xi-eta)
      der(1,4)=-d4*etap*xi
      der(1,5)=etap*(two*xi+eta)
      der(1,6)=d8*etap*etam
      der(1,7)=etam*(two*xi-eta)
      der(1,8)=-d4*etam*xi
      der(2,1)=xim*(xi+two*eta)
      der(2,2)=-d4*xim*eta
      der(2,3)=xim*(two*eta-xi)
      der(2,4)=d8*xim*xip
      der(2,5)=xip*(xi+two*eta)
      der(2,6)=-d4*xip*eta
      der(2,7)=xip*(two*eta-xi)
      der(2,8)=-d8*xim*xip
    CASE(9)
      etam=eta-one
      etap=eta+one
      xim=xi-one
      xip=xi+one
      der(1,1)=pt25*x2m1*eta*etam
      der(1,2)=-pt5*x2m1*etap*etam
      der(1,3)=pt25*x2m1*eta*etap
      der(1,4)=-xi*eta*etap
      der(1,5)=pt25*x2p1*eta*etap
      der(1,6)=-pt5*x2p1*etap*etam
      der(1,7)=pt25*x2p1*eta*etam
      der(1,8)=-xi*eta*etam
      der(1,9)=two*xi*etap*etam
      der(2,1)=pt25*xi*xim*e2m1
      der(2,2)=-xi*xim*eta
      der(2,3)=pt25*xi*xim*e2p1
      der(2,4)=-pt5*xip*xim*e2p1
      der(2,5)=pt25*xi*xip*e2p1
      der(2,6)=-xi*xip*eta
      der(2,7)=pt25*xi*xip*e2m1
      der(2,8)=-pt5*xip*xim*e2m1
      der(2,9)=two*xip*xim*eta
      CASE DEFAULT
      WRITE(*,*)"wrong number of nodes in shape_der"
    END SELECT
  CASE(3)  ! d3 dimensional elements
    xi=points(i,1)
    eta=points(i,2)
    zeta=points(i,3)
    etam=one-eta
    xim=one-xi
    zetam=one-zeta
    etap=eta+one
    xip=xi+one
    zetap=zeta+one
    SELECT CASE(nod)
    CASE(4)
      der(1:3,1:4)=zero
      der(1,1)=one
      der(2,2)=one
      der(3,3)=one
      der(1,4)=-one
      der(2,4)=-one
      der(3,4)=-one
    CASE(8)
      der(1,1)=-pt125*etam*zetam
      der(1,2)=-pt125*etam*zetap
      der(1,3)= pt125*etam*zetap
      der(1,4)= pt125*etam*zetam
      der(1,5)=-pt125*etap*zetam
      der(1,6)=-pt125*etap*zetap
      der(1,7)= pt125*etap*zetap
      der(1,8)= pt125*etap*zetam
      der(2,1)=-pt125*xim*zetam
      der(2,2)=-pt125*xim*zetap
      der(2,3)=-pt125*xip*zetap
      der(2,4)=-pt125*xip*zetam
      der(2,5)= pt125*xim*zetam
      der(2,6)= pt125*xim*zetap
      der(2,7)= pt125*xip*zetap
      der(2,8)= pt125*xip*zetam
      der(3,1)=-pt125*xim*etam
      der(3,2)= pt125*xim*etam
      der(3,3)= pt125*xip*etam
      der(3,4)=-pt125*xip*etam
      der(3,5)=-pt125*xim*etap
      der(3,6)= pt125*xim*etap
      der(3,7)= pt125*xip*etap
      der(3,8)=-pt125*xip*etap
    CASE(10)!ELEMENTO TETRA�DRICO QUADR�TICO--------------------------------
      L1=points(i,1);L2=points(i,2);L3=points(i,3);L4=one-L1-L2-L3
      der(1,1)=two*L1+(two*L1-one)
      der(1,2)=zero
      der(1,3)=zero
      der(1,4)=(-two)*(one-L1-L2-L3)+(two*(one-L1-L2-L3)-one)*(-one)
      der(1,5)=d4*L2
      der(1,6)=zero
      der(1,7)=d4*L3
      der(1,8)=d4*(one-L1-L2-L3)+d4*L1*(-one)
      der(1,9)=d4*L2*(-one)
      der(1,10)=d4*L3*(-one)
      der(2,1)=zero
      der(2,2)=two*L2+(two*L2-one)
      der(2,3)=zero
      der(2,4)=(-two)*(one-L1-L2-L3)+(two*(one-L1-L2-L3)-1)*(-one)
      der(2,5)=d4*L1
      der(2,6)=d4*L3
      der(2,7)=zero
      der(2,8)=d4*L1*(-one)
      der(2,9)=d4*(one-L1-L2-L3)+d4*L2*(-one)
      der(2,10)=d4*L3*(-one)
      der(3,1)=zero
      der(3,2)=zero
      der(3,3)=two*L3+(two*L3-one)
      der(3,4)=(-two)*(one-L1-L2-L3)+(two*(one-L1-L2-L3)-1)*(-one)
      der(3,5)=zero
      der(3,6)=d4*L2
      der(3,7)=d4*L1
      der(3,8)=d4*L1*(-one)
      der(3,9)=d4*L2*(-one)
      der(3,10)=d4*(one-L1-L2-L3)+d4*L3*(-one)
      !Zienkiewicz e Taylor vol.1 pp 188
      !-----------------------------------------------------------------------
    CASE(14) ! type 6 element
      der(1,1)= (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*          &
      (eta-one)*(zeta-one)/d8
      der(1,2)=-(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*          &
      (eta-one)*(zeta+one)/d8
      der(1,3)=-(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*          &
      (eta-one)*(zeta+one)/d8
      der(1,4)= (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*          &
      (eta-one)*(zeta-one)/d8
      der(1,5)= -(eta-one)*(zeta+one)*(zeta-one)*xi
      der(1,6)=-(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
      der(1,7)=  (eta+one)*(eta-one)*(zeta+one)*xi
      der(1,8)= (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
      der(1,9)= -(eta+one)*(eta-one)*(zeta-one)*xi
      der(1,10)= (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*         &
      (eta+one)*(zeta-one)/d8
      der(1,11)=-(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*         &
      (eta+one)*(zeta+one)/d8
      der(1,12)=-(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*         &
      (eta+one)*(zeta+one)/d8
      der(1,13)= (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*         &
      (eta+one)*(zeta-one)/d8
      der(1,14)=  (eta+one)*(zeta+one)*(zeta-one)*xi
      der(2,1)= (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*          &
      (xi-one)*(zeta-one)/d8
      der(2,2)=-(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*          &
      (xi-one)*(zeta+one)/d8
      der(2,3)=-(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*          &
      (xi+one)*(zeta+one)/d8
      der(2,4)= (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*          &
      (xi+one)*(zeta-one)/d8
      der(2,5)=-(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
      der(2,6)= -(xi-one)*(zeta+one)*(zeta-one)*eta
      der(2,7)=  (xi+one)*(xi-one)*(zeta+one)*eta
      der(2,8)=  (xi+one)*(zeta+one)*(zeta-one)*eta
      der(2,9)= -(xi+one)*(xi-one)*(zeta-one)*eta
      der(2,10)= (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*         &
      (xi-one)*(zeta-one)/d8
      der(2,11)=-(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*         &
      (xi-one)*(zeta+one)/d8
      der(2,12)=-(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*         &
      (xi+one)*(zeta+one)/d8
      der(2,13)= (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*         &
      (xi+one)*(zeta-one)/d8
      der(2,14)= (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
      der(3,1)= (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*          &
      (xi-one)*(eta-one)/d8
      der(3,2)=-(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*          &
      (xi-one)*(eta-one)/d8
      der(3,3)=-(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*          &
      (xi+one)*(eta-one)/d8
      der(3,4)= (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*          &
      (xi+one)*(eta-one)/d8
      der(3,5)= -(xi+one)*(xi-one)*(eta-one)*zeta
      der(3,6)= -(xi-one)*(eta+one)*(eta-one)*zeta
      der(3,7)= (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
      der(3,8)=  (xi+one)*(eta+one)*(eta-one)*zeta
      der(3,9)=-(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
      der(3,10)= (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*         &
      (xi-one)*(eta+one)/d8
      der(3,11)=-(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*         &
      (xi-one)*(eta+one)/d8
      der(3,12)=-(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*         &
      (xi+one)*(eta+one)/d8
      der(3,13)= (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*         &
      (xi+one)*(eta+one)/d8
      der(3,14)=  (xi+one)*(xi-one)*(eta+one)*zeta
    CASE(20)
      xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
      etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
      zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
      DO l=1,20
        xi0=xi*xii(l)
        eta0=eta*etai(l)
        zeta0=zeta*zetai(l)
        IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
          der(1,l)=-pt5*xi*(one+eta0)*(one+zeta0)
          der(2,l)=pt25*etai(l)*(one-xi*xi)*(one+zeta0)
          der(3,l)=pt25*zetai(l)*(one-xi*xi)*(one+eta0)
        ELSE IF(l>=9.AND.l<=12)THEN
          der(1,l)=pt25*xii(l)*(one-eta*eta)*(one+zeta0)
          der(2,l)=-pt5*eta*(one+xi0)*(one+zeta0)
          der(3,l)=pt25*zetai(l)*(one+xi0)*(one-eta*eta)
        ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
          der(1,l)=pt25*xii(l)*(one+eta0)*(one-zeta*zeta)
          der(2,l)=pt25*etai(l)*(one+xi0)*(one-zeta*zeta)
          der(3,l)=-pt5*zeta*(one+xi0)*(one+eta0)
        ELSE
          der(1,l)=pt125*xii(l)*(one+eta0)*(one+zeta0)*                    &
          (two*xi0+eta0+zeta0-one)
          der(2,l)=pt125*etai(l)*(one+xi0)*(one+zeta0)*                    &
          (xi0+two*eta0+zeta0-one)
          der(3,l)=pt125*zetai(l)*(one+xi0)*(one+eta0)*                    &
          (xi0+eta0+two*zeta0-one)
        END IF
      END DO
      CASE DEFAULT
      WRITE(*,*)"wrong number of nodes in shape_der"
    END SELECT
    CASE DEFAULT
    WRITE(*,*)"wrong number of dimensions in shape_der"
  END SELECT
  RETURN
  END SUBROUTINE shape_der
  !

  SUBROUTINE checon(loads,oldlds,tol,converged,convergence)
  !
  ! This subroutine sets converged to .FALSE. if relative change in loads
  ! and oldlds is greater than tol and updates oldlds.
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp)::temp1,temp2
  REAL(iwp),INTENT(IN)::loads(0:),tol
  REAL(iwp),INTENT(IN OUT)::convergence
  REAL(iwp),INTENT(IN OUT)::oldlds(0:)
  LOGICAL,INTENT(OUT)::converged
  CONVERGED=.TRUE.
  temp1=MAXVAL(ABS(loads-oldlds))
  temp2=MAXVAL(ABS(loads))
  convergence=temp1/temp2
  CONVERGED=(convergence<=tol)
  oldlds=loads
  RETURN
  END SUBROUTINE checon

  !rotina de Smith e Griffiths para c�lculo das derivadas dos invariantes de tens�o em rela��o �s tens�es.
  SUBROUTINE formm(stress,m1,m2,m3)
  !
  ! This subroutine forms the derivatives of the invariants with respect to
  ! stress in 2- or 3-d.
  !
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp),INTENT(IN)::stress(:)
  REAL(iwp),INTENT(OUT)::m1(:,:),m2(:,:),m3(:,:)
  REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
  two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
  INTEGER::nst,i,j
  nst=UBOUND(stress,1)
  SELECT CASE(nst)
  CASE(4)
    sx=stress(1)
    sy=stress(2)
    txy=stress(3)
    sz=stress(4)
    dx=(two*sx-sy-sz)/three
    dy=(two*sy-sz-sx)/three
    dz=(two*sz-sx-sy)/three
    sigm=(sx+sy+sz)/three
    m1=zero
    m2=zero
    m3=zero
    m1(1,1:2)=one
    m1(2,1:2)=one
    m1(4,1:2)=one
    m1(1,4)=one
    m1(4,4)=one
    m1(2,4)=one
    m1=m1/nine/sigm
    m2(1,1)=two/three
    m2(2,2)=two/three
    m2(4,4)= two/three
    m2(2,4)=-one/three
    m2(4,2)=-one/three
    m2(1,2)=-one/three
    m2(2,1)=-one/three
    m2(1,4)=-one/three
    m2(4,1)=-one/three
    m2(3,3)=two
    m3(3,3)=-dz
    m3(1:2,3)=txy/three
    m3(3,1:2)=txy/three
    m3(3,4)=-two*txy/three
    m3(4,3)=-two*txy/three
    m3(1,1)=dx/three
    m3(2,4)=dx/three
    m3(4,2)=dx/three
    m3(2,2)=dy/three
    m3(1,4)=dy/three
    m3(4,1)=dy/three
    m3(4,4)=dz/three
    m3(1,2)=dz/three
    m3(2,1)=dz/three
  CASE(6)
    sx=stress(1)
    sy=stress(2)
    sz=stress(3)
    txy=stress(4)
    tyz=stress(5)
    tzx=stress(6)
    sigm=(sx+sy+sz)/three
    dx=sx-sigm
    dy=sy-sigm
    dz=sz-sigm
    m1=zero
    m2=zero
    m1(1:3,1:3)=one/(three*sigm)
    DO i=1,3
      m2(i,i)=two
      m2(i+3,i+3)=six
    END DO
    m2(1,2)=-one
    m2(1,3)=-one
    m2(2,3)=-one
    m3(1,1)=dx
    m3(1,2)=dz
    m3(1,3)=dy
    m3(1,4)=txy
    m3(1,5)=-two*tyz
    m3(1,6)=tzx
    m3(2,2)=dy
    m3(2,3)=dx
    m3(2,4)=txy
    m3(2,5)=tyz
    m3(2,6)=-two*tzx
    m3(3,3)=dz
    m3(3,4)=-two*txy
    m3(3,5)=tyz
    m3(3,6)=tzx
    m3(4,4)=-three*dz
    m3(4,5)=three*tzx
    m3(4,6)=three*tyz
    m3(5,5)=-three*dx
    m3(5,6)=three*txy
    m3(6,6)=-three*dy
    DO i=1,6
      DO j=i+1,6
        m1(j,i)=m1(i,j)
        m2(j,i)=m2(i,j)
        m3(j,i)=m3(i,j)
      END DO
    END DO
    m1=m1/three
    m2=m2/three
    m3=m3/three
    CASE DEFAULT
    WRITE(*,*)"nst size not recognised in formm"
  END SELECT
  RETURN
  END SUBROUTINE formm
  !

  END MODULE
