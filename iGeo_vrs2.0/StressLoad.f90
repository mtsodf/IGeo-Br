
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE StressLoad
  
  CONTAINS
  
    !element -> integral(Bt,DeltaSigma)
    SUBROUTINE ElementStressForce(nip,points,weights,coord,deltaStress,unbalForce)
      !Subrotina para calcular a força desbalanceada devido a diferença entre tensões.
      !       fumbal(# nós) = int( transpose(B)* (sigma*-sigma) * dv

      USE MainConstants
      USE NumericalUtilities

      IMPLICIT NONE

      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

      !INPUTS / OUTPUTS.
      INTEGER,INTENT(IN)::nip  !número de pontos de integração
      REAL(iwp),INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integração
                            weights(:),&  !peso dos pontos de integração
                            coord(:,:)  !coordenada dos nós de um elemento
      REAL(iwp),INTENT(IN OUT)::unbalForce(:)  !matriz de rigidez do elemento
      REAL(iwp),INTENT(IN)::deltaStress(:,:)  !delta de tensão (tensão desbalanceada)

      !VARIABLES
      INTEGER::i  !contadores
      REAL(iwp)::det  !determinante do jacobiano

      REAL(iwp),ALLOCATABLE:: bee(:,:),&  !matriz de derivadas das funções de forma
                              jac(:,:),&  !matriz jacobiana das funções de forma
                              fun(:),&  !funções de forma
                              der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                              deriv(:,:)  !derivadas das funções de forma dos deslocamentos em coordenadas globais
      
      ALLOCATE(bee(nst,ndof),jac(ndim,ndim),fun(nod),der(ndim,nod),deriv(ndim,nod))
      fun=0.0_iwp

      !STRESS INTEGRATION.
      !     int( transpose(B)* (sigma*-sigma) * dv ) =
      !                           = SUM(pontos de Gauss, transpose(B) * (sigma*-sigma) * detJacobiano * pesoPtoGauss
      unbalForce = 0.0_iwp
      DO i=1,nip

          CALL shape_fun(fun,points,i)  !funções de forma no ponto de integração i
          CALL shape_der(der,points,i)  !derivada das funções de forma
          jac = MATMUL(der,coord)       !jacobiano do elemento no ponto de integração i
          CALL determinante(jac,det)    !determinante do jacobiano
          CALL invert(jac)              !inverso do jacobiano
          deriv = MATMUL(jac,der)       !derivadas das funções de forma em coordenadas reais
          CALL beemat(bee,deriv)        !matriz de derivadas B

          unbalForce = unbalForce + MATMUL( TRANSPOSE( bee ) , deltaStress(:,i) ) * det * weights(i)

      END DO

    RETURN
    END SUBROUTINE ElementStressForce

    
    SUBROUTINE ElementPoroelasticForce(nip,points,weights,coord,BiotCoef,PorePressure,ElementLoad)

      USE FEMlib
      USE NumericalLib
      USE NumericalUtilities
      USE MainConstants

      IMPLICIT NONE
      INTEGER,INTENT(IN):: nip
      REAL(8),INTENT(IN):: points(:,:),&
                           weights(:),& 
                           coord(:,:),&
                           BiotCoef
      REAL(8),INTENT(IN):: PorePressure(:)
      REAL(8),INTENT(IN OUT):: ElementLoad(:)
      
      INTEGER:: iel,i,j

      REAL(8):: det,&  !MatrixDeterminant do jacobiano
                biot,&  !coeficiente de Biot.
                NTG  !net-to-gross.

      REAL(8),ALLOCATABLE:: fun(:),&  !funções de forma
                            der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                            deriv(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas globais
                            jac(:,:),&  !matriz jacobiana das funções de forma
                            bee(:,:),&  !matriz de derivadas das funções de forma
                            dee(:,:),&  !matriz de relações elásticas
                            pp(:),&  !poro-pressões nodais em um elemenro
                            porolo(:),&  !carga de poro-pressão em um elemento
                            coupl(:,:),&  !matriz de acoplamento poro-elástica de um elemento
                            funf(:),&  !funções de forma de poro-pressão
                            derivf(:,:),&  !derivada das funções de forma de funf em coordenadas locais
                            derf(:,:),& !derivada das funções de forma de funf em coordenadas globais
                            volf(:,:),&  !relacionado às deformações volumétricas
                            vol(:)   !,&  !derivadas correpondentes às deformações volumétricas

      !MEMORY ALLOCATION.
      ALLOCATE( bee(nst,ndof),dee(nst,nst),&
                der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim),&
                porolo(nodf),pp(nodf),coupl(ndof,nodf),&
                funf(nodf),derivf(ndim,nodf),derf(ndim,nodf),volf(ndof,nodf),vol(ndof))
        
      porolo=0.0
      coupl=0.0

      !Numerical integration.
      DO i=1,nip
        !Displacement interpolation.
        CALL shape_function(fun,points,i)
        CALL shape_functionDerivative(der,points,i)
        jac=MATMUL(der,coord)
        CALL MatrixDeterminant(jac,det)
        CALL InvertMatrix(jac)
        deriv=MATMUL(jac,der)
        CALL beemat(bee,deriv)
        vol(:)=bee(1,:)+bee(2,:)+bee(3,:)
        !porepressure interpolation.
        CALL shape_function(funf,points,i)
        CALL shape_functionDerivative(derf,points,i)
        derivf=MATMUL(jac,derf)
        CALL cross_product(vol,funf,volf)
        !coupling matrix u-p.
        coupl(:,:)=coupl(:,:)+volf*det*weights(i)
      END DO

      !Adding element nodal load to global array.
      NTG = 1.0
      pp = PorePressure
      biot = BiotCoef
      ElementLoad = ElementLoad + biot * MATMUL(coupl(:,:),pp) * NTG
        
      RETURN
    END SUBROUTINE
        
  END MODULE