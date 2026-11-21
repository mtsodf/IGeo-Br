
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE PorepressureLoad


  CONTAINS


    SUBROUTINE PoroelasticLoad(GeomechanicalModel,PorePressure,Load)

      USE FEMlib
      USE NumericalLib
      USE DiscretizedModel
      USE NumericalUtilities

      IMPLICIT NONE

      TYPE(Model),INTENT(IN):: GeomechanicalModel     
      REAL(8),INTENT(IN):: PorePressure(:)            !Nodal Porepressure.
      REAL(8),INTENT(IN OUT)::Load(0:)

      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER:: iel,i,j  !contadores.

      INTEGER:: ndim,&  !dimensões do problema
                nst,&  !número de componentes no tensor de tensões
                nn,&  !número de nós
                nels,&  !número de elementos
                nip,&  !número de pontos de integração
                nod,&  !número de nós por elemento
                ndof,&  !número de graus de liberdade de deslocamento por elemento
                nodf,&  !número de nós com grau de liberdade de poro-pressão.
                nthreads  !número de threads ou núcleos de processamento que serão usados
      REAL(iwp):: zero=0.0_iwp,&
                  det,&  !MatrixDeterminant do jacobiano
                  biot,&  !coeficiente de Biot.
                  NTG  !net-to-gross.
      CHARACTER(LEN=320)::element

      !VARIÁVEIS DINÂMICAS
      INTEGER,ALLOCATABLE::num(:) ,&  !conectividade de um elemento.
          g(:)  !vetor com os graus de liberdade de cada elemento.
      REAL(iwp),ALLOCATABLE:: points(:,:),&  !coordenadas dos pontos de integração
                              weights(:),&  !peso dos pontos de integração
                              coord(:,:),&  !coordenada dos nós de um elemento
                              fun(:),&  !funções de forma
                              der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                              deriv(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas globais
                              jac(:,:),&  !matriz jacobiana das funções de forma
                              g_coord(:,:),&  !coordenada dos nós da malha
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
                              !propVetor(:,:),&  !vetor com as propriedades armazenadas (nPropVetor,nelsProp)
                              !propVetorPres(:,:)  !vetor para armazenamento da pressão de fluido (nels,prop)

      element = GeomechanicalModel%element
      nels = GeomechanicalModel%NumElements
      nn = GeomechanicalModel%NumNodes
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = 8
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement
      nthreads = GeomechanicalModel%NumPartitions
      nodf = 8  !porepressure dofs per element for integration by f.e..

      !MEMORY ALLOCATION.
      ALLOCATE(num(nod),g(ndof),points(nip,ndim),weights(nip),bee(nst,ndof),dee(nst,nst),&
                coord(nod,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim),&
                g_coord(ndim,nn),porolo(nodf),pp(nodf),coupl(ndof,nodf),&
                funf(nodf),derivf(ndim,nodf),derf(ndim,nodf),volf(ndof,nodf),vol(ndof))

      !NODAL COORDENATES.
      g_coord=zero
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z

      CALL sample(element,points,weights)
        
      DO iel=1,nels

        IF( GeomechanicalModel%ActiveElements(iel) .AND. GeomechanicalModel%ActiveToFlowElements(iel) ) THEN
          
          porolo=zero
          coupl=zero
          num = GeomechanicalModel%ElementNodes(:,iel)  !g_num(:,iel)
          coord = TRANSPOSE(g_coord(:,num))
          g = GeomechanicalModel%ElementDOFs(:,iel)     !g_g(:,iel)

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
          pp = PorePressure(num)
          biot = GeomechanicalModel%BiotCoef(iel)
          Load(g) = Load(g) + biot * MATMUL(coupl(:,:),pp) * NTG

        END IF
          
      END DO
        
      !Load(0) = zero

        
      RETURN
    END SUBROUTINE PoroelasticLoad


  END MODULE
