
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE GravityLoad
  
  
  CONTAINS
  
       
    SUBROUTINE GravityLoadCalculation(GeomechanicalModel,Load)

        USE FEMlib
        USE NumericalLib
        USE DiscretizedModel
        USE NumericalUtilities

        IMPLICIT NONE

        TYPE(Model),INTENT(IN):: GeomechanicalModel
        REAL(8),INTENT(IN OUT)::Load(0:)

        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER:: iel,i,j  !contadores.

        INTEGER:: ndim,&  !dimensões do problema
                  nst,&  !número de componentes no tensor de tensões
                  nn,&  !número de nós
                  nels,&  !número de elementos
                  nip,&  !número de pontos de integração
                  nod,&  !número de nós por elemento.
                  ndof
        REAL(iwp):: zero=0.0_iwp,&
                    det,&  !MatrixDeterminant do jacobiano
                    biot,&  !coeficiente de Biot.
                    NTG  !net-to-gross.
        CHARACTER(LEN=320)::element
        REAL(iwp):: unitWeight, grav=9.81

        !VARIÁVEIS DINÂMICAS
        INTEGER,ALLOCATABLE:: num(:) ,&  !conectividade de um elemento.
                              g(:)  !vetor com os graus de liberdade de cada elemento.
        REAL(iwp),ALLOCATABLE:: points(:,:),&  !coordenadas dos pontos de integração
                                weights(:),&  !peso dos pontos de integração
                                coord(:,:),&  !coordenada dos nós de um elemento
                                fun(:),&  !funções de forma
                                der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                                deriv(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas globais
                                jac(:,:),&  !matriz jacobiana das funções de forma
                                g_coord(:,:),&  !coordenada dos nós da malha
                                bee(:,:),&  !matriz de derivadas das funções de forma,
                                elemNodalForce(:)

        element = GeomechanicalModel%element
        nels = GeomechanicalModel%NumElements
        nn = GeomechanicalModel%NumNodes
        ndim = GeomechanicalModel%NumDimensions
        nst = GeomechanicalModel%NumTensorComponents
        nip = 8
        nod = GeomechanicalModel%NumNodesPerElement
        ndof = GeomechanicalModel%NumDOFUperElement

        !MEMORY ALLOCATION.
        ALLOCATE(num(nod),g(ndof),points(nip,ndim),weights(nip),bee(nst,ndof),&
                  coord(nod,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim),&
                  g_coord(ndim,nn),elemNodalForce(ndof))

        !NODAL COORDENATES.
        g_coord=zero
        g_coord(1,:) = GeomechanicalModel%X
        g_coord(2,:) = GeomechanicalModel%Y
        g_coord(3,:) = GeomechanicalModel%Z

        CALL sample(element,points,weights)
        
        DO iel=1,nels

          IF( GeomechanicalModel%ActiveElements(iel) .EQV. .TRUE. ) THEN

            num = GeomechanicalModel%ElementNodes(:,iel)  !g_num(:,iel)
            coord = TRANSPOSE(g_coord(:,num))
            g = GeomechanicalModel%ElementDOFs(:,iel)     !g_g(:,iel)

            !Carga gravitacional.
            unitWeight = grav * GeomechanicalModel%Density(iel)
            CALL ElementGravityLoad(nip,points,weights,coord,unitWeight,elemNodalForce)
            Load(g) = Load(g) + elemNodalForce
            
          END IF

        END DO

        RETURN
    END SUBROUTINE
    
    SUBROUTINE ElementGravityLoad(nip,points,weights,coord,unitWeight,elemNodalForce)
      !Subrotina para calcular a força gravitacional.
      !       elemNodalForce(# nós) = unitWeight * int( transpose(N)* dv )
      
      !USE Main_Constants
      !USE iNumericalUtilities
      !USE iMaterials
      USE FEMlib
      USE NumericalLib
      USE NumericalUtilities
      
      IMPLICIT NONE
      
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      
      !VARIÁVEIS INPUT.
      
      INTEGER,INTENT(IN)::nip  !número de pontos de integração
      REAL(iwp),INTENT(IN)::unitWeight  !peso específico (força / volume -> kN/m³)
      REAL(iwp),ALLOCATABLE,INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integração
                                        weights(:),&  !peso dos pontos de integração
                                        coord(:,:)  !coordenada dos nós de um elemento
      REAL(iwp),ALLOCATABLE,INTENT(OUT)::elemNodalForce(:)  !matriz de rigidez do elemento
      
      !VARIÁVEIS ESTÁTICAS.
      
      INTEGER::i  !contadores
      REAL(iwp)::det  !determinante do jacobiano
      
      !VARIÁVEIS DINÂMICAS.
      
      REAL(iwp),ALLOCATABLE::bee(:,:),&  !matriz de derivadas das funções de forma
                             jac(:,:),&  !matriz jacobiana das funções de forma
                             fun(:),&  !funções de forma
                             der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                             deriv(:,:)  !derivadas das funções de forma dos deslocamentos em coordenadas globais

      INTEGER:: ndim = 3
      INTEGER:: nod = 8                       !number of nodes per element.
      INTEGER:: ndof = 24                     !number of displacement dofs per element.
      INTEGER:: nodf = 8                      !number of fem dofs per element.
      INTEGER:: nodof = 3                     !number of displacements per node.
      INTEGER:: nst = 6                       !number of stress tensor.
      
      !ALOCAÇÃO DE MEMÓRIA.
      
      ALLOCATE(bee(nst,ndof),jac(ndim,ndim),fun(nod),der(ndim,nod),deriv(ndim,nod),elemNodalForce(ndof))
      fun=0.0_iwp
      
      !INTEGRAÇÃO NUMÉRICA DA CARGA DE PESO
      !     int( transpose(N)* dv ) = 
      !                           = SOMA(pontos de Gauss, unitWeight * transpose(N) *  detJacobiano * pesoPtoGauss
      
      elemNodalForce = 0.0_iwp
      
      DO i=1,nip
        
        CALL shape_fun(fun,points,i)  !funções de forma no ponto de integração i
        CALL shape_der(der,points,i)  !derivada das funções de forma
        jac = MATMUL(der,coord)       !jacobiano do elemento no ponto de integração i
        CALL determinante(jac,det)    !determinante do jacobiano
        
        elemNodalForce(nodof:ndof:nodof) = elemNodalForce(nodof:ndof:nodof) + &
                                           fun(:) * det * weights(i)
        
      END DO
      
      elemNodalForce = - unitWeight * elemNodalForce
      !(negativo: Z é para cima e a carga gravitacional é para baixo)

    RETURN
    END SUBROUTINE ElementGravityLoad

    
    !SUBROUTINE ElementGravityLoad(nip,points,weights,coord,unitWeight,elemNodalForce)
    !  !Subrotina para calcular a força gravitacional.
    !  !       elemNodalForce(# nós) = unitWeight * int( transpose(N)* dv )
    !  
    !  !USE NumericalLib
    !  USE NumericalUtilities
    !  
    !  IMPLICIT NONE
    !  
    !  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    !  
    !  !VARIÁVEIS INPUT.
    !  INTEGER,INTENT(IN):: nip
    !  REAL(iwp),INTENT(IN)::unitWeight  !peso específico (força / volume -> kN/m³)
    !  REAL(iwp),ALLOCATABLE,INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integração
    !                                    weights(:),&  !peso dos pontos de integração
    !                                    coord(:,:)  !coordenada dos nós de um elemento
    !  REAL(iwp),ALLOCATABLE,INTENT(OUT)::elemNodalForce(:)  !matriz de rigidez do elemento
    !  
    !  !VARIÁVEIS ESTÁTICAS.
    !  
    !  INTEGER::i  !contadores
    !  REAL(iwp)::det  !determinante do jacobiano
    !  
    !  !VARIÁVEIS DINÂMICAS.
    !  
    !  REAL(iwp),ALLOCATABLE::bee(:,:),&  !matriz de derivadas das funções de forma
    !                         jac(:,:),&  !matriz jacobiana das funções de forma
    !                         fun(:),&  !funções de forma
    !                         der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
    !                         deriv(:,:)  !derivadas das funções de forma dos deslocamentos em coordenadas globais
    !    INTEGER:: ndim = 3,&  !dimensões do problema
    !              nst = 6,&   !número de componentes no tensor de tensões
    !              nod = 8,&   !número de nós por elemento.
    !              nodof=3,&   !number of dofs per node.
    !              ndof = 24   !number of element  dofs.
    !    
    !  !ALOCAÇÃO DE MEMÓRIA.
    !  
    !  ALLOCATE(bee(nst,ndof),jac(ndim,ndim),fun(nod),der(ndim,nod),deriv(ndim,nod),elemNodalForce(ndof))
    !  fun=0.0_iwp
    !  
    !  !INTEGRAÇÃO NUMÉRICA DA CARGA DE PESO
    !  !     int( transpose(N)* dv ) = 
    !  !                           = SOMA(pontos de Gauss, unitWeight * transpose(N) *  detJacobiano * pesoPtoGauss
    !  
    !  elemNodalForce = 0.0_iwp
    !  
    !  DO i=1,nip
    !    
    !    CALL shape_fun(fun,points,i)  !funções de forma no ponto de integração i
    !    CALL shape_der(der,points,i)  !derivada das funções de forma
    !    jac = MATMUL(der,coord)       !jacobiano do elemento no ponto de integração i
    !    CALL determinante(jac,det)    !determinante do jacobiano
    !    
    !    elemNodalForce(nodof:ndof:nodof) = elemNodalForce(nodof:ndof:nodof) + fun(:) * det * weights(i)
    !    
    !  END DO
    !  
    !  elemNodalForce = - unitWeight * elemNodalForce
    !  !(negativo: Z é para cima e a carga gravitacional é para baixo)
    !  
    !
    !RETURN
    !END SUBROUTINE ElementGravityLoad

  
  END MODULE