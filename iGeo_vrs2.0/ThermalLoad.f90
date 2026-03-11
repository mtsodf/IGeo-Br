
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE ThermalLoad
  
  CONTAINS
  
    SUBROUTINE ThermalStrainLoad(GeomechanicalModel,Temperature,Load)

      USE FEMlib
      USE NumericalLib
      USE DiscretizedModel
      USE NumericalUtilities
      USE ConstitutiveGeoModels

      IMPLICIT NONE

      TYPE(Model),INTENT(IN):: GeomechanicalModel     
      REAL(8),INTENT(IN):: Temperature(:)            !Nodal Porepressure.
      REAL(8),INTENT(IN OUT):: Load(0:)
    
      !VARIÁVEIS ESTÁTICAS
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER::iel,i,k,&  !contadores
               nlen  !comprimento de palavras
      INTEGER:: ndim,&  !dimensões do problema
                nst,&  !número de componentes no tensor de tensões
                nn,&  !número de nós
                nels,&  !número de elementos
                nip,&  !número de pontos de integração
                nod,&  !número de nós por elemento
                ndof,&  !número de graus de liberdade de deslocamento por elemento
                iy  !estágio de carregamento
      REAL(iwp)::zero=0.0_iwp,one=1.0_iwp,&
                 det,&  !determinantee do jacobiano
                 prop(23),&  !vetor para leitura de propriedades do elemento do arquivo .prp
                 DeltaT,&  !variação de temperatura
                 propDin(8),&  !propriedades dinâmicas
                 Temp0,TempF,&  !temperaturas inicial e final
                 NTG,&  !net-to-gross
                 alfa,&  !coeficiente de dilatação térmica linear
                 e,v,&
                 temp
  
      !VARIÁVEIS DINÂMICAS
  	  INTEGER,ALLOCATABLE::num(:) ,&  !conectividade de um elemento
  	                       g(:)  !vetor com os graus de liberdade de cada elemento
  	  INTEGER,ALLOCATABLE::g_num(:,:),&  !conectividade dos elementos
  	                       g_g(:,:)  !matriz para identificação dos graus de liberdades dos nós
      REAL(iwp),ALLOCATABLE::points(:,:),&  !coordenadas dos pontos de integração
                             weights(:),&  !peso dos pontos de integração
                             coord(:,:),&  !coordenada dos nós de um elemento
                             fun(:),&  !funções de forma
                             der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                             deriv(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas globais
                             jac(:,:),&  !matriz jacobiana das funções de forma
                             g_coord(:,:),&  !coordenada dos nós da malha
                             bee(:,:),&  !matriz de derivadas das funções de forma
                             dee(:,:),&  !matriz de relações elásticas
                             eterm(:),&  !deformação térmica
                             termlo(:)  !carga térmica em um elemento
      CHARACTER(LEN=320)::element
    
      element = GeomechanicalModel%element
      nels = GeomechanicalModel%NumElements
      nn = GeomechanicalModel%NumNodes
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = 8
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement
      
      !ALOCAÇÃO DE MEMÓRIA
      ALLOCATE(num(nod),g(ndof),points(nip,ndim),weights(nip),bee(nst,ndof),dee(nst,nst),&
               coord(nod,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim),&
               g_coord(ndim,nn),termlo(ndof),eterm(nst))
  
      !NODAL COORDENATES.
      g_coord=zero
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
    
      !LOOPING NOS ELEMENTOS PARA CALCULAR VETOR DE CARGAS ESTRUTURAIS
      CALL sample(element,points,weights)
      
      !calculo do acréscimo de temperatura
      DO iel=1,nels
          
        num = GeomechanicalModel%ElementNodes(:,iel)  !g_num(:,iel)
          
        DeltaT = SUM(Temperature(num))/8.0
        
        !calculo carga nodal de poropressão se houver excesso de poro-pressão
      	IF(ABS(DeltaT).GT.(0.01))THEN  !só entra no looping se a variação de temperatura for maior que 0.01oC
      	                                !essa tolerância evita erros de arredondamento em precisão dupla.
                        
          coord = TRANSPOSE(g_coord(:,num))
          g = GeomechanicalModel%ElementDOFs(:,iel)     !g_g(:,iel)
          
      	  !cálculo da matriz elástica [D]
          !e = GeomechanicalModel%YoungMod(iel) * (1.0-GeomechanicalModel%DamageStiffness(iel))
          !v = GeomechanicalModel%Poisson(iel)
          e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
          e = e * (1.0-GeomechanicalModel%DamageStiffness(iel))
          v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
      	  CALL D_elastic_homogenious(dee,e,v)
      	    
      	  alfa=GeomechanicalModel%ThermalExpCoeff(iel)

      	  eterm=zero
          eterm=(/one,one,one,zero,zero,zero/)
          eterm=eterm*alfa*DeltaT  !<<<<<<<falta inserir proporcionalidade ao NTG
          termlo=zero
            
          !calculando o vetor de cargas térmicas do element
          DO i=1,nip  !looping nos pontos de Gauss para integração
            CALL shape_fun(fun,points,i)
            CALL shape_der(der,points,i)
            jac=MATMUL(der,coord)
            CALL determinante(jac,det)
            CALL invert(jac)
            deriv=MATMUL(jac,der)
            CALL beemat(bee,deriv)
            termlo=termlo+MATMUL(MATMUL(TRANSPOSE(bee),dee),eterm)*det*weights(i)  !load=load+BT*D*defTerm*det(J)*peso   / tensões compressivas são negativas
          END DO
            
          !acrescentando o vetor de carga térmica do elemento ao vetor de carga global
          Load(g) = Load(g) + termlo
              
        END IF  !(DeltaT.NE.zero)
          
      END DO  !elementos)
      Load(0)=zero

    RETURN
    END SUBROUTINE ThermalStrainLoad
  
  END MODULE
  