  !*****************************************************************************
  !Geomechanics.
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduaçãoo em Geotecnia
  !Departamento de Engenharia Civil e Ambiental
  !Faculdade de Tecnologia
  !Universidade de Brasília
  !*****************************************************************************
  
  !******************************************
  !STIFFNESS MATRIX BUILD
  !******************************************
  
  MODULE StiffnessMatrix
  
    TYPE ElasticStiffnessMatrix
      
      INTEGER:: NumEquations
      
      !n�mero de elementos n�o nulos na matriz de rigidez.
      INTEGER:: NumValues !nnz
      
      !Vetores do armazenamento CSR.
      REAL(8),ALLOCATABLE:: a(:)
      INTEGER,ALLOCATABLE:: ia(:), ja(:)  !Ponteiros que armazenam as posi��es dos dofs na matriz.
      
    END TYPE
    
    TYPE MassMatrix
      REAL(8),ALLOCATABLE:: Values(:)    !(Number of DOFs)
    END TYPE

    TYPE Damping
      REAL(8),ALLOCATABLE:: Values(:)    !(Number of DOFs)
    END TYPE
    
    TYPE SolverMatricesContainer
      LOGICAL:: Built
      TYPE(ElasticStiffnessMatrix):: Stiff
      TYPE(MassMatrix):: Mass
      TYPE(Damping):: Damp
    END TYPE
  
  CONTAINS 
  
    !Rotina para adicionar elements � matriz de rigidez.
    SUBROUTINE StiffnessConstruction( GeomechanicalModel, a, ia, ja, ElementsToAdd )
      
      !USE IFPORT
      USE DiscretizedModel
      USE NumericalUtilities
      USE FEMlib
      USE ConstitutiveGeoModels
      USE IO
      USE OMP_LIB
      
      IMPLICIT NONE
      
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),ALLOCATABLE:: a(:)
      INTEGER,ALLOCATABLE:: ia(:), ja(:)  !Ponteiros que armazenam as posi��es dos dofs na matriz.
      LOGICAL,INTENT(IN):: ElementsToAdd(:)
      
      !VARI�VEIS EST�TICAS.
      
      INTEGER:: i,j,k,l,iel,&  !contadores.
                posicao,coluna,&  !aux�lio a localiza��o das posi��es no vetor "a".
                nPropVetor,&  !n�mero de propriedades gravadas por elemento.
                no,&  !auxiliar para guardar n�mero do n�.
                GL1,GL2,GL3  !marcadores de grau de liberdade.
      
      INTEGER:: nels,&  !n�mero de elementos.
                nn,&  !n�mero de n�s.
                ndim,&  !n�mero de dimens�es do problema.
                nst,&  !n�mero de componentes no tensor de tens�es.
                nip,&  !n�mero de pontos de integra��o.
                nod,&  !n�mero de n�s por elemento.
                ndof !n�mero de graus de liberdade por elemento.
      
      INTEGER:: nnz  !n�mero de elementos n�o nulos na matriz de rigidez (diagonal + matriz traingular superior).
      INTEGER:: neq  !n�mero de equa��es.

      INTEGER:: nthreads  !n�mero de threads ou n�cleos de processamento que ser�o usados.
      
      REAL(8):: young, poisson
      REAL(8):: Pc,M,kappa,lambda,e0,P0,defVol,evol_vp
      
      INTEGER:: position, matrixPosition
      REAL(8)::zero=0.0  !vari�veis do processo do gradientes conjugados.
      REAL(8)::startTime,endTime  !vari�veis para marca��o de tempo gasto em segundos.
      CHARACTER(LEN=320)::element  !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
      LOGICAL:: propriedadesPvetor(23)  !vetor bin�rio com as propriedades que devem ser gravadas.
      
      INTEGER,ALLOCATABLE:: g(:),&  !vetor com graus de liberdade de um elemento.
                            num(:) ,&  !conectividade de um elemento.
                            degrees(:,:),&  !associa��o graus de liberdade e elementos.
                            posicoes(:),conection(:,:),conectionOrdenado(:),&  !vetores auxiliares para constru��o de conec��o entre graus de liberdade.
                            propVetorEls(:)  !ponteiro os elementos com propriedades armazenadas (em ordem crescente de numera��o).
      INTEGER(4),ALLOCATABLE:: listaAfinidades(:,:)  !lista de afinidades nodais (8*3,nn), 1a posi��o � um ponteiro.

      INTEGER,POINTER:: nf(:,:)  !matriz de graus de liberdade dos n�s.
                        !g_g(:,:),&  !vetor com os graus de liberdade dos n�s de cada elemento.
                        !g_num(:,:),&  !conectividade dos elementos.
      REAL(8),ALLOCATABLE::   km(:,:),&  !matriz de rigidez do elemento,
                              km_fluid(:,:),&
                              g_coord(:,:),&  !coordenadas dos n�s da malha.
                              coord(:,:),&  !coordenada dos n�s de um elemento.
                              points(:,:),&  !coordenadas dos pontos de integra��o.
                              weights(:),&  !peso dos pontos de integra��o.
                              anisotropia_Ec(:),&  !vetor de anisotropia Ec (=Ec/Ea), onde Ea est� no arquivo .prp.
                              anisotropia_niac(:)  !vetor de anisotropia niab (=niab/niac), onde niac est� no arquivo .prp.
      INTEGER:: globalPosition
      REAL(8):: value
      
      REAL(8):: Kf
                
      CHARACTER(LEN=8)::HOUR
      
      !INICIALIZA��O DE VARI�VEIS.
      neq = GeomechanicalModel%NumEquations
      element = GeomechanicalModel%element
      nels = GeomechanicalModel%NumElements
      nn = GeomechanicalModel%NumNodes
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = 8 !GeomechanicalModel%NumIntegrationPoints
              !Os pontos de integra��o para c�lculo de tens�es e deforma��es
              !pode ser adotado como igual a 1 pelo usu�rio.
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement
      nthreads = GeomechanicalModel%NumPartitions
      
      ALLOCATE( km(ndof,ndof), km_fluid(ndof,ndof), g(ndof), num(nod), points(nip,ndim), weights(nip), coord(nod,ndim),&
                g_coord(ndim,nn) )
      
      g_coord=zero
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
      
      CALL sample(element,points,weights)
          
      DO iel=1,nels
            
        IF( ElementsToAdd(iel) )THEN
        
          !Construindo a matriz de rigidez do elemento.
          num = GeomechanicalModel%ElementNodes(:,iel)   !=g_num(:,iel)
          coord = TRANSPOSE(g_coord(:,num))
        
          !Matriz de Rigidez do Elemento Finito s�lido.
          IF(GeomechanicalModel%ConstitutiveGeoModel(iel) == 3)THEN
            Pc = GeomechanicalModel%MCCparameters(1,iel)
            M = GeomechanicalModel%MCCparameters(2,iel)
            kappa = GeomechanicalModel%MCCparameters(3,iel)  !<-kappa,lambda,M,e0,P0
            lambda = GeomechanicalModel%MCCparameters(4,iel)
            e0 = GeomechanicalModel%MCCparameters(5,iel)
            P0 = GeomechanicalModel%MCCparameters(6,iel)
            defVol = SUM(GeomechanicalModel%StrainTensor(1:3,1,iel))
            evol_vp = SUM(GeomechanicalModel%PlasticStrainTensor(1:3,1,iel))
            !poisson = GeomechanicalModel%Poisson(iel)
            poisson = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
            CALL mccdProps(M,e0,P0,lambda,kappa,SUM(GeomechanicalModel%EffectiveStressTensor(:,1,iel))/3.0,&
                            DefVol,evol_vp,Pc,young,poisson) !<-update e and v.
            GeomechanicalModel%TransverIsotropicProps(iel)%E1 = young
            GeomechanicalModel%TransverIsotropicProps(iel)%E2 = young
            CALL ElementoRigidez( ndim,nip,nod,ndof,nst,young,poisson,&
                                points,weights,coord,km)
          ELSE
            !young = GeomechanicalModel%YoungMod(iel)
            !poisson = GeomechanicalModel%Poisson(iel)
            CALL ElementoRigidezTransvIso(ndim,nip,nod,ndof,nst,&
                    GeomechanicalModel%TransverIsotropicProps(iel),points,weights,coord,km)
          END IF
          
          !Matriz de rigidez do fluido (undrained).
          IF(GeomechanicalModel%Undrained(iel))THEN
            Kf = GeomechanicalModel%PoreFluidBulkMod(iel) * GeomechanicalModel%BiotCoef(iel)
            CALL ElementFluidStiffness(ndim,nip,nod,ndof,nst,Kf,points,weights,coord,km_fluid)
            km = km + km_fluid
          END IF
            
          !Armazenando a matriz do Elemento Finito na Matriz de Rigidez Global.
          g = GeomechanicalModel%ElementDOFs(:,iel)  !=g_g(:,iel)
        
          !Transferindo a matriz de rigidez local para a global.
          DO i=1,ndof
              
            IF(g(i).NE.0)THEN
                
              posicao = ia(g(i))
                
              DO j=1,ndof
                IF(g(j).NE.0)THEN  !->X .AND.(g(j)>=g(i))
                    
                  !Encontra a coluna relativa ao dof.
                  coluna=-1
                  DO
                    coluna=coluna+1
                    IF(ja(posicao+coluna)==g(j)) EXIT
                  END DO
                
                  !Armazena o valor do dof na posi��o correta.
                  a(posicao+coluna) = a(posicao+coluna) + km(i,j)
                
                END IF
              END DO
                
            END IF
          
          END DO
           
        END IF
        
      END DO
    
      !Adjusting 0-diagonal terms.
      DO i = 1, neq
        DO j = ia(i), ia(i+1)-1
          IF( ja(j)==i .AND. a(j)==0.0 )THEN
            a(j) = 1.0 
            EXIT
          END IF
        END DO
      END DO
      
    END SUBROUTINE
  
    !Rotina para construir a matriz de rigidez a partir de uma lista de elementos ativos.
    SUBROUTINE BuildStiffness(GeomechanicalModel,Stiff,ActiveElements)
      
      USE DiscretizedModel
      !USE DynamicExplicitSolver
      USE IO
      USE IFPORT
      
      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      LOGICAL,INTENT(IN):: ActiveElements(:)
      TYPE(ElasticStiffnessMatrix):: Stiff
      
      INTEGER:: iel
      
      LOGICAL,ALLOCATABLE:: ElementsToAdd(:)
      
      !*************************************************************************
      !STIFFNESS MATRIX.
      !*************************************************************************
      
      IF( SIZE(Stiff%ia)==0 )THEN
        CALL SparseStiffnessMatrixCSR( GeomechanicalModel, Stiff%a, &
                                     Stiff%ia, Stiff%ja )
      END IF
      
      ALLOCATE(ElementsToAdd(GeomechanicalModel%NumElements))
      ElementsToAdd = .FALSE.
      DO iel=1,GeomechanicalModel%NumElements
        IF(ActiveElements(iel))THEN
          ElementsToAdd(iel) = .TRUE.
        END IF
      END DO
      
      Stiff%a = 0.0
      CALL StiffnessConstruction( GeomechanicalModel, Stiff%a, &
                                  Stiff%ia, Stiff%ja, ElementsToAdd )
      
      RETURN
    END SUBROUTINE
    
    !Remove a rigidez atual do elemento e substitui por elemento anisotr�pico.
    SUBROUTINE ReformElement( GeomechanicalModel, a, ia, ja, ElNumber, E, Enorm, Normal )
      
      !USE IFPORT
      USE DiscretizedModel
      USE NumericalUtilities
      USE FEMlib
      USE IO
      USE ConstitutiveGeoModels
      USE ElementStiffness
      USE OMP_LIB
      
      IMPLICIT NONE
      
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),INTENT(IN OUT):: a(:)
      INTEGER,INTENT(IN OUT):: ia(:), ja(:)  !Ponteiros que armazenam as posi��es dos dofs na matriz.
      INTEGER,INTENT(IN):: ElNumber
      REAL(8):: E,Enorm,Normal(3)
      
      !VARI�VEIS EST�TICAS.
      INTEGER:: i,j,k,l,&  !contadores.
                posicao,coluna,&  !aux�lio a localiza��o das posi��es no vetor "a".
                nPropVetor,&  !n�mero de propriedades gravadas por elemento.
                no,&  !auxiliar para guardar n�mero do n�.
                GL1,GL2,GL3  !marcadores de grau de liberdade.
      
      INTEGER:: nels,&  !n�mero de elementos.
                nn,&  !n�mero de n�s.
                ndim,&  !n�mero de dimens�es do problema.
                nst,&  !n�mero de componentes no tensor de tens�es.
                nip,&  !n�mero de pontos de integra��o.
                nod,&  !n�mero de n�s por elemento.
                ndof !n�mero de graus de liberdade por elemento.
      
      INTEGER:: nnz  !n�mero de elementos n�o nulos na matriz de rigidez (diagonal + matriz traingular superior).
      INTEGER:: neq  !n�mero de equa��es.

      INTEGER:: nthreads  !n�mero de threads ou n�cleos de processamento que ser�o usados.
      
      REAL(8):: young, poisson
      
      INTEGER:: position, matrixPosition
      REAL(8)::zero=0.0  !vari�veis do processo do gradientes conjugados.
      REAL(8)::startTime,endTime  !vari�veis para marca��o de tempo gasto em segundos.
      CHARACTER(LEN=320)::element  !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
      LOGICAL:: propriedadesPvetor(23)  !vetor bin�rio com as propriedades que devem ser gravadas.
      
      INTEGER,ALLOCATABLE:: g(:),&  !vetor com graus de liberdade de um elemento.
                            num(:) ,&  !conectividade de um elemento.
                            degrees(:,:),&  !associa��o graus de liberdade e elementos.
                            posicoes(:),conection(:,:),conectionOrdenado(:),&  !vetores auxiliares para constru��o de conec��o entre graus de liberdade.
                            propVetorEls(:)  !ponteiro os elementos com propriedades armazenadas (em ordem crescente de numera��o).
      INTEGER(4),ALLOCATABLE:: listaAfinidades(:,:)  !lista de afinidades nodais (8*3,nn), 1a posi��o � um ponteiro.

      INTEGER,POINTER:: nf(:,:)  !matriz de graus de liberdade dos n�s.
                        !g_g(:,:),&  !vetor com os graus de liberdade dos n�s de cada elemento.
                        !g_num(:,:),&  !conectividade dos elementos.
      REAL(8),ALLOCATABLE::   km(:,:),&  !matriz de rigidez do elemento,
                              km_fluid(:,:),&
                              g_coord(:,:),&  !coordenadas dos n�s da malha.
                              coord(:,:),&  !coordenada dos n�s de um elemento.
                              points(:,:),&  !coordenadas dos pontos de integra��o.
                              weights(:),&  !peso dos pontos de integra��o.
                              anisotropia_Ec(:),&  !vetor de anisotropia Ec (=Ec/Ea), onde Ea est� no arquivo .prp.
                              anisotropia_niac(:)  !vetor de anisotropia niab (=niab/niac), onde niac est� no arquivo .prp.
      INTEGER:: globalPosition
      REAL(8):: value
      
      REAL(8):: Kf
                
      REAL(8):: Ea, Ec, niac, niab, sf, df, lh, dee(6,6)
      
      REAL(8):: Pc,M,kappa,lambda,e0,P0,defVol,evol_vp
      
      !INICIALIZA��O DE VARI�VEIS.
      neq = GeomechanicalModel% NumEquations
      element = GeomechanicalModel%element
      nels = GeomechanicalModel%NumElements
      nn = GeomechanicalModel%NumNodes
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = 8 !GeomechanicalModel%NumIntegrationPoints
              !Os pontos de integra��o para c�lculo de tens�es e deforma��es
              !pode ser adotado como igual a 1 pelo usu�rio.
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement
      nthreads = GeomechanicalModel%NumPartitions
      
      ALLOCATE( km(ndof,ndof), km_fluid(ndof,ndof), g(ndof), num(nod), points(nip,ndim), weights(nip), coord(nod,ndim),&
                g_coord(ndim,nn) )
      km_fluid = 0.0
      
      g_coord = 0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
      
      CALL sample(element,points,weights)
          
      !Construindo a matriz de rigidez do elemento.
      num = GeomechanicalModel%ElementNodes(:,ElNumber)   !=g_num(:,iel)
      coord = TRANSPOSE(g_coord(:,num))
      g = GeomechanicalModel%ElementDOFs(:,ElNumber)  !=g_g(:,iel)
        
      !REMOVENDO A MATRIZ DE RIGIDEZ HOMOG�NEA DO ELEMENTO----------------------
      
      !Matriz de Rigidez do Elemento Finito s�lido.
      IF(GeomechanicalModel%ConstitutiveGeoModel(ElNumber) == 3)THEN
        Pc = GeomechanicalModel%MCCparameters(1,ElNumber)
        M = GeomechanicalModel%MCCparameters(2,ElNumber)
        kappa = GeomechanicalModel%MCCparameters(3,ElNumber)  !<-kappa,lambda,M,e0,P0
        lambda = GeomechanicalModel%MCCparameters(4,ElNumber)
        e0 = GeomechanicalModel%MCCparameters(5,ElNumber)
        P0 = GeomechanicalModel%MCCparameters(6,ElNumber)
        defVol = SUM(GeomechanicalModel%StrainTensor(1:3,1,ElNumber))
        evol_vp = SUM(GeomechanicalModel%PlasticStrainTensor(1:3,1,ElNumber))
        !poisson = GeomechanicalModel%Poisson(ElNumber)
        poisson = 0.5 * (GeomechanicalModel%TransverIsotropicProps(ElNumber)%n1 + GeomechanicalModel%TransverIsotropicProps(ElNumber)%n2)
        CALL mccdProps(M,e0,P0,lambda,kappa,SUM(GeomechanicalModel%EffectiveStressTensor(:,1,ElNumber))/3.0,&
                        DefVol,evol_vp,Pc,young,poisson) !<-update e and v.
        GeomechanicalModel%TransverIsotropicProps(ElNumber)%E1 = young
        GeomechanicalModel%TransverIsotropicProps(ElNumber)%E2 = young
        CALL ElementoRigidez( ndim,nip,nod,ndof,nst,young,poisson,&
                            points,weights,coord,km)
      ELSE
        !young = GeomechanicalModel%YoungMod(ElNumber)
        !poisson = GeomechanicalModel%Poisson(ElNumber)
        CALL ElementoRigidezTransvIso(ndim,nip,nod,ndof,nst,&
                GeomechanicalModel%TransverIsotropicProps(ElNumber),points,weights,coord,km)
      END IF
      
      !Matriz de rigidez do fluido (undrained).
      IF(GeomechanicalModel%Undrained(ElNumber))THEN
        Kf = GeomechanicalModel%PoreFluidBulkMod(ElNumber) * GeomechanicalModel%BiotCoef(ElNumber)
        CALL ElementFluidStiffness(ndim,nip,nod,ndof,nst,Kf,points,weights,coord,km_fluid)
        km = km + km_fluid
      END IF
      km = -km

      !Transferindo a matriz de rigidez local para a global.
      DO i=1,ndof
              
        IF(g(i).NE.0)THEN
                
          posicao = ia(g(i))
                
          DO j=1,ndof
            IF(g(j).NE.0)THEN  !->X .AND.(g(j)>=g(i))
                    
              !Encontra a coluna relativa ao dof.
              coluna=-1
              DO
                coluna=coluna+1
                IF(ja(posicao+coluna)==g(j)) EXIT
              END DO
                
              !Armazena o valor do dof na posi��o correta.
              a(posicao+coluna) = a(posicao+coluna) + km(i,j)
                
            END IF
          END DO
                
        END IF
          
      END DO
      
      !ADICIONANDO MATRIZ ORTROT�PICA-------------------------------------------
      !   Ea       - M�dulo de Young, no plano de deposi��o ab               !
      !   Ec       - M�dulo de Young, na dire��o de deposi��o c              !
      !   niac     - M�dulo de Poisson para deforma��o no plano de deposi��o !
      !              devido �s tens�es atuantes na dire��o de deposi��o      !
      !   niab     - M�dulo de Poisson para deforma��o no plano de deposi��o !
      !              devido �s tens�es atuantes no mesmo plano               !
      !   sf[DEG]  - �ngulo azimute da estrutura geol�gica planar (strike)   !
      !              em rela��o ao eixo y (norte) e no sentido hor�rio       !
      !	df[DEG]	 - �ngulo de inclina��o da estrutura geol�gica planar (dip)  ! 
      !              em rela��o ao plano horizontal xy e no sentido hor�rio  !
      Ea = Enorm        !normal to fracture plane.
      Ec = E            !tang. to fracture plane.
      niac = poisson
      niab = poisson
      sf = ATAN( Normal(1)/Normal(2) ) * 180.0/3.1415926535897932384626433832795
      lh = SQRT( Normal(1)**2.0 + Normal(2)**2.0 )
      df = ATAN( Normal(3)/lh ) * 180.0/3.1415926535897932384626433832795
      CALL D_HomogeneousElasticMatrixTransvIso_PUCtecgraf( Ea, Ec, niac, niab, sf, df, dee )
      CALL ElementStiffnessCalculator(ndim,nip,nod,ndof,nst,dee,points,weights,coord,km)
      
      !Adicionando rigidez do fluido.
      km = km + km_fluid
      
      !Transferindo a matriz de rigidez local para a global.
      DO i=1,ndof
              
        IF(g(i).NE.0)THEN
                
          posicao = ia(g(i))
                
          DO j=1,ndof
            IF(g(j).NE.0)THEN  !->X .AND.(g(j)>=g(i))
                    
              !Encontra a coluna relativa ao dof.
              coluna=-1
              DO
                coluna=coluna+1
                IF(ja(posicao+coluna)==g(j)) EXIT
              END DO
                
              !Armazena o valor do dof na posi��o correta.
              a(posicao+coluna) = a(posicao+coluna) + km(i,j)
                
            END IF
          END DO
                
        END IF
          
      END DO
      
    RETURN
    END SUBROUTINE
    
    !Armazenamento de matriz e ponteiros para os solvers iterativos intel e rotinas BLAS 2/ BLAS 3.
    !O esquema de armazenamento � o CSR (compressed sparse row).
    SUBROUTINE SparseStiffnessMatrixCSR( GeomechanicalModel, a, ia, ja )

      USE IFPORT
      USE DiscretizedModel
      USE NumericalUtilities
      USE FEMlib
      USE IO
      USE OMP_LIB
      
      IMPLICIT NONE
      
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),ALLOCATABLE,INTENT(IN OUT):: a(:)
      INTEGER,ALLOCATABLE,INTENT(IN OUT):: ia(:), ja(:)  !Ponteiros que armazenam as posi��es dos dofs na matriz.

      INTEGER:: neq  !n�mero de equa��es.
      
      !VARI�VEIS EST�TICAS.
      
      INTEGER:: i,j,k,l,m,iel,&  !contadores.
                posicao,coluna,&  !aux�lio a localiza��o das posi��es no vetor "a".
                nPropVetor,&  !n�mero de propriedades gravadas por elemento.
                no,&  !auxiliar para guardar n�mero do n�.
                GL1,GL2,GL3  !marcadores de grau de liberdade.
      
      INTEGER:: nels,&  !n�mero de elementos.
                nn,&  !n�mero de n�s.
                ndim,&  !n�mero de dimens�es do problema.
                nst,&  !n�mero de componentes no tensor de tens�es.
                nip,&  !n�mero de pontos de integra��o.
                nod,&  !n�mero de n�s por elemento.
                ndof !n�mero de graus de liberdade por elemento.
      
      INTEGER::nnz  !n�mero de elementos n�o nulos na matriz de rigidez (diagonal + matriz traingular superior).

      INTEGER:: nthreads  !n�mero de threads ou n�cleos de processamento que ser�o usados.
      
      REAL(8):: young, poisson
      
      INTEGER:: position, matrixPosition
      REAL(8)::zero=0.0  !vari�veis do processo do gradientes conjugados.
      REAL(8)::startTime,endTime  !vari�veis para marca��o de tempo gasto em segundos.
      CHARACTER(LEN=320)::element  !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
      LOGICAL:: propriedadesPvetor(23)  !vetor bin�rio com as propriedades que devem ser gravadas.
      
      INTEGER,ALLOCATABLE:: g(:),&  !vetor com graus de liberdade de um elemento.
                            num(:) ,&  !conectividade de um elemento.
                            degrees(:,:),&  !associa��o graus de liberdade e elementos.
                            posicoes(:),conection(:,:),conectionOrdenado(:),&  !vetores auxiliares para constru��o de conec��o entre graus de liberdade.
                            propVetorEls(:)  !ponteiro os elementos com propriedades armazenadas (em ordem crescente de numera��o).
      INTEGER(4),ALLOCATABLE:: listaAfinidades(:,:)  !lista de afinidades nodais (8*3,nn), 1a posi��o � um ponteiro.

      INTEGER,POINTER:: nf(:,:)  !matriz de graus de liberdade dos n�s.
                        !g_g(:,:),&  !vetor com os graus de liberdade dos n�s de cada elemento.
                        !g_num(:,:),&  !conectividade dos elementos.
      REAL(8),ALLOCATABLE:: km(:,:),&  !matriz de rigidez do elemento,
                              g_coord(:,:),&  !coordenadas dos n�s da malha.
                              coord(:,:),&  !coordenada dos n�s de um elemento.
                              points(:,:),&  !coordenadas dos pontos de integra��o.
                              weights(:),&  !peso dos pontos de integra��o.
                              anisotropia_Ec(:),&  !vetor de anisotropia Ec (=Ec/Ea), onde Ea est� no arquivo .prp.
                              anisotropia_niac(:)  !vetor de anisotropia niab (=niab/niac), onde niac est� no arquivo .prp.
      INTEGER:: globalPosition
      REAL(8):: value
              
      !integer, external :: compar

      CHARACTER(LEN=8)::HOUR
      
      !INICIALIZA��O DE VARI�VEIS.
      element = GeomechanicalModel%element
      nels = GeomechanicalModel%NumElements
      nn = GeomechanicalModel%NumNodes
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = 8 !GeomechanicalModel%NumIntegrationPoints
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement
      nthreads = GeomechanicalModel%NumPartitions
      
      neq = GeomechanicalModel%NumEquations
      
      ALLOCATE( nf( SIZE(GeomechanicalModel%NodeDofs(:,1)) , SIZE(GeomechanicalModel%NodeDofs(1,:)) ) )
      nf = GeomechanicalModel%NodeDofs
      
      !ALOCA��O DE MEM�RIA.
      ALLOCATE( km(ndof,ndof), g(ndof), num(nod), posicoes(0:neq),&
                conection(27*3+200,1), conectionOrdenado(27*3+199),&
                degrees(8,0:neq), points(nip,ndim), weights(nip), coord(nod,ndim),&
                g_coord(ndim,nn), listaAfinidades(2*(8*3+5),nn) )
      degrees = 0
      km = 0.0
      listaAfinidades = 0
      
      !Marca��o de tempo
      startTime = OMP_get_wtime()
      
      !Montando vetores de armazenamento da matriz de rigidez > a(:),ia(:),ja(:)
      
      DO i=1,nels
        num = GeomechanicalModel%ElementNodes(:,i)   !=g_num(:,i)
        DO j=1,nod
          CALL ListaAd(listaAfinidades,num(1),num(j))
          CALL ListaAd(listaAfinidades,num(2),num(j))
          CALL ListaAd(listaAfinidades,num(3),num(j))
          CALL ListaAd(listaAfinidades,num(4),num(j))
          CALL ListaAd(listaAfinidades,num(5),num(j))
          CALL ListaAd(listaAfinidades,num(6),num(j))
          CALL ListaAd(listaAfinidades,num(7),num(j))
          CALL ListaAd(listaAfinidades,num(8),num(j))
        END DO
      END DO
      
      nnz=0
      DO i=1,nn
        
        !Marcando os n�s ligados ao n� "i".
        k=0
        DO j=2,listaAfinidades(1,i)  !looping na lista de afinidades do n� i
          no = listaAfinidades(j,i)
          
          IF( nf(1,no).NE.0 ) k=k+1
          IF( nf(2,no).NE.0 ) k=k+1
          IF( nf(3,no).NE.0 ) k=k+1
          
        END DO
        
        !graus de liberdade do n� i -> "l"
        l=0
        IF( nf(1,i).NE.0 ) l=l+1
        IF( nf(2,i).NE.0 ) l=l+1
        IF( nf(3,i).NE.0 ) l=l+1
        
        !graus de liberdade ligados ao n� = n�mero de graus de liberdade do n� X graus de liberdade dos n�s ligados.
        nnz = nnz + l * k
        
      END DO
      
      ALLOCATE( a(nnz) )
      a = 0.0
            
      ALLOCATE( ia(neq+1), ja(nnz) )
      ia = 0
      ja = 0
      ia(neq+1) = nnz+1
      
      !Montando matriz de n�mero de posi��es por grau de liberdade.
      !Posicoes(i=0:neq) -> marca quantos graus de liberdade possuem afinidade com os graus de liberdade do n� "i".
      DO i=1,nn
        
        k=0  !marca o n�mero de graus de liberdade associados a "k"
        DO j=2,listaAfinidades(1,i)
          
          no=listaAfinidades(j,i)
          
          GL1=0
          GL2=0
          GL3=0
          GL1=nf(1,no); IF(GL1.NE.0) k=k+1
          GL2=nf(2,no); IF(GL2.NE.0) k=k+1 
          GL3=nf(3,no); IF(GL3.NE.0) k=k+1
          
        END DO
        
        GL1=nf(1,i); posicoes(GL1)=k
        GL2=nf(2,i); posicoes(GL2)=k
        GL3=nf(3,i); posicoes(GL3)=k
        
      END DO
      
      !Acumula posi��es
      ia(1) = 1
      k=1
      DO i=2,neq
        k=k+posicoes(i-1)
        ia(i)=k
      END DO
      
      !Forma��o dos vetores globais ia e ja.
      nodes: DO i=1,nn
        
        !grau de liberdade 1 do n� i.
      
        DO m=1,3
      
          IF(nf(m,i).NE.0)THEN
            conection=0
            conectionOrdenado=0
          
            DO j=2,listaAfinidades(1,i)
            
              no=listaAfinidades(j,i)
            
              GL1=nf(1,no); IF(GL1.NE.0) CALL ListaAd(conection,GL1,1)
              GL2=nf(2,no); IF(GL2.NE.0) CALL ListaAd(conection,GL2,1)
              GL3=nf(3,no); IF(GL3.NE.0) CALL ListaAd(conection,GL3,1)
            
            END DO
          
            !k=27*3+100
            k=27*3+200
            conectionOrdenado=conection(2:k,1)
            CALL SORTQQ(LOC(conectionOrdenado),k-1,SRT$INTEGER4)
            !CALL qsort( conectionOrdenado, K-1, 4, compar )

            GL1=nf(m,i)
            k=ia(GL1)
          
            l=-1
            !DO j=1,(27*3+99)  !dimens�o do array conectionOrdenado(:,1)
            DO j=1,(27*3+199)  !dimens�o do array conectionOrdenado(:,1)
              IF(conectionOrdenado(j).NE.0)THEN
                l=l+1
                coluna=conectionOrdenado(j)
                ja(k+l)=coluna
              END IF
            END DO
          
          END IF
      
        END DO     
        
      END DO nodes
      

      !marca��o de tempo
      endTime = OMP_get_wtime()
        
    RETURN
    END SUBROUTINE SparseStiffnessMatrixCSR
    
    !Subrotina para calcular a matriz de rigidez de um elemento finito 3D hexa�drico de 8 n�s.
    !Matriz el�stica isotr�pica.
    SUBROUTINE ElementoRigidez(ndim,nip,nod,ndof,nst,young,poisson,points,weights,coord,km)
      
      !USE IFPORT
      USE NumericalUtilities
      USE ConstitutiveGeoModels
      
      IMPLICIT NONE
      
      !VARI�VEIS EST�TICAS.
      
      INTEGER::i  !contadores
      INTEGER,INTENT(IN)::ndim,&  !n�mero de dimens�es do problema
                          nip,&  !n�mero de pontos de integra��o
                          nod,&  !n�mero de n�s por elemento
                          ndof,&  !n�mero de graus de liberdade de deslocamento por elemento
                          nst  !n�mero de componentes no tensor de tens�es
      REAL(8)::det                  !determinantee do jacobiano
      REAL(8),INTENT(IN)::young,&   !m�dulo de young
                          poisson   !raz�o de poisson
      
      !VARI�VEIS DIN�MICAS.
      REAL(8),ALLOCATABLE:: dee(:,:),&  !matriz de rela��es el�sticas.
                            bee(:,:),&  !matriz de derivadas das fun��es de forma.
                            jac(:,:),&  !matriz jacobiana das fun��es de forma.
                            fun(:),&  !fun��es de forma.
                            der(:,:),&  !derivadas das fun��es de forma dos deslocamentos em coordenadas naturais.
                            deriv(:,:)  !derivadas das fun��es de forma dos deslocamentos em coordenadas globais.
      REAL(8),INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integra��o.
                          weights(:),&  !peso dos pontos de integra��o.
                          coord(:,:)  !coordenada dos n�s de um elemento.
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::km(:,:)  !matriz de rigidez do elemento.
      
      !ALOCA��O DE MEM�RIA.
      ALLOCATE(dee(nst,nst),bee(nst,ndof),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim))
      fun=0.0
      
      !MATRIZ DE RELA��ES EL�STICAS.
      CALL D_elastic_homogenious(dee,young,poisson)
      
      !INTEGRA��O NUM�RICA DA MATRIZ DE RIGIDEZ.
      km = 0.0
      DO i = 1, nip
        
        CALL shape_fun(fun,points,i)  !fun��es de forma no ponto de integra��o i
        CALL shape_der(der,points,i)  !derivada das fun��es de forma
        jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integra��o i
        CALL determinante(jac,det)  !determinante do jacobiano
        CALL invert(jac)  !inverso do jacobiano
        deriv=MATMUL(jac,der)  !derivadas das fun��es de forma em coordenadas reais
        CALL beemat(bee,deriv)  !Matriz de derivadas B
        
        km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)  !km=km+BT*D*B*det(J)*peso
        
      END DO
    
      RETURN
    END SUBROUTINE ElementoRigidez
    
    SUBROUTINE ElementoRigidezTransvIso(ndim,nip,nod,ndof,nst,TransverIsotropicProps,points,weights,coord,km)
      
      !USE IFPORT
      USE NumericalUtilities
      USE ConstitutiveGeoModels
      
      IMPLICIT NONE
      
      !VARI�VEIS EST�TICAS.
      
      INTEGER::i  !contadores
      INTEGER,INTENT(IN)::ndim,&  !n�mero de dimens�es do problema
                          nip,&  !n�mero de pontos de integra��o
                          nod,&  !n�mero de n�s por elemento
                          ndof,&  !n�mero de graus de liberdade de deslocamento por elemento
                          nst  !n�mero de componentes no tensor de tens�es
      REAL(8)::det                  !determinantee do jacobiano
      
      TYPE(TransverseIsotroplicElasticity):: TransverIsotropicProps
      
      !VARI�VEIS DIN�MICAS.
      REAL(8),ALLOCATABLE:: dee(:,:),&  !matriz de rela��es el�sticas.
                            bee(:,:),&  !matriz de derivadas das fun��es de forma.
                            jac(:,:),&  !matriz jacobiana das fun��es de forma.
                            fun(:),&  !fun��es de forma.
                            der(:,:),&  !derivadas das fun��es de forma dos deslocamentos em coordenadas naturais.
                            deriv(:,:)  !derivadas das fun��es de forma dos deslocamentos em coordenadas globais.
      REAL(8),INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integra��o.
                          weights(:),&  !peso dos pontos de integra��o.
                          coord(:,:)  !coordenada dos n�s de um elemento.
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::km(:,:)  !matriz de rigidez do elemento.
      
      !ALOCA��O DE MEM�RIA.
      ALLOCATE(dee(nst,nst),bee(nst,ndof),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim))
      fun=0.0
      
      !MATRIZ DE RELA��ES EL�STICAS.
      !CALL D_elastic_homogenious(dee,young,poisson)
      dee = D_TransverseIsotropic(TransverIsotropicProps%E1,TransverIsotropicProps%E2,&
          TransverIsotropicProps%N1,TransverIsotropicProps%N2,TransverIsotropicProps%STRIKE,TransverIsotropicProps%DIP)
      
      !INTEGRA��O NUM�RICA DA MATRIZ DE RIGIDEZ.
      km = 0.0
      DO i = 1, nip
        
        CALL shape_fun(fun,points,i)  !fun��es de forma no ponto de integra��o i
        CALL shape_der(der,points,i)  !derivada das fun��es de forma
        jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integra��o i
        CALL determinante(jac,det)  !determinante do jacobiano
        CALL invert(jac)  !inverso do jacobiano
        deriv=MATMUL(jac,der)  !derivadas das fun��es de forma em coordenadas reais
        CALL beemat(bee,deriv)  !Matriz de derivadas B
        
        km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)  !km=km+BT*D*B*det(J)*peso
        
      END DO
    
      RETURN
    END SUBROUTINE
    
    !Subrotina para calcular a matriz de rigidez de um elemento finito 3D hexa�drico de 8 n�s.
    !Matriz de el�stica ortotr�pila.
    SUBROUTINE ElementoRigidezOrtotropico(ndim,nip,nod,ndof,nst,Ea,Ec,niac,niab,sf_gr,df_gr,points,weights,coord,km)
      
      !USE IFPORT
      USE NumericalUtilities
      USE ConstitutiveGeoModels
      
      IMPLICIT NONE
      
      !VARI�VEIS EST�TICAS.
      INTEGER::i  !contadores
      INTEGER,INTENT(IN)::ndim,&  !n�mero de dimens�es do problema
                          nip,&  !n�mero de pontos de integra��o
                          nod,&  !n�mero de n�s por elemento
                          ndof,&  !n�mero de graus de liberdade de deslocamento por elemento
                          nst  !n�mero de componentes no tensor de tens�es
      REAL(8)::det  !determinantee do jacobiano
      REAL(8),INTENT(IN)::Ea,&  !m�dulo de Young, no plano de deposi��o ab  
                            Ec,&  !m�dulo de Young, na dire��o de deposi��o c
                            niac,&  !m�dulo de Poisson para deforma��o no plano de deposi��o
                                    !devido �s tens�es atuantes na dire��o de deposi��o 
                            niab,&  !m�dulo de Poisson para deforma��o no plano de deposi��o
                                    !devido �s tens�es atuantes no mesmo plano
                            sf_gr,&  !�ngulo (graus) azimute da estrutura geol�gica planar (strike)
                                      !em rela��o ao eixo y (norte) e no sentido hor�rio
                            df_gr  !�ngulo (graus) de inclina��o da estrutura geol�gica planar (dip)
                                    !em rela��o ao plano horizontal xy e no sentido hor�rio 
      !VARI�VEIS DIN�MICAS.
      REAL(8),ALLOCATABLE::dee(:,:),&  !matriz de rela��es el�sticas
                              bee(:,:),&  !matriz de derivadas das fun��es de forma
                              jac(:,:),&  !matriz jacobiana das fun��es de forma
                              fun(:),&  !fun��es de forma
                              der(:,:),&  !derivadas das fun��es de forma dos deslocamentos em coordenadas naturais
                              deriv(:,:)  !derivadas das fun��es de forma dos deslocamentos em coordenadas globais
      REAL(8),ALLOCATABLE,INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integra��o
                                        weights(:),&  !peso dos pontos de integra��o
                                        coord(:,:)  !coordenada dos n�s de um elemento
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::km(:,:)  !matriz de rigidez do elemento
      !aloca��o de mem�ria
      ALLOCATE(dee(nst,nst),bee(nst,ndof),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim))
      fun=0.0
      !matriz de rela��es el�sticas
      CALL D_HomogeneousElasticMatrixTransvIso_PUCtecgraf(Ea,Ec,niac,niab,sf_gr,df_gr,dee)
      !integra��o num�rica da matriz de rigidez
      km=0.0
      DO i=1,nip
        CALL shape_fun(fun,points,i)  !fun��es de forma no ponto de integra��o i
        CALL shape_der(der,points,i)  !derivada das fun��es de forma
        jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integra��o i
        CALL determinante(jac,det)  !determinante do jacobiano
        CALL invert(jac)  !inverso do jacobiano
        deriv=MATMUL(jac,der)  !derivadas das fun��es de forma em coordenadas reais
        CALL beemat(bee,deriv)  !Matriz de derivadas B
        !
        km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)  !km=km+BT*D*B*det(J)*peso
        !
      END DO
      
      RETURN
    END SUBROUTINE ElementoRigidezOrtotropico

    !Subrotina para calcular a matriz de rigidez de um elemento finito 3D hexa�drico de 8 n�s.
    !Matriz el�stica de um fluido.
    SUBROUTINE ElementFluidStiffness(ndim,nip,nod,ndof,nst,Kf,points,weights,coord,km)
      
      !USE IFPORT
      USE NumericalUtilities
      USE ConstitutiveGeoModels
      
      IMPLICIT NONE
      
      !VARI�VEIS EST�TICAS.
      
      INTEGER::i  !contadores
      INTEGER,INTENT(IN)::ndim,&  !n�mero de dimens�es do problema
                          nip,&  !n�mero de pontos de integra��o
                          nod,&  !n�mero de n�s por elemento
                          ndof,&  !n�mero de graus de liberdade de deslocamento por elemento
                          nst  !n�mero de componentes no tensor de tens�es
      REAL(8)::det  !determinantee do jacobiano
      REAL(8),INTENT(IN)::Kf  !fluid bulk modulus.
      
      !VARI�VEIS DIN�MICAS.
      REAL(8),ALLOCATABLE:: dee(:,:),&  !matriz de rela��es el�sticas.
                              bee(:,:),&  !matriz de derivadas das fun��es de forma.
                              jac(:,:),&  !matriz jacobiana das fun��es de forma.
                              fun(:),&  !fun��es de forma.
                              der(:,:),&  !derivadas das fun��es de forma dos deslocamentos em coordenadas naturais.
                              deriv(:,:)  !derivadas das fun��es de forma dos deslocamentos em coordenadas globais.
      REAL(8),ALLOCATABLE,INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integra��o.
                                        weights(:),&  !peso dos pontos de integra��o.
                                        coord(:,:)  !coordenada dos n�s de um elemento.
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::km(:,:)  !matriz de rigidez do elemento.
      
      !ALOCA��O DE MEM�RIA.
      ALLOCATE(dee(nst,nst),bee(nst,ndof),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim))
      fun=0.0
      
      !MATRIZ DE RELA��ES EL�STICAS.
      dee = 0.0
      dee(1:3,1:3) = Kf
      
      !INTEGRA��O NUM�RICA DA MATRIZ DE RIGIDEZ.
      km = 0.0
      DO i = 1, nip
        
        CALL shape_fun(fun,points,i)  !fun��es de forma no ponto de integra��o i
        CALL shape_der(der,points,i)  !derivada das fun��es de forma
        jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integra��o i
        CALL determinante(jac,det)  !determinante do jacobiano
        CALL invert(jac)  !inverso do jacobiano
        deriv=MATMUL(jac,der)  !derivadas das fun��es de forma em coordenadas reais
        CALL beemat(bee,deriv)  !Matriz de derivadas B
        
        km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)  !km=km+BT*D*B*det(J)*peso
        
      END DO
    
      RETURN
    END SUBROUTINE ElementFluidStiffness
    
    !!Rotina para adicionar elements � matriz de rigidez.
    !SUBROUTINE AddHomogeniousElementStiffness( GeomechanicalModel,coord,g,points,weights,young,poisson,a,ia,ja )
    !  
    !  USE IFPORT
    !  USE DiscretizedModel
    !  USE NumericalUtilities
    !  USE FEMlib
    !  USE IO
    !  USE OMP_LIB
    !  
    !  IMPLICIT NONE
    !  
    !  TYPE(Model),INTENT(IN):: GeomechanicalModel
    !  REAL(8),INTENT(IN OUT):: a(:)
    !  INTEGER,INTENT(IN):: ia(:), ja(:)  !Ponteiros que armazenam as posi��es dos dofs na matriz.
    !  REAL(8),INTENT(IN):: young, poisson
    !  INTEGER,INTENT(IN):: g(:)           !vetor com graus de liberdade de um elemento.
    !  REAL(8),INTENT(IN):: coord(:,:),&   !coordenada dos n�s de um elemento.
    !                       points(:,:),&  !coordenadas dos pontos de integra��o.
    !                       weights(:)     !peso dos pontos de integra��o
    !  
    !  !VARI�VEIS EST�TICAS.
    !  INTEGER:: i,j,k,l,m,iel,&  !contadores.
    !            posicao,coluna,&  !aux�lio a localiza��o das posi��es no vetor "a".
    !            nPropVetor,&  !n�mero de propriedades gravadas por elemento.
    !            no,&  !auxiliar para guardar n�mero do n�.
    !            GL1,GL2,GL3  !marcadores de grau de liberdade.
    !  
    !  INTEGER:: nels,&  !n�mero de elementos.
    !            nn,&  !n�mero de n�s.
    !            ndim,&  !n�mero de dimens�es do problema.
    !            nst,&  !n�mero de componentes no tensor de tens�es.
    !            nip,&  !n�mero de pontos de integra��o.
    !            nod,&  !n�mero de n�s por elemento.
    !            ndof !n�mero de graus de liberdade por elemento.
    !  
    !  INTEGER:: nnz  !n�mero de elementos n�o nulos na matriz de rigidez (diagonal + matriz traingular superior).
    !  INTEGER:: neq  !n�mero de equa��es.
    !
    !  INTEGER:: position, matrixPosition
    !  REAL(8)::zero=0.0  !vari�veis do processo do gradientes conjugados.
    !  
    !  CHARACTER(LEN=320)::element  !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
    !  
    !  INTEGER,POINTER:: nf(:,:)  !matriz de graus de liberdade dos n�s.
    !  REAL(8),ALLOCATABLE:: km(:,:)  !matriz de rigidez do elemento.
    !  
    !  INTEGER:: globalPosition
    !  REAL(8):: value
    !  
    !  !INICIALIZA��O DE VARI�VEIS.
    !  neq = GeomechanicalModel% NumEquations
    !  element = GeomechanicalModel%element
    !  nels = GeomechanicalModel%NumElements
    !  nn = GeomechanicalModel%NumNodes
    !  ndim = GeomechanicalModel%NumDimensions
    !  nst = GeomechanicalModel%NumTensorComponents
    !  nip = 8 !GeomechanicalModel%NumIntegrationPoints
    !          !Os pontos de integra��o para c�lculo de tens�es e deforma��es
    !          !pode ser adotado como igual a 1 pelo usu�rio.
    !  nod = GeomechanicalModel%NumNodesPerElement
    !  ndof = GeomechanicalModel%NumDOFUperElement
    !  
    !  ALLOCATE( km(ndof,ndof) )
    !
    !  !Matriz de Rigidez do Elemento Finito s�lido.
    !  CALL ElementoRigidez( ndim,nip,nod,ndof,nst,young,poisson,&
    !                        points,weights,coord,km)
    !    
    !  !Transferindo a matriz de rigidez local para a global.
    !  DO i=1,ndof
    !          
    !    IF(g(i).NE.0)THEN
    !            
    !      posicao = ia(g(i))
    !            
    !      DO j=1,ndof
    !        IF(g(j).NE.0)THEN  !->X .AND.(g(j)>=g(i))
    !                
    !          !Encontra a coluna relativa ao dof.
    !          coluna=-1
    !          DO
    !            coluna=coluna+1
    !            IF(ja(posicao+coluna)==g(j)) EXIT
    !          END DO
    !            
    !          !Armazena o valor do dof na posi��o correta.
    !          a(posicao+coluna) = a(posicao+coluna) + km(i,j)
    !            
    !        END IF
    !      END DO
    !            
    !    END IF
    !      
    !  END DO
    !       
    !
    !
    !END SUBROUTINE
    !
            
    !Subrotina para adicionar n�mero a uma lista de conectividade.
    !1a posica de cada vetor da lista (lista(1,posicao)) - cont�m a posicao do �ltimo elemento preenchido do vetor
    SUBROUTINE ListaAd(lista,valor,posicao)

      IMPLICIT NONE
      INTEGER::i,j
      INTEGER,INTENT(IN)::valor,posicao
      INTEGER,ALLOCATABLE,INTENT(IN OUT)::lista(:,:)
      LOGICAL::existe
          
      i=lista(1,posicao)
      IF(i==0) i=1
      
      !Verifica se o valor existe na lista
      existe=.FALSE.
      DO j=2,i
        IF(lista(j,posicao)==valor) existe=.TRUE.
      END DO
      
      !Adiciona � lista se o element ainda n�o existe.
      IF(existe.EQV..FALSE.)THEN
        lista(i+1,posicao)=valor
        lista(1,posicao)=i+1
      END IF
      
      RETURN
    END SUBROUTINE ListaAd
    
    FUNCTION compar(a1, a2)
        
      IMPLICIT NONE
      INTEGER(2) compar
      INTEGER(2):: a1, a2
      compar=a1-a2
      
    END FUNCTION
       
  END MODULE StiffnessMatrix