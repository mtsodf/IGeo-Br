
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE ViscoplasticExplicitySolver                      
  
    TYPE:: TensorContainer
      REAL(8),ALLOCATABLE:: tensor18a(:,:,:),&  !(18)  sigma - incremento de tensão do passo de carga atual.
                            tensor18b(:,:,:),&  !(18)  eps - termos do tensor de deformações totais (elástica + plástica).
                            tensor18c(:,:,:),&  !(18)  evpt + tensor_dpl - deformação plástica acumulada total (estágio de carregamento atual + estágio anterior).
                            tensor19(:,:,:)     !(19)  evpt - deformação plástica total.
    END TYPE
    
    TYPE:: ViscoSolverParameters
        INTEGER:: maxLinIters = 100E3,&
                  maxViscoIters = 50
        REAL(8):: linTol = 1.0E-6, &
                  viscoTol = 1.0E-2
    END TYPE
  
  CONTAINS
  
  
    SUBROUTINE ViscoSolver(GeomechanicalModel,ReportData,Matrix_a,Matrix_ia,Matrix_ja,&
                            RestrictedDOFs,Pnode,DeltaT,&
                            cg_tol,cg_limit,tol,limit,&
                            Displacement)
  
      USE DiscretizedModel
      USE NumericalUtilities
      USE NumericalLib
      USE FEMlib
      USE ConstitutiveGeoModels
      USE HydraulicFractureModel
      USE Damage
      USE IO
      USE OMP_LIB  !bibliotecas OpenMP(R)
  
      IMPLICIT NONE

      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      TYPE(StepReport),INTENT(IN OUT):: ReportData
      REAL(8),INTENT(IN OUT):: Matrix_a(:)
      INTEGER,INTENT(IN OUT):: Matrix_ia(:), Matrix_ja(:) !ponteiros que armazenam as posições dos dofs na matriz.
      INTEGER,INTENT(IN OUT):: RestrictedDOFs(0:)
      REAL(8),INTENT(IN OUT):: Displacement(0:),Pnode(:)
      REAL(8),INTENT(IN):: DeltaT(:)
      REAL(8),INTENT(IN)::cg_tol,&                        !tolerância do solver linear de equações - razão de máximos deslocamentos
                          tol                             !tolerância do solver viscoplástico - razão de máximas forças desbalanceadas
      INTEGER,INTENT(IN)::limit,&                         !limite de iterações para o solver viscoplástico
                          cg_limit                        !limete de iterações para o solver linear de equações
                          
      INTEGER:: iel,&
                iters,&                                   !contador das iterações plásticas
                fracIters                                 !contador das iterações de propagação de fraturas.
      INTEGER:: i,j
      INTEGER::ndim,&                                     !número de dimensões
               nels,&                                     !número de elementos na malha
               nn,&                                       !número de nós na malha
               nod,&                                      !número de nós de deslocamento por elemento (=8 para hexaédro linear)
               nip,&                                      !número de pontos de integração
               ndof,&                                     !número de graus de liberdade de deslocamento por elemento
               nst,&                                      !número de componentes de tensão no tensor de tensões (=6 para 3D)
               neq                                        !número de equações no sistema
      REAL(8):: fmax,&                                    !função de escoamento do modelo de Mohr-Coulomb e CAP
                factor,&                                  !factor para correlação do passo de tempo
                dt=1.0e15,&                               !pseudo-passo de tempo do solver viscoplástico
                convergence                               !convergência
      LOGICAL::converged                                  !flag de convergência do solver viscoplástico
  
      !Variáveis dinâmicas
      REAL(8),ALLOCATABLE:: points(:,:),&                 !coordenadas dos pontos de integração
                            weights(:),&                  !peso dos pontos de integração
                            g_coord(:,:),&                !coordenada dos nós
                            bdylds(:),&                   !cargas nodais desbalanceadas do solver viscoplástico
                            oldis(:),&                    !guarda os deslocamentos da iteração anterior
                            tensor_d(:,:),&               !tensor de deformações totais do estágio de carregamento anterior anterior
                            x(:)                          !solução do sistema de equações
      REAL(8),ALLOCATABLE:: loads(:)                      !carregamento e deslocamentos (OUT)
      
      REAL(8),ALLOCATABLE:: evpt(:,:),&                   !deformação plástica acumulada
                            tensor_dpl(:,:)               !tensor de deformações plásticas do estágio de carregamento anterior anterior
      
      REAL(8),ALLOCATABLE:: diag_precon(:)
      
      CHARACTER(LEN=320):: text
      
      TYPE(TensorContainer):: Container
      
      LOGICAL:: newFractures
      
      INTEGER,ALLOCATABLE:: newRestrictedDOFs(:)
      
      REAL(8):: temp,damageVal,syy,ST
      INTEGER:: modelNumber
      
      REAL(8),ALLOCATABLE:: fmaxHistory(:)
      
      !CALL PrintMessage( .FALSE., .FALSE., '...Structural Viscoplastic Solver' )
            
      reportData%linIters = 0
      reportData%nonLinIters = 0
      
      ALLOCATE(fmaxHistory(0:limit))
      fmaxHistory = 0.0
      fmaxHistory(0)=1.0
      
      neq = GeomechanicalModel%NumEquations
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = GeomechanicalModel%NumIntegrationPoints
      nn = GeomechanicalModel%NumNodes
      
      !ALLOCATE( Matrix_a_reformed(SIZE(Matrix_a)))
      
      ALLOCATE(bdylds(0:neq),oldis(0:neq),g_coord(ndim,nn),&
               loads(0:neq),points(nip,ndim),weights(nip),tensor_d(nst,nip),x(neq))
      bdylds = 0.0
      oldis = Displacement
      Displacement = 0.0
      
      ALLOCATE( newRestrictedDOFs(0:SIZE(RestrictedDOFs)-1) )
      newRestrictedDOFs = 0
      
      g_coord = 0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
      
      ALLOCATE(evpt(nst,nip),tensor_dpl(nst,nip))
      
      ALLOCATE( Container%tensor18a(nst,nip,GeomechanicalModel%NumElements),&
                Container%tensor18b(nst,nip,GeomechanicalModel%NumElements),&
                Container%tensor18c(nst,nip,GeomechanicalModel%NumElements),&
                Container%tensor19(nst,nip,GeomechanicalModel%NumElements))
      Container%tensor18a = 0.0
      Container%tensor18b = 0.0
      Container%tensor18c = 0.0
      Container%tensor19 = 0.0
      
      !DIAGONAL PRE-CONDITIONER.
      diag_precon = 0.0
      CALL CGpreConditioner(Matrix_a,Matrix_ia,Matrix_ja,diag_precon)
      
      CALL sample(GeomechanicalModel%element,points,weights)
    
      !*************************************************************************
      !VISCOPLASTIC STEPS
      !*************************************************************************
      !WRITE(text,'(A)')'Iters  PlasticWork     Converg  MaxYield     MaxDisp'
      !CALL PrintMessage(.FALSE.,.TRUE.,TRIM(text))

      !!Initial viscoplastic time step.
      !CALL ViscoPlasTimeStep(GeomechanicalModel,dt)
      
      fracIters = 0
      iters = 0
      viscoplasticItaration: DO!------------------------------------------------

        iters = iters+1
        fmax = -1.0e15
        
        reportData%nonLinIters = reportData%nonLinIters + 1
      
        !>>>Initial viscoplastic time step.
        CALL ViscoPlasTimeStep(GeomechanicalModel,dt)
      
        !acumulando carga estrutural: carga externa + carga viscoplástica
        loads(1:) = Pnode + bdylds(1:)
    
        !	Solução do Sistema Linear de Equações.(*)
        IF(MAXVAL(ABS(loads(1:)))>0.0)THEN
      
          !LINEAR SOLVER--------------------------------------------------------
          ! [A]*{x}={loads}
          !   [A] -> a,ia,ja :matriz de rigidez no formato CSR
          !   {loads} -> carga estrutural do passo de carregamento (Pnode) mais a carga elastoplástica desbalanceada (bdylds)
          !   {x} -> solução, escrita no vetor "loads"
          !Initial guess -> oldis.
          CALL PCG( cg_limit,&
                    neq,&
                    cg_tol,&
                    diag_precon,&
                    Matrix_a,Matrix_ia,Matrix_ja,&
                    loads,&
                    oldis(1:neq),RestrictedDOFs(1:RestrictedDOFs(0)),&
                    GeomechanicalModel%NumPartitions,reportData%linIters)                          !solver linear iterativo PCG para matriz simétrica.
          loads(0) = 0
        ELSE
          loads = 0.0                                                          !zero displacements.
        END IF

        !	UNBALANCED FORCE------------------------------------------------------
        CALL ViscoplasticMachine( 3,GeomechanicalModel%NumElements,&
                                  GeomechanicalModel%NumNodesPerElement,&
                                  GeomechanicalModel%NumIntegrationPoints,&
                                  GeomechanicalModel%NumDOFUperElement,6,&
                                  iters,limit,&
                                  dt,&
                                  GeomechanicalModel%ElementNodes,&
                                  GeomechanicalModel%ElementDOFs,&
                                  g_coord,&
                                  points,weights,&
                                  loads,&
                                  factor,fmax,&
                                  bdylds,DeltaT,GeomechanicalModel,Container)   
        fmaxHistory(iters) = fmax
        
        ! UPDATING TIME STEP----------------------------------------------------
        !IF( (0.1 * factor) > dt .AND. &
        !  fmaxHistory(iters)<fmaxHistory(iters-1) )THEN
        !  dt = 1.1 * dt
        !ELSEIF(fmaxHistory(iters)>fmaxHistory(iters-1))THEN
        !  !dt = dt * fmaxHistory(iters-1)/fmaxHistory(iters)
        !ELSE 
        !END IF
          
        
        !SAVING TENSORS---------------------------------------------------------
        DO iel=1,GeomechanicalModel%NumElements
          !READ(17)tensor_dpl;READ(18)sigma; READ(18)eps; READ(18)evpt
          !WRITE(19)evpt-tensor_dpl  !evpt - deformação plástica acumulada nas iterações atuais / tensor_dpl - deformação plástica acumulada total (todos os passos de carga)
          tensor_dpl = GeomechanicalModel%PlasticStrainTensor(:,:,iel)
          evpt = Container%tensor18c(:,:,iel)
          Container%tensor19(:,:,iel) = evpt - tensor_dpl
        END DO
        
        ! Damage----------------------------------------------------------------
        CALL DamageUpdate(GeomechanicalModel,Container%tensor18a,Container%tensor18c)
                  
        !PRINT STUFF------------------------------------------------------------
        !Iters  plasticWork  Converg  MaxYield'
        !WRITE(text,'(I5,E13.2,E12.4,E10.2,E12.3)')iters,SUM(GeomechanicalModel%PlasticStrainWork),&
        !                                            (1-convergence),fmax,MAXVAL(ABS(loads(1:)))
        !CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
                
        reportData%maxVAL = MAXVAL(loads(1:))
        reportData%stepResult = 'P'
        CALL PrintStepReport(reportData%time,reportData%globalStep,reportData%flowStep,&
                                reportData%structuralStep,reportData%fractureStep,&
                                reportData%linIters,reportData%nonLinIters,&
                                reportData%maxVAL,reportData%stepResult) 
        
        !	STOP TEST-------------------------------------------------------------
        convergence = 1.0
        
        !Discplacement criteria.
        CALL checon(loads(1:),oldis(1:),tol,converged,convergence)                  !"loads" -> "oldis".
        
        converged = convergence<=tol
        
        IF(iters==1)THEN
          convergence = 1.0
          converged=.FALSE.
        END IF
        IF(converged.OR.iters==limit) bdylds=0.0                                     !zero unbalanced force at final (!)
        IF(converged.OR.iters==limit) EXIT
        IF(fmax.LE.0.0) EXIT
        
        oldis = loads
        
      END DO viscoplasticItaration!---------------------------------------------
      !CALL PrintMessage(.FALSE.,.TRUE.,'...Structural Viscoplastic Solver Finished.')
      !CALL PrintMessage(.FALSE.,.FALSE.,'----------------------------------------------------')
      
      Displacement = Displacement + loads
      
      CALL UpdatePorepressure(GeomechanicalModel%ElementNodes,GeomechanicalModel%ElementDOFs,&
                              g_coord,Displacement,points,weights,GeomechanicalModel)
      
      !*************************************************************************
      ! SAVING STRESS AND STRAIN.
      !*************************************************************************
      GeomechanicalModel%Displacement = GeomechanicalModel%Displacement + Displacement(1:)
      GeomechanicalModel%EffectiveStressTensor = Container%tensor18a
      GeomechanicalModel%StrainTensor = Container%tensor18b
      GeomechanicalModel%PlasticStrainTensor = Container%tensor18c
    
    END SUBROUTINE
  
    
    ! Unbalanced force calculator.
    SUBROUTINE ViscoplasticMachine(ndim,nels,nod,nip,ndof,nst,iters,limit,dt,g_num,&
                                   g_g,g_coord,points,weights,loads,factor,fmax,bdylds,&
                                   DeltaT,GeomechanicalModel,Container)
  
      USE DiscretizedModel
      USE ConstitutiveGeoModels
      USE SolidMechanicsLib
      USE FEMlib
      USE NumericalUtilities
      USE OMP_LIB  !bibliotecas OpenMP(R)
  
      IMPLICIT NONE
  
      !******************************************
      !VARIÁVEIS
      !******************************************
      
      !Input.
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      INTEGER,INTENT(IN)::ndim,&  !número de dimensões
                          nels,&  !número de elementos na malha
                          nod,&  !número de nós de deslocamento por elemento (=8 para hexaédro linear)
                          nip,&  !número de pontos de integração
                          ndof,&  !número de graus de liberdade de deslocamento por elemento
                          nst,&  !número de componentes de tensão no tensor de tensões (=6 para 3D)
                          iters,limit  !contador das iterações plásticase e limite de iterações para o solver viscoplástico  
      REAL(8),INTENT(IN OUT):: dt  !pseudo-passo de tempo do solver viscoplástico
      REAL(8),INTENT(OUT):: factor  !factor para correlação do passo de tempo
      REAL(8),INTENT(IN OUT):: fmax  !máxima função de escoamento do modelo de Mohr-Coulomb
      INTEGER,ALLOCATABLE,INTENT(IN):: g_num(:,:),&  !conectividade dos elementos
                                       g_g(:,:)  !matriz para identificação dos graus de liberdades dos nós
      REAL(8),ALLOCATABLE,INTENT(IN):: points(:,:),&  !coordenadas dos pontos de integração
                                       weights(:),&  !peso dos pontos de integração
                                       g_coord(:,:),&  !coordenada dos nós
                                       loads(:)  !deslocamentos nodais
      REAL(8),INTENT(IN):: DeltaT(:)
      REAL(8),ALLOCATABLE,INTENT(IN OUT):: bdylds(:)  !cargas nodais viscoplástiacas (desbalanceadas) do solver viscoplástico.
      TYPE(TensorContainer),INTENT(IN OUT):: Container
      
      !Variáveis estáticas.
      INTEGER:: iel,i  !contadores
      REAL(8):: f,&  !função de escoamento do modelo de Mohr-Coulomb e CAP
                defVOL,&  !deformação volumétrica
                evol_vp,& !viscoplastic volumetric strain.
                e,v,&  !young e poisson
                ST,&  !traction strength
                eterm(6),&  !deformação de origem térmica
                det,&  !determinante de matriz
                prop(23),&  !vetor para leitura de propriedades do elemento do arquivo .prp
                Ea,Ec,niac,niab,&  !parâmetros do modelo ortotrópico (ver interface no módulo MODULE iMaterials)
                ddt,new_dt,&   !passo de tempo.
                zero=0.0,&
                pi=ACOS(-1.)
      REAL(8):: biotCoef, Kf
      REAL(8):: sigm, dsbar, theta  !invariantes de tensão.
      REAL(8):: f_cc,&    !função de escoamento do Cam-Clay.
                fator_cc  !factor para retorno de tensões para o envelope do CAM-CLAY.
      LOGICAL:: converged,&  !flag de convergência do solver viscoplástico
                termo,&  !identifica se é uma simulação com esforço térmico (cartão #TERMO no arquivo base)
                anisotropia,&  !identificador para modelo transversalmente isotrópico
                found,&  !auxiliar
                update_dt !flag para atualização do passo de tempo crítico. 
  
      !Variáveis dinâmicas.
      INTEGER,ALLOCATABLE:: num(:) ,&                     !conectividade de um elemento
                            g(:)                          !vetor com os graus de liberdade de cada elemento
      REAL(8),ALLOCATABLE:: coord(:,:),&                  !coordenada dos nós de um elemento
                            jac(:,:),&                    !matriz jacobiana das funções de forma.
                            der(:,:),&                    !derivadas das funções de forma dos deslocamentos em coordenadas naturais.                        
                            deriv(:,:),&                  !derivadas das funções de forma com relação às coordenadas locais e globais
                            bload(:),&                    !cargas nodais desbalanceadas do solver viscoplástico para um elemento
                            bee(:,:),&                    !matriz de derivadas das funções de forma
                            dee(:,:),&                    !matriz de relações elásticas
                            devp(:),&                     !pseudo-tensões elásticas devido a deformação plástica (=D*delta_evp)
                            deltaSyy(:),&                 !unbalanced YY stress due traction criteria
                            eld(:),&                      !deslocamentos nodais de um elemento
                            eload(:),&                    !contribuição de um ponto de integração para o vetor global bload
                            erate(:),&                    !taxa de deformação plástiva
                            evp(:),&                      !incremento de deformação plástica
                            evpt(:,:),&                   !deformação plástica acumulada
                            eps(:,:),&                    !termos do tensor de deformações totais (elástica + plástica)
                            sigma(:,:),&                  !incremento do tensor de tensões
                            tensor_s(:,:),&               !tensor de tensões do estágio de carregamento anterior
                            tensor_dpl(:,:),&             !tensor de deformações plásticas do estágio de carregamento anterior anterior
                            tensor_d(:,:),&               !tensor de deformações totais do estágio de carregamento anterior anterior
                            anisotropia_Ec(:),&           !vetor de anisotropia Ec (=Ec/Ea), onde Ea está no arquivo .prp.
                            anisotropia_niac(:),&         !vetor de anisotropia niab (=niab/niac), onde niac está no arquivo .prp.
                            shearDamageArray(:,:),&       !matriz de 2 coluna para dano da resistência ao cisalhamento nos elementos.
                            elasticDamageArray(:)         !vetor de dano ao módulo de elasticidade.
      REAL(8),ALLOCATABLE:: viscoLoad(:)                  !carga viscoplástica (prévio a multiplicação pelo passo de tempo).
      
      REAL(8):: currentStress(6,1)
      
      REAL(8):: dilation !plasticDilation
      REAL(8),ALLOCATABLE:: damageTraction(:)
      INTEGER:: modelNumber
      REAL(8):: temp
      
      REAL(8):: deltaTemp
      
      REAL(8):: Pc,M,kappa,lambda,e0,P0
      
      REAL(8):: I1,I2,I3
      
      !Alocação de memória.
      ALLOCATE( g(ndof),num(nod),bee(nst,ndof),eld(ndof),bload(ndof),                   &
                coord(nod,ndim),dee(nst,nst),deriv(ndim,nod),devp(nst),deltaSyy(nst),   &
                eload(ndof),erate(nst),evp(nst),sigma(nst,nip),evpt(nst,nip),           &
                eps(nst,nip),tensor_s(nst,nip),tensor_d(nst,nip),tensor_dpl(nst,nip),   &
                damageTraction(nip))
      ALLOCATE( viscoLoad(0:GeomechanicalModel%NumEquations) )
      viscoLoad = 0.0
      new_dt = 1.0e15
      update_dt = .FALSE.
      
      ALLOCATE(der(ndim,nod),jac(ndim,ndim))
  
      !******************************************
      !CÁLCULO DE DEFORMAÇÕES, TENSÕES E CARGA
      !VISCOPLÁSTICA DESBALANCEADA.
      !******************************************
      evpt=zero
      tensor_dpl=zero
      bee=zero
      f=zero
      f_cc=zero
      evp=zero
      factor=1.0E15
  
      !Element looping for stress calculation.
      DO iel=1,nels
    
        !Properties.
        prop(2) = 0
        IF( GeomechanicalModel%ConstitutiveGeoModel(iel)==2 ) prop(2)=1
        !prop(12) = GeomechanicalModel%YoungMod(iel)
        prop(12) = GeomechanicalModel%TransverIsotropicProps(iel)%minE
        !prop(13) = GeomechanicalModel%Poisson(iel)
        prop(13) = GeomechanicalModel%TransverIsotropicProps(iel)%minNi
        prop(14) = GeomechanicalModel%BiotCoef(iel)
        prop(15) = GeomechanicalModel%ThermalExpCoeff(iel)
        prop(16) = GeomechanicalModel%Cohesion(iel) * (1.-GeomechanicalModel%DamageShear(iel))
        prop(17) = TAN(GeomechanicalModel%FrictionAngle(iel) * pi/180) * (1.-GeomechanicalModel%DamageShear(iel))
        prop(17) = ATAN(prop(17)) * 180./pi 
        prop(18) = GeomechanicalModel%DilationAngle(iel)
        
        IF(GeomechanicalModel%ConstitutiveGeoModel(iel) == 3)THEN
          Pc = GeomechanicalModel%MCCparameters(1,iel)
          M = GeomechanicalModel%MCCparameters(2,iel)
          kappa = GeomechanicalModel%MCCparameters(3,iel)  !<-kappa,lambda,M,e0,P0
          lambda = GeomechanicalModel%MCCparameters(4,iel)
          e0 = GeomechanicalModel%MCCparameters(5,iel)
          P0 = GeomechanicalModel%MCCparameters(6,iel)
          defVol = SUM(GeomechanicalModel%StrainTensor(1:3,1,iel))
          evol_vp = SUM(GeomechanicalModel%PlasticStrainTensor(1:3,1,iel))
          !v = GeomechanicalModel%Poisson(iel)
          v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
          CALL mccdProps(M,e0,P0,lambda,kappa,SUM(GeomechanicalModel%EffectiveStressTensor(:,1,iel))/3.0,&
                          DefVol,evol_vp,Pc,e,v) !<-update e and v.
          GeomechanicalModel%TransverIsotropicProps(iel)%E1 = e
          GeomechanicalModel%TransverIsotropicProps(iel)%E2 = e
          CALL D_elastic_homogenious(dee,e,v)
        ELSE
          !e = GeomechanicalModel%YoungMod(iel)
          !v = GeomechanicalModel%Poisson(iel)
          dee = D_TransverseIsotropic(GeomechanicalModel%TransverIsotropicProps(iel)%E1,&
            GeomechanicalModel%TransverIsotropicProps(iel)%E2,&
            GeomechanicalModel%TransverIsotropicProps(iel)%n1,&
            GeomechanicalModel%TransverIsotropicProps(iel)%n2,&
            GeomechanicalModel%TransverIsotropicProps(iel)%strike,&
            GeomechanicalModel%TransverIsotropicProps(iel)%dip)
        END IF
        
        biotCoef = GeomechanicalModel%BiotCoef(iel)
        Kf = GeomechanicalModel%PoreFluidBulkMod(iel)
    
        !FE incidence.
        num = g_num(:,iel)
        coord = TRANSPOSE(g_coord(:,num))
        g = g_g(:,iel)
    
        !Stress-Strain relation.
        !CALL D_HomogeneousElasticMatrix(dee,e,v)
    
        !Previous stress and strain.
        tensor_s = GeomechanicalModel%EffectiveStressTensor(:,:,iel)  !READ(15)tensor_s
        tensor_d = GeomechanicalModel%StrainTensor(:,:,iel)           !READ(16)tensor_d
        tensor_dpl = GeomechanicalModel%PlasticStrainTensor(:,:,iel)  !READ(17)tensor_dpl
        evpt =  Container%tensor19(:,:,iel)                           !READ(19)evpt  !deformação plástica acumulada da última iteração do passo atual
    
        !Thermal strain.
        deltaTemp = SUM(DeltaT(num)) / 8.0
        CALL TensorDeformacaoTermica(GeomechanicalModel%ThermalExpCoeff(iel),deltaTemp,eterm)
    
        !Gauss Points Looping.
        eld = loads(g)
        bload = zero
        eps = zero
        !damageTraction = 0.0
        
        DO i=1,nip
          f = zero
          
          !Derivative of shape functions [B].(*)
          CALL shape_der(der,points,i)  !derivada das funções de forma
          jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integração i
          CALL determinante(jac,det)  !determinante do jacobiano
          CALL invert(jac)  !inverso do jacobiano
          deriv=MATMUL(jac,der)  !derivadas das funções de forma em coordenadas reais
          CALL beemat(bee,deriv)  !Matriz de derivadas B
        
          !Strain ( {def}=[B]*{disp} ).
          eps(:,i) = MATMUL(bee,eld)
          eps(:,i) = eps(:,i) - evpt(:,i)   !retirando a deformação viscoplástica de eps()
          eps(:,i) = eps(:,i) - GeomechanicalModel%CreepStrainIncrementTensor(:,i,iel)  !retirando a deformação de fluência
      
          !Elastic stress ( {st}=[D]*{defTerm} ).
          sigma(:,i) = zero
          sigma(:,i) = MATMUL(dee,eps(:,i))  !D*(e_total-e_visco)+tensões do step anterior
      
          !Thermal stress.
          sigma(:,i) = sigma(:,i) - MATMUL(dee,eterm)
      
          !delta viscoplastic strain "devp" -> "evpt".
          IF( converged.OR.(iters==limit) )THEN
            !final da iteração. Não roda plasticidade.
            devp = sigma(:,i)  !passing the current stress increment.
          ELSE
          
            !************************************************************************************
            !VISCOPLASTIC CALCULATION.
            f=-1.0E6
            erate = 0.0
          
            !MOHR-COULOMB------------------------------------------------------------------------
            !Calculating yield function (f) and strain rate (erate).
            IF( GeomechanicalModel%ConstitutiveGeoModel(iel)==2 .OR. GeomechanicalModel%ConstitutiveGeoModel(iel)==4  ) THEN
              CALL MCDefVispl( i, prop, nst, sigma, tensor_s, f, erate)
            END IF
            
            !CAM-CLAY----------------------------------------------------------------------------
            IF( GeomechanicalModel%ConstitutiveGeoModel(iel)==3  ) THEN
              GeomechanicalModel%MCCparameters(8,iel) = SUM(evpt(1:3,i) + tensor_dpl(1:3,i))
              CALL CCDefVispl(tensor_s(:,i)+sigma(:,i),GeomechanicalModel%MCCparameters,f,erate)
            END IF
              
            !YIELD FUNCTION----------------------------------------------------------------------
            IF( GeomechanicalModel%ConstitutiveGeoModel(iel)==2 ) f = f - GeomechanicalModel%YieldEnvelopeTOL(iel)
            IF( ( f>fmax ) ) fmax=f
          
            !VISCOPLASTIC STRAIN AND UNBALANCED STRESS-------------------------------------------
            IF( (f>=zero).OR.(converged.OR.iters==limit) )THEN
              evp = erate * dt * GeomechanicalModel%Gamma(iel)          !{DeltaDeformacaoPlastica}={TaxaDeformacaoViscoplastica}*Dt .
              devp = MATMUL(dee,evp)                                    !DeltaTensão devido o DeltaDeformacaoPlastica (TENSÃO DESBALANCEADA).
            ELSE
              evp = zero
              devp = zero
            END IF
            
            !Acumulando Def.Plas. total
            evpt(:,i) = evpt(:,i) + evp
            
            !INCREMENT OF PLASTIC STRAIN ENERGY--------------------------------------------------
            !currentStress(:,1) = tensor_s(:,i)+sigma(:,i)-devp
            GeomechanicalModel%PlasticStrainWork(iel) = GeomechanicalModel%PlasticStrainWork(iel) + &
                                                        0.5*DOT_PRODUCT(evp,sigma(:,i)-devp+tensor_s(:,i)) * det * weights(i)

            !FIM DO CÁLCULO VISCOPLÁSTICO.
            !************************************************************************************
          END IF
      
          !Integration of umbalanced load "eload" ( {b} = devp*[B] ).
          IF(f>=zero.OR.(converged.OR.iters==limit)) eload = MATMUL(devp,bee)               !transposto[B] * {DeltaTensão}
          IF(f>=zero.OR.(converged.OR.iters==limit)) bload = bload + eload*det*weights(i)   !integral(BT*Delta.Tensão)

          GeomechanicalModel%YieldFunc(i,iel) = f
      
          !Time step increment ("factor")--------------------------------------------------------
          !volumetric strain criteria.
          !IF( (erate(1)+erate(2)+erate(3)).NE.0.0 )THEN
          !  defVOL = eps(1,i) + eps(2,i) + eps(3,i) + evp(1) + evp(2) + evp(3)
          !  temp = ABS( defVOL / (erate(1)+erate(2)+erate(3)) )
          !  IF(temp<factor) factor = temp
          !END IF
          !!shear strain criteria.
          !CALL invar(erate,I1,I2,I3)
          !IF( I2>1.0e-3 )THEN
          !  temp = I2
          !  CALL invar(eps(:,i),I1,I2,I3)
          !  temp = 2/3 * I2/temp
          !  IF(TEMP>0.0.AND.(temp<factor)) factor = temp
          !END IF
          
        END DO  !Gauss points
        
        !DAMAGE TRACTION-------------------------------------------------------------------------
        damageTraction(1) = MAXVAL(damageTraction)
        IF(damageTraction(1)>GeomechanicalModel%DamageTraction(iel)) GeomechanicalModel%DamageTraction(iel) = damageTraction(1)
        
        !Atualização das forças desbalanceadas globais com a contribuição do elemento------------
        !bdylds(g) = bdylds(g) + bload
        !bdylds(0) = zero
        viscoLoad(g) = viscoLoad(g) + bload
    
        !Saving stress and strain from current step----------------------------------------------
        Container%tensor18a(:,:,iel) = sigma + tensor_s       !WRITE(18) sigma + tensor_s
        Container%tensor18b(:,:,iel) = eps + evpt + tensor_d  !WRITE(18) eps + evpt + tensor_d
        Container%tensor18c(:,:,iel) = evpt + tensor_dpl      !WRITE(18) evpt + tensor_dpl
    
      END DO !elements***************************************************************************
    
      bdylds = bdylds + viscoLoad
      bdylds(0) = zero
    
      IF(factor<9.0E14)THEN
        factor = SQRT(factor)
      ELSE
        factor = dt
      END IF
                                   
    RETURN
    END SUBROUTINE                                
                        
                                   
    ! CRITICAL TIME STEP.                                 
    SUBROUTINE ViscoPlasTimeStep(GeomechanicalModel,dt)
  
      USE DiscretizedModel
      USE ConstitutiveGeoModels
      USE SolidMechanicsLib
  
      IMPLICIT NONE
      REAL(8),INTENT(IN OUT)::dt  !pseudo-passo de tempo do solver viscoplástico
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      
      INTEGER::iel  !contador
      REAL(8):: ddt,start_dt=1.e15,&                          !pseudo-passo de tempo do solver viscoplástico
                !prop(23),&                                   !vetor para leitura de propriedades do elemento do arquivo .prp
                phi,snph,&                                    !ângulo de atrito e seno
                psi,snps,&                                    !ângulo de dilatação e seno
                e,v,&                                         !módulo de young e coficiente de poisson
                pi,one=1.0,two=2.0,d180=180.0
      
      REAL(8):: dee(6,6)
      
      REAL(8):: c1,a,b,gamma,He,Hp,p,q,lode_theta,f0,f
      REAL(8):: Pc,M,kappa,lambda,e0,P0,defVol,evol_vp,G,HeHp,Pc_int
      
      INTEGER:: i,nip
      
      LOGICAL::found
      
      nip = GeomechanicalModel%NumIntegrationPoints
      
      GeomechanicalModel%Gamma = 1.
      
      !Cálculo do passo de tempo crítico.
      pi=ACOS(-one)
      dt=start_dt
      DO iel=1,GeomechanicalModel%NumElements
        
        IF(GeomechanicalModel%ConstitutiveGeoModel(iel)==2)THEN       !Mohr-Coulomb
          phi=TAN(GeomechanicalModel%FrictionAngle(iel) * pi/180) * (1.-GeomechanicalModel%DamageShear(iel))
          phi=ATAN(phi) * 180./pi
          snph=SIN(phi*pi/d180)
          psi=GeomechanicalModel%DilationAngle(iel)
          snps=SIN(psi*pi/d180)
          !e=GeomechanicalModel%YoungMod(iel)
          !v=GeomechanicalModel%Poisson(iel)
          e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
          v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)

          ddt=two*(one+v)*(one-two*v)/(e*(one-two*v+snph*snps))     !critério 1 de Stolle 1991, Eq. (32)
          IF(ddt<dt)dt=ddt
          
        ELSEIF(GeomechanicalModel%ConstitutiveGeoModel(iel)==3)THEN   !Cam-Clay.
          Pc = GeomechanicalModel%MCCparameters(1,iel)
          M = GeomechanicalModel%MCCparameters(2,iel)
          kappa = GeomechanicalModel%MCCparameters(3,iel)             !<-kappa,lambda,M,e0,P0
          lambda = GeomechanicalModel%MCCparameters(4,iel)
          e0 = GeomechanicalModel%MCCparameters(5,iel)
          P0 = GeomechanicalModel%MCCparameters(6,iel)
          defVol = SUM(GeomechanicalModel%StrainTensor(1:3,1,iel))
          evol_vp = SUM(GeomechanicalModel%PlasticStrainTensor(1:3,1,iel))
          !v = GeomechanicalModel%Poisson(iel)
          v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
          CALL mccdProps(M,e0,P0,lambda,kappa,SUM(GeomechanicalModel%EffectiveStressTensor(:,1,iel))/3.0,&
                            DefVol,evol_vp,Pc,e,v)                    !<-update e and v.
          GeomechanicalModel%TransverIsotropicProps(iel)%E1 = e
          GeomechanicalModel%TransverIsotropicProps(iel)%E2 = e
          gamma = 1.0
          
          !TRESCA
          !ddt = two*(one+v)*(one-two*v)/(e*(one-two*v)) !Treca
          !ddt = ddt**3.0
          !IF(ddt<dt)dt=ddt
          
          !STOLLE.
          CALL D_elastic_homogenious(dee,e,v)
          
          DO i=1,nip
            CALL invar(GeomechanicalModel%EffectiveStressTensor(:,i,iel),&
                        p, q, lode_theta)
            Pc_int = Pc_intercept(M,p,q)
            f = mccf(M,Pc_int,p,q)
            HeHP = mcc_HeHp(M,Pc_int,e0,P0,lambda,kappa,dee,&
                            GeomechanicalModel%EffectiveStressTensor(:,i,iel))
            ddt = 1.0 / gamma / HeHp * 0.25 / f
            IF(ddt<dt)dt=ddt
          END DO
          
          !PUC
          !c1 = e/(1.0+v)/(1.0-2.0*v)
          !a = v*c1
          !b = (1.0-2.0*v)*c1
          !gamma = 1.0
          !DO i=1,nip
          !  CALL invar(GeomechanicalModel%EffectiveStressTensor(:,i,iel),&
          !              p, q, lode_theta) 
          !  f0 = 1.0  !mccf(m,pc,p,q)
          !  He = 2*b*q**2.0/((M*p)**4.0)+1.0/3.0*(3.0*a+b)*(-2.0*q**2.0/M**2.0/p**3.0+pc/p**2.0)**2.0
          !  Hp = pc/p * v / (lambda-kappa)*(-2.0*q**2.0/m**2.0*p**3.0+pc/p**2.0)
          !  ddt = 2.0 * f0 / gamma / (He + Hp)
            IF(ddt<dt)dt=ddt
          !END DO
          
        ELSE
        END IF
        IF(GeomechanicalModel%FracturedElementModel(iel) )THEN
          !e=GeomechanicalModel%YoungMod(iel)
          !v=GeomechanicalModel%Poisson(iel)
          e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
          v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
          ddt=two*(one+v)*(one-two*v)/(e*(one-two*v+one))
          IF(GeomechanicalModel%IsFractured(iel))THEN
            ddt=0.50*ddt
          END IF
          IF(ddt<dt)dt=ddt
        END IF
        
        GeomechanicalModel%Gamma(iel) = ddt
        
      END DO
      
      !CALIBRATING VISCOSITY PARAMETER (GAMMA)
      DO iel=1,GeomechanicalModel%NumElements
        IF(GeomechanicalModel%ConstitutiveGeoModel(iel)==2.OR.&
             GeomechanicalModel%ConstitutiveGeoModel(iel)==3)THEN
            GeomechanicalModel%Gamma(iel) = GeomechanicalModel%Gamma(iel) / dt
        END IF
      END DO
  
      RETURN
    END SUBROUTINE
                     
    ! Calcula tensor de deformação viscoplástica para o critério Mohr-Coulomb.
    SUBROUTINE MCDefVispl(i,prop,nst,sigma,tensor_s,f,erate)
      !Subtorina para verificar critério de ruptura de Mohr-Coulomb e 
      !calcular taxa de deformação viscoplástica.
  
      USE ConstitutiveGeoModels
      USE SolidMechanicsLib
  
      !VARIÁVEIS.
      IMPLICIT NONE
      INTEGER,INTENT(IN)::  i,&  !número do ponto de Gauss.
                            nst  !número de componentes de tensão no tensor de tensões (=6 para 3D).
      REAL(8):: sigm,dsbar,lode_theta,&  !invariantes de tensão.
                  dq1,dq2,dq3,&  !derivadas do potencial de escoamento de Mohr-Coulomb com relação aos três invariantes de tensão.
                  phi,c,&  !ângulo de atrito (graus) e coesão.
                  psi,&  !ângulo de dilatação.
                  zero=0.0
      REAL(8),INTENT(IN)::prop(23)  !vetor com propriedades do elemento (lido do arquivo .prp)
      REAL(8),INTENT(IN OUT)::f  !função de escoamento
      REAL(8),ALLOCATABLE:: m1(:,:),m2(:,:),m3(:,:),&  !derivada dos inveriantes com relação as tensões
                              flow(:,:)  !derivada do potencial plástico com relação as tensões
      REAL(8),ALLOCATABLE,INTENT(IN)::  sigma(:,:),&  !incremento do tensor de tensões
                                          tensor_s(:,:)  !tensor de tensões do estágio de carregamento anterior
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::  erate(:)  !taxa de deformação plástiva 
  
      !ALOCAÇÃO DE MEMÓRIA.
      ALLOCATE( m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst) )
      !>>>Verificar se a alocação acima não está piorando a performance.
      !>>>Caso contrário, declarar na rotina de equilíbrio estrutural.
      
      !PROPRIEDADES MECÂNICAS.
      c = prop(16)    !coesão (kPa)
      phi = prop(17)  !ângulo de atrito (graus)
      psi = prop(18)  !ângulo de dilatação (graus)
  
      !CÁLCULO DA FUNÇÃO DE ESCOAMENTO.
      CALL invar(sigma(:,i)+tensor_s(:,i),sigm,dsbar,lode_theta)
      CALL mocouf(phi,c,sigm,dsbar,lode_theta,f)
  
      !CÁLCULO DA TAXA DE ESCOAMENTO.
      !dQ/dsigma = dQ/dsigm*dsigm/ds + dQ/dJ2*dJ2/ds + dQ/dJ3*dJ3/ds
      !devp/td = f * dQ/dsigma
      IF(f>=zero)THEN
        CALL mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)  !dq1=dQ/dsigm, dq2=dQ/dJ2, dq3=dQ/dJ3
        CALL formm(sigma(:,i)+tensor_s(:,i),m1,m2,m3)  !m1=dsigm/ds, m2=dJ2/ds, m3=dJ3/ds
        flow = f * ( m1 * dq1 + m2 * dq2 + m3 * dq3 )
        erate = erate + MATMUL(flow,sigma(:,i)+tensor_s(:,i))  !f*dQ/dsigma
      END IF
  
    RETURN
    END SUBROUTINE MCDefVispl

                         
    ! Calcula tensor de deformação viscoplástica para o critério Cam-Ccay modificado.
    SUBROUTINE CCDefVispl(Sigma,MCCparameters,f,erate)
      !Subtorina para verificar critério de ruptura de Mohr-Coulomb e 
      !calcular taxa de deformação viscoplástica.
  
      USE ConstitutiveGeoModels
      USE SolidMechanicsLib
      USE DiscretizedModel
  
      !VARIÁVEIS.
      IMPLICIT NONE
      REAL(8),INTENT(IN):: Sigma(6)                 
      REAL(4),INTENT(IN OUT):: MCCparameters(8)       !Pc(kPa)  M  kappa  lambda   e0   P0(kPa)
      REAL(8),INTENT(IN OUT):: f                      !função de escoamento
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::  erate(:)  !taxa de deformação plástiva 
      
      REAL(8):: sigm,dsbar,lode_theta,&   !invariantes de tensão.
                dq1,dq2,dq3               !derivadas do potencial de escoamento de Mohr-Coulomb com relação aos três invariantes de tensão.
      REAL(8):: Pc, M, kappa, lambda, e0, P0, ec, evol_vp, ni
      REAL(8):: Pc_int, por
      INTEGER,PARAMETER:: nst = 6
      REAL(8):: m1(nst,nst),m2(nst,nst),m3(nst,nst),&   !derivada dos inveriantes com relação as tensões
                flow(nst,nst)                           !derivada do potencial plástico com relação as tensões
      
      !STRESS INCARIANTS
      CALL invar(Sigma,sigm,dsbar,lode_theta)
      
      !PROPRIEDADES MECÂNICAS.
      Pc = MCCparameters(1)
      M = MCCparameters(2)
      kappa = MCCparameters(3)
      lambda = MCCparameters(4)
      e0 = MCCparameters(5)
      P0 = MCCparameters(6)
      ec = MCCparameters(7)
      evol_vp = MCCparameters(8)
      
      !HARDENING.
      ni = 1+ec
      Pc = P0 * EXP( -evol_vp*ni/(lambda-kappa) )
      MCCparameters(1) = Pc
      ni = 1+e0
      ni = ni - lambda*LOG(Pc/P0)
      MCCparameters(7) = ec
      
      !YIELD.
      f = mccf(M,Pc,sigm,dsbar)
  
      !VISCOPLASTIC RATE.
      !dQ/dsigma = dQ/dsigm*dsigm/ds + dQ/dJ2*dJ2/ds + dQ/dJ3*dJ3/ds
      !devp/td = f * dQ/dsigma
      IF(f>0.0)THEN
        Pc_int = Pc_intercept(M,sigm,dsbar)
        CALL mcc_q(Pc_int,M,sigm,dsbar,dq1,dq2,dq3)             !dq1=dQ/dsigm, dq2=dQ/dJ2, dq3=dQ/dJ3
        CALL formm(Sigma,m1,m2,m3)                              !m1=dsigm/ds, m2=dJ2/ds, m3=dJ3/ds
        flow = f * ( m1 * dq1 + m2 * dq2 + m3 * dq3 )
        erate = erate + MATMUL(flow,Sigma)                      !f*dQ/dsigma
      END IF
      
write(16,*)sigm/1e3,dsbar/1e3,Pc/1e3
  
    RETURN
    END SUBROUTINE CCDefVispl
       
    ! Calcula tensor de deformação viscoplástica para o critério Traction.
    SUBROUTINE TracDefVispl(i,ST,nst,sigma,tensor_s,f,erate)
      !Subtorina para verificar critério de ruptura de Mohr-Coulomb e 
      !calcular taxa de deformação viscoplástica.
  
      USE ConstitutiveGeoModels
      USE SolidMechanicsLib
  
      !VARIÁVEIS.
      IMPLICIT NONE
      INTEGER,INTENT(IN)::  i,&  !número do ponto de Gauss.
                            nst  !número de componentes de tensão no tensor de tensões (=6 para 3D).
      REAL(8),INTENT(IN OUT):: f                        !função de escoamento
      REAL(8),INTENT(IN):: ST                           !traction strength
      REAL(8),ALLOCATABLE,INTENT(IN):: sigma(:,:),&     !incremento do tensor de tensões
                                       tensor_s(:,:)    !tensor de tensões do estágio de carregamento anterior
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::  erate(:)    !taxa de deformação plástiva 
      
      REAL(8):: sigm,dsbar,lode_theta,S3,&              !invariantes de tensão
                pi=ACOS(-1.0)
      REAL(8),ALLOCATABLE:: DQDS(:)                     !dQdStress.
      
      !FAILURE ENVELOPE.
      !CALL invar(sigma(:,i)+tensor_s(:,i),sigm,dsbar,lode_theta)
      S3 = sigma(2,i)+tensor_s(2,i) !sigm + 2.0/3.0 * dsbar * SIN( lode_theta + 2.0 * pi / 3.0 )      
      f = S3-ST
  
      !CÁLCULO DA TAXA DE ESCOAMENTO.
      !dQ/dsigma = dQ/dsigm*dsigm/ds + dQ/dJ2*dJ2/ds + dQ/dJ3*dJ3/ds
      !devp/td = f * dQ/dsigma
      IF( f>=0.0 )THEN
        ALLOCATE( DQDS(nst) )
        CALL tractionDQDS(sigma(:,i)+tensor_s(:,i),ST,DQDS)
        !DQDS = 0.0
        !DQDS(2) = 1.0
        erate = erate + f * DQDS
      END IF
  
    RETURN
    END SUBROUTINE TracDefVispl
    
    
    SUBROUTINE UpdatePorepressure(g_num,g_g,g_coord,Displacement,points,weights,GeomechanicalModel)
  
      USE DiscretizedModel
      USE ConstitutiveGeoModels
      USE SolidMechanicsLib
      USE FEMlib
      USE NumericalUtilities
      USE OMP_LIB  !bibliotecas OpenMP(R)
  
      IMPLICIT NONE
  
      !******************************************
      !VARIÁVEIS
      !******************************************
      
      !Input.
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel 
      INTEGER,INTENT(IN):: g_num(:,:),&   !conectividade dos elementos
                           g_g(:,:)       !matriz para identificação dos graus de liberdades dos nós
      REAL(8),INTENT(IN):: Displacement(0:),&
                           points(:,:),&  !coordenadas dos pontos de integração
                                       weights(:),&   !peso dos pontos de integração
                                       g_coord(:,:)   !coordenada dos nós
      
      REAL(8),ALLOCATABLE:: deltaPorepressure(:),Pp(:)
      INTEGER,ALLOCATABLE:: nsum(:)
      
      INTEGER:: ndim,&    !número de dimensões
                nels,&    !número de elementos na malha
                nod,&     !número de nós de deslocamento por elemento (=8 para hexaédro linear)
                nip,&     !número de pontos de integração
                ndof,&    !número de graus de liberdade de deslocamento por elemento
                nst       !número de componentes de tensão no tensor de tensões (=6 para 3D) 
            
      INTEGER:: iel,i
      REAL(8):: det,&     !determinante de matriz
                biotCoef, Kf
      
      INTEGER,ALLOCATABLE:: num(:) ,&                     !conectividade de um elemento
                            g(:)                          !vetor com os graus de liberdade de cada elemento
      REAL(8),ALLOCATABLE:: coord(:,:),&                  !coordenada dos nós de um elemento
                            jac(:,:),&                    !matriz jacobiana das funções de forma.
                            der(:,:),&                    !derivadas das funções de forma dos deslocamentos em coordenadas naturais.                        
                            deriv(:,:),&                  !derivadas das funções de forma com relação às coordenadas locais e globais
                            bee(:,:),&                    !matriz de derivadas das funções de forma
                            dee(:,:),&                    !matriz de relações elásticas
                            eld(:),&                      !deslocamentos nodais de um elemento
                            strain(:),&                   !termos do tensor de deformações
                            stress(:)                     !incremento do tensor de tensões
      
      nels = GeomechanicalModel%NumElements
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = GeomechanicalModel%NumIntegrationPoints
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement
      
      !Alocação de memória.
      ALLOCATE( g(ndof),num(nod),bee(nst,ndof),eld(ndof),der(ndim,nod),jac(ndim,ndim),&
                coord(nod,ndim),dee(nst,nst),deriv(ndim,nod),stress(nst),strain(nst) )
      
      ALLOCATE( deltaPorepressure(GeomechanicalModel%NumElements) )
      deltaPorepressure = 0.0
          
      !*************************************************************************
      ! UNDRAINED ELEMENTS - INCREMENT OF POREPRESSURE.
      !*************************************************************************
      
      ! Pore fluid stress (undrained solid elements).
      DO iel=1,GeomechanicalModel%NumElements
        IF( GeomechanicalModel%Undrained(iel) )THEN
        
          num = GeomechanicalModel%ElementNodes(:,iel)
          coord=TRANSPOSE(g_coord(:,num))
          g = GeomechanicalModel%ElementDOFs(:,iel)
          eld = Displacement(g)                                     !TOTAL DISPLACEMENT.
              
          biotCoef = GeomechanicalModel%BiotCoef(iel)
          Kf = GeomechanicalModel%PoreFluidBulkMod(iel)

          deltaPorepressure(iel) = 0.0
              
          DO i=1,nip
              
            ! FE matrices
            CALL shape_der(der,points,i)    !-> der
            jac=MATMUL(der,coord)           !-> jac
            CALL determinante(jac,det)      !-> det
            jac=MATMUL(der,coord)           !-> jac
            CALL invert(jac)                !-> inv(jac)
            deriv=MATMUL(jac,der)           !-> deriv
            CALL beemat(bee,deriv)          !-> bee
              
            ! Strain
            strain = MATMUL(bee,eld)                                !strain due TOTAL DISPLACEMENT.
              
            ! Constitutive model
            stress = FluidHydrostaticStress(Kf,biotCoef,strain)     !stress due TOTAL DISPLACEMENT.
            
            ! Undrained porepressure increment.
            deltaPorepressure(iel) = deltaPorepressure(iel) + SUM(stress(1:3))/3.0
                
          END DO
              
          deltaPorepressure(iel) = 1.0/nip * (  deltaPorepressure(iel) )
           
        END IF
      END DO
             
      ALLOCATE( nsum(GeomechanicalModel%NumNodes), Pp(GeomechanicalModel%NumNodes) )
      nsum = 0
      Pp = 0.0
      DO iel=1,GeomechanicalModel%NumElements
        IF( GeomechanicalModel%Undrained(iel) )THEN 
          num = GeomechanicalModel%ElementNodes(:,iel)
          Pp(num) = Pp(num) - deltaPorepressure(iel)
          nsum(num) = nsum(num) + 1
        END IF
      END DO      
      DO i=1,SIZE(Pp)
        IF( nsum(i)>0 ) Pp(i) = Pp(i) / nsum(i)
      END DO
      GeomechanicalModel%Porepressure = GeomechanicalModel%Porepressure + Pp
      
    RETURN
    END SUBROUTINE UpdatePorepressure
    
  
    SUBROUTINE PCG(cg_limit,neq,cg_tol,diag_precon,a,ia,ja,loads,initialGuess,RestrictedDOFs,nthreads,linSolverIters)
      !---------------------------------------------------------------------------------------------------------------
      !solver de sistema linear de equações com o método dos GRADIENTES CONJUGADOS pré-condionado diagonalmente (PCG).
      !rotina paralelizada para o OpenMP vrs. 3.0 .
      !uso de rotinas seriais do Intel Fortran Math Kernel Library - MKL
      !armazenamento de matrizes formato CSR Intel
      !
      !Referência: Ian M. Smith, D. V. Griffiths.
      !Professor Smith is currently the Head of the Geotechnics Division at the Manchester Centre for Civil and Construction Engineering, UMIST
      !Professor Griffiths is Professor of Civil Engineering, Division of Engineering at the Colorado School of Mines
      !---------------------------------------------------------------------------------------------------------------
      
      USE OMP_LIB  !bibliotecas OpenMP(R)
      USE ParallelRoutines
      USE IO
      
      IMPLICIT NONE
      
      REAL(8),INTENT(IN):: diag_precon(:),& !pré-condicionador do solver
                            a(:)  !vetor de elementos não zero da matriz A do sistema A*x=b
      REAL(8),INTENT(IN OUT)::loads(:)  !vetor de carregamento estrutural (IN) e resultados (OUT)
      REAL(8),INTENT(IN)::initialGuess(:)  !vetor com a solução inicial de entrada no método
      INTEGER,INTENT(IN)::neq,&  !número de equações
                          cg_limit  !número máximo de iterações do solver      
      REAL(8),INTENT(IN)::cg_tol  !tolerância para convergência do solver linear
      INTEGER,INTENT(IN)::ia(:),&  !identificação da posição de início das linhas no vetor a(:)
                          ja(:)  !identificação das colunas dos elementos não-zero da matriz
      INTEGER,INTENT(IN):: RestrictedDOFs(:)
      INTEGER,INTENT(IN):: nthreads
      INTEGER,INTENT(IN OUT):: linSolverIters
      
      !VARIÁVEIS ESTÁTICAS
      INTEGER::i,&  !contadores
                nParticoes,j,k,&  !variáveis auxilares para esquema de partição
                cg_iters,&  !contador de iterações do solver linear
                PCG_tempo(3,2)  !marcação de tempo do solver linear

      REAL(8)::alpha,beta,up,dot_p_u,dot_loads_d,&  !variáveis do processo do gradientes conjugados
                  convergence,&         !convergência
                  temp1,temp2,zero=0.0  !auxiliares
      REAL(8)::startTime,endTime  !variáveis para marcação de tempo gasto em segundos
      LOGICAL::cg_converged  !flag para convergência do solver
      
      !VARIÁVEIS DINÂMICAS
      INTEGER,ALLOCATABLE::particoes(:,:)  !guarda informações da partição dos vetores
      REAL(8),ALLOCATABLE::time(:)  !vetor que guarda o tempo computacional do produto de cada linha 
                                      !da matriz de rigidez (a) com o vetor "p". Para balanceamento de cargas do 
                                      !processamento paralelo.
      REAL(8),ALLOCATABLE::d(:),p(:),u(:),x(:),xnew(:),&  !variáveis dos gradientes conjudados
                              auxRealDb(:)  !auxiliar
      
      CHARACTER(LEN=320):: text
      
      !VARIÁVEIS DE TELEMETRIA
      real(8)::t1,t2,dt(7)
      
      ALLOCATE(d(neq),p(neq),u(neq),x(neq),xnew(neq),time(neq),auxRealDb(neq))
      
      !SOLUÇÃO DO SISTEMA DE EQUAÇÕES
      
      !marcação de tempo
      startTime = OMP_get_wtime()
      
      !INICIALIZAÇÃO
      CALL itime(PCG_tempo(:,1))
      
      !1a aproximação da solução
      x = initialGuess
      
      !Resíduo inicial da equação de equilíbrio: p = loads - [K]{u}
      auxRealDb=zero
      CALL mkl_dcsrgemv('N', neq, a, ia, ja, initialGuess, auxRealDb)  !MKL - rotina blas2 para multiplicação matriz esparsa por vetor.
      !CALL mkl_dcsrsymv('U', neq, a, ia, ja, initialGuess, auxRealDb)  !matriz simétrica.
      auxRealDb = loads(2:neq+1) - auxRealDb
      loads(2:neq+1) = auxRealDb
      loads(1+RestrictedDOFs) = 0.0
      
      !Aplicando pré-condicionador diagonal no resíduo da equação de equilíbrio
      d = zero
      CALL vdmul( neq, diag_precon(2:neq+1), loads(2:neq+1), d )  !MKL VMF
      d(RestrictedDOFs) = 0.0
      p = d
      
      !PARTICIONANDO A MATRIZ DE RIGIDEZ
      nParticoes = nthreads
      ALLOCATE(particoes(2,nthreads))
      particoes(2,nthreads)=neq
      !
      !particionando os graus de liberdade igualmente (vetor unitários).
      xnew = 1.0
      CALL PARTICIONAMENTO(xnew,particoes)
      xnew = 0.0    
      
      !LOOPING ITERATIVO++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      auxRealDb=zero
      cg_iters=0
      cg_converged=.TRUE.
      conjug_grad: DO
        cg_iters = cg_iters+1
        
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(i,j,k) SCHEDULE(STATIC)
        DO i=1,nParticoes
          !marcação dos graus de liberdade da partição
          j = particoes(1,i)
          k = particoes(2,i)  !grau de liberdade de início da próxima partição
          !u = k * p
          CALL SpAx( a, ia, ja, j, k, p, u )
        END DO
  !$OMP END DO
  !$OMP END PARALLEL
        u(RestrictedDOFs) = 0.0

        !up = DOT( loads, d )
        up = DOT_PRODUCT_pp( loads(2:neq+1), d, particoes)
        !
        !dot_p_u = DOT( p, u)
        dot_p_u = DOT_PRODUCT_pp( p, u, particoes)
        !
        !alpha
        alpha = up / dot_p_u
        !xnew = alpha * p
        CALL MULTaeb_pp(xnew,p,alpha,particoes)
        !xnew = xnew + x
        CALL ADICAOaaeb_pp(xnew,x,DBLE(1.0),particoes)
        !loads(2:neq+1) = loads(2:neq+1) - u * alpha
        CALL ADICAOaaeb_pp(loads(2:neq+1),u,DBLE(-alpha),particoes)
        !d = diag_precon(2:neq+1) * loads(2:neq+1)
        CALL MULTab_pp(d, diag_precon(2:neq+1), loads(2:neq+1),particoes)
        d(RestrictedDOFs) = 0.0
        !beta = DOT( dot_loads, d ) / up
        dot_loads_d = DOT_PRODUCT_pp( loads(2:neq+1), d, particoes)
        beta = dot_loads_d / up
        !p = d + p * beta
        CALL ADICAOabea_pp(p,d,beta,particoes)
        !verificando convergência
        CALL MAX_pp( xnew, x, particoes, temp2, temp1)
        convergence = temp1 / temp2
        cg_converged = (convergence<=cg_tol)
        !critério de parada
        IF( ( (cg_converged) .OR. (cg_iters==cg_limit)) .AND. (cg_iters>1) ) EXIT
        !x = xnew
        CALL dswap(neq, x, 1, xnew, 1)  !mkl  blas1 - troca os elementos entre os vetores / objetivo: x=xnew / apresentou melhor desempenho do que dcopy(neq, xnew, 1, x, 1)
      
      END DO conjug_grad
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !GUARDANDO RESULTADO
      loads(1)=0.0
      loads(2:(neq+1))=xnew
      
      !PRINT STUFF.
      !convergence = 1.0 - convergence
      !temp2 = MAXVAL( ABS(loads) )
      !WRITE(text,'(A,I6,A,E10.2,A,E10.2)')'...num. linear iter. :',cg_iters,' converg.:',convergence,' max.disp(m)',ABS(temp2)
      !CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
      
      !FINALIZANDO E ESCREVENDO RESULTADOS DA TELEMETRIA DO SOLVER
      IF (cg_iters==cg_limit)THEN
        CALL PrintMessage(.TRUE.,.TRUE.,"LINEAR SOLVER DID NOT CONVERGED...")
        CALL PrintMessage(.TRUE.,.TRUE.,"")
        STOP
      END IF
      
      !DESALOCANDO MEMÓRIA
      !DEALLOCATE(d,p,u,x,xnew)
      
      linSolverIters = cg_iters
      
    RETURN
    END SUBROUTINE

    
    SUBROUTINE CGpreConditioner(a,ia,ja, diag_precon)
      !-----------------------------------------------------------------------------
      !pre-condicionador para o solver iterativo gradientes conjugados
      !rotina paralelizada para o OpenMP vrs. 3.0.
      !Pré-condicionador tradicional do Smith e Griffiths
      !-----------------------------------------------------------------------------
  
      USE OMP_LIB  !bibliotecas OpenMP(R)
    
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::ia(:),&   !identificação da posição de início das linhas no vetor a(:)
                          ja(:)     !identificação das colunas dos elementos não-zero da matriz
      REAL(8),INTENT(IN)::a(:)      !vetor de elementos não zero da matriz A do sistema A*x=b
      REAL(8),ALLOCATABLE,INTENT(IN OUT)::diag_precon(:)  !pré-dondicionador.
    
      INTEGER:: i,&  !contador
                posicaoInicioLinha,posicao  !auxilares para posicionamento na matriz de rigidez
      INTEGER::neq                          !número de graus de liberdade na matriz de rigidez
      REAL(8):: one=1.0, zero=0.0

      neq = SIZE(ia)-1
                            
      ALLOCATE(diag_precon(0:neq))
      diag_precon = 0.0
    
      DO i=1,neq
  	    posicaoInicioLinha = ia(i)
  	    posicao = -1
  	    DO
  	      posicao = posicao + 1
  	      IF( ja(posicaoInicioLinha+posicao) == i ) EXIT
  	    END DO
  	    diag_precon(i) = one/( a(posicaoInicioLinha+posicao) )
  	  END DO
  	  
      diag_precon(0) = zero
  
    RETURN 
    END SUBROUTINE CGpreConditioner

    
    SUBROUTINE checon(loads,oldlds,tol,converged,convergence)
      ! This subroutine sets converged to .FALSE. if relative change in loads
      ! and oldlds is greater than tol and updates oldlds.
      
      IMPLICIT NONE
      REAL(8)::temp1,temp2
      REAL(8),INTENT(IN)::loads(0:),tol
      REAL(8),INTENT(IN OUT)::convergence
      REAL(8),INTENT(IN OUT)::oldlds(0:)
      LOGICAL,INTENT(OUT)::converged
      
      CONVERGED = .TRUE.
      temp1 = MAXVAL(ABS(loads-oldlds))
      temp2 = MAXVAL(ABS(loads))
      convergence = temp1/temp2
      CONVERGED = (convergence<=tol)
      oldlds = loads
      
    RETURN
    END SUBROUTINE checon
    
    
    SUBROUTINE PrepareViscoplasticSolver(Chanel_base,GeomechanicalModel,ActiveElements,&
                                            SolverParameters,Stiff)
      
      USE DiscretizedModel
      USE FEMlib
      USE StiffnessMatrix
      USE IO
      
      IMPLICIT NONE
      INTEGER,INTENT(IN):: Chanel_base
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      LOGICAL,INTENT(IN):: ActiveElements(:)
      TYPE(ViscoSolverParameters):: SolverParameters
      TYPE(ElasticStiffnessMatrix):: Stiff
      
      INTEGER:: iel
      LOGICAL:: found
      LOGICAL,ALLOCATABLE:: ElementsToAdd(:)
      
      !*************************************************************************
      !SOLVER PARAMETERS.
      !*************************************************************************
      
      CALL SearchInFile(Chanel_base,'#STRUCT_ViscoplasticSolver',found)
      IF(found) THEN
        READ(Chanel_base,*)SolverParameters%maxLinIters
        READ(Chanel_base,*)SolverParameters%linTol
        READ(Chanel_base,*)SolverParameters%maxViscoIters
        READ(Chanel_base,*)SolverParameters%viscoTol
      END IF
      
      !*************************************************************************
      !STIFFNESS MATRIX.
      !*************************************************************************
      
      IF( .NOT.ALLOCATED(Stiff%ia) )THEN
        CALL SparseStiffnessMatrixCSR(GeomechanicalModel,Stiff%a,Stiff%ia,Stiff%ja)
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
    

  END MODULE
  