
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE CreepMod
    
    CONTAINS
    
    !Creep solver
    SUBROUTINE CreepSolver( GeomechanicalModel,Matrices,&
                            StructuralBC,SolverParameters,ReportData,DeltaT,Load,TotalTime)
        
        USE OMP_LIB
        USE IO
        USE DiscretizedModel
        USE InterfaceMechanics
        USE ViscoplasticExplicitySolver
        USE StiffnessMatrix
        USE SolidMechanicsLib

        IMPLICIT NONE

        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE (SolverMatricesContainer),INTENT(IN OUT):: Matrices
        TYPE(InterfaceStiffnessVector),INTENT(IN OUT):: StructuralBC
        TYPE (ViscoSolverParameters),INTENT(IN):: SolverParameters
        TYPE(StepReport),INTENT(IN OUT):: ReportData
        REAL(8),INTENT(IN):: DeltaT(:),Load(:)
        REAL(8),INTENT(IN):: TotalTime

        REAL(8),ALLOCATABLE:: displacement(:)
        REAL(8),ALLOCATABLE:: CreepDevStrainRate(:,:)
        REAL(8),ALLOCATABLE:: NlinDevStrainRatio(:,:)        !non-linear deviatoric strain ratio: rate/cumulative maximum
        
        INTEGER:: i,iel,n
        INTEGER:: numSteps, subStep

        REAL(8):: cg_tol,&                         !tolerância do solver linear de equações - razão de máximos deslocamentos
                  tol                              !tolerância do solver viscoplástico - razão de máximas forças desbalanceadas
        INTEGER:: limit,&                          !limite de iterações para o solver viscoplástico
                  cg_limit                         !limete de iterações para o solver linear de equações
        
        REAL(8):: t,dt,dt_err,dt_const,elapsedTime,tal
        REAL(8):: maxYieldFunction
        REAL(8):: sigm,devStrain,lode_theta
        LOGICAL:: found,&
                  test=.TRUE.
        
        REAL(8),ALLOCATABLE:: RHS_load(:)          !RHS for structural equilibrium
        
        ALLOCATE( displacement(0:GeomechanicalModel%NumNodes*3),&
                  CreepDevStrainRate(GeomechanicalModel%NumIntegrationPoints,GeomechanicalModel%NumElements),&
                  NlinDevStrainRatio(GeomechanicalModel%NumIntegrationPoints,GeomechanicalModel%NumElements))
        displacement = 0.0
        CreepDevStrainRate = 0.0
        NlinDevStrainRatio = 0.0
        
        ALLOCATE( RHS_load(0:GeomechanicalModel%NumEquations) )
        
        !TRANSIENT SOLVER - POROMECHANICS + CREEP
        !Explicit Coupling
        dt = totalTime
        subStep = 0
        t = 0
        DO WHILE(t<totalTime)
            subStep = subStep + 1
            reportData%structuralStep = subStep
            RHS_load = 0.0
            
            !TIME STEP
            IF( subStep>1 )THEN
                !Error control criteria
                dt_err = MAXVAL(NlinDevStrainRatio)   !e_dev_rate / e_dev
                tal = dt_err * dt
                IF( tal>0.0 .AND. tal<0.10 )THEN  !.AND. GeomechanicalModel%CreepWork(3)>GeomechanicalModel%CreepWork(2)) THEN
                    dt = dt * 1.02
                    CALL CreepCriticalTimeStep(GeomechanicalModel,dt_const)
                    IF(dt_const>dt) dt = dt_const
                ELSE
                    CALL CreepCriticalTimeStep(GeomechanicalModel,dt_const)
                    dt = 0.50 * dt_const
                END IF
            ELSE
                dt = 0.0
            END IF
            IF( t+dt>totalTime ) dt = totalTime - t
            t = t + dt
            elapsedTime = elapsedTime + dt
            
            !CREEP LOAD
            IF( subStep>1 ) CALL SaltCreepLoad(GeomechanicalModel,CreepDevStrainRate,Dt,RHS_load,maxYieldFunction)
            
            !GENERAL LOADS
            IF( subStep==1) RHS_load = RHS_load + Load
            
            !STRUCTURAL SOLVER
            displacement = 0.0
            cg_tol = SolverParameters%linTol
            cg_limit = SolverParameters%maxLinIters
            tol = SolverParameters%viscoTol
            limit = SolverParameters%maxViscoIters
            CALL ViscoSolver(GeomechanicalModel,reportData,Matrices%Stiff%a, Matrices%Stiff%ia, &
                                Matrices%Stiff%ja,&
                                StructuralBC%ListOfRestrictedDOFs,RHS_load(1:),DeltaT,&
                                cg_tol,cg_limit,tol,limit,&
                                displacement)
            
            !Non-linear desviatory increment
            DO iel=1,GeomechanicalModel%NumElements
                DO i=1,GeomechanicalModel%NumIntegrationPoints
                    CALL invar(GeomechanicalModel%StrainTensor(:,i,iel),sigm,devStrain,lode_theta)
                    IF(devStrain>0) &
                        NlinDevStrainRatio(i,iel) = CreepDevStrainRate(i,iel) / devStrain     !tal/dt
                END DO
            END DO
            
        END DO
    
    END SUBROUTINE
    
    !Salt creep load.
    SUBROUTINE SaltCreepLoad(GeomechanicalModel,CreepStrainRate,Dt,Loads,MaxYieldFunction)
  
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
      REAL(8),INTENT(IN):: Dt
      REAL(8),INTENT(IN OUT):: Loads(0:),&  !Nodal Load
                               CreepStrainRate(:,:)
      REAL(8),INTENT(IN OUT):: MaxYieldFunction
      
      INTEGER,ALLOCATABLE:: g_num(:,:),&  !conectividade dos elementos
                             g_g(:,:)  !matriz para identificação dos graus de liberdades dos nós
      REAL(8),ALLOCATABLE:: points(:,:),&  !coordenadas dos pontos de integração
                            weights(:),&  !peso dos pontos de integração
                            g_coord(:,:)  !coordenada dos nós
      
      INTEGER:: ndim,&  !número de dimensões
                nels,&  !número de elementos na malha
                nod,&  !número de nós de deslocamento por elemento (=8 para hexaédro linear)
                nip,&  !número de pontos de integração
                ndof,&  !número de graus de liberdade de deslocamento por elemento
                nst,&  !número de componentes de tensão no tensor de tensões (=6 para 3D)
                neq,nn
            
      !Variáveis estáticas.
      INTEGER:: iel,i  !contadores
      REAL(8):: det,&  !determinante de matriz
                Ea,Ec,niac,niab,&  !parâmetros do modelo ortotrópico (ver interface no módulo MODULE iMaterials)
                zero=0.0,&
                pi=ACOS(-1.)
      REAL(8):: sigm, dsbar, theta  !invariantes
      REAL(8):: f
  
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
                            eld(:),&                      !deslocamentos nodais de um elemento
                            eload(:),&                    !contribuição de um ponto de integração para o vetor global bload
                            erate(:),&                    !taxa de deformação de fluência
                            erate_primary(:),&            !taxa de deformação de fluência primária
                            evp(:),&                      !incremento de deformação plástica
                            evpt(:,:),&                   !deformação plástica acumulada
                            eps(:,:),&                    !termos do tensor de deformações totais (elástica + plástica)
                            sigma(:,:),&                  !incremento do tensor de tensões
                            tensor_s(:,:),&               !tensor de tensões do estágio de carregamento anterior
                            anisotropia_Ec(:),&           !vetor de anisotropia Ec (=Ec/Ea), onde Ea está no arquivo .prp.
                            anisotropia_niac(:)           !vetor de anisotropia niab (=niab/niac), onde niac está no arquivo .prp.

      INTEGER:: modelNumber
      REAL(8):: temp
      
      neq = GeomechanicalModel%NumEquations
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = GeomechanicalModel%NumIntegrationPoints
      nn = GeomechanicalModel%NumNodes
      ndof = GeomechanicalModel%NumDOFUperElement
      nels = GeomechanicalModel%NumElements
      nod = GeomechanicalModel%NumNodesPerElement
          
      !Alocação de memória.
      ALLOCATE( g(ndof),num(nod),bee(nst,ndof),eld(ndof),bload(ndof),                   &
                coord(nod,ndim),dee(nst,nst),deriv(ndim,nod),devp(nst),                 &
                eload(ndof),erate(nst),evp(nst),sigma(nst,nip),eps(nst,nip),            &
                tensor_s(nst,nip),evpt(nst,nip),erate_primary(nst) )
      
      ALLOCATE(der(ndim,nod),jac(ndim,ndim))
  
      ALLOCATE(g_coord(ndim,nn),points(nip,ndim),weights(nip))
      
      g_coord = 0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
      
      CALL sample(GeomechanicalModel%element,points,weights)
      
      !******************************************
      !CÁLCULO DE DEFORMAÇÕES, TENSÕES E CARGA
      !VISCOPLÁSTICA DESBALANCEADA.
      !******************************************
      GeomechanicalModel%CreepWork(1) = zero
      GeomechanicalModel%CreepStrainIncrementTensor = 0.0
      bee=zero
      evp=zero
  
      !Element looping for stress calculation.
      DO iel=1,nels
    
          IF( GeomechanicalModel%ActiveElements(iel) )THEN
              
              !FE incidence.
              num = GeomechanicalModel%ElementNodes(:,iel)    !g_num(:,iel)
              coord = TRANSPOSE(g_coord(:,num))
              g = GeomechanicalModel%ElementDOFs(:,iel)       !g_g(:,iel)
                      
              dee = D_TransverseIsotropic(GeomechanicalModel%TransverIsotropicProps(iel)%E1,&
                                              GeomechanicalModel%TransverIsotropicProps(iel)%E2,&
                                              GeomechanicalModel%TransverIsotropicProps(iel)%n1,&
                                              GeomechanicalModel%TransverIsotropicProps(iel)%n2,&
                                              GeomechanicalModel%TransverIsotropicProps(iel)%strike,&
                                              GeomechanicalModel%TransverIsotropicProps(iel)%dip)
    
              tensor_s = GeomechanicalModel%EffectiveStressTensor(:,:,iel)    !efetive stress tensor
    
              !Gauss Points Looping.
              bload = zero
              eps = zero
              
              DO i=1,nip
                
                !Derivative of shape functions [B].(*)
                CALL shape_der(der,points,i)  !derivada das funções de forma
                jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integração i
                CALL determinante(jac,det)  !determinante do jacobiano
                CALL invert(jac)  !inverso do jacobiano
                deriv=MATMUL(jac,der)  !derivadas das funções de forma em coordenadas reais
                CALL beemat(bee,deriv)  !Matriz de derivadas B

                !************************************************************************************
                !VISCOPLASTIC CALCULATION - WIPP MODEL
                erate = 0.0
                CALL SaltCreepDefVispl( GeomechanicalModel, iel, GeomechanicalModel%Temperature(iel), &
                                          tensor_s(:,i), F,  erate)
                
                !Acumulando Def.Plas. total
                evp = erate * dt                                          !creep strain increment {DeltaDeformacaoPlastica} = {TaxaDeformacaoViscoplastica} * Dt
                GeomechanicalModel%CreepStrainIncrementTensor(:,i,iel) = evp
                GeomechanicalModel%CreepStrainTensor(:,i,iel) = GeomechanicalModel%CreepStrainTensor(:,i,iel) + evp
                GeomechanicalModel%PlasticStrainTensor(:,i,iel) = GeomechanicalModel%PlasticStrainTensor(:,i,iel) + evp
                
                !Umbalanced stress due evp
                devp = MATMUL(dee,evp)                                    !DeltaTensão devido o DeltaDeformacaoPlastica (TENSÃO DESBALANCEADA).
                
                !Saving creep strain rate
                CALL invar(erate,sigm,CreepStrainRate(i,iel),theta)

                !Integration of umbalanced load "eload" ( {b} = devp*[B] ).
                eload = MATMUL(devp,bee)               !transposto[B] * {DeltaTensão}
                bload = bload + eload*det*weights(i)   !integral(BT*Delta.Tensão)
                
                !Creep work
                GeomechanicalModel%CreepWork(1) = GeomechanicalModel%CreepWork(1) + &
                    CreepStrainRate(i,iel)*Dt * det * weights(i)      !int(creep_strainDev)
                
              END DO  !Gauss points

              !Atualização das forças desbalanceadas globais com a contribuição do elemento------------
              Loads(g) = Loads(g) + bload
    
          END IF !GeomechanicalModel%ActiveElements(iel)
      END DO !elements***************************************************************************
      
      GeomechanicalModel%CreepWork(2) = GeomechanicalModel%CreepWork(3)
      GeomechanicalModel%CreepWork(3) = GeomechanicalModel%CreepWork(1)
    
      Loads(0) = zero
                                   
    RETURN
    END SUBROUTINE

    
    SUBROUTINE SaltCreepDefVispl( GeomechanicalModel, ElementNumber, Temperature, tensor_s, &
        F,  erate)
  
        USE DiscretizedModel
        USE ConstitutiveGeoModels
        USE SolidMechanicsLib
  
        !VARIÁVEIS.
        IMPLICIT NONE
        TYPE(Model),INTENT(IN):: GeomechanicalModel
        INTEGER,INTENT(IN):: ElementNumber
        REAL(8),INTENT(IN)::  tensor_s(:)                       !tensor de tensões do estágio de carregamento anterior
        REAL(8),INTENT(IN OUT):: Temperature,F                  !função de escoamento
        REAL(8),INTENT(IN OUT)::  erate(:)                      !taxa de deformação de fluência total
        
        REAL(8):: sigm,dsbar,lode_theta,&       !invariantes de tensão.
                  dq1,dq2,dq3,&                 !derivadas do potencial de escoamento de Mohr-Coulomb com relação aos três invariantes de tensão.
                  phi,c,&                       !ângulo de atrito (graus) e coesão.
                  psi,&                         !ângulo de dilatação.
                  zero=0.0
        REAL(8),ALLOCATABLE:: m1(:,:),m2(:,:),m3(:,:),&     !derivada dos inveriantes com relação as tensões
                              flow(:,:)                     !derivada do potencial plástico com relação as tensões
        INTEGER::  i,&          !número do ponto de Gauss.
                   nst          !número de componentes de tensão no tensor de tensões (=6 para 3D).
                
        REAL(8):: MU0,Qc,N2,T_Celsius,MU2
        REAL(8):: ep_rate
        REAL(8):: R=8.34145
        
        !ALOCAÇÃO DE MEMÓRIA.
        nst = GeomechanicalModel%NumTensorComponents
        ALLOCATE( m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst) )
        
        !PROPRIEDADES MECÂNICAS.
        c = 0.0
        phi = 0.0
        psi = 0.0
  
        !CÁLCULO DA FUNÇÃO DE ESCOAMENTO.
        CALL invar(tensor_s,sigm,dsbar,lode_theta)
        F = dsbar

        !CÁLCULO DA TAXA DE ESCOAMENTO.
        !dQ/dsigma = dQ/dsigm*dsigm/ds + dQ/dJ2*dJ2/ds + dQ/dJ3*dJ3/ds
        !devp/td = f * dQ/dsigma
        IF(f>zero)THEN
          CALL mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)  !dq1=dQ/dsigm, dq2=dQ/dJ2, dq3=dQ/dJ3
          CALL formm(tensor_s,m1,m2,m3)  !m1=dsigm/ds, m2=dJ2/ds, m3=dJ3/ds
          
          !Flow function
          i = GeomechanicalModel%ElementCreep(ElementNumber)
          IF(i>0)THEN
              
              !secondary or steady state
              MU0 = GeomechanicalModel%CreepData(i)%MU0
              Qc = GeomechanicalModel%CreepData(i)%Qc
              N2 = GeomechanicalModel%CreepData(i)%N2
              T_Celsius = Temperature
              ep_rate = MU0 * F**N2 * EXP(-Qc/R/Temp_Kelvin(T_Celsius))
              
              !3D stress/deformation space flow
              !total creep rate
              flow = ep_rate * ( m1 * dq1 + m2 * dq2 + m3 * dq3 )
              erate = erate + MATMUL(flow,tensor_s)                         !f*dQ/dsigma
          ELSE
              erate = 0.0
          END IF
          
        END IF
  
    RETURN
    END SUBROUTINE

    
    !Critical time step.    
    SUBROUTINE CreepCriticalTimeStep(GeomechanicalModel,TimeStep)
        
        USE DiscretizedModel
    
        IMPLICIT NONE
        TYPE(Model),INTENT(IN):: GeomechanicalModel
        REAL(8),INTENT(IN OUT):: TimeStep
        
        INTEGER:: iel,ip,i,nip
        REAL(8):: dt,&
                  phi,snph,&                                    !ângulo de atrito e seno
                  psi,snps,&                                    !ângulo de dilatação e seno
                  e,v,&                                         !módulo de young e coficiente de poisson
                  pi,one=1.0,two=2.0,d180=180.0
        REAL(8):: f,fcr,g
        REAL(8):: a(2)
        REAL(8):: MU0,Qc,N2,T_Celsius,MU2
        REAL(8):: R=8.34145
        
        nip = GeomechanicalModel%NumIntegrationPoints
        
        TimeStep = 1.0E30
        pi=ACOS(-one)
        DO iel=1,GeomechanicalModel%NumElements
            i=GeomechanicalModel%ElementCreep(iel)
            IF( GeomechanicalModel%ActiveElements(iel) .AND. i>0 )THEN
                !MC with c=0, phi=0 and qsi=0
                phi = TAN(0.0)
                phi = ATAN(phi) * 180./pi
                snph = SIN(phi*pi/d180)
                psi = 0.0
                snps = SIN(psi*pi/d180)
                a(1) = GeomechanicalModel%TransverIsotropicProps(iel)%E1
                a(2) = GeomechanicalModel%TransverIsotropicProps(iel)%E2
                e = MAXVAL(a)
                a(1) = GeomechanicalModel%TransverIsotropicProps(iel)%n1
                a(2) = GeomechanicalModel%TransverIsotropicProps(iel)%n2
                v = MINVAL(a)
                g = e / 2.0 / (1.0+v)
                
                !Yield function
                MU0 = GeomechanicalModel%CreepData(i)%MU0
                Qc = GeomechanicalModel%CreepData(i)%Qc
                N2 = GeomechanicalModel%CreepData(i)%N2
                T_Celsius = GeomechanicalModel%Temperature(iel)
                
                DO ip=1,nip
                    fcr = CreepFunction(GeomechanicalModel%EffectiveStressTensor(ip,ip,iel))
                    IF( fcr<=0.0 ) fcr = 1.0
                    !dt = 4.0 * ( 1.0 + v)/( 3.0 * e * MU0 * EXP(-Qc/R/Temp_Kelvin(T_Celsius)) * fcr**(N2-1) ) * N2      !critério 1 de Stolle 1991, Eq. (32)
                    dt = EXP(Qc/R/Temp_Kelvin(T_Celsius))/(g*MU0*fcr**(N2-1))                                            !FLAC
                    IF( dt < TimeStep ) TimeStep = dt
                END DO
            END IF
        END DO
        
    RETURN
    END SUBROUTINE
    
    SUBROUTINE Lamb(MU0,Qc,T_kelvin,MU2)
        
        IMPLICIT NONE
        REAL(8),INTENT(IN):: MU0,Qc,T_kelvin
        REAL(8),INTENT(IN OUT):: MU2
    
        REAL(8):: R=8.34145
        
        MU2 = MU0 * EXP( -Qc / R / T_kelvin )
        
    END SUBROUTINE
    
    FUNCTION Temp_Kelvin(Temp_Celsius)
    
        IMPLICIT NONE
        REAL(8):: Temp_Kelvin
        REAL(8),INTENT(IN):: Temp_Celsius
        
        Temp_Kelvin = Temp_Celsius + 273.15
        
    END FUNCTION
    
    FUNCTION CreepFunction(EffectiveStress)
    
        USE SolidMechanicsLib
        
        IMPLICIT NONE
        REAL(8):: CreepFunction
        REAL(8),INTENT(IN):: EffectiveStress(6)
        
        REAL(8):: sigm,dsbar,lode_theta
        
        CALL invar(EffectiveStress,sigm,dsbar,lode_theta)
        
        CreepFunction = dsbar
        
    END FUNCTION
    
  END MODULE
  