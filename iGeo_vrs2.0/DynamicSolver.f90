!*****************************************************************************
!COMPUTATIONAL GEOMECHANICS
!Eng. Luis Carlos de Sousa Junior
!Computational Geomechanics
!*****************************************************************************

!-----------------------------------------------------------------------------
! #MEDSolver - Multiphysics Explicit Dynamic Explicit Solver
!-----------------------------------------------------------------------------

MODULE DynamicSolver

    TYPE DynamicSolverParameters
        REAL(8):: tol
        REAL(8):: dampFractor
        REAL(8):: timeStep
        INTEGER:: nsteps
        REAL(8):: MaxDisplacementControl = 10.0
        INTEGER:: SolverConvergence                     !0 - success
                                                        !1 - not converged
                                                        !2 - failed
                                                        !3 - return due excess of displacement.
                                                        !4 - return due excess of plastic shear strain.
    END TYPE

    TYPE PoreFlowModel
      ! ELEMENTS DATA
      REAL(8),ALLOCATABLE:: FlowElementPoreVolume(:),&           !current pore volume
                            FlowElementTransmissibilities(:,:),& !I,J,K transmissibilities
                            FlowElementKrel(:,:)                 !phases relative perms
      ! FLOW INTERFACES
      INTEGER,ALLOCATABLE:: FlowElementIndex(:),&        !index for each element in the grid
                            FlowInterfaceNodes(:,:),&     !four nodes interface
                            FlowInterfaceElements(:,:)    !adjacent flow elements
      REAL(8),ALLOCATABLE:: FlowInterfaceArea(:),&        !area normal to flow
                            FlowInterfaceNormal(:,:)      !normal vector to area
    END TYPE
    
    TYPE DisplacementBC
        INTEGER,ALLOCATABLE:: DOFs(:)                     !DOFs with displacement function
        REAL(8),ALLOCATABLE:: Displacement(:)             !displacement at dofs listed at DOFs(:)
        REAL(8),ALLOCATABLE:: RampFunction(:)             !ramp function for each time step
        INTEGER,ALLOCATABLE:: ListYnodalRestriction(:)    !list of nodes to be restricted
                                                          !for negative Y displacement
    END TYPE

CONTAINS
  
    SUBROUTINE DynamicSolverHead(LogChanel)
        USE IO

        IMPLICIT NONE
        INTEGER,INTENT(IN):: LogChanel

        CALL PrintMessage(.TRUE.,.TRUE.,'...EXPLICIT NONLINEAR DYNAMIC SOLVER' )
        
    END SUBROUTINE

    SUBROUTINE RunDynamicSolver(GeomechanicalModel,Matrices,SolverParameters,&
                                RestrictedDOFs,DispBC,Fext,x1,&
                                OutputFreq,ElementOutput,NodeOutput,&
                                LogChanel,ActiveElements)
  
        USE DiscretizedModel
        USE NumericalUtilities
        USE FEMlib
        USE ConstitutiveGeoModels
        USE HydraulicFractureModel
        USE Damage
        USE IO
        USE SolidMechanicsLib
        USE StiffnessMatrix
        USE OMP_LIB  !bibliotecas OpenMP(R)
  
        IMPLICIT NONE
  
        !Input.
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE(SolverMatricesContainer),INTENT(IN OUT):: Matrices
        TYPE(DynamicSolverParameters):: SolverParameters
        INTEGER,INTENT(IN OUT):: RestrictedDOFs(0:)
        TYPE(DisplacementBC):: DispBC
        REAL(8),INTENT(IN):: Fext(:)
        LOGICAL,INTENT(IN OUT):: ActiveElements(:)
        REAL(8),INTENT(IN OUT):: x1(0:)                               !initial displacement.
        
        INTEGER:: OutputFreq,&
                  NodeOutput, &
                  ElementOutput
        INTEGER,INTENT(IN):: LogChanel
        
        TYPE(Damping):: Damp
        REAL(8),ALLOCATABLE:: dampMultiplier(:)
  
        INTEGER:: i,j,k,n,q,w,counter,i1,i2,dof,step,iel,numElements,numEquations,numNodes,  &
                  numElasticElements,&          
                  numDamageElements,&                                 !number of damage elements.
                  numMohrCoulombElements,&
                  numCAMCLAYElements,&
                  numHyperbolicElements,&
                  numSandModel,&
                  numFracSpringDOFs,&
                  convergCounter,convergFreq,&                        !convergence check in number of steps.
                  fracCounter,fracFreq,&                              !fracture propagation counter.
                  nsteps                                              !total number of transient steps.
        REAL(8):: det,dtim,dt0,one=1.0,pt5=0.5,sigm,time,zero=0.0
        REAL(8):: e,v,biotCoef,Kf
        
        REAL(8):: max,maxOld
        
        REAL(8):: factor = 0.100
        
        CHARACTER(LEN=15)::element='hexahedron'
        
        INTEGER,ALLOCATABLE:: g(:),g_g(:,:),g_num(:,:),nf(:,:),num(:),freeDOFs(:)
        
        REAL(8),ALLOCATABLE:: bdylds(:),bee(:,:),bload(:),coord(:,:),dee(:,:),           &
                              der(:,:),deriv(:,:),d1x1(:),d2x1(:),eld(:),eload(:),       &
                              diag(:),g_coord(:,:),jac(:,:),                             &
                              points(:,:),weights(:),xold(:),xnew(:)
        
        REAL(8),ALLOCATABLE:: currentStress(:,:,:),&       
                              currentStrain(:,:,:),&
                              currentPlasticStrain(:,:,:),&
                              stress(:),stressIncrement(:),&            !stress.
                              stressPP(:),&                             !fluid stress. 
                              strEP(:),strain(:),strainIncrement(:),&   !strain.
                              str(:,:,:)
        
        REAL(8),ALLOCATABLE:: deltaPorepressure(:)
        
        REAL(8),ALLOCATABLE:: fint(:)                                    !internal force.
          
        REAL(8),ALLOCATABLE:: detJac(:,:)                                !FE element matrix.
        !REAL(8),ALLOCATABLE:: B(:,:,:,:)                                !FE element matrix.
                  
        INTEGER,ALLOCATABLE:: ListOfElasticElements(:),&           
                              ListOfMohrCoulombElements(:),&           
                              ListOfDamageElements(:),&          
                              ListOfCAMCLAYElements(:),&           
                              ListOfUndrainedElements(:),&           
                              ListOfHyperbolicElements(:),&          
                              ListOfSandModel(:)           
                  
        INTEGER:: ndim = 3                                                !number of dimensions: 3D -> 3.
        INTEGER:: nod = 8                                                 !number of nodes per element.
        INTEGER:: ndof = 24                                               !number of displacement dofs per element.
        INTEGER:: nodf = 8                                                !number of fem dofs per element.
        INTEGER:: nodof = 3                                               !number of displacements per node.
        INTEGER:: nst = 6                                                 !number of stress tensor.
        INTEGER:: nip                                                     !number of integration points.
        
        CHARACTER(LEN=320):: text
        
        INTEGER:: dofMonitor,&                                            !DOF to be printed into the screem.
                  outputCounter
        REAL(8),ALLOCATABLE:: nodeDispOutput(:,:),&
                              elementStressOutput(:,:),&
                              elementPorepressureOutput(:),&
                              elementStrainOutput(:,:)
        
        INTEGER:: numDispDofs
        
        REAL(8),ALLOCATABLE:: rampFunction(:)
        
        REAL(8):: p, area, syy,  fy, deltaForce, DamageFracCutoff=0.80
        
        REAL(8):: startTime,endTime
        
        LOGICAL:: converged
        
        REAL(8):: c,phi,dil
        REAL(8):: temp
        
        REAL(8):: dampFractor
        
        REAL(8):: M,Pc,Pt,A,P0,e0,lambda,kappa,defVol
        
        INTEGER,ALLOCATABLE:: nsum(:)
        REAL(8),ALLOCATABLE:: Pp(:)
        
        LOGICAL:: doNotStop
        
        LOGICAL:: updateGeometry
        
        REAL(8):: xtol
        
        REAL(8):: Strength(3)
        INTEGER:: modelNumber
        
        REAL(8):: T0,MaxPlasticDilation,ST
        
        REAL(8):: syy_avg
        
        LOGICAL:: ElastifStiffness
        
        REAL(8):: e_ref, p_ref
        
        REAL(8):: eps
        REAL(8):: Cp
        REAL(8):: Cr
        !REAL(8):: Phi
        REAL(8):: Phir
        REAL(8):: DilAngle
        
        REAL(8):: rho
        
        REAL(8):: kinecticEnergy(10),maxKinecticEnergy,kDer(10)

        REAL(8):: dsbar, lode_theta, f
        
        !*************************************************************************
        ! INITIALIZING.
        !*************************************************************************
        CALL DynamicSolverHead(LogChanel)
        
        numElements = GeomechanicalModel%NumElements
        numNodes = GeomechanicalModel%NumNodes
        nip = GeomechanicalModel%NumIntegrationPoints
        nsteps = SolverParameters%nsteps
        
        numEquations = GeomechanicalModel%NumEquations
        
        ALLOCATE( nf(nodof,numNodes),points(nip,ndim),weights(nip),g_coord(ndim,numNodes),      &
                  num(nod),dee(nst,nst),currentStress(nst,nip,numElements),coord(nod,ndim),     &
                  jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),                                 &
                  g_num(nod,numElements),bee(nst,ndof),eld(ndof),strain(nst),                   &
                  bload(ndof),eload(ndof),g(ndof),stress(nst),stressPP(nst),                    &
                  stressIncrement(nst),strainIncrement(nst),                                    &
                  g_g(ndof,numElements),currentStrain(nst,nip,numElements),                     &
                  currentPlasticStrain(nst,nip,numElements),strEP(nst),                         &
                  str(nst,nip,numElements))
        currentPlasticStrain = 0.0
        g_coord = 0.0
        g_coord(1,:) = GeomechanicalModel%X
        g_coord(2,:) = GeomechanicalModel%Y
        g_coord(3,:) = GeomechanicalModel%Z
        
        !LIST OF FREE STRUCTURAL DOFS
        j = 0
        DO i=1,RestrictedDOFs(0)
          IF(RestrictedDOFs(i).NE.0) j=j+1
        END DO
        ALLOCATE( freeDOFs(numEquations-j) )
        freeDOFs = 0
        j=0
        DO i=1,numEquations
          IF(.NOT.ANY(RestrictedDOFs(1:RestrictedDOFs(0))==i))THEN
            j=j+1
            freeDOFs(j) = i
          END IF
        END DO
        
        !INITIALIZING VARIABLES
        ALLOCATE(bdylds(0:numEquations),d1x1(0:numEquations),&
                  d2x1(0:numEquations),diag(0:numEquations),fint(0:numEquations),&
                  xold(0:numEquations),xnew(0:numEquations))
        diag = 0.0
        diag = Matrices%Mass%Values
        nf = GeomechanicalModel%NodeDofs
        fint = 0.0
        xold = 0.0
        xnew = 0.0
        
        CALL sample(element,points,weights)
        
        currentStress = 0.0
        currentPlasticStrain = GeomechanicalModel%PlasticStrainTensor
        currentStrain = 0.0
        
        ALLOCATE(nodeDispOutput(3,INT(nsteps/outputFreq)+1),&
                  elementStressOutput(9,INT(nsteps/outputFreq)+1),&
                  elementPorepressureOutput(INT(nsteps/outputFreq)+1),&
                  elementStrainOutput(6,INT(nsteps/outputFreq)+1))
        nodeDispOutput = 0.0
        elementStressOutput = 0.0
        elementPorepressureOutput = 0.0
        elementStrainOutput = 0.0
        
        numDispDofs = SIZE(DispBC%Dofs)
        
        ALLOCATE(deltaPorepressure(numElements))
        deltaPorepressure = 0.0

        kinecticEnergy=0.0
        kDer=0.0

        !*************************************************************************
        ! STRUCTURAL MODEL
        ! list of elements for each model type
        !*************************************************************************
        
        !-CONSTITUTIVE MODELS-----------------------------------------------------
        numElasticElements = COUNT(GeomechanicalModel%ConstitutiveGeoModel(1:)==0,1)
        numDamageElements = COUNT(GeomechanicalModel%ConstitutiveGeoModel(1:)==1,1)
        numMohrCoulombElements = COUNT(GeomechanicalModel%ConstitutiveGeoModel(1:)==2,1)
        numCAMCLAYElements = COUNT(GeomechanicalModel%ConstitutiveGeoModel(1:)==3,1)
        numHyperbolicElements = COUNT(GeomechanicalModel%ConstitutiveGeoModel(1:)==4,1)
        numSandModel = COUNT(GeomechanicalModel%ConstitutiveGeoModel(1:)==5,1)
        ALLOCATE( ListOfElasticElements(0:numElasticElements),&
                  ListOfDamageElements(0:numDamageElements),&
                  ListOfMohrCoulombElements(0:numMohrCoulombElements),&
                  ListOfCAMCLAYElements(0:numCAMCLAYElements),&
                  ListOfHyperbolicElements(0:numHyperbolicElements),&
                  ListOfSandModel(0:numSandModel))
        ListOfElasticElements = 0
        ListOfDamageElements = 0
        ListOfElasticElements = 0
        ListOfCAMCLAYElements = 0
        ListOfHyperbolicElements = 0
        ListOfSandModel = 0
        i = 0
        j = 0
        k = 0
        n = 0
        q = 0
        w=0
        DO iel=1,NumElements
          IF(ActiveElements(iel))THEN
            IF( GeomechanicalModel%ConstitutiveGeoModel(iel)==0 )THEN
              i=i+1
              ListOfElasticElements(i) = iel
            ELSEIF( GeomechanicalModel%ConstitutiveGeoModel(iel)==1 )THEN
              j=j+1
              ListOfDamageElements(j) = iel
            ELSEIF( GeomechanicalModel%ConstitutiveGeoModel(iel)==2 )THEN
              k=k+1
              ListOfMohrCoulombElements(k) = iel
            ELSEIF( GeomechanicalModel%ConstitutiveGeoModel(iel)==3 )THEN
              n=n+1
              ListOfCAMCLAYElements(n) = iel
            ELSEIF( GeomechanicalModel%ConstitutiveGeoModel(iel)==5 )THEN
              q=q+1
              ListOfSandModel(q) = iel
            ELSEIF( GeomechanicalModel%ConstitutiveGeoModel(iel)==4 )THEN
              w=w+1
              ListOfHyperbolicElements(w) = iel
            ELSE
            END IF
          END IF
        END DO
        
        ListOfElasticElements(0) = i
        ListOfDamageElements(0) = j
        ListOfMohrCoulombElements(0) = k
        ListOfCAMCLAYElements(0) = n
        ListOfSandModel(0) = q
        ListOfHyperbolicElements(0) = w
        
        ElastifStiffness = .FALSE.
        IF(i>0) ElastifStiffness = .TRUE.
        
        !-UNDRAINED ELEMENTS------------------------------------------------------

        i = COUNT(GeomechanicalModel%Undrained,1) + q
        ALLOCATE( ListOfUndrainedElements(0:i) )
        i = 0
        DO iel=1,NumElements
          IF( GeomechanicalModel%Undrained(iel) )THEN
            i = i + 1
            ListOfUndrainedElements(i) = iel
          END IF
        END DO
        ListOfUndrainedElements(0) = i

        !-MATRICES----------------------------------------------------------------

        ALLOCATE(Damp%Values(numEquations),dampMultiplier(0:numEquations))
        Damp%Values = 0.0
        DO i=1,numEquations
          Damp%Values(i) = Matrices%Mass%Values(i) * SolverParameters%dampFractor
        END DO
        dampMultiplier = 1.0

        !FE MATRICES FOR LARGE DISPLACEMENT APPLICATIONS
        !>>>ALLOCATE( B(nst,ndof,nip,numElements) )
        ALLOCATE( detJac(nip,numElements) )
        DO iel=1,numElements
          num = GeomechanicalModel%ElementNodes(:,iel)
          coord=TRANSPOSE(g_coord(:,num))
          
          DO i=1,nip
            
            CALL shape_der(der,points,i)      !-> der
            jac=MATMUL(der,coord)             !-> jac
            CALL determinante(jac,det)        !-> det 
            detJac(i,iel) = det               !<- det
            
            !CALL invert(jac)
            !deriv=MATMUL(jac,der)
            !CALL beemat(bee,deriv)
            !>>>B(:,:,i,iel) = bee
      
          END DO
          
        END DO

        !*************************************************************************
        ! FLOW MODEL
        ! list of elements for each model type
        !*************************************************************************

        
        !*************************************************************************
        ! THERMAL MODEL
        ! list of elements for each model type
        !*************************************************************************


        !*************************************************************************
        ! LOAD RAMP FUNCTION.
        !*************************************************************************
        
        ALLOCATE(rampFunction(SolverParameters%nsteps))
        rampFunction = 0.0
        
        DO step=1,SolverParameters%nsteps
          IF(step<factor*nsteps)THEN
            rampFunction(step) = step / (factor*nsteps)
          ELSE
            rampFunction(step) = 1.0
          END IF
        END DO
        
        !*************************************************************************
        ! COUPLED DYNAMIC SOLUTTION
        !*************************************************************************
        
        SolverParameters%SolverConvergence = 0
        d1x1=zero
        d2x1=zero  
        time=zero
        outputCounter = 0
        fracCounter = 0
        step = 0
        max = 1.0E29
        maxOld = 1.0E30
        dampFractor = 1.0
        
        ! Tolerance for maximum displacement
        xtol = 5.0 * MAXVAL( ABS(GeomechanicalModel%Displacement) )
        IF(xtol==0.0) xtol = 1.0E10
        
        ! Initial time step OR critical time step
        str = GeomechanicalModel%EffectiveStressTensor
        dt0 = Criticaldt(GeomechanicalModel,Matrices%Mass,str)*0.10
        IF( dt0 > SolverParameters%timeStep ) dt0 = SolverParameters%timeStep
        dtim = dt0
        
        !PREPARING LOG TABLE
        WRITE(text,'(A)')"   step timeStep    time       acel.    veloc.    disp.      tol"
        CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
        WRITE(text,'(I7,5E10.2,E10.2)')0,dtim,time,0.0,0.0,0.0,0.0
        CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
        
        startTime = OMP_get_wtime()
        
        !DEC$ NOPARALLEL
        DO !LOOPING---------------------------------------------------------------
            !Time stepping

            !str = currentStress+GeomechanicalModel%EffectiveStressTensor
            !dt0 = Criticaldt(GeomechanicalModel,Matrices%Mass,str)*0.10
            IF(dt0>SolverParameters%timeStep) dt0 = SolverParameters%timeStep
            dtim = dt0
          
            !DEC$ NOPARALLEL
            DO counter=1,OutputFreq !SUBSTEPING-------------------------------------
            !no stop criteria analysis in this stepping sequency of OutputFreq steps

                step = step + 1
                time = time + dtim

                ! KINECTIC ENERGY DERIVATIVE__________________________________________
                !verifies if a stable state was achieved
                IF( Kder(10)*Kder(9)<0 )THEN
                    d1x1 = 0.0
                    d2x1 = 0.0
                END IF

                ! DISPLACEMENT UPDATE_________________________________________________
                !$OMP WORKSHARE
                xold = x1
                !$OMP END WORKSHARE            
                !$OMP WORKSHARE
                x1 = x1 + dtim * d1x1 + pt5 * dtim**2 * d2x1
                !$OMP END WORKSHARE
          
                ! ENFORCE DISPLACEMENT BOUNDARY CONDITIONS____________________________
                ! Zero displacement B.C.
                !$OMP WORKSHARE
                x1(RestrictedDOFs(1:RestrictedDOFs(0))) = 0.0
                !$OMP END WORKSHARE

                ! Forced displacement B.C.
                !DEC$ PARALLEL
                DO i=1,numDispDofs
                    x1(DispBC%Dofs(i)) = DispBC%Displacement(i) * DispBC%RampFunction(step)
                END DO
            
                ! Negative Y displacement restriction.
                !DEC$ PARALLEL
                DO i=1,DispBC%ListYnodalRestriction(0)
                    g(1) = GeomechanicalModel%NodeDOFs(2,i)
                    IF(x1(g(1))<0.0)THEN
                        x1(g(1)) = 0.0
                    END IF
                END DO
          
                ! CONVERGENCE STTUFF__________________________________________________
                maxOld = max
                IF(step>1)max = MAXVAL(ABS(x1(1:)-xold(1:)))/MAXVAL(ABS(xold(1:)))

                ! INTERNAL FORCE______________________________________________________
                fint = zero

                ! Linear Elastic elements.
                ! {fint} = [K] * {x1}
                ! MKL blas2:mkl_dcsrgemv
                IF(ElastifStiffness)THEN
                    CALL mkl_dcsrgemv('N', numEquations, Matrices%Stiff%a, Matrices%Stiff%ia, &
                                        Matrices%Stiff%ja, x1(1:), fint(1:) )
                    fint = -fint
                END IF

                ! Traction Damage solid elements + shear strength constitutive model.
                !$OMP PARALLEL DO DEFAULT(SHARED) &
                !$OMP & PRIVATE(i1,iel,num,coord,g,bload,eld,e,v,i,det,der,jac,deriv,bee,&
                !$OMP &         strain,strainIncrement,Strength,stress,strEP,eload,syy_avg)
                DO i1=1,ListOfDamageElements(0)
              
                    iel = ListOfDamageElements(i1)
                    num = GeomechanicalModel%ElementNodes(:,iel)
                    coord=TRANSPOSE(g_coord(:,num))
                    g = GeomechanicalModel%ElementDOFs(:,iel)
                    bload = 0.0
                    eld = x1(g)                                     !USES TOTAL DISPLACEMENT.

                    modelNumber = GeomechanicalModel%CDMmodelNum(iel)

                    DO i=1,nip

                        ! FE matrices
                        det = detJac(i,iel)
                        CALL shape_der(der,points,i)    !-> der
                        jac=MATMUL(der,coord)           !-> jac
                        CALL invert(jac)                !-> inv(jac)
                        deriv=MATMUL(jac,der)           !-> deriv
                        CALL beemat(bee,deriv)          !-> bee

                        ! Strain
                        strain = MATMUL(bee,eld)                      !strain due TOTAL DISPLACEMENT.

                        ! Constitutive model
                        !e = GeomechanicalModel%YoungMod(iel) * (1.0-GeomechanicalModel%DamageStiffness(iel))
                        e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
                        e = e  * (1.0-GeomechanicalModel%DamageStiffness(iel))
                        !v = GeomechanicalModel%Poisson(iel)
                        v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
                        strainIncrement = strain - currentStrain(:,i,iel)                                                         !STRAIN INCREMENT for this time step.
                        Strength(1) = GeomechanicalModel%CDM(modelNumber)%ST        !1 - traction stregth.
                        Strength(2) = 0.0                                           !2 - traction strain at 1.
                        Strength(3) = GeomechanicalModel%CDM(modelNumber)%EyyMax    !3 - maximum traction strain.
                        CALL StressDamage( currentStress(:,i,iel)+GeomechanicalModel%EffectiveStressTensor(:,i,iel),&             !NEW STRESS STATE: stress -> current total stress.
                                           currentStrain(:,i,iel)+GeomechanicalModel%StrainTensor(:,i,iel),&
                                           strainIncrement,e,v,&
                                           Strength,&
                                           GeomechanicalModel%DamageTraction(iel),stress)                                                
                        stress = stress - GeomechanicalModel%EffectiveStressTensor(:,i,iel)                                       !REMOVING INITIAL STRESS.                
                        !stress = StressElastic(e,v,strain)  !ELASTIC BEHAVIOR
                        IF(ANY(ISNAN(stress)))THEN
                          !constitutive model failed to converged -> hydrostatic compressible fluid.
                          SolverParameters%SolverConvergence = 2  !<- structural solver failed.
                        END IF


                        ! Numerical integration
                        eload = MATMUL(stress,bee)                    !int(Bt*stress*dV)
                        bload = bload + eload * det * weights(i)

                        ! Plastic Strain.
                        strEP = 0.0
                        currentPlasticStrain(:,i,iel) = currentPlasticStrain(:,i,iel) + strEP

                        ! Saving strain.
                        currentStrain(:,i,iel) = strain

                        ! Y J Integral.
                        syy_avg = GeomechanicalModel%EffectiveStressTensor(2,i,iel) + 0.5 * (stress(2) + currentStress(2,i,iel))
                        IF(syy_avg>0.0)THEN
                          eld = x1(g) - xold(g)     !increment of displacement.
                          strain = MATMUL(bee,eld)  !increment of strain.
                          GeomechanicalModel%YYstrainDesityEnergy(iel) = GeomechanicalModel%YYstrainDesityEnergy(iel) + syy_avg * strain(2) / nip
                        END IF

                        ! Saving stress.
                        currentStress(:,i,iel) = stress

                    END DO

                    DO i=1,ndof
                        fint(g(i)) = fint(g(i)) - bload(i)
                    END DO
            
                    !if(iel==6) write(15,'(E12.4,F12.2)')currentStrain(2,1,iel)+GeomechanicalModel%StrainTensor(2,1,iel),&
                    !                                    stress(2)+GeomechanicalModel%EffectiveStressTensor(2,1,iel)
              
                END DO
                !$OMP END PARALLEL DO

                ! Mohr-Coulomb solid elements.
                !$OMP PARALLEL DO DEFAULT(SHARED) &
                !$OMP & PRIVATE(i1,iel,num,coord,g,bload,eld,e,v,i,det,der,jac,deriv,bee, &
                !$OMP &         strain,stress,strEP,eload,   &
                !$OMP &         c,phi,dil)
                DO i1=1,ListOfMohrCoulombElements(0)

                    iel = ListOfMohrCoulombElements(i1)
                    num = GeomechanicalModel%ElementNodes(:,iel)
                    coord=TRANSPOSE(g_coord(:,num))
                    g = GeomechanicalModel%ElementDOFs(:,iel)
                    bload = 0.0
                    eld = x1(g) - xold(g)                         !INCREMENT OF DISPLACEMENT.

                    DO i=1,nip

                        ! FE matrices
                        det = detJac(i,iel)
                        CALL shape_der(der,points,i)    !-> der
                        jac=MATMUL(der,coord)           !-> jac
                        CALL invert(jac)                !-> inv(jac)
                        deriv=MATMUL(jac,der)           !-> deriv
                        CALL beemat(bee,deriv)          !-> bee

                        ! Strain
                        strain = MATMUL(bee,eld)                    !strain due INCREMENT OF DISPLACEMENT.

                        ! Constitutive model
                        !e = GeomechanicalModel%YoungMod(iel) * (1.0-GeomechanicalModel%DamageStiffness(iel))
                        e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
                        e = e * (1.0-GeomechanicalModel%DamageStiffness(iel))
                        !v = GeomechanicalModel%Poisson(iel)
                        v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
                        c = GeomechanicalModel%Cohesion(iel) * (1.0-GeomechanicalModel%DamageShear(iel))
                        !phi = GeomechanicalModel%FrictionAngle(iel) * ( 3.1415926535897932384626433832795 / 180.0 )
                        !phi = TAN(phi) * (1.0-GeomechanicalModel%DamageShear(iel))
                        !phi = ATAN(phi) * (180.0 / 3.1415926535897932384626433832795)
                        phi = GeomechanicalModel%FrictionAngle(iel) * (1.0-GeomechanicalModel%DamageShear(iel))
                        dil = GeomechanicalModel%DilationAngle(iel) 
                        stress = StressElastoplasticMC( e,v,&
                                                        currentStress(:,i,iel)+GeomechanicalModel%EffectiveStressTensor(:,i,iel),&
                                                        c,phi,dil,strain,strEP,GeomechanicalModel%YieldFunc(i,iel))
                        IF(ANY(ISNAN(stress)))THEN
                          !constitutive model failed to converged -> hydrostatic compressible fluid.
                          SolverParameters%SolverConvergence = 2  !<- structural solver failed.
                        END IF

                        ! Numerical integration: int(Bt*stress*dV)
                        stress = stress + currentStress(:,i,iel)  !stress increment + cummulated stress increment.
                        eload = MATMUL(stress,bee)
                        bload = bload + eload*det*weights(i)

                        ! Saving stress and strain
                        currentStress(:,i,iel) = stress
                        currentStrain(:,i,iel) = currentStrain(:,i,iel) + strain
                        currentPlasticStrain(:,i,iel) = currentPlasticStrain(:,i,iel) + strEP

                    END DO

                    DO i=1,ndof
                        fint(g(i)) = fint(g(i)) - bload(i)
                    END DO

                END DO
                !$OMP END PARALLEL DO

                ! Hyperbolic solid elements.
                !!$OMP PARALLEL DO DEFAULT(SHARED) &
                !!$OMP & PRIVATE(i1,iel,num,coord,g,bload,eld,e,v,i,det,der,jac,deriv,bee, &
                !!$OMP &         strain,stress,strEP,eload,   &
                !!$OMP &         c,phi,dil)
                DO i1=1,ListOfSandModel(0)

                    iel = ListOfSandModel(i1)
                    num = GeomechanicalModel%ElementNodes(:,iel)
                    coord=TRANSPOSE(g_coord(:,num))
                    g = GeomechanicalModel%ElementDOFs(:,iel)
                    bload = 0.0
                    eld = x1(g) - xold(g)                         !INCREMENT OF DISPLACEMENT.

                    DO i=1,nip

                        ! FE matrices
                        det = detJac(i,iel)
                        CALL shape_der(der,points,i)    !-> der
                        jac=MATMUL(der,coord)           !-> jac
                        CALL invert(jac)                !-> inv(jac)
                        deriv=MATMUL(jac,der)           !-> deriv
                        CALL beemat(bee,deriv)          !-> bee

                        ! Strain
                        strain = MATMUL(bee,eld)                    !strain due INCREMENT OF DISPLACEMENT.

                        ! Constitutive model
                        !v = GeomechanicalModel%Poisson(iel)
                        v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
                        c = GeomechanicalModel%Cohesion(iel)
                        phi = GeomechanicalModel%FrictionAngle(iel)
                        !e_ref = GeomechanicalModel%YoungMod(iel)
                        e_ref = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
                        p_ref = GeomechanicalModel%SandModel(i1)%p_ref
                        m = GeomechanicalModel%SandModel(i1)%m
                        CALL invar(currentPlasticStrain(:,i,iel),sigm, eps, lode_theta)
                        eps = eps + GeomechanicalModel%SandModel(i1)%eq0
                        CALL InterpolateSandModelTable(GeomechanicalModel,eps,Cp,Cr,Phi,Phir,DilAngle)
                        CALL StressHYP( v,e_ref,p_ref,m,Cp,Cr,Phi,Phir,DilAngle,&
                                        currentStress(:,i,iel)+GeomechanicalModel%EffectiveStressTensor(:,i,iel),&
                                        strain,stress,strEP,GeomechanicalModel%YieldFunc(i,iel))
                        IF(ANY(ISNAN(stress)))THEN
                          !constitutive model failed to converged -> hydrostatic compressible fluid.
                          SolverParameters%SolverConvergence = 2  !<- structural solver failed.
                        END IF

                        ! Numerical integration: int(Bt*stress*dV)
                        stress = stress + currentStress(:,i,iel)  !stress increment + cummulated stress increment.
                        eload = MATMUL(stress,bee)
                        bload = bload + eload*det*weights(i)

                        ! Saving stress and strain
                        currentStress(:,i,iel) = stress
                        currentStrain(:,i,iel) = currentStrain(:,i,iel) + strain
                        currentPlasticStrain(:,i,iel) = currentPlasticStrain(:,i,iel) + strEP

                    END DO

                    DO i=1,ndof
                        fint(g(i)) = fint(g(i)) - bload(i)
                    END DO

                END DO
                !!$OMP END PARALLEL DO
            
                ! Cam-clay solid elements.
                !!$OMP PARALLEL DO DEFAULT(SHARED) &
                !!$OMP & PRIVATE(i1,iel,num,coord,g,bload,eld,e,v,i,det,der,points,jac,deriv,bee, &
                !!$OMP &         strain,stress,strEP,weights,eload,   &
                !!$OMP &         kappa,lambda,M,e0,P0)
                ! Modified Cam-clay solid elements.
                DO i1=1,ListOfCAMCLAYElements(0)
            
                    iel = ListOfCAMCLAYElements(i1)
                    num = GeomechanicalModel%ElementNodes(:,iel)
                    coord=TRANSPOSE(g_coord(:,num))
                    g = GeomechanicalModel%ElementDOFs(:,iel)
                    bload = 0.0
                    eld = x1(g) - xold(g)                         !current increment of displacement.

                    !Properties.
                    !e = GeomechanicalModel%YoungMod(iel) * (1.0-GeomechanicalModel%DamageStiffness(iel))
                    e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
                    e = e * (1.0-GeomechanicalModel%DamageStiffness(iel))
                    !v = GeomechanicalModel%Poisson(iel)
                    v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
                    !Pc(kPa)  M  kappa  lambda   e0   P0(kPa)
                    Pc = GeomechanicalModel%MCCparameters(1,iel)
                    M = GeomechanicalModel%MCCparameters(2,iel)
                    kappa = GeomechanicalModel%MCCparameters(3,iel)  !<-kappa,lambda,M,e0,P0
                    lambda = GeomechanicalModel%MCCparameters(4,iel)
                    e0 = GeomechanicalModel%MCCparameters(5,iel)
                    P0 = GeomechanicalModel%MCCparameters(6,iel)

                    DO i=1,nip

                        ! FE matrices
                        det = detJac(i,iel)
                        CALL shape_der(der,points,i)    !-> der
                        jac=MATMUL(der,coord)           !-> jac
                        CALL invert(jac)                !-> inv(jac)
                        deriv=MATMUL(jac,der)           !-> deriv
                        CALL beemat(bee,deriv)          !-> bee

                        ! Strain
                        strain = MATMUL(bee,eld)                    !strain due total displacement.

                        ! Constitutive model
                        strEP = currentPlasticStrain(:,i,iel) + GeomechanicalModel%PlasticStrainTensor(:,i,iel)
                        stress = StressMCC(strain,e,v,currentStress(:,i,iel)+GeomechanicalModel%EffectiveStressTensor(:,i,iel),&
                                           Pc,M,kappa,lambda,e0,P0,strEP)
                        IF(ANY(ISNAN(stress)))THEN
                          !constitutive model failed to converged -> hydrostatic compressible fluid.
                          SolverParameters%SolverConvergence = 2  !<- structural solver failed.
                        END IF

                        ! Strain.
                        currentPlasticStrain(:,i,iel) = strEP - GeomechanicalModel%PlasticStrainTensor(:,i,iel)

                        ! Numerical integration: int(Bt*stress*dV)
                        stress = stress + currentStress(:,i,iel)
                        eload = MATMUL(stress,bee)
                        bload = bload + eload*det*weights(i)

                        ! Saving stress and strain
                        currentStress(:,i,iel) = stress
                        currentStrain(:,i,iel) = currentStrain(:,i,iel) + strain


                        !Time step update.
                        rho = GeomechanicalModel%Density(iel)
                        c = SQRT( e*(1.0-v)/(1.0+v)/(1.0-2.0*v)/rho )
                        c = c * SQRT( (1.0-v)/(1.0+v)/(1.0-2.0*v))
                        dt0 = GeomechanicalModel%ElementMinLength(iel) / c
                        IF( 0.80*dt0 < dtim ) dtim = 0.10*dt0

                    END DO

                    DO i=1,ndof
                        fint(g(i)) = fint(g(i)) - bload(i)
                    END DO
          
                END DO

                fint(0) = zero

                ! TOTAL FORCE_________________________________________________________
                bdylds = 0.0
            
                !$OMP WORKSHARE
                bdylds(1:) = bdylds(1:) + fint(1:) + rampFunction(step) * Fext
                !$OMP END WORKSHARE

                !$OMP WORKSHARE
                bdylds(1:) = bdylds(1:) - dampMultiplier(1:) * Damp%Values * d1x1(1:)
                !$OMP END WORKSHARE

                !$OMP WORKSHARE
                bdylds(RestrictedDOFs(1:RestrictedDOFs(0))) = 0.0
                !$OMP END WORKSHARE
          
                ! VELOCITY AND ACCELERTION UPDATE_____________________________________
                !$OMP WORKSHARE
                bdylds(freeDOFs) = bdylds(freeDOFs) / diag(freeDOFs)
                !$OMP END WORKSHARE

                !$OMP WORKSHARE
                d1x1 = d1x1 + (d2x1+bdylds ) * pt5 * dtim
                d2x1 = bdylds
                !$OMP END WORKSHARE
            
                ! DAMAGE______________________________________________________________
                CALL DamageUpdate(GeomechanicalModel,currentStress,currentPlasticStrain )

                ! UPDATE GEOMETRY_____________________________________________________
                IF(GeomechanicalModel%UpdateGeometry)THEN
                    g_coord(1,:) = GeomechanicalModel%X - x1(1:3*numNodes-2:3)
                    g_coord(2,:) = GeomechanicalModel%Y - x1(2:3*numNodes-1:3)
                    g_coord(3,:) = GeomechanicalModel%Z - x1(3:3*numNodes:3)
                    !$OMP PARALLEL DO DEFAULT(SHARED) &
                    !$OMP & PRIVATE(iel,num,coord,i,der,jac,det)
                    DO iel=1,numElements
                        num = GeomechanicalModel%ElementNodes(:,iel)
                        coord=TRANSPOSE(g_coord(:,num))
                        DO i=1,nip
                            CALL shape_der(der,points,i)      !-> der
                            jac=MATMUL(der,coord)             !-> jac
                            CALL determinante(jac,det)        !-> det 
                            detJac(i,iel) = det               !<- det
                        END DO
                    END DO
                    !$OMP END PARALLEL DO
                END IF

                ! KINECT ENERGY____________________________________________________
                kinecticEnergy(1:9)=kinecticEnergy(2:10)
                kinecticEnergy(10)=0.0
                !$OMP PARALLEL DO DEFAULT(SHARED) &
                !$OMP & PRIVATE(iel,num,dof,i)
                DO iel=1,numElements
                    DO i=1,8
                      num = GeomechanicalModel%ElementNodes(:,iel)
                      dof=num(i)*3
                      kinecticEnergy(10)=kinecticEnergy(10)+diag(dof)*d1x1(dof)*d1x1(dof)
                      kinecticEnergy(10)=kinecticEnergy(10)+diag(dof-1)*d1x1(dof-1)*d1x1(dof-1)
                      kinecticEnergy(10)=kinecticEnergy(10)+diag(dof-2)*d1x1(dof-2)*d1x1(dof-2)
                    END DO
                END DO
                !$OMP END PARALLEL DO
                IF(kinecticEnergy(10)>maxKinecticEnergy) maxKinecticEnergy=kinecticEnergy(10)
                kDer(1:9)=kDer(2:10)
                kDer(10) = (kinecticEnergy(10)-kinecticEnergy(9)) / dtim
                !write(101,'(2F16.5)')kinecticEnergy(10),kDer(10)

            END DO!SUBSTEPING-------------------------------------------------------
          
            outputCounter = outputCounter + 1
        
            !next velocities.
            xold = xnew
            xnew = x1

            dofMonitor = MAXLOC(ABS(x1(1:)),1)

            !parameter for solver convergence
            max = MAXVAL(ABS(xnew(1:)-xold(1:)))/MAXVAL(ABS(xold(1:)))

            !console+log file
            WRITE(text,'(I7,5E10.2,E10.2)')step,dtim,time,d1x1(dofMonitor),d2x1(dofMonitor),x1(dofMonitor),max
            CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))

            !file output.
            IF(nodeOutput>0)THEN
                nodeDispOutput(1,outputCounter) = GeomechanicalModel%Displacement(nodeOutput*3-2) + x1(nodeOutput*3-2)
                nodeDispOutput(2,outputCounter) = GeomechanicalModel%Displacement(nodeOutput*3-2) + x1(nodeOutput*3-1)
                nodeDispOutput(3,outputCounter) = GeomechanicalModel%Displacement(nodeOutput*3-2) + x1(nodeOutput*3)
                elementStressOutput(1:6,outputCounter) = GeomechanicalModel%EffectiveStressTensor(:,1,elementOutput) + &
                                                        currentStress(:,1,elementOutput)
                elementStrainOutput(1:6,outputCounter) = GeomechanicalModel%StrainTensor(:,1,elementOutput) + &
                                                        currentStrain(:,1,elementOutput)
                CALL invar(elementStressOutput(1:6,outputCounter), &
                           elementStressOutput(7,outputCounter), &
                           elementStressOutput(8,outputCounter), &
                           elementStressOutput(9,outputCounter))
            END IF
            ! SOLID DEGRADATION ANALYSIS--------------------------------------------
            !SolverParameters%SolverConvergence  !0 - success
                                                 !1 - not converged
                                                 !2 - failed
            DO iel=1,NumElements
              IF(ISNAN(max))THEN
                SolverParameters%SolverConvergence = 2
              END IF
            END DO
          
            ! STOP CRITERIA---------------------------------------------------------
            ! solver fail test
            IF(ISNAN(max))THEN
              SolverParameters%SolverConvergence = 2
            END IF
            ! excess of displacement test
            IF(MAXVAL(ABS(x1(1:)))>SolverParameters%MaxDisplacementControl)THEN
              SolverParameters%SolverConvergence = 2
            END IF
            ! update kinect energy

            ! convergence analysis
            IF(SolverParameters%SolverConvergence==2)RETURN
            IF(max<=SolverParameters%tol .OR.&
              ( KinecticEnergy(10)/maxKinecticEnergy < SolverParameters%tol )  ) THEN
                converged=.TRUE.
            END IF
            IF(step>=SolverParameters%nsteps) THEN
              SolverParameters%SolverConvergence = 1
              EXIT
            END IF
            IF(  (converged).AND.(step>1).AND.&
                 (.NOT.doNotStop).AND.(rampFunction(step)>=0.99) )THEN
              SolverParameters%SolverConvergence = 0
              EXIT
            END IF
            SolverParameters%timeStep = dtim
          
        END DO !-END STEP LOOP----------------------------------------------------
        
        endTime = OMP_get_wtime()

        WRITE(text,'(A,F16.5,A)')'...tempo de solver (wallclock): ',(endTime-startTime),'sec.'
        CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))

        !*************************************************************************
        ! CALCULATING STRESS AT ELASTIC ELEMENTS.
        !*************************************************************************
        !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:fint) &
        !$OMP & PRIVATE(i1,iel,num,coord,g,bload,eld,e,v,i,det,der,jac,deriv,bee,&
        !$OMP &         strain,stress,strEP,eload)
        DO i1=1,ListOfElasticElements(0)
            
            iel = ListOfElasticElements(i1)
            num = GeomechanicalModel%ElementNodes(:,iel)
            coord=TRANSPOSE(g_coord(:,num))
            g = GeomechanicalModel%ElementDOFs(:,iel)
            bload = 0.0
            eld = x1(g)                                     !USES TOTAL DISPLACEMENT.

            !e = GeomechanicalModel%YoungMod(iel) !* (1.0-GeomechanicalModel%DamageStiffness(iel))
            e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
            !e = e * (1.0-GeomechanicalModel%DamageStiffness(iel))
            !v = GeomechanicalModel%Poisson(iel)
            v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)

            DO i=1,nip

                ! FE matrices.
                CALL shape_der(der,points,i)    !-> der
                jac=MATMUL(der,coord)           !-> jac
                CALL invert(jac)                !-> inv(jac)
                deriv=MATMUL(jac,der)           !-> deriv
                CALL beemat(bee,deriv)          !-> bee

                ! Strain
                strain = MATMUL(bee,eld)                      !strain due TOTAL DISPLACEMENT.

                ! Constitutive model
                stress = StressElastic(e,v,strain)            !stress due TOTAL DISPLACEMENT.

                ! Plastic Strain.
                ! NO PLASTIC STRAIN.

                ! Saving stress and strain
                currentStress(:,i,iel) = stress
                currentStrain(:,i,iel) = strain

            END DO
              
        END DO
        !$OMP END PARALLEL DO

        DO iel=1,numElements
            DO i=1,nip
              det = detJac(i,iel)
              GeomechanicalModel%PlasticStrainWork(iel) = GeomechanicalModel%PlasticStrainWork(iel) + &
                    0.5*DOT_PRODUCT(currentPlasticStrain(:,i,iel),&
                                    GeomechanicalModel%EffectiveStressTensor(:,i,iel) + currentStress(:,i,iel)) &
                                   * det * weights(i)
            END DO
              
        END DO
        
        WRITE(text,'(A,F16.5,A)')'...Total Plastic Strain: ',SUM(GeomechanicalModel%PlasticStrainWork),' J'
        CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
        
        !*************************************************************************
        ! SAVING STRESS AND STRAIN.
        !*************************************************************************
        
        GeomechanicalModel%Displacement = GeomechanicalModel%Displacement + x1(1:)        
        GeomechanicalModel%EffectiveStressTensor = GeomechanicalModel%EffectiveStressTensor + currentStress        
        GeomechanicalModel%StrainTensor = GeomechanicalModel%StrainTensor + currentStrain
        GeomechanicalModel%PlasticStrainTensor = GeomechanicalModel%PlasticStrainTensor + &
                                                                       currentPlasticStrain
        !*************************************************************************
        !UPDATING YIELD FUNCTION
        !*************************************************************************
        DO i1=1,ListOfCAMCLAYElements(0)
            
            iel = ListOfCAMCLAYElements(i1)

            !Properties.
            !Pc(kPa)  M  kappa  lambda   e0   P0(kPa)
            Pc = GeomechanicalModel%MCCparameters(1,iel)
            M = GeomechanicalModel%MCCparameters(2,iel)

            DO i=1,nip
              CALL invar(GeomechanicalModel%EffectiveStressTensor(:,i,iel),&
                         sigm, dsbar, lode_theta)
              GeomechanicalModel%YieldFunc(i,iel) = mccf(M,Pc,sigm,dsbar)
            END DO
          
        END DO
        
        !*************************************************************************
        ! INCREMENT OF POREPRESSURE DUE UNDRAINED BEHAVIOR.
        !*************************************************************************
        
        ALLOCATE( nsum(GeomechanicalModel%NumNodes), Pp(GeomechanicalModel%NumNodes) )
        nsum = 0
        Pp = 0.0
        DO i1=1,ListOfUndrainedElements(0)
            
            iel = ListOfUndrainedElements(i1)
            num = GeomechanicalModel%ElementNodes(:,iel)
            coord=TRANSPOSE(g_coord(:,num))
            g = GeomechanicalModel%ElementDOFs(:,iel)

            eld = x1(g)                                             !TOTAL DISPLACEMENT.

            biotCoef = GeomechanicalModel%BiotCoef(iel)
            Kf = GeomechanicalModel%PoreFluidBulkMod(iel)

            deltaPorepressure(iel) = 0.0

            DO i=1,nip

                ! FE matrices
                det = detJac(i,iel)
                CALL shape_der(der,points,i)    !-> der
                jac=MATMUL(der,coord)           !-> jac
                CALL invert(jac)                !-> inv(jac)
                deriv=MATMUL(jac,der)           !-> deriv
                CALL beemat(bee,deriv)          !-> bee

                ! Strain
                strain = MATMUL(bee,eld)                              !strain due TOTAL DISPLACEMENT.

                ! Constitutive model
                stress = FluidHydrostaticStress(Kf,biotCoef,strain)    !stress due TOTAL DISPLACEMENT.

                ! Udrained porepressure increment.
                deltaPorepressure(iel) = deltaPorepressure(iel) + SUM(stress(1:3))/3.0

            END DO

            deltaPorepressure(iel) = 1.0/nip * (  deltaPorepressure(iel) )

            Pp(num) = Pp(num) - deltaPorepressure(iel)
            nsum(num) = nsum(num) + 1
          
        END DO
        
        DO i=1,SIZE(Pp)
            IF( nsum(i)>0 ) Pp(i) = Pp(i) / nsum(i)
        END DO
        
        GeomechanicalModel%Porepressure = GeomechanicalModel%Porepressure + Pp
        
        !*************************************************************************
        ! OUTPUT MONITOR NODE/ELEMENT.
        !*************************************************************************
        IF(nodeOutput>0)THEN
          
            WRITE(LogChanel,'(A,I5)')'MONITOR NODE - DISPLACEMENTS',nodeOutput
            WRITE(LogChanel,*)'    Ux      Uy      Uz'
            DO i=1,outputCounter
                WRITE(LogChanel,'(3E15.5)') nodeDispOutput(1,i),&
                                            nodeDispOutput(2,i),&
                                            nodeDispOutput(3,i)
            END DO

            WRITE(LogChanel,'(A,I5)')'MONITOR ELEMENT - STRESS (kPa)',elementOutput
            DO i=1,outputCounter
                WRITE(LogChanel,'(8F10.1,E14.5)')elementStressOutput(:,i)/1000.0
            END DO

            WRITE(LogChanel,'(A,I5)')'MONITOR ELEMENT - POREPRESSURE (kPa)',elementOutput
            DO i=1,outputCounter
                 WRITE(LogChanel,'(8F10.1,E14.5)')elementStressOutput(:,i)/1000.0
            END DO

            WRITE(LogChanel,'(A,I5)')'MONITOR ELEMENT - STRAIN',elementOutput
            DO i=1,outputCounter
                WRITE(LogChanel,'(6E12.3)')elementStrainOutput(:,i)
            END DO
          
            !print stuff-------------------------------------------------------------------
            !lambda = GeomechanicalModel%MCCparameters(4,1)
            !e0 = GeomechanicalModel%MCCparameters(5,1)
            !P0 = GeomechanicalModel%MCCparameters(6,1)
            !e = e0 - lambda * log(elementStressOutput(7,outputCounter)/P0)  !void ratio
            !write(15,*) elementStressOutput(7,outputCounter)/1000.0,&
            !            elementStressOutput(8,outputCounter)/1000.0,&   !p,q,
            !            e
            !------------------------------------------------------------------------------
          
        END IF
        
        CALL PrintMessage(.FALSE.,.FALSE.,'...Explicit structural finite element solver FINISHED')
        
        RETURN
    END SUBROUTINE
                              
    FUNCTION Criticaldt(GeomechanicalModel,Mass,stress)
    
        USE DiscretizedModel
        USE HydraulicFractureModel
        USE ConstitutiveGeoModels
        USE StiffnessMatrix
        IMPLICIT NONE

        REAL(8):: Criticaldt
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE(MassMatrix),INTENT(IN):: Mass
        REAL(8),INTENT(IN):: stress(:,:,:)

        INTEGER:: NumElements
        INTEGER:: num(8)

        INTEGER:: iel,dof

        REAL(8):: dx,dy,dz,Lmin,c,e,v,k,rho,timeStep,defVol,evol_vp

        REAL(8):: M,Pc,Pt,A,P0,e0,lambda,kappa

        NumElements = GeomechanicalModel%NumElements

        Criticaldt = 1.0E30

        DO iel=1,NumElements

          IF(GeomechanicalModel%ActiveElements(iel))THEN
            num = GeomechanicalModel%ElementNodes(:,iel)

            dx = ABS(GeomechanicalModel%X(num(1)) - GeomechanicalModel%X(num(4)) )
            dy = ABS(GeomechanicalModel%Y(num(1)) - GeomechanicalModel%Y(num(5)) )
            dz = ABS(GeomechanicalModel%Z(num(1)) - GeomechanicalModel%Z(num(2)) )

            Lmin = dz
            IF(dx<Lmin) Lmin = dx
            IF(dy<Lmin) Lmin = dy

            !e = GeomechanicalModel%YoungMod(iel) * (1.0-GeomechanicalModel%DamageStiffness(iel))
            e = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%E1 + GeomechanicalModel%TransverIsotropicProps(iel)%E2)
            e = e * (1.0-GeomechanicalModel%DamageStiffness(iel))
            !v = GeomechanicalModel%Poisson(iel)
            v = 0.5 * (GeomechanicalModel%TransverIsotropicProps(iel)%n1 + GeomechanicalModel%TransverIsotropicProps(iel)%n2)
            IF(GeomechanicalModel%Undrained(iel))THEN
              k = e/3.0/(1.0-2.0*v) + GeomechanicalModel%PoreFluidBulkMod(iel)
              e = k*3.0*(1.0-2.0*v)
            END IF

            ! CAM-CLAY
            IF(GeomechanicalModel%ConstitutiveGeoModel(iel)==3)THEN
              !Pc(kPa)  M  kappa  lambda   e0   P0(kPa)
              Pc = GeomechanicalModel%MCCparameters(1,iel)
              M = GeomechanicalModel%MCCparameters(2,iel)
              kappa = GeomechanicalModel%MCCparameters(3,iel)  !<-kappa,lambda,M,e0,P0
              lambda = GeomechanicalModel%MCCparameters(4,iel)
              e0 = GeomechanicalModel%MCCparameters(5,iel)
              P0 = GeomechanicalModel%MCCparameters(6,iel)
              defVol = SUM(GeomechanicalModel%StrainTensor(1:3,1,iel))
              evol_vp = SUM(GeomechanicalModel%PlasticStrainTensor(1:3,1,iel))
              CALL mccdProps(M,e0,P0,lambda,kappa,SUM(GeomechanicalModel%EffectiveStressTensor(:,1,iel))/3.0,&
                            DefVol,evol_vp,Pc,e,v)!<-update e and v.
            END IF

            rho = GeomechanicalModel%Density(iel)
            c = SQRT( e*(1.0-v)/(1.0+v)/(1.0-2.0*v)/rho )
            c = c * SQRT( (1.0-v)/(1.0+v)/(1.0-2.0*v))

            timeStep = Lmin / c

            IF(timeStep<Criticaldt) Criticaldt = timeStep

            GeomechanicalModel%ElementMinLength(iel) = Lmin
          END IF

        END DO

        !Spring-mass elements.
        !DO iel=1,SIZE(HFmodel%SpringDOF)
        !  dof = HFmodel%SpringDOF(iel)
        !  IF(dof>0) timeStep = SQRT(HFmodel%SpringKy(iel)/Mass%Values(dof))
        !  IF(timeStep<Criticaldt) Criticaldt = timeStep
        !END DO
    
    END FUNCTION
    
    
    SUBROUTINE PrepareDynamicExplicitSolver(GeomechanicalModel,Matrices,DynamicParameters,&
                                            RestrictedDOFs,Chanel_base)
      
        USE DiscretizedModel
        USE FEMlib
        USE StiffnessMatrix
        USE IO

        IMPLICIT NONE
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE(SolverMatricesContainer),INTENT(IN OUT):: Matrices
        TYPE(DynamicSolverParameters):: DynamicParameters
        INTEGER,INTENT(IN):: RestrictedDOFs(0:)
        INTEGER,INTENT(IN):: Chanel_base

        INTEGER:: NumDOFs,NumElements
        INTEGER,ALLOCATABLE:: num(:),&                    !element nodes
                              g(:)                        !element dofs
        INTEGER:: i,iel,n,node,dir
        REAL(8):: vol, m
        INTEGER:: nod = 8                                 !number of nodes per element
        INTEGER:: ndof = 24                               !number of displacement dofs per element
        INTEGER:: nfreedofs
        LOGICAL:: found
        CHARACTER(LEN=320):: text
        REAL(8),ALLOCATABLE:: g_coord(:,:),coord(:,:)
        CHARACTER(LEN=15)::element='hexahedron'
        REAL(8),ALLOCATABLE:: points(:,:),weights(:)
        INTEGER:: nip=8

        NumDOFs = GeomechanicalModel%NumEquations
        NumElements = GeomechanicalModel%NumElements

        ! MASS MATRIX
        ALLOCATE(Matrices%Mass%Values(0:NumDOFs),g(ndof),num(nod))
        Matrices%Mass%Values = 0.0
        ALLOCATE( g_coord(3,GeomechanicalModel%NumNodes),coord(nod,3))
        g_coord = 0.0
        g_coord(1,:) = GeomechanicalModel%X
        g_coord(2,:) = GeomechanicalModel%Y
        g_coord(3,:) = GeomechanicalModel%Z
        ALLOCATe( points(nip,3),weights(nip) )
        CALL sample(element,points,weights)
        DO iel=1,NumElements

          num = GeomechanicalModel%ElementNodes(:,iel)
          g = GeomechanicalModel%ElementDOFs(:,iel)
          coord = TRANSPOSE(g_coord(:,num))

          GeomechanicalModel%Volume(iel) = ElementVolume(num, coord, nip, points, weights)
          vol = GeomechanicalModel%Volume(iel)
          m = vol * GeomechanicalModel%Density(iel)

          nfreedofs = ndof
 
          Matrices%Mass%Values(g) = Matrices%Mass%Values(g) + m / nfreedofs

        END DO

        ! SOLVER PARAMETERS
        CALL SearchInFile(Chanel_base,'#STRUCT_DynaSOLVER',found)
        IF(found) THEN
          READ(Chanel_base,*)DynamicParameters%tol
          READ(Chanel_base,*)DynamicParameters%nsteps
          READ(Chanel_base,*)DynamicParameters%dampFractor
          READ(Chanel_base,*)DynamicParameters%timeStep
        ELSE
          DynamicParameters%tol=1.0e-6
          DynamicParameters%dampFractor = 1.0E3
          DynamicParameters%nsteps = 50000
          DynamicParameters%timeStep = 1.0e-6
        END IF

        RETURN
    END SUBROUTINE
    
    ! Builds an CSR stiffness matrix for linear elastic elements
    ! Porpouse: speedup the dynamic solver
    SUBROUTINE BuildElasticStiffness(GeomechanicalModel,Matrices,ActiveElements)
      
      USE DiscretizedModel
      USE StiffnessMatrix
      USE IO

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      LOGICAL,INTENT(IN):: ActiveElements(:)
      TYPE(SolverMatricesContainer),INTENT(IN OUT):: Matrices

      !   TYPE SolverMatricesContainer
      !     TYPE(ElasticStiffnessMatrix):: Stiff
      !     TYPE(MassMatrix):: Mass
      !     TYPE(Damping):: Damp
      !   END TYPE

      INTEGER:: iel
      LOGICAL,ALLOCATABLE:: ElementsToAdd(:)

      ! Allocating full stiffness matrix
      IF( .NOT.ALLOCATED(Matrices%Stiff%ia) )THEN
        CALL SparseStiffnessMatrixCSR( GeomechanicalModel, Matrices%Stiff%a, &
                                        Matrices%Stiff%ia, Matrices%Stiff%ja )
      END IF

      !Identifying the elastic elements to built the stiffness
      ALLOCATE(ElementsToAdd(GeomechanicalModel%NumElements))
      ElementsToAdd = .FALSE.
      DO iel=1,GeomechanicalModel%NumElements
        IF(ActiveElements(iel))THEN
          IF( GeomechanicalModel%ConstitutiveGeoModel(iel)==0)THEN
            ElementsToAdd(iel) = .TRUE.
          END IF
        END IF
      END DO

      ! Builds an stiffness matrix with only elastic matrix
      Matrices%Stiff%a = 0.0
      CALL StiffnessConstruction( GeomechanicalModel, Matrices%Stiff%a, &
                                  Matrices%Stiff%ia, Matrices%Stiff%ja, ElementsToAdd )

      RETURN
    END SUBROUTINE
    
END MODULE