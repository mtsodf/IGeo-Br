
    !*****************************************************************************
    !COMPUTATIONAL GEOMECHANICS
    !Eng. Luis Carlos de Sousa Junior
    !Programa de Pos-graduação em Geotecnia.
    !University of Brasilia
    !*****************************************************************************

    MODULE FluidFlowSolver

    !PIPE MODEL.
    TYPE PipeModel
        !Inside cell (meadle)
        CHARACTER(LEN=10),ALLOCATABLE:: Path(:,:)   !I,J,K,DIR PATH
                                                    !DIR: 1 for horizontal, 2 for vertical
        REAL(8):: Diameter
        INTEGER,ALLOCATABLE:: GlobalNodesNumber(:)  !nodes ordered
        REAL(8),ALLOCATABLE:: PipeTrans(:)          !transmissibility for each pipe element
    END TYPE
    
    !POROUS FLOW MODEL.
    ! NODE BASED.
    ! FINITE DIFFERENCE TRANSMISSIBILITY.
    TYPE, PUBLIC:: PorousFlowModel

        REAL(8):: CriticalTimeStep = 1.0          !1 second as initial value

        !Fluid data.
        REAL(8):: Cw,Bwref,Pwref,Rhoref,Viscw     !water.

        REAL(8),ALLOCATABLE:: NodalPoreCompressbility(:)

        LOGICAL:: RAD = .FALSE.,&
                  HorWellFlag = .FALSE.                     !BOX MODEL AND HORIZONTAL WELL (X DIR)
        
        !Active to flow list.
        INTEGER,ALLOCATABLE:: ActiveToFlowConections(:)

        !Depth.
        REAL(8),ALLOCATABLE:: NodeDepth(:)

        REAL(8),ALLOCATABLE:: NodeVolume(:),NodeFracVolume(:),ElementVolume(:)

        !Volume da fratura já contabilizado no balanço de massa do fluxo, por
        !nó (veja "BALANÇO DE MASSA DA FRATURA" em FlowStepSolution).
        REAL(8),ALLOCATABLE:: FracVolumeAccounted(:)

        !Node pair trasmissibility.
        REAL(8),ALLOCATABLE:: TransArrayINI(:)    !INITIAL transmissibility array (all nodes transmissibilities).
        REAL(8),ALLOCATABLE:: TransArray(:)       !CURRENT transmissibility array
        INTEGER,ALLOCATABLE:: TransPointer(:,:)   !node pair transmissibilities for each TransArray.
        REAL(8),ALLOCATABLE:: TransAreas(:,:)     !(1,element) - AX; (2,element) - AY; (3,element) - AZ.
                                                  !(4,element) - TX; (5,element) - TY; (6,element) - TZ.
        INTEGER,ALLOCATABLE:: Trans2Element(:)
        LOGICAL,ALLOCATABLE:: IsLiquefied(:)
        REAL(8),ALLOCATABLE:: ConectionFlowRates(:)   !stores internodes flow rates.

        !Static Pressure
        REAL(8):: StaticPressure, StaticPressureZref, GradStaticPressure
        
        !Outflow
        !Q = -J * (Pi - Pest)
        REAL(8),ALLOCATABLE:: Jfactor(:)            !nodal value
        REAL(8),ALLOCATABLE:: Pest(:)               !nodal value
        
        !Accumulation.
        REAL(8),ALLOCATABLE:: Accum(:)            !nodal accumulation term.

        !Boundary condtitions.
        REAL(8),ALLOCATABLE:: Rate(:)             !Instantaneous node rate.
        REAL(8):: AverageTotalRate
        !REAL(8),ALLOCATABLE:: PressureBC(:)      !Node pressure boundary condition.
        LOGICAL,ALLOCATABLE:: PressureBC_flag(:)  !flag for node pressure BC.
        REAL(8):: TotalFlowRate

        !WELL BORE DATA.. ( #DirectionalWell )
        REAL(8),ALLOCATABLE:: ListOfWellNodes(:)
        REAL(8),ALLOCATABLE:: ListOfWellElements(:)

        !Flows.
        LOGICAL,ALLOCATABLE:: InjNode(:)          !Nodes with fluid injection.
        REAL(8),ALLOCATABLE:: NodeFlows(:),&      !porous + FRACTURE injection
                              ElementFlows(:),RateFactor(:)
        REAL(8):: TotalInjectionRate              !porous + FRACTURE injection
        REAL(8):: WellIndex

        !Matrices
        REAL(8),ALLOCATABLE:: a(:)
        INTEGER,ALLOCATABLE:: ia(:),ja(:)
        
        !Fracture Unbalanced Volume
        LOGICAL:: FracUnbalVolume = .FALSE.
        
        !Permeability damage/stimulation
        REAL(8),ALLOCATABLE:: OriginalPerm(:,:)             !(1,:) - KH  , (2,:) - KH
        REAL(8),ALLOCATABLE:: TimeOfInjection(:)            !saves the cumulative time of exposition to fluid injection done for well or fracture.
        REAL(8),ALLOCATABLE:: InjectionDamageModel(:,:)     !for each layer, contains time and perm. multiplier for injection damage model.

        !Pipe model
        LOGICAL:: ContainsPipeModel = .FALSE.
        TYPE(PipeModel):: Pipe
        
    END TYPE
    
    CONTAINS


    ! TRANSIENT FLOW SOLVER++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE FlowStepSolution( nthreads,FlowModel,HFmodel,GeomechanicalModel,WellData,&
                                 DeltaTime,LimitRate,reportData,Porepressure,WaterSaturation )
        ![A]*{p} + [C]*(d/dt){p} + {Q_std} = {RHS}
        !USE IFPORT
        USE IO
        USE ProblemDefinition
        USE DiscretizedModel
        USE HydraulicFractureModel
        USE PCG_Flow
        USE OMP_LIB  !bibliotecas OpenMP(R)

        IMPLICIT NONE
        INTEGER,INTENT(IN):: nthreads
        TYPE(Model),INTENT(IN):: GeomechanicalModel
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
        REAL(8),INTENT(IN):: DeltaTime,LimitRate
        TYPE(StepReport),INTENT(IN OUT):: reportData
        REAL(8),INTENT(IN OUT):: Porepressure(:),WaterSaturation(:)

        LOGICAL,ALLOCATABLE:: RestrictedDOFs(:),presureBC_flag(:)

        REAL(8):: dt,t(0:5)

        CHARACTER(LEN=320):: text
        INTEGER:: iter,maxIters=1,linSolverIters
        INTEGER:: nTrans
        LOGICAL:: converged

        INTEGER:: step,i,j,k,pos,iel,counter,nCells,numNodes,neq,ix,iz,irad,dof,dof1,dof2,cell,n,nnz
        REAL(8),ALLOCATABLE:: A(:),Afrac(:),RHS(:),C(:),diag(:)
        REAL(8),ALLOCATABLE:: csrA_c(:)

        REAL(8),ALLOCATABLE:: P0(:),solution(:)
        REAL(8),ALLOCATABLE:: HydraulicHead(:)
        REAL(8):: maxPressureVariation,maxAllowedPressureVariation

        INTEGER:: num(8),g(4)

        REAL(8),ALLOCATABLE:: fracTrans(:,:)
        INTEGER,ALLOCATABLE:: fracTransJ(:,:)

        REAL(8):: V, V0, V1, phi, Cw, Bwref, Pref, Viscw, mob, grav, deltaP
        REAL(8):: totalCompressibility
        REAL(8):: HydraulicPotential
        REAL(8),ALLOCATABLE:: fractureUnbalVolume(:)        !fracture dilation due its displacement on mechanical step.

        !Balanço de massa da fratura.
        REAL(8),ALLOCATABLE:: fracCompliance(:),&           !dV/dp estimado de cada nó da face aberta (m3/Pa).
                              fracMismatch(:),&             !volume da mecânica ainda não contabilizado (m3).
                              Pstart(:)                     !pressão no início da chamada.
        INTEGER,ALLOCATABLE:: faceElements(:),openFaceElements(:)
        LOGICAL:: mismatchPending
        REAL(8):: zLow,zHigh,xLow,xHigh,charLength,planeStrainModulus
        INTEGER:: nOpen,fracEl,solidEl

        REAL(8):: min,temp

        REAL(8),ALLOCATABLE:: Bw(:), ViscosityW(:), Rhow(:)

        REAL(8),ALLOCATABLE:: Rates(:)

        !Fracture Matrix.
        INTEGER:: nFracTrans
        REAL(8),ALLOCATABLE:: HFa(:)    !non-zero terms of transmissibility matrix (with BC).

        !Solver parameters.
        INTEGER:: cg_limit = 1.0E6
        REAL(8):: cg_tol = 1.0E-6

        REAL(8),ALLOCATABLE::diag_precon(:)  !pré-dondicionador.

        REAL(8),ALLOCATABLE:: initialGuess(:)

        REAL(8):: Pinj(2),Pest,II,Zref,Qinj(2)
        INTEGER:: IIcounter

        REAL(8):: limitRateDeviation = 0.01

        REAL(8):: maxRate

        LOGICAL:: rateControled

        REAL(8):: penaltyValue = 1.0E30

        !************************************************************************
        !POROUS MEDIA
        !************************************************************************
        nCells = GeomechanicalModel%NumElements
        numNodes = GeomechanicalModel%NumNodes
        neq = numNodes
        nTrans = SIZE(FlowModel%TransArray)
        nFracTrans = SIZE(HFmodel%Trans)

        ALLOCATE( RHS(0:neq), initialGuess(neq), solution(0:neq) )
        ALLOCATE( P0(0:neq),diag(0:neq) )
        !!dir$ loop count min(512)
        RHS = 0.0
        !!dir$ loop count min(512)
        diag = 0.0
        !!dir$ loop count min(512)
        solution = 0.0

        !Transmissibility / Mobility matrix.
        ALLOCATE( A(nTrans), Afrac(nFracTrans),C(0:numNodes) )
        !!dir$ loop count min(512)
        A = 0.0
        !!dir$ loop count min(512)
        C = 0.0

        ALLOCATE(RestrictedDOFs(neq))
        RestrictedDOFs = .FALSE.
        counter = 0
        DO i=1,numNodes
            IF(FlowModel%PressureBC_flag(i)) THEN
                counter = counter + 1
                RestrictedDOFs(i) = .TRUE.
            END IF
        END DO

        ALLOCATE(presureBC_flag(size(FlowModel%PressureBC_flag)))
        presureBC_flag = FlowModel%PressureBC_flag

        ALLOCATE(Rates(neq))
        Rates = 0.0

        ALLOCATE( csrA_c(SIZE(FlowModel%a)) )

        Pinj = 0.0
        Qinj = 0.0
        maxRate = LimitRate
        IF( (LimitRate>0.9E30) ) THEN
            rateControled = .FALSE.
        ELSE
            rateControled = .TRUE.
        END IF

        !************************************************************************
        !NODE MOBIITY POTENTIAL.
        !************************************************************************
        ALLOCATE( Bw(numNodes), ViscosityW(numNodes), Rhow(numNodes), HydraulicHead(numNodes) )
        Bwref = FlowModel%Bwref
        Pref =  FlowModel%Pwref
        Cw = FlowModel%Cw
        ViscosityW = FlowModel%Viscw
        Rhow = 0.0
        Bw = 0.0
        HydraulicHead = 0.0

        !************************************************************************
        !FRACTURE POROSITY.
        !************************************************************************
        !COMPUTE TOTAL FLUID VOUME (POROUS + FRACTURE)
        !FRACTURE APERTURE UPDATE AND COMPUTE PRESSURE VARIATION
        !-> needed to add storage for fracking non-reservoir cells
        ALLOCATE( fractureUnbalVolume(numNodes) )
        fractureUnbalVolume = 0.0
        FlowModel%NodeFracVolume = 0.0
        i = 0
        DO iz=1,GeomechanicalModel%GridDivision(2)+1
            DO ix=1,GeomechanicalModel%GridDivision(1)+1
                i = i + 1

                V0 = FlowModel%NodeFracVolume(i) !+ FlowModel%NodeVolume(i)   !previous volume

                dof = GeomechanicalModel%NodeDofs(2,i)
                FlowModel%NodeFracVolume(i) = HFmodel%NodeArea(i) * (HFmodel%NodeAperture0(i) + GeomechanicalModel%Displacement(dof) )

                fractureUnbalVolume(i) = FlowModel%NodeFracVolume(i) - V0

                !IF(.NOT.FlowModel%PressureBC_flag(i).AND. .NOT.FlowModel%InjNode(i) )THEN
                !    V1 = FlowModel%NodeVolume(i) + FlowModel%NodeFracVolume(i)  !current volume
                !    !Pressure drop due fracture dilation
                !    deltaP = 0.0
                !    IF(V0>0 .AND. Cw>0) deltaP = -1.0/V0 * (V1-V0) / Cw
                !    Porepressure(i) = Porepressure(i) + deltaP
                !    IF( Porepressure(i)<0 ) Porepressure(i) = 0.0
                !END IF

            END DO
        END DO

        !************************************************************************
        !BALANÇO DE MASSA DA FRATURA.
        !************************************************************************
        !O fluido que entra na fratura é o aumento do seu volume, dV/dt. Como a
        !mecânica é resolvida depois do fluxo (acoplamento sequencial), o volume
        !é escrito como
        !    V(p) = V_mec + C_f * (p - p_inicio)
        !com V_mec = volume da fratura da última solução mecânica (calculado
        !acima, em NodeFracVolume) e C_f = complacência estimada de cada nó da
        !face aberta. O termo C_f entra no armazenamento, implícito na pressão.
        !A diferença entre o volume da mecânica e o volume já contabilizado nas
        !chamadas anteriores (fracMismatch) entra como fonte no primeiro subpasso
        !aceito, de modo que o volume de fluido na fratura acompanha o volume da
        !mecânica, com um passo de atraso, e a massa é conservada.
        !Sem esse termo, a injeção a vazão controlada numa rocha pouco permeável
        !não tinha para onde ir além da compressibilidade do fluido, e a pressão
        !divergia.
        !C_f vem da solução de deformação plana, u_max = (1-nu^2) p Lc / E, com
        !Lc = MIN(altura, 2 x comprimento) da região aberta. O esquema é estável
        !para C_f maior que metade da complacência real; ele só afeta a
        !velocidade com que a diferença entre os volumes é corrigida.
        ALLOCATE( fracCompliance(numNodes),fracMismatch(numNodes),Pstart(numNodes) )
        fracCompliance = 0.0
        fracMismatch = 0.0
        Pstart = Porepressure
        mismatchPending = .FALSE.
        IF( HFmodel%ContainsFracModel )THEN
            ALLOCATE( faceElements(numNodes),openFaceElements(numNodes) )
            faceElements = 0
            openFaceElements = 0
            nOpen = 0
            planeStrainModulus = 0.0
            zLow = 1.0E30
            zHigh = -1.0E30
            xLow = 1.0E30
            xHigh = -1.0E30
            DO fracEl=1,HFmodel%NumElements
                num(1:4) = HFmodel%ListOfyyNodes(:,fracEl)
                faceElements(num(1:4)) = faceElements(num(1:4)) + 1
                IF( HFmodel%OpenFracElements(fracEl) )THEN
                    openFaceElements(num(1:4)) = openFaceElements(num(1:4)) + 1
                    nOpen = nOpen + 1
                    solidEl = HFmodel%FracToSolidElement(fracEl)
                    planeStrainModulus = planeStrainModulus + &
                        GeomechanicalModel%TransverIsotropicProps(solidEl)%E1 / &
                        (1.0 - GeomechanicalModel%TransverIsotropicProps(solidEl)%n1**2)
                    DO j=1,4
                        IF( GeomechanicalModel%Z(num(j)) < zLow ) zLow = GeomechanicalModel%Z(num(j))
                        IF( GeomechanicalModel%Z(num(j)) > zHigh ) zHigh = GeomechanicalModel%Z(num(j))
                        IF( GeomechanicalModel%X(num(j)) < xLow ) xLow = GeomechanicalModel%X(num(j))
                        IF( GeomechanicalModel%X(num(j)) > xHigh ) xHigh = GeomechanicalModel%X(num(j))
                    END DO
                END IF
            END DO
            IF( nOpen > 0 )THEN
                planeStrainModulus = planeStrainModulus / nOpen
                charLength = zHigh - zLow
                IF( 2.0*(xHigh - xLow) < charLength ) charLength = 2.0*(xHigh - xLow)
                DO i=1,numNodes
                    !nó livre: todos os elementos da face que o contêm estão abertos.
                    IF( faceElements(i) > 0 .AND. openFaceElements(i) == faceElements(i) ) &
                        fracCompliance(i) = HFmodel%NodeArea(i) * charLength / planeStrainModulus
                END DO
            END IF
            fracMismatch = FlowModel%NodeFracVolume - FlowModel%FracVolumeAccounted
            mismatchPending = .TRUE.
        END IF

        !************************************************************************
        !TIME STEP SOLVER.
        !************************************************************************
        t(0) = 0.00
        !t(1) = 0.10
        !t(2) = 0.20
        !t(3) = 0.50
        !t(4) = 0.75
        !t(5) = 1.00

        dt = DeltaTime !/ 5.0

        maxAllowedPressureVariation = 1.0E30
        !IF(WellData%InjectionPressure)THEN
        !    !maxAllowedPressureVariation = 1.0E30
        !    !maxAllowedPressureVariation = 1.05 * ABS( MAXVAL(Porepressure) - &
        !    !                                          MINVAL(Porepressure) )
        !    maxAllowedPressureVariation = 1.05 * ABS( MAXVAL(Porepressure) )
        !ELSE
        !    maxAllowedPressureVariation = 0.10 * MAXVAL(Porepressure)
        !END IF

        FlowModel%AverageTotalRate = 0.0

        IIcounter = 0
        iter = 0
        step = 0
        DO
            step = step + 1
            !t(step) = t(step) * DeltaTime
            !dt = t(step) - t(step-1)
            t(0) = t(0) + dt

            reportData%time = reportData%time + dt

            P0(1:) = Porepressure

            A = 0.0
            RHS = 0.0
            Afrac = 0.0

            !FLUID PROPS AND FLUID HEAD--------------------------------------------
            !!dir$ loop count min(16)
            DO i=1,numNodes
                Bw(i) = CalcBW(FlowModel,P0(i))
                Rhow(i) = CalcDENSITY(Bw(i),FlowModel)
                HydraulicHead(i) = CalcHYDRAULICPOTENTIAL(P0(i),Rhow(i),FlowModel%NodeDepth(i))
            END DO

            !UPDATE TRANSMISSIBILITY MATRIX AND RHS--------------------------------
            counter = 0
            DO i=1,FlowModel%ActiveToFlowConections(0)  !nTrans
                iel = FlowModel%ActiveToFlowConections(i)                           !<- active transmissibilility (iel -> transmissibility element index)
                dof1 = FlowModel%TransPointer(1,iel)
                dof2 = FlowModel%TransPointer(2,iel)
                IF(HydraulicHead(dof1)>HydraulicHead(dof2))THEN
                    dof = dof1
                ELSE
                    dof = dof2
                END IF
                mob = 1.0 / ( ViscosityW(dof) * Bw(dof) )
                A(iel) = mob * FlowModel%TransArray(iel)
                !gravitational term.
                RHS(dof1) = RHS(dof1) - A(iel) * Rhow(dof1) * ( FlowModel%NodeDepth(dof2) - FlowModel%NodeDepth(dof1) )
                RHS(dof2) = RHS(dof2) - A(iel) * Rhow(dof2) * ( FlowModel%NodeDepth(dof1) - FlowModel%NodeDepth(dof2) )
            END DO
            !Fracture transmissibility.
            DO i=1,HFmodel%ListOfActiveElements(0)!nFracTrans
                iel = HFmodel%ListOfActiveElements(i)
                dof1 = HFmodel%TransDOFs(1,iel)
                dof2 = HFmodel%TransDOFs(2,iel)
                IF(HydraulicHead(dof1)>HydraulicHead(dof2))THEN
                    dof = dof1
                ELSE
                    dof = dof2
                END IF
                mob = 1.0 / ( ViscosityW(dof) * Bw(dof) )
                Afrac(iel) = mob * HFmodel%Trans(iel)
                !gravitational term.
                RHS(dof1) = RHS(dof1) - Afrac(iel) * Rhow(dof1) * ( FlowModel%NodeDepth(dof2) - FlowModel%NodeDepth(dof1) )
                RHS(dof2) = RHS(dof2) - Afrac(iel) * Rhow(dof2) * ( FlowModel%NodeDepth(dof1) - FlowModel%NodeDepth(dof2) )
            END DO

            !CUMULATIVE TERM-------------------------------------------------------
            ! C = V * phi * cp
            DO i=1,numNodes
                totalCompressibility = Cw + FlowModel%NodalPoreCompressbility(i)
                V0 = FlowModel%NodeVolume(i)                !porous fluid volume
                V1 = FlowModel%NodeFracVolume(i)            !fracture fluid volume

                !C(i) = - V * totalCompressibility / Bw(i)
                !fracCompliance: dV/dp da fratura (balanço de massa da fratura).
                C(i) = -(V0*totalCompressibility + V1*Cw + fracCompliance(i)) / Bw(i)

                RHS(i) = RHS(i) + C(i) * P0(i) / dt
            END DO

            !Volume da fratura aberto pela mecânica e ainda não contabilizado:
            !sai do fluido no primeiro subpasso aceito.
            IF( mismatchPending )THEN
                DO i=1,numNodes
                    RHS(i) = RHS(i) + fracMismatch(i) / (Bw(i) * dt)
                END DO
            END IF

            !RATES-----------------------------------------------------------------
            temp = 0.0
            DO i=1,numNodes
                RHS(i) = RHS(i) - FlowModel%Rate(i)
            END DO

            !FAR FIELD (OUTFLOW)---------------------------------------------------
            !Option activated with #OUTFLOW
            IF(ALLOCATED(FlowModel%Jfactor))THEN
                DO i=1,numNodes
                    C(i) = C(i) - FlowModel%Jfactor(i)*dt
                    RHS(i) = RHS(i) - FlowModel%Pest(i)
                END DO
            END IF

            !FRACTURE DILATION-----------------------------------------------------
            !Option activated with #FRAC_UNBALVOLUME
            IF( FlowModel%FracUnbalVolume )THEN
                IF(step==1)THEN
                    V0 = SUM( fractureUnbalVolume )
                    IF( V0 > SUM(FlowModel%Rate) ) V0 = SUM(FlowModel%Rate)
                    DO i=1,numNodes
                        IF (FlowModel%InjNode(i)) THEN
                            RHS(i) = RHS(i) + FlowModel%RateFactor(i) * V0 / dt
                        END IF
                    END DO
                END IF
            END IF

            !MATRIX----------------------------------------------------------------
            FlowModel%a = 0.0

            DO i=1,FlowModel%ActiveToFlowConections(0)
                iel = FlowModel%ActiveToFlowConections(i)
                dof1 = FlowModel%TransPointer(1,iel)
                dof2 = FlowModel%TransPointer(2,iel)

                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof2,A(iel))
                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof1,-A(iel))

                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof1,A(iel))
                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof2,-A(iel))

            END DO

            DO i=1,HFmodel%ListOfActiveElements(0)
                iel = HFmodel%ListOfActiveElements(i)
                dof1 = HFmodel%TransDOFs(1,iel)
                dof2 = HFmodel%TransDOFs(2,iel)

                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof2,Afrac(iel))
                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof1,-Afrac(iel))

                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof1,Afrac(iel))
                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof2,-Afrac(iel))
            END DO

            !Add storage (C)
            DO i=1,numNodes
                dof1 = i
                CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof1,+C(dof1)/dt)
            END DO

            !PENALTY
            csrA_c = FlowModel%a
            CALL AddPenalty(FlowModel%a,FlowModel%ia,FlowModel%ja,RHS(1:),FlowModel%PressureBC_flag,penaltyValue,P0(1:))

            CALL CGpreCond(FlowModel%a,FlowModel%ia,FlowModel%ja,diag_precon)

            !LINEAR SOLVER---------------------------------------------------------
            ! ([A]+[C])*{p} = {RHS}
            !CALL PrintMessage(.FALSE.,.FALSE.,'...linear solver')
            solution(1:) = P0(1:)
            CALL PCG_CSR( cg_limit,neq,cg_tol,diag_precon,&
                            FlowModel%a,FlowModel%ia,FlowModel%ja,&
                            RHS,solution,nthreads,reportData%linIters)
            maxPressureVariation = MAXVAL( ABS( P0(1:) - solution(1:) ) )

            !COMPUTE RATES---------------------------------------------------------
            CALL mkl_dcsrgemv('N', numNodes, csrA_c,FlowModel%ia,FlowModel%ja, solution(1:), FlowModel%NodeFlows)
            FlowModel%NodeFlows = FlowModel%NodeFlows - C(1:)*P0(1:)/dt
            IF( mismatchPending ) FlowModel%NodeFlows = FlowModel%NodeFlows - fracMismatch/(Bw*dt)
            Rates = 0.0
            DO i=1,numNodes
                IF( FlowModel%InjNode(i) ) Rates(i) = FlowModel%NodeFlows(i)
            END DO
            FlowModel%TotalInjectionRate = SUM(Rates)
            !CALL ComputeRates(FlowModel,HFmodel,A,Afrac,Rhow,C(1:),dt,P0(1:),solution(1:),Rates)

            !PROCEED STEP DECISION-------------------------------------------------
            Qinj(2) = FlowModel%TotalInjectionRate
            IF( rateControled .AND.  IIcounter < 5 .AND.  &
                    ( ABS(Qinj(2)) >  ABS(maxRate*(1.0+limitRateDeviation)) .OR. &
                    ABS(Qinj(2)) <  ABS(maxRate*(1.0-limitRateDeviation)) ) )THEN

                    reportData%stepResult = 'R'
                    !rateControled = .TRUE.
                    IIcounter = IIcounter + 1
                    CALL InjectivityIndex(FlowModel,maxRate,Rhow,Porepressure,II,Pinj,Pest,Zref,Qinj)
                    t(0) = t(0) - dt

            ELSEIF( maxPressureVariation <= maxAllowedPressureVariation )THEN
                    reportData%stepResult = 'P'
                    Porepressure = solution(1:)
                    mismatchPending = .FALSE.

                    FlowModel%AverageTotalRate = FlowModel%AverageTotalRate + FlowModel%TotalInjectionRate * dt

                    dt = dt * 2.0
                    !update porous fluid volume
                    !c = +1/V * dV/dP
                    DO i=1,numNodes
                        deltaP = solution(i) - P0(i)
                        FlowModel%NodeVolume(i) = FlowModel%NodeVolume(i) * (1.0 + deltaP * FlowModel%NodalPoreCompressbility(i)) 
                    END DO
                    IIcounter = 0

            ELSE
                    reportData%stepResult = 'R'
                    IF(t(0)>0.0) t(0) = t(0) - dt
                    dt = dt / 2
                    !rateControled = .FALSE.
                    IIcounter = 0
            END IF

            !OUTPUT----------------------------------------------------------------
            reportData%structuralStep = 0
            reportData%fractureStep = 0
            reportData%nonLinIters = 0
            reportData%maxVAL = MAXVAL(solution(1:))/1.0E3
            CALL PrintStepReport(reportData%time,reportData%globalStep,reportData%flowStep,&
                                    reportData%structuralStep,reportData%fractureStep,&
                                    reportData%linIters,reportData%nonLinIters,&
                                    reportData%maxVAL,reportData%stepResult)
            IF(reportData%stepResult == 'R') reportData%time = reportData%time - 2.0 * dt

            !STOP CRITERIA---------------------------------------------------------
            !IF(t(step)>=DeltaTime) EXIT
            IF(t(0)>=DeltaTime) EXIT
            IF(t(0)+dt>DeltaTime) dt = DeltaTime - t(0)

        END DO!******************************************************************

        FlowModel%AverageTotalRate = FlowModel%AverageTotalRate / t(0)

        !Volume da fratura contabilizado até o fim desta chamada.
        IF( HFmodel%ContainsFracModel ) FlowModel%FracVolumeAccounted = &
            FlowModel%NodeFracVolume + fracCompliance * (Porepressure - Pstart)

        RETURN
        END SUBROUTINE

        ! TRANSIENT FLOW SOLVER++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SteadyStateFlowSolution( nthreads,FlowModel,HFmodel,GeomechanicalModel,WellData,&
                                            Porepressure )
        ! [A]*{p} = {RHS}

        !USE IFPORT
        USE IO
        USE ProblemDefinition
        USE DiscretizedModel
        USE HydraulicFractureModel
        USE PCG_Flow
        USE OMP_LIB  !bibliotecas OpenMP(R)

        IMPLICIT NONE
        INTEGER,INTENT(IN):: nthreads
        TYPE(Model),INTENT(IN):: GeomechanicalModel
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
        REAL(8),INTENT(IN OUT):: Porepressure(:)

        LOGICAL,ALLOCATABLE:: RestrictedDOFs(:),presureBC_flag(:)

        CHARACTER(LEN=320):: text
        INTEGER:: iter,maxIters=1
        INTEGER:: nTrans
        LOGICAL:: converged

        INTEGER:: i,j,k,pos,iel,counter,nCells,numNodes,neq,iz,irad,dof,dof1,dof2,cell,n,nnz
        REAL(8),ALLOCATABLE:: A(:),Afrac(:),RHS(:),diag(:),C(:)

        REAL(8),ALLOCATABLE:: P0(:),solution(:)
        REAL(8),ALLOCATABLE:: HydraulicHead(:)

        INTEGER:: num(8),g(4)

        REAL(8),ALLOCATABLE:: fracTrans(:,:)
        INTEGER,ALLOCATABLE:: fracTransJ(:,:)

        REAL(8):: V, phi, Cw, Bwref, Pref, Viscw, mob, grav
        REAL(8):: totalCompressibility
        REAL(8):: HydraulicPotential

        REAL(8):: min,temp

        REAL(8),ALLOCATABLE:: Bw(:), ViscosityW(:), Rhow(:)

        !Fracture Matrix.
        INTEGER:: nFracTrans
        REAL(8),ALLOCATABLE:: HFa(:)    !non-zero terms of transmissibility matrix (with BC).

        !Solver parameters.
        INTEGER:: cg_limit = 1.0E6, linIters
        REAL(8):: cg_tol = 1.0E-6

        REAL(8),ALLOCATABLE::diag_precon(:)  !pré-dondicionador.

        REAL(8),ALLOCATABLE:: initialGuess(:)

        !************************************************************************
        !POROUS MEDIA
        !************************************************************************
        nCells = GeomechanicalModel%NumElements
        numNodes = GeomechanicalModel%NumNodes
        neq = numNodes
        nTrans = SIZE(FlowModel%TransArray)
        nFracTrans = SIZE(HFmodel%Trans)

        ALLOCATE( RHS(neq), initialGuess(neq), solution(0:neq) )
        ALLOCATE( P0(0:neq),diag(0:neq) )
        !!dir$ loop count min(512)
        RHS = 0.0
        !!dir$ loop count min(512)
        diag = 0.0
        !!dir$ loop count min(512)
        solution = 0.0

        !Transmissibility / Mobility matrix.
        ALLOCATE( A(nTrans), Afrac(nFracTrans), C(0:numNodes) )
        !!dir$ loop count min(512)
        A = 0.
        Afrac = 0.
        C = 0.

        ALLOCATE(RestrictedDOFs(neq))
        RestrictedDOFs = .FALSE.
        DO i=1,numNodes
            IF(FlowModel%PressureBC_flag(i)) RestrictedDOFs(i) = .TRUE.
        END DO

        ALLOCATE(presureBC_flag(size(FlowModel%PressureBC_flag)))
        presureBC_flag = FlowModel%PressureBC_flag

        !************************************************************************
        !NODE MOBIITY POTENTIAL.
        !************************************************************************
        ALLOCATE( Bw(numNodes), ViscosityW(numNodes), Rhow(numNodes), HydraulicHead(numNodes) )
        Bwref = FlowModel%Bwref
        Pref =  FlowModel%Pwref
        Cw = FlowModel%Cw
        ViscosityW = FlowModel%Viscw
        Rhow = 0.0
        Bw = 0.0
        HydraulicHead = 0.0

        !************************************************************************
        !TIME STEP SOLVER.
        !************************************************************************

        maxIters = 1
        iter = 0
        DO
            iter = iter + 1

            P0(1:) = Porepressure

            !FLUID PROPS-----------------------------------------------------------
            !!dir$ loop count min(16)
            DO i=1,numNodes
                Bw(i) = CalcBW(FlowModel,P0(i))
                Rhow(i) = CalcDENSITY(Bw(i),FlowModel)
                HydraulicHead(i) = CalcHYDRAULICPOTENTIAL(P0(i),Rhow(i),FlowModel%NodeDepth(i))
            END DO

            !RHS-------------------------------------------------------------------
            !!dir$ loop count min(512)
            RHS = 0.0
            !!dir$ loop count min(256)
            DO i=1,numNodes
                RHS(i) = RHS(i) - FlowModel%Rate(i)
            END DO

            !UPDATE TRANSMISSIBILITY MATRIX AND RHS--------------------------------
            counter = 0
            A = 0.0
            DO i=1,FlowModel%ActiveToFlowConections(0)  !nTrans
                !IF(FlowModel%ActiveToFlowConections(i))THEN
                    iel = FlowModel%ActiveToFlowConections(i)                           !<- active transmissibilility
                    dof1 = FlowModel%TransPointer(1,iel)
                    dof2 = FlowModel%TransPointer(2,iel)
                    IF(HydraulicHead(dof1)>HydraulicHead(dof2))THEN
                        dof = dof1
                    ELSE
                        dof = dof2
                    END IF
                    mob = 1.0 / ( ViscosityW(dof) * Bw(dof) )
                    A(iel) = mob * FlowModel%TransArray(iel)
                    !gravitational term.
                    RHS(dof1) = RHS(dof1) - A(iel) * Rhow(dof1) * ( FlowModel%NodeDepth(dof2) - FlowModel%NodeDepth(dof1) )
                    RHS(dof2) = RHS(dof2) - A(iel) * Rhow(dof2) * ( FlowModel%NodeDepth(dof1) - FlowModel%NodeDepth(dof2) )
                !END IF
            END DO
            !Fracture transmissibility.
            Afrac = 0.0
            DO i=1,HFmodel%ListOfActiveElements(0)!nFracTrans
                iel = HFmodel%ListOfActiveElements(i)
                dof1 = HFmodel%TransDOFs(1,iel)
                dof2 = HFmodel%TransDOFs(2,iel)
                IF(HydraulicHead(dof1)>HydraulicHead(dof2))THEN
                    dof = dof1
                ELSE
                    dof = dof2
                END IF
                mob = 1.0 / ( ViscosityW(dof) * Bw(dof) )
                Afrac(iel) = mob * HFmodel%Trans(iel)
                !gravitational term.
                RHS(dof1) = RHS(dof1) - Afrac(iel) * Rhow(dof1) * ( FlowModel%NodeDepth(dof2) - FlowModel%NodeDepth(dof1) )
                RHS(dof2) = RHS(dof2) - Afrac(iel) * Rhow(dof2) * ( FlowModel%NodeDepth(dof1) - FlowModel%NodeDepth(dof2) )
            END DO

            !MATRIX----------------------------------------------------------------
            FlowModel%a = 0.0

            DO i=1,FlowModel%ActiveToFlowConections(0)
                iel = FlowModel%ActiveToFlowConections(i)
                dof1 = FlowModel%TransPointer(1,iel)
                dof2 = FlowModel%TransPointer(2,iel)
                IF(.NOT.RestrictedDOFs(dof1))THEN
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof2,A(iel))
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof1,-A(iel))
                END IF
                IF(.NOT.RestrictedDOFs(dof2))THEN
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof1,A(iel))
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof2,-A(iel))
                END IF
            END DO

            DO i=1,HFmodel%ListOfActiveElements(0)
                iel = HFmodel%ListOfActiveElements(i)
                dof1 = HFmodel%TransDOFs(1,iel)
                dof2 = HFmodel%TransDOFs(2,iel)
                IF(.NOT.RestrictedDOFs(dof1))THEN
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof2,Afrac(iel))
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof1,dof1,-Afrac(iel))
                END IF
                IF(.NOT.RestrictedDOFs(dof2))THEN
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof1,Afrac(iel))
                    CALL CSRAddValue(FlowModel%a,FlowModel%ia,FlowModel%ja,dof2,dof2,-Afrac(iel))
                END IF
            END DO

            CALL CSRdiaginaCheck(FlowModel%a,FlowModel%ia,FlowModel%ja,presureBC_flag)

            CALL CGpreCond(FlowModel%a,FlowModel%ia,FlowModel%ja,diag_precon)

            !LINEAR SOLVER---------------------------------------------------------
            ! [A]{p} = {RHS}
            !CALL PrintMessage(.FALSE.,.FALSE.,'...linear solver')
            solution(1:) = 0. !P0(1:)
            C(1:) = RHS
            CALL PCG_CSR( cg_limit,neq,cg_tol,diag_precon,&
                FlowModel%a,FlowModel%ia,FlowModel%ja,&
                C,solution,nthreads,linIters)
            Porepressure = Porepressure + solution(1:)                             !-> P0(1:) + solution

            !STOP CRITERIA---------------------------------------------------------
            IF(iter > maxIters) EXIT

        END DO!******************************************************************

        RETURN
        END SUBROUTINE

        SUBROUTINE InjectivityIndex(FlowModel,TargetRate,&
                                    Rho,Porepressure,II,Pinj,Pest,Zref,Qinj)
            USE DiscretizedModel

            IMPLICIT NONE
            TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
            REAL(8),INTENT(IN OUT):: Rho(:),Porepressure(:)
            REAL(8),INTENT(IN):: TargetRate
            REAL(8),INTENT(IN OUT):: II,Pinj(2),Zref,Pest,Qinj(2)

            INTEGER:: numNodes,i
            REAL(8):: newPinj,z,tang

            numNodes = SIZE(FlowModel%NodeDepth)

            DO i=1,numNodes
                IF(FlowModel%InjNode(i))THEN
                    Pinj(2) = Porepressure(i)
                    Zref = -FlowModel%NodeDepth(i)
                    EXIT
                END IF
            END DO

            IF( Pinj(1)==0.0) THEN
                !FIRST TRIAL
                Pinj(1) = Pinj(2)
                Qinj(1) = FlowModel%TotalInjectionRate
                Pest = FlowModel%StaticPressure + FlowModel%GradStaticPressure * ( FlowModel%StaticPressureZref - Zref )
                II = Qinj(2)/(Pinj(2)-Pest)
                Pinj(2) = Pest + TargetRate / II
            ELSE
                !NEWTON-RAPHSON
                tang = ( Qinj(2)-Qinj(1) ) / ( Pinj(2)-Pinj(1) )
                Pinj(2) = Pinj(1) + ( TargetRate - Qinj(1) ) / tang
            END IF

            DO i=1,numNodes
                IF(FlowModel%InjNode(i))THEN
                    z = FlowModel%NodeDepth(i)
                    newPinj = Pinj(2) + Rho(i) * 9.81 * ( Zref + z )
                    Porepressure(i) = newPinj
                END IF
            END DO

        RETURN
    END SUBROUTINE
    
    ! PRESSURE BOUNDARY CONDITIONS+++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE ExternalPressureBC(GeomechanicalModel,Xl,Xp,Yl,Yp,Zl,Zp,FlowModel)
        !For Radial:
        !   X -> RAD
        !   Y -> THETA

        USE DiscretizedModel

        IMPLICIT NONE
        LOGICAL,INTENT(IN):: Xl,Xp,Yl,Yp,Zl,Zp
        TYPE(Model),INTENT(IN):: GeomechanicalModel
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel

        INTEGER:: i,Nodes,nNodes
        REAL(8):: maxY,maxX,maxZ,minZ
        REAL(8):: tol=0.001

        INTEGER:: irad, NRAD, iz, NZ, itheta, NTHETA

        FlowModel%PressureBC_flag = .FALSE.

        nNodes = GeomechanicalModel%NumNodes

        !BOUNDARY DIFUSION OR DAMPING----------------------------------

        IF(FlowModel%RAD .EQV. .FALSE.)THEN

            maxX = MAXVAL(GeomechanicalModel%X)
            maxY = MAXVAL(GeomechanicalModel%Y)
            maxZ = MAXVAL(GeomechanicalModel%Z)
            minZ = MINVAL(GeomechanicalModel%Z)
            DO i=1,nNodes
                IF( ABS(GeomechanicalModel%X(i))<tol.AND.Xl ) FlowModel%PressureBC_flag(i) = .TRUE.
                IF( ABS(GeomechanicalModel%X(i)-maxX)<tol.AND.Xp ) FlowModel%PressureBC_flag(i) = .TRUE.
                IF( ABS(GeomechanicalModel%Y(i))<tol.AND.Yl ) FlowModel%PressureBC_flag(i) = .TRUE.
                IF( ABS(GeomechanicalModel%Y(i)-maxY)<tol.AND.Yp ) FlowModel%PressureBC_flag(i) = .TRUE.
                IF( ABS(GeomechanicalModel%Z(i)-maxZ)<tol.AND.Zp ) FlowModel%PressureBC_flag(i) = .TRUE.    !top
                IF( ABS(GeomechanicalModel%Z(i)-minZ)<tol.AND.Zl ) FlowModel%PressureBC_flag(i) = .TRUE.    !base
            END DO

        ELSE

            NRAD = GeomechanicalModel%GridDivision(1)
            NTHETA = GeomechanicalModel%GridDivision(3)
            NZ = GeomechanicalModel%GridDivision(2)

            i = 0
            DO itheta = 1, NTHETA+1
                DO iz = 1, NZ+1
                    DO irad=1, NRAD+1
                        i = i + 1
                        IF( irad==1.AND.Xl ) FlowModel%PressureBC_flag(i) = .TRUE.
                        IF( irad==(NRAD+1).AND.Xp ) FlowModel%PressureBC_flag(i) = .TRUE.
                        IF( itheta==1.AND.Yl ) FlowModel%PressureBC_flag(i) = .TRUE.
                        IF( itheta==(NTHETA+1).AND.Yp ) FlowModel%PressureBC_flag(i) = .TRUE.
                        IF( iz==1.AND.Zp ) FlowModel%PressureBC_flag(i) = .TRUE.    !top
                        IF( iz==(NZ+1).AND.Zl ) FlowModel%PressureBC_flag(i) = .TRUE.    !base
                    END DO
                END DO
            END DO

        END IF

        RETURN
    END SUBROUTINE

    !Well wall pressure b.c. for open hole well.
    SUBROUTINE InternalPressureBC(FlowModel,WellData,GeomechanicalModel,Step)

        !USE omp_lib
        USE ProblemDefinition
        USE IO
        USE HydraulicFractureModel
        USE DiscretizedModel

        IMPLICIT NONE
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        INTEGER,INTENT(IN):: Step

        INTEGER:: i,numNodes
        REAL(8):: tol = 0.001, radius, r0

        REAL(8):: Bw,Rhow,hydraulicHead_ref
        REAL(8):: Pinj,Zinj

        LOGICAL:: isTable

        numNodes = SIZE(FlowModel%NodeDepth(:))

        r0 = GeomechanicalModel%X(1)

        isTable = .FALSE.
        IF( SIZE(WellData%InjPressureTable(1,:))>0 ) isTable = .TRUE.

        IF( FlowModel%RAD .AND. Step>0 )THEN
            !RADIAL MODEL / VERTICAL WELL -------------------------------------------------------------
            DO i=1,numNodes
                radius = SQRT(GeomechanicalModel%X(i)**2.0+GeomechanicalModel%Y(i)**2.0)
                IF( radius<=r0+WellData%InjX                            .AND.&
                        (GeomechanicalModel%Z(i)<=WellData%InjZ1)       .AND.&
                        (GeomechanicalModel%Z(i)>=WellData%InjZ2) )THEN
                    Bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
                    Rhow = CalcDENSITY(Bw,FlowModel)
                    IF(isTable)THEN
                        Pinj = WellData%InjPressureTable(2,step)
                        Zinj = WellData%InjPressureTable(3,step)
                    ELSE
                        Pinj = WellData%Pinj
                        Zinj = WellData%Zinj
                    END IF
                    IF(Pinj>0.0)THEN
                        FlowModel%PressureBC_flag(i) = .TRUE.
                        FlowModel%InjNode(i) = .TRUE.
                        GeomechanicalModel%Porepressure(i) = Pinj + Rhow * 9.81 * (Zinj-GeomechanicalModel%Z(i))
                    ELSE
                        FlowModel%PressureBC_flag(i) = .FALSE.
                        FlowModel%InjNode(i) = .FALSE.
                    END IF
                END IF
            END DO
        ELSEIF( FlowModel%HorWellFlag)THEN
            !BOX MODEL / HORIZONTAL WELL --------------------------------------------------------------
            DO i = 1,numNodes
                IF(WellData%WellBoreWallFlag(i))THEN
                    Bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
                    Rhow = CalcDENSITY(Bw,FlowModel)
                    IF(isTable)THEN
                        Pinj = WellData%InjPressureTable(2,step)
                        Zinj = WellData%InjPressureTable(3,step)
                    ELSE
                        Pinj = WellData%Pinj
                        Zinj = WellData%Zinj
                    END IF
                    IF(Pinj>0.0)THEN
                        FlowModel%PressureBC_flag(i) = .TRUE.
                        FlowModel%InjNode(i) = .TRUE.
                        GeomechanicalModel%Porepressure(i) = Pinj + Rhow * 9.81 * (Zinj-GeomechanicalModel%Z(i))
                    ELSE
                        FlowModel%PressureBC_flag(i) = .FALSE.
                        FlowModel%InjNode(i) = .FALSE.
                    END IF
                END IF
            END DO
        ELSE
            IF( Step>0 )THEN
                !BOX MODEL / VERTICAL WELL ----------------------------------------------------------------
                DO i=1,numNodes
                    IF( (GeomechanicalModel%Y(i)<=WellData%InjX)       .AND.&
                        (GeomechanicalModel%X(i)<=WellData%InjX)       .AND.&
                        (GeomechanicalModel%Z(i)<=WellData%InjZ1)      .AND.&
                        (GeomechanicalModel%Z(i)>=WellData%InjZ2) )THEN
                    Bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
                    Rhow = CalcDENSITY(Bw,FlowModel)
                    IF(isTable)THEN
                        Pinj = WellData%InjPressureTable(2,step)
                        Zinj = WellData%InjPressureTable(3,step)
                    ELSE
                        Pinj = WellData%Pinj
                        Zinj = WellData%Zinj
                    END IF
                    IF(Pinj>0.0)THEN
                        FlowModel%PressureBC_flag(i) = .TRUE.
                        FlowModel%InjNode(i) = .TRUE.
                        GeomechanicalModel%Porepressure(i) = Pinj + Rhow * 9.81 * (Zinj-GeomechanicalModel%Z(i))
                    ELSE
                        FlowModel%PressureBC_flag(i) = .FALSE.
                        FlowModel%InjNode(i) = .FALSE.
                    END IF
                    END IF
                END DO
            END IF
        END IF

        RETURN
    END SUBROUTINE

    !Well wall flow b.c. for open hole well.
    SUBROUTINE InternalRateBC(FlowModel,WellData,GeomechanicalModel,Step)

        !USE omp_lib
        USE ProblemDefinition
        USE IO
        USE HydraulicFractureModel
        USE DiscretizedModel

        IMPLICIT NONE
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        INTEGER,INTENT(IN):: Step

        INTEGER:: i,numNodes

        REAL(8):: Bw,rate_std

        !REAL(8):: pi=3.14,k,h,mu,re,rw

        FlowModel%Rate = 0.0

        !SS RADIAL FLOW
        !q=2pikh/mu*(pe-pwf)/ln(re/rw)
        !FlowModel%WellIndex = 0.0
        !DO i=1,GeomechanicalModel%numElements
        !    
        !    IF (FlowModel%InjNode(i)) THEN
        !        k=
        !        h=
        !        mu=
        !        re=
        !        rw=
        !        FlowModel%WellIndex = FlowModel%WellIndex + 2.0*pi*k*h/mu/ln(re/rw)
        !    END IF
        !END DO

        !FlowModel%Rate(nodeList(1:n)) = WellData%Rate/n
        !FlowModel%InjNode(nodeList(1:n)) = .TRUE.

        IF( SIZE(WellData%InjRateTable(1,:))>0 ) THEN
            numNodes = GeomechanicalModel%NumNodes
            rate_std = WellData%InjRateTable(2,Step)

            DO i=1,numNodes
                IF (FlowModel%InjNode(i)) THEN
                    FlowModel%Rate(i) = rate_std * FlowModel%RateFactor(i)
                END IF
            END DO

        END IF
    RETURN
    END SUBROUTINE

    ! Builds sink/source based on nodes.
    SUBROUTINE InjectionRate(FlowModel,WellData,GeomechanicalModel)

        !USE omp_lib
        USE ProblemDefinition
        USE IO
        USE HydraulicFractureModel
        USE DiscretizedModel

        IMPLICIT NONE
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel

        INTEGER:: i,ix,iz,iy,numNodes,NX,NY,NZ,node,iel,num(8),n
        REAL(8):: tol = 0.001, radius, r0
        REAL(8):: sumK,z_avg,h

        REAL(8),ALLOCATABLE:: nodePerm(:)

        numNodes = GeomechanicalModel%NumNodes

        !r0 = GeomechanicalModel%X(1)
        IF( FlowModel%RAD )THEN
            !RADIAL MODEL / VERTICAL WELL -------------------------------------------------------------
            DO i=1,numNodes
                radius = SQRT(GeomechanicalModel%X(i)**2.0+GeomechanicalModel%Y(i)**2.0)
                !IF( radius<=r0+WellData%InjX                       .AND.&
                IF( radius<=WellData%InjX                          .AND.&
                    (GeomechanicalModel%Z(i)<=WellData%InjZ1)      .AND.&
                    (GeomechanicalModel%Z(i)>=WellData%InjZ2) )THEN
                    FlowModel%InjNode(i) = .TRUE.
                END IF
            END DO
        ELSEIF( FlowModel%HorWellFlag)THEN
            !BOX MODEL / HORIZONTAL WELL --------------------------------------------------------------
            DO i = 1,numNodes
                IF(WellData%WellBoreWallFlag(i))THEN
                    FlowModel%InjNode(i) = .TRUE.
                END IF
            END DO
        ELSE
            !BOX MODEL / VERTICAL WELL ----------------------------------------------------------------
            DO i=1,numNodes
                IF( (GeomechanicalModel%Y(i)<=WellData%InjX)       .AND.&
                    (GeomechanicalModel%X(i)<=WellData%InjX)       .AND.&
                    (GeomechanicalModel%Z(i)<=WellData%InjZ1)      .AND.&
                    (GeomechanicalModel%Z(i)>=WellData%InjZ2) )THEN
                    FlowModel%InjNode(i) = .TRUE.
                END IF
            END DO
        END IF

        n = 0
        DO node=1,numNodes
            IF(FlowModel%InjNode(node)) n = n + 1
        END DO
        IF(n==0)THEN
            CALL PrintMessage( .TRUE., .TRUE., 'No intervall for fluid injection (review data)' )
            STOP
        END IF

        !CALCULATION OF RATE FACTORS ASSOCIATED TO INJECTION NODES ()
        FlowModel%RateFactor = 0.0
        ALLOCATE( nodePerm(numNodes) )
        nodePerm = 0.0

        NX = GeomechanicalModel%GridDivision(1)
        NZ = GeomechanicalModel%GridDivision(2)
        NY = GeomechanicalModel%GridDivision(3)

        sumK = 0.0
        iel = 0
        DO iy=1,NY
            DO iz=1,NZ
                DO ix=1,NX

                    iel = iel + 1
                    num = GeomechanicalModel%ElementNodes(:,iel)
                    z_avg = SUM( GeomechanicalModel%Z(num) ) / 8.

                    IF( ix==1 .AND. &
                        (z_avg <=WellData%InjZ1) .AND.&
                        (z_avg >=WellData%InjZ2) )THEN

                        h = GeomechanicalModel%Z(num(2))-GeomechanicalModel%Z(num(1))
                        DO i=1,8
                            IF(FlowModel%InjNode(num(i)))THEN
                                nodePerm(num(i)) = nodePerm(num(i)) + GeomechanicalModel%Kh(iel) * h
                                sumK = sumK + GeomechanicalModel%Kh(iel) * h
                            END IF
                        END DO

                    END IF
                END DO
            END DO
        END DO

        IF(sumK==0.0)THEN
            sumK = n
            nodePerm = 1.0
        END IF
        DO node=1,numNodes
            IF(FlowModel%InjNode(node))THEN
                FlowModel%RateFactor(node) = nodePerm(node) / sumK
            END IF
        END DO

        RETURN
    END SUBROUTINE

    !Well wall flow b.c. for open hole well.
    !Stages of injection pressure.
    SUBROUTINE InternalPressureBC_STAGES(FlowModel,WellData,GeomechanicalModel,Stage)

        !USE omp_lib
        USE ProblemDefinition
        USE IO
        USE HydraulicFractureModel
        USE DiscretizedModel

        IMPLICIT NONE
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(WellboreData),INTENT(IN):: WellData
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        INTEGER,INTENT(IN):: Stage

        INTEGER:: i,numNodes
        REAL(8):: tol = 0.001, radius

        REAL(8):: Bw,Rhow,hydraulicHead_ref

        numNodes = SIZE(FlowModel%NodeDepth(:))

        !LIST OF NODES AND ELEMENTS.
        IF( FlowModel%RAD )THEN
            !RADIAL MODEL / VERTICAL WELL -------------------------------------------------------------
            DO i=1,numNodes
                radius = SQRT(GeomechanicalModel%X(i)**2.0+GeomechanicalModel%Y(i)**2.0)
                IF( radius<=WellData%InjX                          .AND.&
                    (GeomechanicalModel%Z(i)<=WellData%InjZ1)      .AND.&
                    (GeomechanicalModel%Z(i)>=WellData%InjZ2) )THEN
                FlowModel%PressureBC_flag(i) = .TRUE.
                FlowModel%InjNode(i) = .TRUE.
                Bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
                Rhow = CalcDENSITY(Bw,FlowModel)
                GeomechanicalModel%Porepressure(i) = WellData%InjPressureTable(2,Stage) + Rhow * 9.81 * &
                    (WellData%InjPressureTable(3,Stage)-GeomechanicalModel%Z(i))
                END IF
            END DO
        ELSEIF( FlowModel%HorWellFlag)THEN
            !BOX MODEL / HORIZONTAL WELL --------------------------------------------------------------
            DO i = 1,numNodes
                IF(WellData%WellBoreWallFlag(i))THEN
                    FlowModel%PressureBC_flag(i) = .TRUE.
                    FlowModel%InjNode(i) = .TRUE.
                    Bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
                    Rhow = CalcDENSITY(Bw,FlowModel)
                    GeomechanicalModel%Porepressure(i) = WellData%InjPressureTable(2,Stage) + Rhow * 9.81 * &
                        (WellData%InjPressureTable(3,Stage)-GeomechanicalModel%Z(i))
                END IF
            END DO
        ELSE
            !BOX MODEL / VERTICAL WELL ----------------------------------------------------------------
            DO i=1,numNodes
                IF(    (GeomechanicalModel%Y(i)<=WellData%InjX)       .AND.&
                    (GeomechanicalModel%X(i)<=WellData%InjX)       .AND.&
                    (GeomechanicalModel%Z(i)<=WellData%InjZ1)      .AND.&
                    (GeomechanicalModel%Z(i)>=WellData%InjZ2) )THEN
                FlowModel%PressureBC_flag(i) = .TRUE.
                FlowModel%InjNode(i) = .TRUE.
                Bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
                Rhow = CalcDENSITY(Bw,FlowModel)
                GeomechanicalModel%Porepressure(i) = WellData%InjPressureTable(2,Stage) + Rhow * 9.81 * &
                    (WellData%InjPressureTable(3,Stage)-GeomechanicalModel%Z(i))
                END IF
            END DO
        END IF

        !IF(n>0) FlowModel%Rate(nodeList(1:n)) = WellData%Rate/n

        RETURN
    END SUBROUTINE

    SUBROUTINE WellWallPressureBC(FlowModel,MudWeight,GeomechanicalModel,Z1,Z2)

        !USE omp_lib
        USE ProblemDefinition
        USE IO
        USE HydraulicFractureModel
        USE DiscretizedModel

        IMPLICIT NONE
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        REAL(8),INTENT(IN):: MudWeight
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        REAL(8),INTENT(IN):: Z1,Z2

        INTEGER:: i,numNodes
        REAL(8):: tol = 0.001, radius(2)

        REAL(8):: Bw,Rhow,hydraulicHead_ref

        numNodes = SIZE(FlowModel%NodeDepth(:))

        !RADIAL MODEL / VERTICAL WELL -------------------------------------------------------------
        radius(1) = 0.5* (GeomechanicalModel%X(1)+GeomechanicalModel%X(2) )
        DO i=1,numNodes
            radius(2) = SQRT(GeomechanicalModel%X(i)**2.0+GeomechanicalModel%Y(i)**2.0)
            IF( (GeomechanicalModel%Z(i)<=Z1) .AND. (GeomechanicalModel%Z(i)>=Z2) .AND. &
                (radius(2)<radius(1)) )THEN
            FlowModel%PressureBC_flag(i) = .TRUE.
            FlowModel%InjNode(i) = .TRUE.
            Rhow = MudWeight
            GeomechanicalModel%Porepressure(i) =  -Rhow * GeomechanicalModel%Z(i)
            END IF
        END DO

        RETURN
    END SUBROUTINE

    
    ! UPDATE CUMULATIVE TERM AND MOBILITY++++++++++++++++++++++++++++++++++++++
    ! Account permeability changes due stress or strain and sand production.
    SUBROUTINE UpdateTransmissibility(Problem,TimeStep,FlowModel,GeomechanicalModel)

        USE DiscretizedModel
        USE SolidMechanicsLib
        USE ProblemDefinition

        IMPLICIT NONE
        TYPE(PD),INTENT(IN):: Problem
        REAL(8),INTENT(IN):: TimeStep
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel

        INTEGER:: i,ip,iel,j,k,nCells,node1,node2,num(8),nip,cell,counter,NX,NZ

        REAL(8),ALLOCATABLE:: permMult(:),permDamageMult(:)
        REAL(8):: deltaPerm
        REAL(8):: lowPerm = 1.E-6 * 1.0*9.869D-16

        REAL(8):: qStrain, meanStrain, lode_theta

        LOGICAL,ALLOCATABLE:: liquefied(:)

        REAL(8):: der

        nip = GeomechanicalModel%NumIntegrationPoints
        nCells = GeomechanicalModel%NumElements
        ALLOCATE( liquefied(SIZE(FlowModel%IsLiquefied)) )
        liquefied = .FALSE.

        ALLOCATE( permMult(nCells), permDamageMult(nCells) )
        permMult = 1.0
        permDamageMult = 1.0

        ! INJECTION DAMAGE-----------------------------------------------

        NX = GeomechanicalModel%GridDivision(1)
        NZ = GeomechanicalModel%GridDivision(2)
        iel = 0.0
        DO k=1,NZ
            DO i=1,NX
                iel = iel + 1

                IF(GeomechanicalModel%IsFractured(iel))THEN

                    der = ( 1.0 - FlowModel%InjectionDamageModel(2,k) ) / FlowModel%InjectionDamageModel(1,k)
                    FlowModel%TimeOfInjection(iel) = FlowModel%TimeOfInjection(iel) + TimeStep
                    permDamageMult(iel) = 1.0 - der * FlowModel%TimeOfInjection(iel)
                    IF( permDamageMult(iel)<FlowModel%InjectionDamageModel(2,k) )&
                        permDamageMult(iel) = FlowModel%InjectionDamageModel(2,k)
                    GeomechanicalModel%DamagePerm(iel) = 1.0 - permDamageMult(iel)

                    !for output
                    GeomechanicalModel%KH(iel) = permDamageMult(iel) * FlowModel%OriginalPerm(1,iel)
                    GeomechanicalModel%KV(iel) = permDamageMult(iel) * FlowModel%OriginalPerm(2,iel)
                END IF

            END DO
        END DO

        ! SHEAR ENHANCED PERMEABILITY------------------------------------
        DO i=1,SIZE(GeomechanicalModel%PermModel)
            DO k=1,GeomechanicalModel%PermModel(i)%ElementList(0)
                iel = GeomechanicalModel%PermModel(i)%ElementList(k)

                !Perm multiplier as a function of plastic shear strain.
                deltaPerm = 0.0
                DO ip=1,nip
                    CALL invar(GeomechanicalModel%PlasticStrainTensor(:,ip,iel), meanStrain, qStrain, lode_theta)
                    deltaPerm = deltaPerm + GeomechanicalModel%PermModel(i)%PermIncrementByPlasticShear * qStrain
                END DO
                deltaPerm = deltaPerm / nip

                IF(GeomechanicalModel%KH(iel)>0.0 .AND. deltaPerm>0 .AND.&
                    GeomechanicalModel%KH(iel) .NE. lowPerm)THEN
                    !active do flow elements
                    permMult(iel) = deltaPerm + GeomechanicalModel%KH(iel)
                    IF(permMult(iel)>10000.*9.869D-16) permMult(iel) = 10000.*9.869D-16  !<<<LIMIT TO 10 Darcy)
                    permMult(iel) = permMult(iel) / GeomechanicalModel%KH(iel)

                ELSEIF(deltaPerm>0)THEN
                    !non-active do flow elements
                    GeomechanicalModel%KH(iel) = 2.0*lowPerm
                    GeomechanicalModel%KV(iel) = 2.0*lowPerm
                    IF(deltaPerm>10000.*9.869D-16) deltaPerm = 10000.*9.869D-16  !<<<LIMIT TO 10 Darcy)
                    permMult(iel) = deltaPerm / lowPerm
                ELSE
                END IF

                IF( permMult(iel)<1.0 ) permMult(iel) = 1.0

            END DO
        END DO

        ! TRANSMISSIBILITY-----------------------------------------------
        FlowModel%TransArray = 0.0
        DO i=1,SIZE(FlowModel%TransArray)
            iel = FlowModel%Trans2Element(i)
            FlowModel%TransArray(i) = permDamageMult(iel) * permMult(iel) * FlowModel%TransArrayINI(i)
        END DO

        !UPDATE ACTIVE TO FLOW CONECTIONS
        FlowModel%ActiveToFlowConections = 0
        counter = 0
        DO iel=1,SIZE(FlowModel%TransArray)
            cell = FlowModel%Trans2Element(iel)
            IF( FlowModel%TransArray(iel).GT.0.0 .AND. &
                GeomechanicalModel%KH(cell) .GT. lowPerm .AND. &
                GeomechanicalModel%KH(cell)*permMult(cell) .GT. lowPerm )THEN   !<<<LIMIT TO 100 mD)
            counter = counter + 1
            FlowModel%ActiveToFlowConections(counter) = iel
            END IF
        END DO
        FlowModel%ActiveToFlowConections(0) = counter

        RETURN
    END SUBROUTINE

    ! CALCULATES FLUID PRESURE FROM HYDRAULIC HEAD+++++++++++++++++++++++++++++
    SUBROUTINE RecoverPressureFromHydraulicHead(HydraulicHead,FlowModel,Bw,Rhow,&
                                                Zdown,Porepressure)

        IMPLICIT NONE
        REAL(8),INTENT(IN):: HydraulicHead(:)
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        REAL(8),INTENT(IN OUT):: Bw(:),Rhow(:),Zdown(:)
        REAL(8),INTENT(IN OUT):: Porepressure(:)

        INTEGER:: i,numNodes,iter

        REAL(8):: Cw, Bwref, Pref

        numNodes = SIZE(Porepressure)

        !Hydraulic Head.
        Bwref = FlowModel%Bwref
        Pref =  FlowModel%Pwref
        Cw = FlowModel%Cw
        Rhow = 0.0
        Bw = 0.0
        DO i=1,numNodes
            DO iter=1,5
                Bw(i) = CalcBW(FlowModel,Porepressure(i))
                Rhow(i) = CalcDENSITY(Bw(i),FlowModel)
                Porepressure(i) = HydraulicHead(i) + Rhow(i) * 9.81 * Zdown(i)
            END DO
        END DO

        RETURN
    END SUBROUTINE

    ! CALCULATES FORMATION VOLUME FACTOR+++++++++++++++++++++++++++++++++++++++
    FUNCTION CalcBW(FlowModel,Porepressure)

        IMPLICIT NONE
        REAL(8):: CalcBW
        TYPE(PorousFlowModel),INTENT(IN):: FlowModel
        REAL(8),INTENT(IN):: Porepressure

        REAL(8):: Cw, Bwref, Pref

        Pref =  FlowModel%Pwref
        Bwref = FlowModel%Bwref
        Cw = FlowModel%Cw

        CalcBW = ( 1.0/Bwref + Cw/Bwref * (Porepressure-Pref) )**(-1.0)

    END FUNCTION

    ! CALCULATES FLUID DENSITY+++++++++++++++++++++++++++++++++++++++++++++++++
    FUNCTION CalcDENSITY(Bw,FlowModel)

        IMPLICIT NONE
        REAL(8):: CalcDENSITY
        TYPE(PorousFlowModel),INTENT(IN):: FlowModel
        REAL(8),INTENT(IN):: Bw

        CalcDENSITY = 1.0/Bw * FlowModel%Rhoref

    END FUNCTION

    ! CALCULATES HYDRAULIC POTENTIAL+++++++++++++++++++++++++++++++++++++++++++
    FUNCTION CalcHYDRAULICPOTENTIAL(Porepressure,Rho,Zdown)

        IMPLICIT NONE
        REAL(8):: CalcHYDRAULICPOTENTIAL
        REAL(8),INTENT(IN):: Porepressure,Rho,Zdown

        CalcHYDRAULICPOTENTIAL = Porepressure - Rho * 9.81 * Zdown !at zero datum.

    END FUNCTION

    ! INITIALIZE FLOW MODEl++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE InitializeFlowModel(GeomechanicalModel,Problem,FlowModel)

        USE omp_lib
        USE IO
        USE NumericalUtilities
        USE FEMlib
        USE DiscretizedModel
        USE ProblemDefinition

        IMPLICIT NONE
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE (PD), INTENT(IN):: Problem
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel

        INTEGER:: i,numNodes,iter

        REAL(8):: Cw, Bwref, Pref, POTref
        REAL(8),ALLOCATABLE:: Bw(:),Rhow(:)

        numNodes = GeomechanicalModel%NumNodes

        ALLOCATE( Bw(numNodes), Rhow(numNodes) )

        !Hydraulic Head.
        !Cw = FlowModel%Cw
        !Rhow = 0.0
        !Bw = 0.0
        !POTref = CalcHYDRAULICPOTENTIAL(Problem%OWpressure,Problem%nwGrad,&
        !                                -Problem%OWdepth)
        DO i=1,numNodes
            !DO iter=1,5
            !    Bw(i) = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
            !    Rhow(i) = CalcDENSITY(Bw(i),FlowModel)
            !    GeomechanicalModel%Porepressure(i) = POTref + (Rhow(i)*9.81) * FlowModel%NodeDepth(i)
            !END DO
            GeomechanicalModel%Porepressure(i) = Problem%OWpressure*1000.0 + (Problem%nwGrad*1000.0) * ( Problem%OWdepth + FlowModel%NodeDepth(i) )
        END DO

        FlowModel%StaticPressure = Problem%OWpressure*1000.0
        FlowModel%StaticPressureZref = Problem%OWdepth
        FlowModel%GradStaticPressure = Problem%nwGrad*1000.0

        RETURN
    END SUBROUTINE

    ! BUILDS POROUS FLOW MODEL+++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE BuildPorousFlowModel(GeomechanicalModel,Problem,FlowModel)

        USE IO
        USE NumericalUtilities
        USE FEMlib
        USE DiscretizedModel
        USE ProblemDefinition

        IMPLICIT NONE
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE (PD), INTENT(IN):: Problem
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel

        INTEGER,ALLOCATABLE:: DofNodes(:,:)

        INTEGER:: Ncells,Nnodes,iel,i,j,counter,cell,n,ind(3),nodes(4),NZ

        CHARACTER(LEN=320)::element             !'hexahedron' ou 'tetrahedron'.
        INTEGER,ALLOCATABLE:: num(:)            !element incidence.
        REAL(8):: det
        REAL(8),ALLOCATABLE:: g_coord(:,:),&    !nodal coordenates.
        coord(:,:)        !element nodal coordenates.
        REAL(8),ALLOCATABLE:: points(:,:),&     !coordenadas dos pontos de integração.
        weights(:)        !peso dos pontos de integração.
        REAL(8),ALLOCATABLE:: jac(:,:),der(:,:),deriv(:,:)

        REAL(8):: elemVOL

        REAL(8):: PermTensor(3,3)
        REAL(8),ALLOCATABLE:: elTransm(:,:)
        INTEGER:: dof
        REAL(8):: lowPerm = 1.E-6 * 1.0*9.869D-16

        INTEGER:: ndim = 3                      !number of dimensions: 3D -> 3.
        INTEGER:: nod = 8                       !number of nodes per element.
        INTEGER:: nip = 8                       !number of integration points.

        REAL(8):: LX,LY,LZ,AX,AY,AZ,K,L,area,T,TX,TY,TZ
        REAL(8):: coordX1_avg,coordX2_avg,coordY1_avg,coordY2_avg
        INTEGER:: pointers(12,2)

        INTEGER,ALLOCATABLE:: numPos(:)

        REAL(8):: Pest, Jfactor

        INTEGER:: k1,k2
        REAL(8):: time,damage

        IF(TRIM(Problem%MeshType)=='RAD')FlowModel%RAD = .TRUE.
        IF(Problem%WellData%HorWellFlag)FlowModel%HorWellFlag = .TRUE.

        Ncells = GeomechanicalModel%NumElements
        Nnodes = GeomechanicalModel%NumNodes
        nod = GeomechanicalModel%NumNodesPerElement
        ndim = GeomechanicalModel%NumDimensions

        ALLOCATE( FlowModel%Accum(Nnodes),FlowModel%NodeDepth(Nnodes),&
            FlowModel%NodeVolume(Nnodes),FlowModel%NodeFracVolume(Nnodes),&
            FlowModel%FracVolumeAccounted(Nnodes),&
            FlowModel%ElementVolume(Ncells),FlowModel%NodalPoreCompressbility(Nnodes),&
            FlowModel%Rate(Nnodes),FlowModel%PressureBC_flag(Nnodes) )
        FlowModel%FracVolumeAccounted = 0.0
        !, FlowModel%PressureBC(Nnodes)
        FlowModel%Rate = 0.0
        !FlowModel%PressureBC = 0.0
        FlowModel%PressureBC_flag = .FALSE.
        FlowModel%Accum = 0.0
        FlowModel%NodeVolume = 0.0
        FlowModel%NodeDepth = ABS(GeomechanicalModel%Z)
        FlowModel%ElementVolume = 0.0
        FlowModel%NodalPoreCompressbility = 0.0
        FlowModel%NodeFracVolume = 0.0

        !FLUID DATA.
        FlowModel%Bwref = Problem%water_bwI
        FlowModel%Pwref = Problem%water_pref
        FlowModel%Cw = Problem%water_cw
        FlowModel%Rhoref = Problem%density_pref
        FlowModel%Viscw = Problem%water_viscI

        !NODAL COORDENATES.
        ALLOCATE(g_coord(3,Nnodes))
        g_coord = 0.0
        g_coord(1,:) = GeomechanicalModel%X
        g_coord(2,:) = GeomechanicalModel%Y
        g_coord(3,:) = GeomechanicalModel%Z

        FlowModel%TotalInjectionRate = 0.0

        !ELEMENT INTEGRATION-----------------------------------------------------
        ALLOCATE(points(nip,ndim),coord(nod,ndim),  &
            jac(ndim,ndim),weights(nip),num(nod),der(ndim,nod),    &
            deriv(ndim,nod),elTransm(nod,nod) )
        element = GeomechanicalModel%element
        CALL sample(element,points,weights)
        counter = 0 !active to flow element counter.
        DO iel=1,Ncells

            PermTensor = 0.0
            PermTensor(1,1) = GeomechanicalModel%KH(iel)
            PermTensor(2,2) = PermTensor(1,1)
            PermTensor(3,3) = GeomechanicalModel%KV(iel)

            IF( (GeomechanicalModel%KH(iel)>0.0).OR.(GeomechanicalModel%KH(iel)>0.0) ) counter = counter + 1

            num = GeomechanicalModel%ElementNodes(:,iel)
            coord=TRANSPOSE(g_coord(:,num))

            !FlowModel%TransDOFs(:,iel) = num    !node=pressure dof.

            !ELEMENT VOLUME (FINITE ELEMENT)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            !elTransm = 0.0
            elemVOL = 0.0
            DO i=1,nip
                CALL shape_der(der,points,i)
                jac = MATMUL(der,coord)
                CALL determinante(jac,det)
                !CALL invert(jac)
                !deriv = MATMUL(jac,der)
                elemVOL = elemVOL + det * weights(i)
                !elTransm = elTransm + MATMUL(MATMUL(TRANSPOSE(deriv),PermTensor),deriv)*det*weights(i)
            END DO
            !FlowModel%Trans(:,:,iel) = elTransm
            !Bw @ P0 -> FlowModel%Bwref e FlowModel%Pwref
            FlowModel%ElementVolume(iel) = elemVol
            FlowModel%Accum(num) = FlowModel%Accum(num) + elemVol * GeomechanicalModel%POR(iel) / nod
            FlowModel%NodeVolume(num) = FlowModel%NodeVolume(num) + GeomechanicalModel%POR(iel) * elemVol/8.0

            !COMPRESSIBILITY
            FlowModel%NodalPoreCompressbility(num) = FlowModel%NodalPoreCompressbility(num) + &
                                                        GeomechanicalModel%PoreCompressibility(iel) * elemVol / nod
        END DO

        DO i=1,Nnodes
            IF(FlowModel%NodeVolume(i) > 0.0)THEN
                FlowModel%NodalPoreCompressbility(i) = FlowModel%NodalPoreCompressbility(i) / FlowModel%NodeVolume(i)
            ELSE
                FlowModel%NodalPoreCompressbility(i) = 0.0
            END IF
        END DO

        !NODE TO NODE TRANSMISSIBILITIES-----------------------------------------
        !REAL(8),ALLOCATABLE:: TransArray(:)        !transmissibility array (all nodes transmissibilities).
        !INTEGER,ALLOCATABLE:: TransPointer(:,:)    !node pair transmissibilities for each TransArray.

        ALLOCATE( FlowModel%TransArray(Ncells*12), FlowModel%TransPointer(2,Ncells*12),&
            FlowModel%IsLiquefied(Ncells*12), FlowModel%Trans2Element(Ncells*12),&
            FlowModel%TransAreas(6,Ncells))
        ALLOCATE( FlowModel%TransArrayINI(SIZE(FlowModel%TransArray)))
        FlowModel%TransArray = 0.0
        FlowModel%TransPointer = 0.0
        FlowModel%IsLiquefied = .FALSE.
        FlowModel%Trans2Element = 0
        FlowModel%TransAreas = 0.0

        pointers(1,1)=2  ; pointers(1,2)=3
        pointers(2,1)=2  ; pointers(2,2)=1
        pointers(3,1)=2  ; pointers(3,2)=6
        pointers(4,1)=3  ; pointers(4,2)=7
        pointers(5,1)=3  ; pointers(5,2)=4
        pointers(6,1)=1  ; pointers(6,2)=5
        pointers(7,1)=1  ; pointers(7,2)=4
        pointers(8,1)=4  ; pointers(8,2)=8
        pointers(9,1)=6  ; pointers(9,2)=7
        pointers(10,1)=6 ; pointers(10,2)=5
        pointers(11,1)=7 ; pointers(11,2)=8
        pointers(12,1)=5 ; pointers(12,2)=8

        counter = 0
        DO iel=1,Ncells

            num = GeomechanicalModel%ElementNodes(:,iel)
            coord=TRANSPOSE(g_coord(:,num))

            !IF(GeomechanicalModel%KH(iel)==0.)THEN
            !  GeomechanicalModel%KH(iel) = lowPerm
            !  GeomechanicalModel%KV(iel) = lowPerm
            !END IF
            PermTensor = 0.0
            PermTensor(1,1) = GeomechanicalModel%KH(iel)
            PermTensor(2,2) = PermTensor(1,1)
            PermTensor(3,3) = GeomechanicalModel%KV(iel)

            LZ = ABS( GeomechanicalModel%Z(num(1))-GeomechanicalModel%Z(num(2)) )
            !LX = ABS( GeomechanicalModel%X(num(1))-GeomechanicalModel%X(num(4)) )
            LX = SQRT( ( GeomechanicalModel%X(num(1))-GeomechanicalModel%X(num(4)) )**2.0 + &
                       ( GeomechanicalModel%Y(num(1))-GeomechanicalModel%Y(num(4)) )**2.0  )
            !LY = ABS( GeomechanicalModel%Y(num(1))-GeomechanicalModel%Y(num(5)) )
            !LY = SQRT( (GeomechanicalModel%X(num(1))-GeomechanicalModel%X(num(5)) )**2.0 + &
            !           (GeomechanicalModel%Y(num(1))-GeomechanicalModel%Y(num(5)) )**2.0  )
            coordX1_avg = 0.5 * ( GeomechanicalModel%X(num(1))+GeomechanicalModel%X(num(4)) )
            coordX2_avg = 0.5 * ( GeomechanicalModel%X(num(5))+GeomechanicalModel%X(num(8)) )
            coordY1_avg = 0.5 * ( GeomechanicalModel%Y(num(1))+GeomechanicalModel%Y(num(4)) )
            coordY2_avg = 0.5 * ( GeomechanicalModel%Y(num(5))+GeomechanicalModel%Y(num(8)) )
            LY = SQRT( ( coordX2_avg - coordX1_avg )**2.0 + ( coordY2_avg - coordY1_avg )**2.0  )
            AX = LY * LZ
            AY = LX * LZ
            AZ = LX * LY

            FlowModel%TransAreas(1,iel) = AX
            FlowModel%TransAreas(2,iel) = AY
            FlowModel%TransAreas(3,iel) = AZ
            FlowModel%TransAreas(4,iel) = LX
            FlowModel%TransAreas(5,iel) = LY
            FlowModel%TransAreas(6,iel) = LZ

            ! Nodal trasnmissibilities (1/4)
            TX = PermTensor(1,1) * AX/4.0 / LX
            TY = PermTensor(2,2) * AY/4.0 / LY
            TZ = PermTensor(3,3) * AZ/4.0 / LZ

            ! conection 1.
            counter = counter + 1
            FlowModel%TransArray(counter) = TX
            FlowModel%TransPointer(1,counter) = num( pointers(1,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(1,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 2.
            counter = counter + 1
            FlowModel%TransArray(counter) = TZ
            FlowModel%TransPointer(1,counter) = num( pointers(2,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(2,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 3.
            counter = counter + 1
            FlowModel%TransArray(counter) = TY
            FlowModel%TransPointer(1,counter) = num( pointers(3,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(3,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 4.
            counter = counter + 1
            FlowModel%TransArray(counter) = TY
            FlowModel%TransPointer(1,counter) = num( pointers(4,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(4,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 5.
            counter = counter + 1
            FlowModel%TransArray(counter) = TZ
            FlowModel%TransPointer(1,counter) = num( pointers(5,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(5,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 6.
            counter = counter + 1
            FlowModel%TransArray(counter) = TY
            FlowModel%TransPointer(1,counter) = num( pointers(6,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(6,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 7.
            counter = counter + 1
            FlowModel%TransArray(counter) = TX
            FlowModel%TransPointer(1,counter) = num( pointers(7,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(7,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 8.
            counter = counter + 1
            FlowModel%TransArray(counter) = TY
            FlowModel%TransPointer(1,counter) = num( pointers(8,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(8,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 9.
            counter = counter + 1
            FlowModel%TransArray(counter) = TX
            FlowModel%TransPointer(1,counter) = num( pointers(9,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(9,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 10.
            counter = counter + 1
            FlowModel%TransArray(counter) = TZ
            FlowModel%TransPointer(1,counter) = num( pointers(10,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(10,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 11.
            counter = counter + 1
            FlowModel%TransArray(counter) = TZ
            FlowModel%TransPointer(1,counter) = num( pointers(11,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(11,2) )
            FlowModel%Trans2Element(counter) = iel
            ! conection 12.
            counter = counter + 1
            FlowModel%TransArray(counter) = TX
            FlowModel%TransPointer(1,counter) = num( pointers(12,1) )
            FlowModel%TransPointer(2,counter) = num( pointers(12,2) )
            FlowModel%Trans2Element(counter) = iel

        END DO

        !Saving initial transmissibility.
        FlowModel%TransArrayINI = FlowModel%TransArray

        !FLAG FO ACTIVE TO FLOW CONECTIONS---------------------------------------
        ALLOCATE(FlowModel%ActiveToFlowConections(0:SIZE(FlowModel%TransArray)))
        FlowModel%ActiveToFlowConections = 0
        counter = 0
        DO iel=1,SIZE(FlowModel%TransArray)
            cell = FlowModel%Trans2Element(iel)
            IF( FlowModel%TransArray(iel).GT.0.0 .AND. &
                GeomechanicalModel%KH(cell) .NE. lowPerm )THEN
            counter = counter + 1
            FlowModel%ActiveToFlowConections(counter) = iel
            END IF
        END DO
        FlowModel%ActiveToFlowConections(0) = counter

        !INITIALIZING NODAL RATES------------------------------------------------
        ALLOCATE( FlowModel%ConectionFlowRates(SIZE(FlowModel%TransArray)) )
        FlowModel%ConectionFlowRates = 0.0
        ALLOCATE( FlowModel%NodeFlows(GeomechanicalModel%NumNodes), &
                  FlowModel%InjNode(GeomechanicalModel%NumNodes),&
                  FlowModel%RateFactor(GeomechanicalModel%NumNodes) )
        FlowModel%NodeFlows = 0.0
        FlowModel%InjNode = .FALSE.
        FlowModel%RateFactor = 0.0

        !FRACTURE UMBALANCED VOLUME OPTION---------------------------------------
        FlowModel%FracUnbalVolume = Problem%FracUnbalVolume

        !OUTFLOW ARRAY-----------------------------------------------------------
        IF( ALLOCATED(Problem%FarFieldFlow) )THEN
            ALLOCATE(FlowModel%Jfactor(GeomechanicalModel%NumNodes),&
                FlowModel%Pest(GeomechanicalModel%NumNodes))
            FlowModel%Jfactor = 0.0
            FlowModel%Pest = 0.0
            DO i=1,SIZE(Problem%FarFieldFlow)
                DO iel=1,GeomechanicalModel%NumElements
                    ind = GeomechanicalModel%ElementIJK(:,iel)

                    IF(ind(1)>=Problem%FarFieldFlow(i)%I1 .AND.&
                       ind(1)<=Problem%FarFieldFlow(i)%I2 .AND.&
                       ind(3)>=Problem%FarFieldFlow(i)%J1 .AND.&
                       ind(3)<=Problem%FarFieldFlow(i)%J2 .AND.&
                       ind(2)>=Problem%FarFieldFlow(i)%K1 .AND.&
                       ind(2)<=Problem%FarFieldFlow(i)%K2)THEN

                        num = GeomechanicalModel%ElementNodes(:,iel)
                        IF(Problem%FarFieldFlow(i)%DIR==1)THEN
                            nodes(1) = num(3)
                            nodes(2) = num(4)
                            nodes(3) = num(7)
                            nodes(4) = num(8)
                            area = FlowModel%TransAreas(1,iel)
                        ELSEIF(Problem%FarFieldFlow(i)%DIR==2)THEN
                            nodes(1) = num(5)
                            nodes(2) = num(6)
                            nodes(3) = num(7)
                            nodes(4) = num(8)
                            area = FlowModel%TransAreas(2,iel)
                        ELSE
                        END IF

                        DO j=1,4
                            Jfactor = Problem%FarFieldFlow(i)%J / (24*60*60) / 1000.0       !m³/d /kPa /m² -> m³/sec / Pa / m²
                            FlowModel%Jfactor(nodes(j)) = FlowModel%Jfactor(nodes(j)) + Jfactor * area / 4.0
                            Pest = Problem%FarFieldFlow(i)%Pest*1000.0 + (Problem%nwGrad*1000.0) * &
                                    ( Problem%FarFieldFlow(i)%Zref_Pest + FlowModel%nodeDepth(nodes(j)) )
                            FlowModel%Pest(nodes(j)) =FlowModel%Pest(nodes(j)) + Pest * ( Jfactor * area / 4.0 )

                        END DO

                    END IF
                END DO
            END DO
        END IF

        !MATRICES----------------------------------------------------------------
        !REAL(8),ALLOCATABLE:: a(:)
        !INTEGER,ALLOCATABLE:: ia(:),ja(:)
        ALLOCATE( DofNodes(1,Nnodes) )
        DO i=1,Nnodes
            DofNodes(1,i) = i
        END DO

        CALL SparseMatrixCSR( FlowModel%TransPointer, DofNodes, &
            FlowModel%a, FlowModel%ia, FlowModel%ja )

        !PERM DAMAGE--------------------------------------------------------------
        NZ = GeomechanicalModel%GridDivision(2)
        ALLOCATE( FlowModel%TimeOfInjection(Ncells), FlowModel%InjectionDamageModel(2,NZ),&
                    FlowModel%OriginalPerm(2,Ncells))
        FlowModel%TimeOfInjection = 0.0
        FlowModel%InjectionDamageModel(1,:) = 1.0
        FlowModel%InjectionDamageModel(2,:) = 1.0
        FlowModel%OriginalPerm = 0.0

        DO i=1,SIZE( Problem%FluidPermDamage(1,:) )
            READ(Problem%FluidPermDamage(1,i),*) K1
            READ(Problem%FluidPermDamage(2,i),*) K2
            READ(Problem%FluidPermDamage(3,i),*) time
            READ(Problem%FluidPermDamage(4,i),*) damage
            FlowModel%InjectionDamageModel(1,K1:K2) = time*60
            FlowModel%InjectionDamageModel(2,K1:K2) = damage
        END DO
        DO iel=1,Ncells
            FlowModel%OriginalPerm(1,iel) = GeomechanicalModel%KH(iel)
            FlowModel%OriginalPerm(2,iel) = GeomechanicalModel%KV(iel)
        END DO

        RETURN
    END SUBROUTINE
    
    ! BUILDS POROUS FLOW MODEL+++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE BuildPipeModel(GeomechanicalModel,Problem,FlowModel)
        !PipeModel
        !    REAL(8):: Diameter
        !    INTEGER,ALLOCATABLE:: GlobalNodesNumber(:)  !nodes ordered
        !    REAL(8),ALLOCATABLE:: PipeTrans(:)          !transmissibility for each pipe element
        
        USE DiscretizedModel
        USE ProblemDefinition

        IMPLICIT NONE
        TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
        TYPE (PD), INTENT(IN):: Problem
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel

        INTEGER:: i,n,iel,II,JJ,KK,NI,NJ,NK
        INTEGER:: num(8)
        
        NI=GeomechanicalModel%GridDivision(1)
        NK=GeomechanicalModel%GridDivision(2)
        NJ=GeomechanicalModel%GridDivision(3)
        
        FlowModel%ContainsPipeModel = .TRUE.
        
        FlowModel%Pipe%Diameter = Problem%Pipe%Diameter
        
        n = SIZE( Problem%Pipe%WellPath(1,:) )
        
        ALLOCATE( FlowModel%Pipe%GlobalNodesNumber(n+1) )
        FlowModel%Pipe%GlobalNodesNumber = 0.0
        
        DO i=1,n
            II=Problem%Pipe%WellPath(2,i)
            JJ=Problem%Pipe%WellPath(3,i)
            KK=Problem%Pipe%WellPath(4,i)
            
            iel = II + (KK-1)*NI + (JJ-1)*(NI*NK)
            num = GeomechanicalModel%ElementNodes(:,iel)
            
            IF( Problem%Pipe%WellPath(1,i)==1 )THEN         !X
            ELSEIF( Problem%Pipe%WellPath(1,i)==2 )THEN     !Y
            ELSEIF( Problem%Pipe%WellPath(1,i)==3 )THEN     !Z: vertical direction
                FlowModel%Pipe%GlobalNodesNumber(i) = num(2)
            ELSE
            END IF
            
            IF(i==n) FlowModel%Pipe%GlobalNodesNumber(i+1) = num(1)
            
        END DO
        
        RETURN
    END SUBROUTINE

    SUBROUTINE ComputeRates(FlowModel,HFmodel,A,Afrac,Rhow,C,dt,Press0,Press,Rates)
    
        USE HydraulicFractureModel
        
        IMPLICIT NONE
        TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
        TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
        REAL(8),INTENT(IN):: Rhow(:),Press0(:),Press(:),A(:),C(:),Afrac(:)
        REAL(8),INTENT(IN):: dt
        REAL(8),INTENT(IN OUT):: Rates(:)
        
        INTEGER:: i,iel,dof1,dof2,numNodes
        REAL(8):: rate,Pot1,Pot2
        
        Rates = 0.0
        
        numNodes = SIZE(FlowModel%NodeDepth)
        
        !mob = 1.0 / ( ViscosityW(dof) * Bw(dof) )
        !A(iel) = mob * FlowModel%TransArray(iel)
        !Q = A * (Pot2 - Pot1)
        
        FlowModel%ConectionFlowRates = 0.0
        
        !Matrix
        DO i=1,FlowModel%ActiveToFlowConections(0)  !nTrans
            iel = FlowModel%ActiveToFlowConections(i)                           !<- active transmissibilility (iel -> transmissibility element index)
            dof1 = FlowModel%TransPointer(1,iel)
            dof2 = FlowModel%TransPointer(2,iel)
            
            IF( .NOT.FlowModel%InjNode(dof1) .OR. .NOT.FlowModel%InjNode(dof2) )THEN
                
                Pot1 = Press(dof1) - Rhow(dof1) * FlowModel%NodeDepth(dof1)
                Pot2 = Press(dof2) - Rhow(dof2) * FlowModel%NodeDepth(dof2)
                rate = A(iel) * (Pot2 - Pot1)
                    
                IF( FlowModel%InjNode(dof1) ) Rates(dof1) = Rates(dof1) + rate
                IF( FlowModel%InjNode(dof2) ) Rates(dof2) = Rates(dof2) - rate
            
                FlowModel%ConectionFlowRates(iel) = FlowModel%ConectionFlowRates(iel) + rate
                                
                FlowModel%NodeFlows(dof1) = FlowModel%NodeFlows(dof1) + rate
                FlowModel%NodeFlows(dof2) = FlowModel%NodeFlows(dof2) - rate   
                
                IF( FlowModel%InjNode(dof1) ) Rates(dof1) = Rates(dof1) + rate
                IF( FlowModel%InjNode(dof2) ) Rates(dof2) = Rates(dof2) - rate
                
            END IF
        
        END DO
        
        !Fracture
        DO i=1,HFmodel%ListOfActiveElements(0)!nFracTrans
            iel = HFmodel%ListOfActiveElements(i)
            dof1 = HFmodel%TransDOFs(1,iel)
            dof2 = HFmodel%TransDOFs(2,iel)
            
            IF( .NOT.FlowModel%InjNode(dof1) .OR. .NOT.FlowModel%InjNode(dof2) )THEN

                Pot1 = Press(dof1) - Rhow(dof1) * FlowModel%NodeDepth(dof1)
                Pot2 = Press(dof2) - Rhow(dof2) * FlowModel%NodeDepth(dof2)
                rate = Afrac(iel) * (Pot2 - Pot1)
                
                HFmodel%ConectionFlowRates(iel) = HFmodel%ConectionFlowRates(iel) + rate
                
                FlowModel%NodeFlows(dof1) = FlowModel%NodeFlows(dof1) + rate
                FlowModel%NodeFlows(dof2) = FlowModel%NodeFlows(dof2) - rate   
                    
                IF( FlowModel%InjNode(dof1) ) Rates(dof1) = Rates(dof1) + rate
                IF( FlowModel%InjNode(dof2) ) Rates(dof2) = Rates(dof2) - rate
                
            END IF
            
        END DO
        
        FlowModel%TotalInjectionRate = SUM(Rates)
    
    RETURN
    END SUBROUTINE

    !Linear Solver for fluid flow.
    !LINE ELEMENT (2 NODES) BASED VERSION
    ![LHS]*{x}={RHS}
    !PARALLEL VERSION.
    SUBROUTINE PCG_CSR(cg_limit,neq,cg_tol,diag_precon,a,ia,ja,loads,initialGuess,nthreads,linSolverIters)
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
        REAL(8),INTENT(IN OUT)::initialGuess(:)  !vetor com a solução inicial de entrada no método
        INTEGER,INTENT(IN)::neq,&  !número de equações
        cg_limit  !número máximo de iterações do solver
        REAL(8),INTENT(IN)::cg_tol  !tolerância para convergência do solver linear
        INTEGER,INTENT(IN)::ia(:),&  !identificação da posição de início das linhas no vetor a(:)
        ja(:)  !identificação das colunas dos elementos não-zero da matriz
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
        !loads(1+RestrictedDOFs) = 0.0

        !Aplicando pré-condicionador diagonal no resíduo da equação de equilíbrio
        d = zero
        CALL vdmul( neq, diag_precon(2:neq+1), loads(2:neq+1), d )  !MKL VMF
        !d(RestrictedDOFs) = 0.0
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
            !u(RestrictedDOFs) = 0.0

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
            !d(RestrictedDOFs) = 0.0
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
        initialGuess = loads

        !PRINT STUFF.
        !WRITE(text,'(A,I6,A,E10.3,A,E10.3)')'...linear iters.:',cg_iters,&
        !    ' pres.(kPa):',MINVAL( ABS(loads(2:(neq+1))) )/1000.0,&
        !    ' TO ',MAXVAL( ABS(loads(2:(neq+1))) )/1000.0
        !CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))

        !FINALIZANDO E ESCREVENDO RESULTADOS DA TELEMETRIA DO SOLVER
        IF (cg_iters==cg_limit)THEN
            CALL PrintMessage(.TRUE.,.TRUE.,"LINEAR SOLVER DID NOT CONVERGE...")
            CALL PrintMessage(.TRUE.,.TRUE.,"")
            STOP
        END IF

        !DESALOCANDO MEMÓRIA
        !DEALLOCATE(d,p,u,x,xnew)

        linSolverIters = cg_iters

        RETURN
    END SUBROUTINE

    SUBROUTINE CGpreCond(a,ia,ja, diag_precon)
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

        IF(.NOT.ALLOCATED(diag_precon))ALLOCATE(diag_precon(0:neq))
        !!dir$ loop count min(512)
        diag_precon = 0.0

        !!dir$ loop count min(16)
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
    END SUBROUTINE CGpreCond

    !Armazenamento de matriz e ponteiros para os solvers iterativos intel e rotinas BLAS 2/ BLAS 3.
    !O esquema de armazenamento é o CSR (compressed sparse row).
    SUBROUTINE SparseMatrixCSR( ConectionList, DofNodes, a, ia, ja )

        USE IFPORT
        USE DiscretizedModel
        USE NumericalUtilities
        USE FEMlib
        USE StiffnessMatrix
        USE IO
        USE OMP_LIB

        IMPLICIT NONE

        INTEGER,INTENT(IN):: ConectionList(:,:)
        INTEGER,INTENT(IN):: DofNodes(:,:)
        REAL(8),ALLOCATABLE:: a(:)
        INTEGER,ALLOCATABLE:: ia(:), ja(:)  !Ponteiros que armazenam as posições dos dofs na matriz.

        INTEGER:: neq  !número de equações.

        !VARIÁVEIS ESTÁTICAS.

        INTEGER:: i,j,k,l,m,iel,&  !contadores.
        posicao,coluna,&  !auxílio a localização das posições no vetor "a".
        no,&  !auxiliar para guardar número do nó.
        GL1,GL2,GL3  !marcadores de grau de liberdade.

        INTEGER:: nels,&  !número de elementos.
        nn,&  !número de nós.
        nod,&  !número de nós por elemento.
        ndof !número de graus de liberdade por elemento.

        INTEGER::nnz  !número de elementos não nulos na matriz de rigidez (diagonal + matriz traingular superior).

        INTEGER:: position, matrixPosition

        INTEGER,ALLOCATABLE:: g(:),&  !vetor com graus de liberdade de um elemento.
        num(:) ,&  !conectividade de um elemento.
        degrees(:,:),&  !associação graus de liberdade e elementos.
        posicoes(:),conection(:,:),conectionOrdenado(:)  !vetores auxiliares para construção de conecção entre graus de liberdade.
        INTEGER(4),ALLOCATABLE:: listaAfinidades(:,:)  !lista de afinidades nodais (8*3,nn), 1a posição é um ponteiro.

        INTEGER:: globalPosition
        REAL(8):: value

        !INTEGER, external :: compar

        !INICIALIZAÇÃO DE VARIÁVEIS.
        nels = SIZE(ConectionList(1,:))
        nn = SIZE(DofNodes(1,:))
        nod = SIZE(ConectionList(:,1))
        ndof = SIZE(DofNodes(:,1))

        neq = nn*ndof

        !ALOCAÇÃO DE MEMÓRIA.
        ALLOCATE( num(nod), posicoes(0:neq),&
            conection(27+200,1), conectionOrdenado(27+199),&
            degrees(8,0:neq), listaAfinidades(2*(8+5),nn) )
        degrees = 0
        listaAfinidades = 0

        !Montando vetores de armazenamento da matriz de rigidez > a(:),ia(:),ja(:)

        DO i=1,nels
            num = ConectionList(:,i)
            DO j=1,nod
                DO k=1,nod
                    CALL ListaAd(listaAfinidades,num(k),num(j))
                    CALL ListaAd(listaAfinidades,num(k),num(j))
                END DO
            END DO
        END DO

        nnz=0
        DO i=1,nn

            !Marcando os nós ligados ao nó "i".
            k=0
            DO j=2,listaAfinidades(1,i)  !looping na lista de afinidades do nó i
                no = listaAfinidades(j,i)

                DO m=1,ndof
                    IF( DofNodes(m,no).NE.0 ) k=k+1
                END DO

            END DO

            !graus de liberdade do nó i -> "l"
            l=0
            DO m=1,ndof
                IF( DofNodes(m,i).NE.0 ) l=l+1
            END DO

            !graus de liberdade ligados ao nó = número de graus de liberdade do nó X graus de liberdade dos nós ligados.
            nnz = nnz + l * k

        END DO

        ALLOCATE( a(nnz) )
        a = 0.0

        ALLOCATE( ia(neq+1), ja(nnz) )
        ia = 0
        ja = 0
        ia(neq+1) = nnz+1

        !Montando matriz de número de posições por grau de liberdade.
        !Posicoes(i=0:neq) -> marca quantos graus de liberdade possuem afinidade com os graus de liberdade do nó "i".
        DO i=1,nn

            k=0  !marca o número de graus de liberdade associados a "k"
            DO j=2,listaAfinidades(1,i)

                no=listaAfinidades(j,i)

                GL1=0
                GL2=0
                GL3=0
                GL1=DofNodes(1,no); IF(GL1.NE.0) k=k+1
                !GL2=nf(2,no); IF(GL2.NE.0) k=k+1
                !GL3=nf(3,no); IF(GL3.NE.0) k=k+1

            END DO

            GL1=DofNodes(1,i); posicoes(GL1)=k
            !GL2=nf(2,i); posicoes(GL2)=k
            !GL3=nf(3,i); posicoes(GL3)=k

        END DO

        !Acumula posições
        ia(1) = 1
        k=1
        DO i=2,neq
            k=k+posicoes(i-1)
            ia(i)=k
        END DO

        !Formação dos vetores globais ia e ja.
        DO i=1,nn

            !grau de liberdade 1 do nó i.

            DO m=1,ndof

                IF(DofNodes(m,i).NE.0)THEN
                    conection=0
                    conectionOrdenado=0

                    DO j=2,listaAfinidades(1,i)

                        no=listaAfinidades(j,i)

                        GL1=DofNodes(1,no); IF(GL1.NE.0) CALL ListaAd(conection,GL1,1)
                        !GL2=nf(2,no); IF(GL2.NE.0) CALL ListaAd(conection,GL2,1)
                        !GL3=nf(3,no); IF(GL3.NE.0) CALL ListaAd(conection,GL3,1)

                    END DO

                    k=27+200
                    conectionOrdenado=conection(2:k,1)
                    CALL SORTQQ(LOC(conectionOrdenado),k-1,SRT$INTEGER4)
                    !CALL qsort( conectionOrdenado, K-1, 4, compar )

                    GL1=DofNodes(m,i)
                    k=ia(GL1)

                    l=-1
                    !DO j=1,(27*3+99)  !dimensão do array conectionOrdenado(:,1)
                    DO j=1,(27+199)  !dimensão do array conectionOrdenado(:,1)
                        IF(conectionOrdenado(j).NE.0)THEN
                            l=l+1
                            coluna=conectionOrdenado(j)
                            ja(k+l)=coluna
                        END IF
                    END DO

                END IF

            END DO

        END DO

        RETURN
    END SUBROUTINE SparseMatrixCSR

    SUBROUTINE CSRAddValue(a,ia,ja,dof1,dof2,Value)

        IMPLICIT NONE
        REAL(8),INTENT(IN OUT):: a(:)
        INTEGER,INTENT(IN):: ia(:), ja(:)  !Ponteiros que armazenam as posições dos dofs na matriz.
        INTEGER,INTENT(IN):: dof1,dof2
        REAL(8),INTENT(IN):: Value

        INTEGER:: i,j,posicao,coluna

        !DOF1--------------------------------------------------------------------
        !Transferindo a matriz de rigidez local para a global.
        posicao = ia(dof1)
        !Encontra a coluna relativa ao dof.
        coluna=-1
        DO
            coluna=coluna+1
            IF(ja(posicao+coluna)==dof2) EXIT
        END DO
        !Armazena o valor do dof na posição correta.
        a(posicao+coluna) = a(posicao+coluna) + Value

        !!DOF2--------------------------------------------------------------------
        !!Transferindo a matriz de rigidez local para a global.
        !posicao = ia(dof2)
        !!Encontra a coluna relativa ao dof.
        !coluna=-1
        !DO
        !  coluna=coluna+1
        !  IF(ja(posicao+coluna)==dof1) EXIT
        !END DO
        !!Armazena o valor do dof na posição correta.
        !a(posicao+coluna) = a(posicao+coluna) + Value

    END SUBROUTINE

    SUBROUTINE AddPenalty(a,ia,ja,RHS,PenaltyFlag,PenaltyValue,PrescribedValues)

        IMPLICIT NONE
        REAL(8),INTENT(IN OUT):: a(:),RHS(:),PrescribedValues(:)
        INTEGER,INTENT(IN):: ia(:), ja(:)  !Ponteiros que armazenam as posições dos dofs na matriz.
        LOGICAL,INTENT(IN):: PenaltyFlag(:)
        REAL(8),INTENT(IN):: PenaltyValue

        INTEGER:: i,neq,posicao,coluna

        neq = SIZE(ia)-1

        DO i=1,neq

            posicao = ia(i)

            !Find collumn of dof i
            coluna=-1
            DO
                coluna=coluna+1
                IF(ja(posicao+coluna)==i) EXIT
            END DO

            !Applying penalty.
            IF(PenaltyFlag(i))THEN
                a(posicao+coluna) = a(posicao+coluna) + PenaltyValue
                RHS(i) = a(posicao+coluna) * PrescribedValues(i)
            END IF

        END DO

    END SUBROUTINE

    SUBROUTINE CSRdiaginaCheck(a,ia,ja,PressureBC_flag)

        IMPLICIT NONE
        REAL(8),INTENT(IN OUT):: a(:)
        INTEGER,INTENT(IN):: ia(:), ja(:)  !Ponteiros que armazenam as posições dos dofs na matriz.
        LOGICAL,INTENT(IN OUT):: PressureBC_flag(:)

        INTEGER:: i,neq,posicao,coluna

        neq = SIZE(ia)-1

        DO i=1,neq

            posicao = ia(i)

            !Encontra a coluna relativa ao dof.
            coluna=-1
            DO
                coluna=coluna+1
                IF(ja(posicao+coluna)==i) EXIT
            END DO

            !Armazena o valor do dof na posição correta.
            IF(a(posicao+coluna)==0.0)THEN
                a(posicao+coluna) = 1.0
                PressureBC_flag(i) = .TRUE.
            END IF

        END DO

    END SUBROUTINE

    END MODULE
