
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  !******************************************
  !THERMAL MODEL.
  !******************************************

  
  MODULE Thermal
  
    TYPE, PUBLIC:: ThermalModel
      
      LOGICAL:: Numerical                           !.TRUE. for Finite Volume.
                                                    !.FALSE. for analytical (PG).
      
      REAL(8):: Tinj,Zinj,GradTinj                  !input data
      
      REAL(8),ALLOCATABLE:: Thickness(:)            !layer thickness.
      
      !DYNAMIC DATA.
      REAL(8),ALLOCATABLE:: RockTemperature(:)      !nodal temperature.
      REAL(8),ALLOCATABLE:: FluidTemperature(:)     !nodal temperature.
      
      LOGICAL,ALLOCATABLE:: BC_flag(:)
      
      !TRANSMISSIBILITY.
      REAL(8),ALLOCATABLE:: TransArray(:)           !CURRENT transmissibility array 
      INTEGER,ALLOCATABLE:: TransPointer(:,:)       !node pair transmissibilities for each TransArray.
      REAL(8),ALLOCATABLE:: TransPOR(:)
      
      !PROPERTIES.
      REAL(8),ALLOCATABLE:: FluidCP(:)
      REAL(8),ALLOCATABLE:: RockCP(:)
      REAL(8),ALLOCATABLE:: h(:)                    !rock-fluid heat transfer function.
      REAL(8),ALLOCATABLE:: NodeVolume(:),NodeRockDensity(:),NodePOR(:)
    
    END TYPE
  
    CONTAINS
  
      
    SUBROUTINE ThermalStepSolution( Nthreads,GeomechanicalModel,FlowModel,TempModel,&
                                    DeltaTime,Temperature,reportData )

      USE IO
      USE ProblemDefinition
      USE DiscretizedModel
      USE FluidFlowSolver
      USE PCG_Flow

      IMPLICIT NONE
      INTEGER,INTENT(IN):: Nthreads
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      TYPE(PorousFlowModel),INTENT(IN):: FlowModel
      TYPE(ThermalModel),INTENT(IN OUT):: TempModel
      REAL(8),INTENT(IN):: DeltaTime                                               
      REAL(8),INTENT(IN OUT):: Temperature(:)
      TYPE(StepReport),INTENT(IN OUT):: reportData
    
      REAL(8):: dt,time

      CHARACTER(LEN=320):: text
      INTEGER:: iter,maxIters=1
      INTEGER:: nTrans
      LOGICAL:: converged

      INTEGER:: i,j,k,iel,counter,nCells,numNodes,neq,&
                iz,irad,dof,dof1,dof2,cell,n,num(8)
      REAL(8),ALLOCATABLE:: A(:),RHS(:),diag(:)
      REAL(8),ALLOCATABLE:: T0(:),solution(:)
        
      INTEGER:: node1,node2,numConections

      REAL(8):: V,cr,phi,rho,qfluid,bw,cp,deltaTemp,z
      
      REAL(8):: min,temp
      
      !Solver parameters.
      INTEGER:: cg_limit = 1.0E7
      REAL(8):: cg_tol = 1.0E-8
      REAL(8):: one=1.0
      
      REAL(8),ALLOCATABLE:: initialGuess(:)
      REAL(8):: maxTemperatureVariation
      REAL(8):: maxAllowedTemperatureVariation
      
      REAL(8),ALLOCATABLE:: csr_a(:)
      INTEGER,ALLOCATABLE:: csr_ia(:),csr_ja(:),DofNodes(:,:)
      
      REAL(8),ALLOCATABLE:: diag_precon(:)  !pré-dondicionador.
      INTEGER,ALLOCATABLE:: dofs(:),bcDofs(:)

      CALL PrintMessage( .FALSE., .FALSE., '...Thermal Solver' )
      
      nCells = GeomechanicalModel%NumElements
      numNodes = GeomechanicalModel%NumNodes
      neq = numNodes
      nTrans = SIZE(TempModel%TransArray)
      numConections = SIZE(FlowModel%TransArray)
      
      !CSR matrix
      ALLOCATE( DofNodes(1,numNodes) )
      DO i=1,numNodes
          DofNodes(1,i) = i
      END DO
      CALL SparseMatrixCSR( FlowModel%TransPointer, DofNodes, csr_a, csr_ia, csr_ja )
      
      !************************************************************************
      !POROUS MEDIA
      !************************************************************************
          
      ALLOCATE( RHS(0:neq), initialGuess(neq), solution(0:neq) )
      ALLOCATE( T0(0:neq),diag(0:neq) )
      RHS = 0.0
      diag = 0.0
      solution = 0.0
        
      !Transmissibility / Mobility matrix.
      ALLOCATE( A(numNodes) )
      A = 0.0
      
      DO iel=1,GeomechanicalModel%NumElements
          IF( GeomechanicalModel%IsFractured(iel) )THEN
              num = GeomechanicalModel%ElementNodes(:,iel)
              DO i=1,4
                  dof = num(i)
                  z = GeomechanicalModel%Z(dof) - TempModel%Zinj
                  TempModel%FluidTemperature(dof) = TempModel%Tinj - z * TempModel%GradTinj
              END DO
          END IF
      END DO
      
      !************************************************************************
      !DEGREES OF FREEDON
      !************************************************************************
      ALLOCATE( dofs(0:numNodes), bcDofs(0:numNodes) )
      dofs = 0
      j = 0
      k = 0
      DO i=1,numNodes
          IF( .NOT.TempModel%BC_flag(i) )THEN
              j = j + 1
              dofs(j) = i
          ELSE
              k = k + 1
              bcDofs(k) = i
          END IF
      END DO
      dofs(0) = j
      bcDofs(0) = k
      
      !************************************************************************
      !TIME STEP SOLVER.
      !************************************************************************
      
      maxAllowedTemperatureVariation = 0.50 * ABS( MAXVAL(Temperature) - &
                                                   MINVAL(Temperature) )
      
      dt = 0.25 * DeltaTime
      iter = 0.0
      time = 0.0
      counter = 0
      DO
        time = time + dt
        csr_a = 0.0
        A = 0.0
        RHS = 0.0
        
        T0(1:) = Temperature
                
        !CUMULATIVE TERM (C)---------------------------------------------------
        DO i=1,numNodes
            bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(i))
            rho = CalcDENSITY(bw,FlowModel)
            V = TempModel%NodeVolume(i)
            cp = TempModel%FluidCP(i)
            A(i) = A(i) + V * cp * rho
        END DO
        A( bcDofs(1:k) ) = one
        !DO i=1,numNodes
        !    dof1 = i
        !    IF(.NOT.TempModel%BC_flag(dof1))THEN
        !        CALL CSRAddValue(csr_a,csr_ia,csr_ja,dof1,dof1,A(dof1))
        !    ELSE
        !        CALL CSRAddValue(csr_a,csr_ia,csr_ja,dof1,dof1,one)
        !    END IF
        !END DO
        
        !CONDUTION-------------------------------------------------------------
        ! -> already done in TempModel%TransArray and TempModel%TransPointer
        !DO iel=1,numConections
        !    dof1 = TempModel%TransPointer(1,iel)
        !    dof2 = TempModel%TransPointer(2,iel)
        !    IF(.NOT.TempModel%BC_flag(dof1))THEN
        !        CALL CSRAddValue(csr_a,csr_ia,csr_ja,dof1,dof2,TempModel%TransArray(iel))
        !        CALL CSRAddValue(csr_a,csr_ia,csr_ja,dof1,dof1,-TempModel%TransArray(iel))
        !    END IF
        !    IF(.NOT.TempModel%BC_flag(dof2))THEN
        !        CALL CSRAddValue(csr_a,csr_ia,csr_ja,dof2,dof1,TempModel%TransArray(iel))
        !        CALL CSRAddValue(csr_a,csr_ia,csr_ja,dof2,dof2,-TempModel%TransArray(iel))
        !    END IF
        !END DO
        
        !ADVECTION (RHS)-------------------------------------------------------
        DO i=1,FlowModel%ActiveToFlowConections(0)
            iel = FlowModel%ActiveToFlowConections(i)
            node1 = FlowModel%TransPointer(1,iel)
            node2 = FlowModel%TransPointer(2,iel)
            qfluid = FlowModel%ConectionFlowRates(iel)                              !POROUS flow from node 1 to node 2
            
            ! qfluid = rho * cp * q * T_movingFluid
            deltaTemp = T0(node1) - T0(node2)
            IF( qfluid<0.0 ) THEN
                bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(node1))
                rho = CalcDENSITY(bw,FlowModel)
                qfluid = -rho * TempModel%FluidCP(node1) * qfluid * deltaTemp
                RHS(node2) = RHS(node2) + qfluid
            ELSE
                bw = CalcBW(FlowModel,GeomechanicalModel%Porepressure(node2))
                rho = CalcDENSITY(bw,FlowModel)
                qfluid = rho * TempModel%FluidCP(node2) * qfluid * deltaTemp
                RHS(node1) = RHS(node1) + qfluid
            END IF
            
        END DO
        
        !FLUID-ROCK HEAT EXCHANGE (RHS)----------------------------------------
        !DO i=1,numNodes
        !    RHS(i) = RHS(i) - TempModel%H(i) * ( TempModel%RockTemperature(i)-&
        !                                         TempModel%FluidTemperature(i) )
        !END DO
        
        !BOUNDARY CONDITION----------------------------------------------------
        i = bcDofs(0)
        RHS( bcDofs(1:k) ) = 0.0
        !DO i=1,numNodes
        !    IF(TempModel%BC_flag(i)) RHS(i) = 0.0
        !END DO
        
        !RHS INTEGRATION-------------------------------------------------------
        RHS = RHS * dt
        
        !CSR MATRIX------------------------------------------------------------
        CALL CGpreCond(csr_a,csr_ia,csr_ja,diag_precon)
        
        !LINEAR SOLVER---------------------------------------------------------
        ![A+TempModel%TransPointer]{x} = {RHS}
        initialGuess = 0.0
        solution = 0.0
        !CALL PCG_FlowLineElementBased_noFrac(cg_limit,cg_tol,TempModel%TransArray,&
        !                                        TempModel%TransPointer,A,TempModel%BC_flag,&
        !                                        solution,initialGuess)
        
        !CALL PCG_CSR( cg_limit,numNodes,cg_tol,diag_precon,&
        !              csr_a,csr_ia,csr_ja,RHS,solution,Nthreads,reportData%linIters)
        
        DO i=1,dofs(0)
            !IF( .NOT.TempModel%BC_flag(i) .AND. A(i).NE.0.0 )THEN
            j = dofs(i)
            solution(j) = RHS(j) / A(j)
            !END IF
        END DO
        
        !solution <- delta temperature
        solution(1:) = Temperature + solution(1:)

        !STEPPING AND STOP CRITERIA--------------------------------------------
        maxTemperatureVariation = MAXVAL(ABS(Temperature - solution(1:)))
        
        IF( maxTemperatureVariation < maxAllowedTemperatureVariation )THEN
            counter = 0
            Temperature = solution(1:)
            TempModel%FluidTemperature = Temperature
            TempModel%RockTemperature = TempModel%FluidTemperature
            IF( time>=DeltaTime ) EXIT
            dt = 1.10 * dt
        ELSE
            counter = counter + 1
            time = time - dt
            dt = 0.8 * dt
            IF( (time + dt) > DeltaTime) dt = DeltaTime - time
            !IF(counter>=5) CALL ExpandBC(solution,Temperature,maxAllowedTemperatureVariation,&
            !                            bcDofs,FlowModel%ActiveToFlowConections,FlowModel%TransPointer)
        END IF

      END DO!****************************************************************** 
      
      CALL PrintMessage( .FALSE., .FALSE., '...End Thermal Solver' )
      
    RETURN
    END SUBROUTINE                                 

    !SUBROUTINE ExpandBC(solution,Temperature,maxAllowedTemperatureVariation,bcDofs,&
    !                    ActiveToFlowConections,TransPointer)
    !
    !    IMPLICIT NONE
    !    REAL(8),INTENT(IN):: solution(:)
    !    REAL(8),INTENT(IN OUT):: Temperature(:)
    !    REAL(8),INTENT(IN):: maxAllowedTemperatureVariation
    !    INTEGER,INTENT(IN OUT):: bcDofs(0:)
    !    INTEGER,INTENT(IN):: ActiveToFlowConections(0:),TransPointer(:,:)    
    !    
    !    INTEGER:: i,numDofs,iel,dof1,dof2
    !    LOGICAL,ALLOCATABLE:: criticalDofs(:)
    !    
    !    numDofs = SIZE(Temperature)
    !    
    !    ALLOCATE( criticalDofs(numDofs) )
    !    criticalDofs = .FALSE.
    !    DO i=1,numDofs
    !        IF( ABS(Temperature(i) - solution(i)) > maxAllowedTemperatureVariation ) criticalDofs(i) = .TRUE.
    !    END DO
    !    
    !    DO i=1,ActiveToFlowConections(0)
    !        iel = ActiveToFlowConections(i)
    !        dof1 = TransPointer(1,iel)
    !        dof2 = TransPointer(2,iel)
    !        IF( criticalDofs(dof1) )THEN
    !            IF( ANY(bcDofs(1:)==dof2) ) Temperature(dof1) = Temperature(dof2)
    !        ELSEIF( criticalDofs(dof2) )THEN
    !            IF( ANY(bcDofs(1:)==dof1) ) Temperature(dof2) = Temperature(dof1)
    !        ELSE
    !        END IF
    !    END DO
    !
    !RETURN
    !END SUBROUTINE

    SUBROUTINE UpdateTemperature( TempModel,GeomechanicalModel,FlowModel,&
                                  DeltaTime )

      USE IO
      USE ProblemDefinition
      USE DiscretizedModel
      USE FluidFlowSolver

      IMPLICIT NONE
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      TYPE(PorousFlowModel),INTENT(IN):: FlowModel
      TYPE(ThermalModel),INTENT(IN OUT):: TempModel
      REAL(8),INTENT(IN):: DeltaTime     
      
      INTEGER:: i,numConections,node1,node2,num(8)
      REAL(8):: qfluid,bw,p,rhor
      REAL(8),ALLOCATABLE:: flow(:),storage(:),Vbphi(:),Vb1_phi(:),rhow(:)
      
      numConections = SIZE(FlowModel%TransArray)
      
      ALLOCATE( storage(GeomechanicalModel%NumNodes),rhow(GeomechanicalModel%NumNodes),&
                Vbphi(GeomechanicalModel%NumNodes),Vb1_phi(GeomechanicalModel%NumNodes) )
      storage = 0      
      Vbphi = 0.0
      Vb1_phi = 0.0
      rhow = 0.0
      
      DO i=1,GeomechanicalModel%NumElements
        num = GeomechanicalModel%ElementNodes(:,i)
        
        Vbphi(num) = Vbphi(num) + FlowModel%ElementVolume(i) / 8.0 * GeomechanicalModel%POR(i)
        
        Vb1_phi(num) = Vb1_phi(num) + FlowModel%ElementVolume(i) / 8.0 * (1.0-GeomechanicalModel%POR(i))
        
      END DO
      
      DO i=1,GeomechanicalModel%NumNodes
        p = GeomechanicalModel%Porepressure(i)
        bw = CalcBW(FlowModel,p)
        rhow(i) = CalcDENSITY(bw,FlowModel)
        rhor = TempModel%NodeRockDensity(i)
        
        storage(i) = Vbphi(i) * rhow(i) * TempModel%FluidCP(i) + Vb1_phi(i) * rhor * TempModel%RockCP(i)
        
      END DO

      ALLOCATE(flow(GeomechanicalModel%NumNodes))
      flow = 0.0
      
      DO i=1,numConections
        node1 = TempModel%TransPointer(1,i)
        node2 = TempModel%TransPointer(2,i)
        
        qfluid = FlowModel%ConectionFlowRates(i)
      
        IF( qfluid<0.0 ) THEN
          qfluid = rhow(node1) * TempModel%FluidCP(node1) * qfluid * TempModel%FluidTemperature(node1)
        ELSE
          qfluid = rhow(node2) * TempModel%FluidCP(node2) * qfluid * TempModel%FluidTemperature(node2)
        END IF
        
        flow(node1) = flow(node1) + qfluid * TempModel%TransPOR(i)
        flow(node2) = flow(node2) - qfluid * TempModel%TransPOR(i)
        
      END DO
      
      DO i=1,GeomechanicalModel%NumNodes
        IF(TempModel%BC_flag(i)==.FALSE.) &
          
          TempModel%FluidTemperature(i) = TempModel%FluidTemperature(i) - flow(i) * DeltaTime / storage(i)
          
      END DO
    
      TempModel%RockTemperature = TempModel%FluidTemperature
      
    RETURN
    END SUBROUTINE
    
                        
    SUBROUTINE BuildThermalModel(GeomechanicalModel,Problem,FlowModel,TempModel)
      
      USE IO
      USE NumericalUtilities
      USE FEMlib
      USE DiscretizedModel
      USE FluidFlowSolver
      USE ProblemDefinition

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      TYPE (PD), INTENT(IN):: Problem
      TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
      TYPE(ThermalModel),INTENT(IN OUT):: TempModel

      INTEGER:: Ncells,Nnodes,iel,i,counter,n

      CHARACTER(LEN=320)::element             !'hexahedron' ou 'tetrahedron'.
      INTEGER,ALLOCATABLE:: num(:)            !element incidence.
      REAL(8):: det
      REAL(8),ALLOCATABLE:: g_coord(:,:),&    !nodal coordenates.
                            coord(:,:)        !element nodal coordenates.
      REAL(8):: elemVOL
      
      REAL(8):: CondTensor(3,3)
      REAL(8),ALLOCATABLE:: elTransm(:,:)
      INTEGER:: dof

      INTEGER:: ndim = 3                      !number of dimensions: 3D -> 3.
      INTEGER:: nod = 8                       !number of nodes per element.
      INTEGER:: nip = 8                       !number of integration points.
      
      REAL(8):: LX,LY,LZ,AX,AY,AZ,K,L,area,T,TX,TY,TZ,z,h
      INTEGER:: pointers(12,2)
      
      REAL(8),ALLOCATABLE:: points(:,:),&     !coordenadas dos pontos de integra��o.
                            weights(:)        !peso dos pontos de integra��o.
      REAL(8),ALLOCATABLE:: jac(:,:),der(:,:),deriv(:,:)
      
      REAL(8),ALLOCATABLE:: density(:)     !node void volume.
      
      INTEGER:: Nz,Nrad,iz,irad
      
      IF(Problem%PerkinsGonzales==.FALSE.) TempModel%Numerical = .TRUE.
      
      TempModel%Tinj = Problem%Tinj
      TempModel%Zinj = Problem%TinjZ
      TempModel%GradTinj = Problem%TinjGrad
      
      Ncells = GeomechanicalModel%NumElements
      Nnodes = GeomechanicalModel%NumNodes
      nod = GeomechanicalModel%NumNodesPerElement
      ndim = GeomechanicalModel%NumDimensions

      ALLOCATE( TempModel%RockTemperature(Nnodes),TempModel%FluidTemperature(Nnodes),&
                TempModel%FluidCP(Nnodes),TempModel%RockCP(Nnodes),num(8),&
                TempModel%BC_flag(Nnodes),TempModel%H(Nnodes) )
      TempModel%RockTemperature = 0.0
      TempModel%FluidTemperature = 0.0
      TempModel%FluidCP = 0.0
      TempModel%RockCP = 0.0
      TempModel%H = 0.0
      TempModel%BC_flag = .FALSE.
            
      Nz = GeomechanicalModel%GridDivision(2)
      Nrad = GeomechanicalModel%GridDivision(1)
      ALLOCATE( TempModel%Thickness(GeomechanicalModel%GridDivision(2)+1) )
      TempModel%Thickness = 0.0
      iel = 0
      DO iz=1, Nz
        DO irad=1, Nrad
          iel = iel + 1
          num = GeomechanicalModel%ElementNodes(:,iel)
          h = GeomechanicalModel%Z(num(2)) - GeomechanicalModel%Z(num(1))
          TempModel%Thickness(iz) = TempModel%Thickness(iz) + h/2.0
          TempModel%Thickness(iz+1) = TempModel%Thickness(iz+1) + h/2.0
        END DO
      END DO
      
      ALLOCATE(g_coord(3,Nnodes))
      g_coord = 0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z
      
      !PROPERTIES--------------------------------------------------------------
      
      DO iel=1,Ncells
        num = GeomechanicalModel%ElementNodes(:,iel)
        TempModel%RockCP(num) = GeomechanicalModel%RockSpecificHeatCapacity(iel)
        TempModel%FluidCP(num) = GeomechanicalModel%FluidSpecificHeatCapacity(iel)
        TempModel%H(num) = GeomechanicalModel%RockFluidHeatTransferCoef(iel)
      END DO
      
      DO i=1,Nnodes
        z = GeomechanicalModel%Z(i) - Problem%TopTemperatureZ
        TempModel%RockTemperature(i) = Problem%TopTemperature - z * Problem%TemperatureGrad
        TempModel%FluidTemperature(i) = TempModel%RockTemperature(i)
      END DO
      GeomechanicalModel%Temperature = TempModel%RockTemperature

      !NODE TO NODE TRANSMISSIBILITIES-----------------------------------------
      !REAL(8),ALLOCATABLE:: TransArray(:)        !transmissibility array (all nodes transmissibilities).
      !INTEGER,ALLOCATABLE:: TransPointer(:,:)    !node pair transmissibilities for each TransArray.
      
      ALLOCATE( TempModel%TransArray(Ncells*12), &
                TempModel%TransPointer(2,Ncells*12),TempModel%TransPOR(Ncells*12) )
      TempModel%TransArray = 0.0
      TempModel%TransPointer = 0.0
      TempModel%TransPOR = 0.0
      
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
        
        CondTensor = 0.0
        CondTensor(1,1) = GeomechanicalModel%ThermalConductivity(iel)
        CondTensor(2,2) = CondTensor(1,1)
        CondTensor(3,3) = CondTensor(1,1)
        
        !LX = ABS( GeomechanicalModel%X(num(1))-GeomechanicalModel%X(num(4)) )
        LX = SQRT( ( GeomechanicalModel%X(num(1))-GeomechanicalModel%X(num(4)) )**2.0 + &
                   ( GeomechanicalModel%Y(num(1))-GeomechanicalModel%Y(num(4)) )**2.0  )
        !LY = ABS( GeomechanicalModel%Y(num(1))-GeomechanicalModel%Y(num(5)) )
        LY = SQRT( (GeomechanicalModel%X(num(1))-GeomechanicalModel%X(num(5)) )**2.0 + &
                   (GeomechanicalModel%Y(num(1))-GeomechanicalModel%Y(num(5)) )**2.0  )
        LZ = ABS( GeomechanicalModel%Z(num(1))-GeomechanicalModel%Z(num(2)) )
        AX = LY*LZ
        AY = LX*LZ
        AZ = LX*LY
        
        ! Nodal trasnmissibilities (1/4)
        TX = CondTensor(1,1) * AX/4.0 / LX
        TY = CondTensor(2,2) * AY/4.0 / LY
        TZ = CondTensor(3,3) * AZ/4.0 / LZ
        
        ! conection 1.
        counter = counter + 1
        TempModel%TransArray(counter) = TX
        TempModel%TransPointer(1,counter) = num( pointers(1,1) )
        TempModel%TransPointer(2,counter) = num( pointers(1,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 2.
        counter = counter + 1
        TempModel%TransArray(counter) = TZ
        TempModel%TransPointer(1,counter) = num( pointers(2,1) )
        TempModel%TransPointer(2,counter) = num( pointers(2,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 3.
        counter = counter + 1
        TempModel%TransArray(counter) = TY
        TempModel%TransPointer(1,counter) = num( pointers(3,1) )
        TempModel%TransPointer(2,counter) = num( pointers(3,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 4.
        counter = counter + 1
        TempModel%TransArray(counter) = TY
        TempModel%TransPointer(1,counter) = num( pointers(4,1) )
        TempModel%TransPointer(2,counter) = num( pointers(4,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 5.
        counter = counter + 1
        TempModel%TransArray(counter) = TZ
        TempModel%TransPointer(1,counter) = num( pointers(5,1) )
        TempModel%TransPointer(2,counter) = num( pointers(5,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 6.
        counter = counter + 1
        TempModel%TransArray(counter) = TY
        TempModel%TransPointer(1,counter) = num( pointers(6,1) )
        TempModel%TransPointer(2,counter) = num( pointers(6,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 7.
        counter = counter + 1
        TempModel%TransArray(counter) = TX
        TempModel%TransPointer(1,counter) = num( pointers(7,1) )
        TempModel%TransPointer(2,counter) = num( pointers(7,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 8.
        counter = counter + 1
        TempModel%TransArray(counter) = TY
        TempModel%TransPointer(1,counter) = num( pointers(8,1) )
        TempModel%TransPointer(2,counter) = num( pointers(8,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 9.
        counter = counter + 1
        TempModel%TransArray(counter) = TX
        TempModel%TransPointer(1,counter) = num( pointers(9,1) )
        TempModel%TransPointer(2,counter) = num( pointers(9,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 10.
        counter = counter + 1
        TempModel%TransArray(counter) = TZ
        TempModel%TransPointer(1,counter) = num( pointers(10,1) )
        TempModel%TransPointer(2,counter) = num( pointers(10,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 11.
        counter = counter + 1
        TempModel%TransArray(counter) = TZ
        TempModel%TransPointer(1,counter) = num( pointers(11,1) )
        TempModel%TransPointer(2,counter) = num( pointers(11,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        ! conection 12.
        counter = counter + 1
        TempModel%TransArray(counter) = TX
        TempModel%TransPointer(1,counter) = num( pointers(12,1) )
        TempModel%TransPointer(2,counter) = num( pointers(12,2) )
        TempModel%TransPOR(counter) = GeomechanicalModel%POR(iel)
        
      END DO
      
      !NODE ROCK VOLUME AND DENSITY--------------------------------------------
      ALLOCATE( TempModel%NodeVolume(Nnodes),TempModel%NodeRockDensity(Nnodes),&
                TempModel%NodePOR(Nnodes),density(Nnodes) )
      TempModel%NodeVolume = 0.0
      TempModel%NodePOR = 0.0
      density = 0.0
      
      ALLOCATE( points(nip,ndim),coord(nod,ndim),  &
                jac(ndim,ndim),weights(nip),der(ndim,nod) )
      element = GeomechanicalModel%element
      CALL sample(element,points,weights)
      
      DO iel=1,Ncells

        num = GeomechanicalModel%ElementNodes(:,iel)
        coord=TRANSPOSE(g_coord(:,num))
        
        elemVOL = 0.0
        DO i=1,nip
          CALL shape_der(der,points,i)
          jac = MATMUL(der,coord)
          CALL determinante(jac,det)
          elemVOL = elemVOL + det * weights(i)
        END DO
        
        TempModel%NodeVolume(num) = TempModel%NodeVolume(num) + elemVol/8.0
        TempModel%NodePOR(num) = TempModel%NodePOR(num) + GeomechanicalModel%POR(iel) * elemVol/8.0
        density(num) = density(num) + elemVol/8.0 * GeomechanicalModel%Density(iel)
        
      END DO
      
      DO i=1,Nnodes
        TempModel%NodeRockDensity(i) = density(i) / TempModel%NodeVolume(i)
        TempModel%NodePOR(i) = TempModel%NodePOR(i) / TempModel%NodeVolume(i)
      END DO
      
      RETURN
    END SUBROUTINE
    
        
    SUBROUTINE TemperatureBC(TempModel,Problem,FlowModel,GeomechanicalModel)
    
      USE ProblemDefinition
      USE IO
      USE FluidFlowSolver
      USE DiscretizedModel
    
      IMPLICIT NONE
      TYPE(ThermalModel),INTENT(IN OUT):: TempModel
      TYPE (PD), INTENT(IN):: Problem
      TYPE(PorousFlowModel),INTENT(IN):: FlowModel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
    
      INTEGER:: i,NRAD,NTHETA,NZ,itheta,iz,irad
      REAL(8):: tol = 0.001, radius
      REAL(8):: x,y,z
      
      NRAD = GeomechanicalModel%GridDivision(1)
      NTHETA = GeomechanicalModel%GridDivision(3)
      NZ = GeomechanicalModel%GridDivision(2)
      
      IF( TRIM(Problem%WellData%MeshType) == 'RAD' )THEN
        !RADIAL MODEL / VERTICAL WELL -------------------------------------------------------------
        
        i = 0
        DO itheta = 1, NTHETA+1
          DO iz = 1, NZ+1
            DO irad=1, NRAD+1
              i = i + 1
              !y = GeomechanicalModel%Y(i)
              !x = SQRT( GeomechanicalModel%X(i)**2.0 + y**2 )
              z = GeomechanicalModel%Z(i)
              !IF( z<=Problem%WellData%InjZ1 .AND. z>=Problem%WellData%InjZ2 .AND. &
              !    x<=Problem%WellData%InjX  )THEN
              IF( z<=Problem%WellData%InjZ1 .AND. z>=Problem%WellData%InjZ2 .AND. &
                  FlowModel%InjNode(i)==.TRUE. )THEN
                TempModel%BC_flag(i) = .TRUE.
                z = GeomechanicalModel%Z(i) - Problem%TinjZ
                TempModel%FluidTemperature(i) = Problem%Tinj - z * Problem%TinjGrad
              !ELSEIF( Problem%WellData%WellBoreWallFlag(i) )THEN
              !  TempModel%BC_flag(i) = .TRUE.
              !ELSE
              END IF
            END DO
          END DO
        END DO
        
      ELSEIF( Problem%WellData%HorWellFlag == .TRUE.)THEN
        
        !BOX MODEL / HORIZONTAL WELL --------------------------------------------------------------
        i = 0
        DO itheta = 1, NTHETA+1
          DO iz = 1, NZ+1
            DO irad=1, NRAD+1
              i = i + 1
              IF( iz==1 )THEN
                TempModel%BC_flag(i) = .TRUE.
              ELSEIF(Problem%WellData%WellBoreWallFlag(i))THEN
                TempModel%BC_flag(i) = .TRUE.
                z = GeomechanicalModel%Z(i) - Problem%TinjZ
                TempModel%FluidTemperature(i) = Problem%Tinj - z * Problem%TinjGrad
              END IF
            END DO
          END DO
        END DO
        
      ELSE
        
        !BOX MODEL / VERTICAL WELL ----------------------------------------------------------------
        i = 0
        DO itheta = 1, NTHETA+1
          DO iz = 1, NZ+1
            DO irad=1, NRAD+1
              i = i + 1
              y = GeomechanicalModel%Y(i)
              x = SQRT( GeomechanicalModel%X(i)**2.0 + y**2.0 )
              z = GeomechanicalModel%Z(i)
              IF( z<=Problem%WellData%InjZ1 .AND. z>=Problem%WellData%InjZ2 .AND. &
                  x<=Problem%WellData%InjX  )THEN
                TempModel%BC_flag(i) = .TRUE.
                z = GeomechanicalModel%Z(i) - Problem%TinjZ
                TempModel%FluidTemperature(i) = Problem%Tinj - z * Problem%TinjGrad
              !ELSEIF( itheta==NTHETA+1 )THEN
              !  TempModel%BC_flag(i) = .TRUE.
              !  TempModel%RockTemperature(i) = 100.0
              !ELSE
              END IF
              !IF( iz==1 ) TempModel%BC_flag(i) = .TRUE.
            END DO
          END DO
        END DO
        
      END IF
      
      TempModel%RockTemperature = TempModel%FluidTemperature
      
    RETURN
    END SUBROUTINE
  

    SUBROUTINE PerkinsGonzalesThermalZone( GeomechanicalModel, TempModel, FlowModel, &
                                           DeltaTime, CumFlowVolume )
      
      USE IO
      USE DiscretizedModel
      USE FluidFlowSolver

      IMPLICIT NONE
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      TYPE(ThermalModel),INTENT(IN OUT):: TempModel  
      TYPE(PorousFlowModel),INTENT(IN OUT):: FlowModel
      REAL(8),INTENT(IN):: DeltaTime
      REAL(8),INTENT(IN OUT):: CumFlowVolume(:)
      
      REAL(8),ALLOCATABLE:: coldRadius(:), injTemp(:)
      
      INTEGER:: Ntheta,Nz,Nrad,itheta,i,iz,irad,num(8)
      REAL(8):: dh,radius,Wi
      REAL(8):: p,bw,rhor,rhof
      REAL(8):: Cf,Cr,por
      REAL(8):: r0,r,x,y
      
      Ntheta = GeomechanicalModel%GridDivision(3)
      Nz = GeomechanicalModel%GridDivision(2)
      Nrad = GeomechanicalModel%GridDivision(1)
      
      ALLOCATE( coldRadius(nz+1), injTemp(nz+1) )
      coldRadius = 0.0
      injTemp = 0.0
      
      i = 0
      DO iz=1, nz+1
        
        injTemp(iz) = TempModel%FluidTemperature(i+1)
        
        DO irad=1, nrad+1
          i = i + 1
          
          IF( FlowModel%InjNode(i)==.TRUE. )THEN
            
            !CumFlowVolume(iz) = CumFlowVolume(iz) + FlowModel%NodeFlows(i) * DeltaTime
            CumFlowVolume(iz) = CumFlowVolume(iz) - FlowModel%NodeFlows(i) * DeltaTime
            
            dh = TempModel%Thickness(iz)
            Wi = CumFlowVolume(iz)
            
            p = GeomechanicalModel%Porepressure(i)
            bw = CalcBW(FlowModel,p)
            rhof = CalcDENSITY(bw,FlowModel)
          
            rhor = TempModel%NodeRockDensity(i)
        
            Cf = TempModel%FluidCP(i)
            Cr = TempModel%RockCP(i)
            por = TempModel%NodePOR(i)
            
            coldRadius(iz) = SQRT( 4.0/3.14/dh * (rhof*Cf*Wi)/(rhor*Cr*(1.0-por)+rhof*Cf*por) )
          
          END IF
          
        END DO
      END DO

      r0 = GeomechanicalModel%X(1)
      i = 0
      DO itheta=1, ntheta+1
        DO iz=1, nz+1
          DO irad=1, nrad+1
            i = i + 1
            
            x = GeomechanicalModel%X(i)
            y = GeomechanicalModel%Y(i)
            r = SQRT( x**2.0 + y**2.0 )
            
            IF( (r-r0) < coldRadius(iz) )THEN
              
              TempModel%FluidTemperature(i) = InjTemp(iz)
              
            END IF
            
          END DO
        END DO
      END DO
      
      TempModel%RockTemperature = TempModel%FluidTemperature
      
    END SUBROUTINE
    

  
  END MODULE
  