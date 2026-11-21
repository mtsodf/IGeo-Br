  
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE PCG_Flow

  CONTAINS


    !Linear Solver for fluid flow (no frac version).
    !LINE ELEMENT (2 NODES) BASED VERSION
    ![LHS]*{x}={RHS}
    SUBROUTINE PCG_FlowLineElementBased_noFrac(cg_limit,cg_tol,Trans,TransDOFs,Trans_vector,RestrictedDOFs,&
                                                RHS,initialGuess)
    
      !---------------------------------------------------------------------------------------------------------------
      !solver de sistema linear de equações com o método dos GRADIENTES CONJUGADOS pré-condionado diagonalmente (PCG).
      !rotina paralelizada para o OpenMP vrs. 3.0 .
      !uso de rotinas seriais do Intel Fortran Math Kernel Library - MKL
      !armazenamento de matrizes formato CSR Intel
      !
      !Referência: Ian M. Smith, D. V. Griffiths.
      !Professor Smith is currently the Head of the Geotechnics Division at the Manchester Centre for Civil and Construction
      !Engineering, UMIST
      !Professor Griffiths is Professor of Civil Engineering, Division of Engineering at the Colorado School of Mines
      !---------------------------------------------------------------------------------------------------------------
    
      USE OMP_LIB  !OpenMP(R) lib.
      USE IO
    
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    
      INTEGER,INTENT(IN)::cg_limit    !maximum number of iterations.
      REAL(iwp),INTENT(IN)::cg_tol    !tolerance.
      REAL(iwp),INTENT(IN OUT)::Trans(:),Trans_vector(:)
      INTEGER,INTENT(IN)::TransDOFs(:,:)
      LOGICAL,INTENT(IN):: RestrictedDOFs(:)  !flag for prescribed pressure b.c..
      REAL(iwp),INTENT(IN OUT)::RHS(:)  !vetor de carregamento estrutural (IN) e resultados (OUT)
      REAL(iwp),INTENT(IN)::initialGuess(:)  !vetor com LHSa solução inicial de entrada no método
    
      INTEGER,ALLOCATABLE:: activeDOFs(:),vectorActiveDOFs(:)
      REAL(iwp),ALLOCATABLE:: diag_precon(:),& !pre-conditioner.
                              store(:)
                                
      INTEGER::i,iel,k,counter,cg_iters,nPorousDOFs
      INTEGER:: iz,irad
      REAL(iwp):: alpha,beta,up,&
                  convergence,&
                  temp1,temp2
      INTEGER:: i1,i2
      INTEGER:: neq,nels,numActiveElements
        
      INTEGER:: nod = 2   !number of nodes per element.
      INTEGER,ALLOCATABLE:: g(:),node(:)
      REAL(iwp):: elTransm
    
      REAL(iwp),ALLOCATABLE:: d(:),p(:),u(:),x(:),xnew(:),auxRealDb(:)
        
      REAL(iwp):: penalty=1.0e30
    
      CHARACTER(LEN=320):: text
        
      LOGICAL:: cg_converged
        
      INTEGER:: nFracElements
    
      !CALL PrintMessage(.FALSE.,.FALSE.,'...fluid flow linear solver.')
        
      neq = SIZE(RHS) !number of DOFs.
      numActiveElements = SIZE(Trans)
    
      !ACTIVE DOFS.
      !Inactive Elements.
      ALLOCATE( vectorActiveDOFs(neq), g(nod) )
      vectorActiveDOFs = 0
      DO iel=1,numActiveElements
        g = TransDOFs(:,iel)
        vectorActiveDOFs(g) = 1
      END DO
      !Pressure B.C.
      DO i=1,SIZE(RHS)
        IF( RestrictedDOFs(i) ) vectorActiveDOFs(i) = 0
      END DO
      counter = SUM(vectorActiveDOFs)
        
      ALLOCATE( activeDOFs(counter) )
      ALLOCATE( node(neq-counter), store(neq-counter) )
      node = 0
      counter = 0
      k = 0
      DO i=1,SIZE(RHS)
        IF(vectorActiveDOFs(i) > 0  )THEN
          counter = counter + 1
          activeDOFs(counter) = i
        ELSE
          k = k + 1
          node(k) = i
        END IF
      END DO
    
      !MEMORY ALLOCATION.
      ALLOCATE(d(0:neq),p(0:neq),u(0:neq),x(0:neq),xnew(0:neq),auxRealDb(0:neq))
      d=0.0
      p=0.0
      u=0.0
      x=0.0
      xnew=0.0
      auxRealDb=0.0
    
      !--------------------------------------------------------------------------------------
      !PRE-CONDITIONER.
    
      nels = SIZE(Trans)
        
      ALLOCATE(diag_precon(0:neq))
      diag_precon = 0.0
    
      !Porous.
      DO iel=1,numActiveElements
        g = TransDOFs(:,iel)
        DO k=1,nod
          diag_precon(g(k)) = diag_precon(g(k)) - Trans(iel)
        END DO
      END DO
      diag_precon(1:neq) = diag_precon(1:neq) + Trans_vector
        
      !Boundary conditions.
      DO i=1,neq
        IF(RestrictedDOFs(i))THEN
          diag_precon(i) = diag_precon(i) + penalty
          RHS(i) = diag_precon(i) * RHS(i)
        END IF
      END DO
      store = diag_precon(node)
        
      !Inverting.
      DO i=1,neq
        IF(diag_precon(i).NE.0.0)THEN
          diag_precon(i) = 1.0 / diag_precon(i)
        ELSE
          diag_precon(i) = 0.0
        END IF
      END DO
    
      !--------------------------------------------------------------------------------------
      d(1:neq) = diag_precon(1:neq) * RHS
      p = d
      x = 0.0
    
      !LOOPING ITERATIVO++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      cg_iters=0
      cg_converged=.TRUE.
    
      conjug_grad: DO
        cg_iters = cg_iters+1
    
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !u = k * p
        u = 0.0
        !CALL mkl_dcsrgemv('N', neq, LHSa, LHSia, LHSja, p, u)  !MKL - rotina blas2 para multiplicao matriz esparsa por vetor.
        !Porous.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i,iel,g,elTransm) REDUCTION(+: u)
        DO iel=1,numActiveElements
          g = TransDOFs(:,iel)
          elTransm = Trans(iel)
          u(g(1)) = u(g(1)) - elTransm * p(g(1)) + elTransm * p(g(2))
          u(g(2)) = u(g(2)) - elTransm * p(g(2)) + elTransm * p(g(1))
        END DO
!$OMP END PARALLEL DO
!$OMP WORKSHARE
        u(1:neq) = u(1:neq) + Trans_vector * p(1:neq)
!$OMP END WORKSHARE
        !BC.
!$OMP WORKSHARE
        u(node) = p(node) * store
!$OMP END WORKSHARE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
        up=DOT_PRODUCT(RHS,d(1:neq))
          
        alpha=up/DOT_PRODUCT(p,u)
!$OMP WORKSHARE
        xnew=x+p*alpha
!$OMP END WORKSHARE
!$OMP WORKSHARE
        RHS=RHS-u(1:neq)*alpha
!$OMP END WORKSHARE

        d(1:neq)=diag_precon(1:neq)*RHS
        beta=DOT_PRODUCT(RHS,d(1:neq))/up
        
!$OMP WORKSHARE
        p=d+p*beta
!$OMP END WORKSHARE
        
        temp1=MAXVAL(ABS(xnew(activeDOFs)-x(activeDOFs))) !MAXVAL(ABS(xnew-x))
        temp2=MAXVAL(ABS(xnew))
        convergence=temp1/temp2
        cg_converged=(convergence<=cg_tol)
        x=xnew
    
        IF(cg_converged.OR.cg_iters==cg_limit)EXIT
    
      END DO conjug_grad
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !GUARDANDO RESULTADO
      RHS = xnew(1:neq)
    
      !FINALIZANDO E ESCREVENDO RESULTADOS DA TELEMETRIA DO SOLVER
      IF (cg_iters==cg_limit)THEN
        WRITE(*,*)" "
        WRITE(*,*)"...PCG flow did not converge... finalizing"
        READ(*,*)
        STOP
      END IF
    
      RETURN
    END SUBROUTINE

                                        
    
                                        
                                        
                                        
                                        
    !Linear Solver for fluid flow.
    !ELEMENT BASED VERSION.
    ![LHS]*{x}={RHS}
    SUBROUTINE PCG_FlowLinearSolver(cg_limit,cg_tol,Trans,Trans_vector,TransDOFs,ActiveToFlowElements,RestrictedDOFs,&
                                    FracTrans,FracTransDOFs,ListOfActiveFracElements,RHS,initialGuess)

        !---------------------------------------------------------------------------------------------------------------
        !solver de sistema linear de equações com o método dos GRADIENTES CONJUGADOS pré-condionado diagonalmente (PCG).
        !rotina paralelizada para o OpenMP vrs. 3.0 .
        !uso de rotinas seriais do Intel Fortran Math Kernel Library - MKL
        !armazenamento de matrizes formato CSR Intel
        !
        !Referência: Ian M. Smith, D. V. Griffiths.
        !Professor Smith is currently the Head of the Geotechnics Division at the Manchester Centre for Civil and Construction
        !Engineering, UMIST
        !Professor Griffiths is Professor of Civil Engineering, Division of Engineering at the Colorado School of Mines
        !---------------------------------------------------------------------------------------------------------------

        USE OMP_LIB  !OpenMP(R) lib.
        USE IO

        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

        INTEGER,INTENT(IN)::cg_limit    !maximum number of iterations.
        REAL(iwp),INTENT(IN)::cg_tol    !tolerance.
        REAL(iwp),INTENT(IN OUT)::Trans(:,:,:),Trans_vector(:)
        INTEGER,INTENT(IN)::TransDOFs(:,:)
        INTEGER,INTENT(IN):: ActiveToFlowElements(:)
        LOGICAL,INTENT(IN):: RestrictedDOFs(:)  !flag for prescribed pressure b.c..
        REAL(iwp),INTENT(IN):: FracTrans(:,:,:)
        INTEGER,INTENT(IN):: FracTransDOFs(:,:),ListOfActiveFracElements(:)
        REAL(iwp),INTENT(IN OUT)::RHS(:)  !vetor de carregamento estrutural (IN) e resultados (OUT)
        REAL(iwp),INTENT(IN)::initialGuess(:)  !vetor com LHSa solução inicial de entrada no método

        INTEGER,ALLOCATABLE:: activeDOFs(:),vectorActiveDOFs(:)
        REAL(iwp),ALLOCATABLE:: diag_precon(:),& !pre-conditioner.
                                store(:)
                                
        INTEGER::i,iel,k,counter,cg_iters,nPorousDOFs
        INTEGER:: iz,irad
        REAL(iwp):: alpha,beta,up,&
                    convergence,&
                    temp1,temp2
        INTEGER:: i1,i2
        INTEGER:: neq,nels,numActiveElements
        
        INTEGER:: nod = 8   !number of nodes per element.
        INTEGER,ALLOCATABLE:: g(:),node(:)
        REAL(iwp),ALLOCATABLE:: elTransm(:,:)

        REAL(iwp),ALLOCATABLE:: d(:),p(:),u(:),x(:),xnew(:),auxRealDb(:)
        
        REAL(iwp):: penalty=1.0e30

        CHARACTER(LEN=320):: text
        
        LOGICAL:: cg_converged
        
        INTEGER:: nFracElements

        !CALL PrintMessage(.FALSE.,.FALSE.,'...fluid flow linear solver.')
        
        neq = SIZE(RHS) !number of DOFs.
        numActiveElements = SIZE(ActiveToFlowElements)
        nFracElements = SIZE(ListOfActiveFracElements)

        !ACTIVE DOFS.
        !Inactive Elements.
        ALLOCATE( vectorActiveDOFs(neq), g(nod) )
        vectorActiveDOFs = 0
        DO i=1,numActiveElements
          iel = ActiveToFlowElements(i)
          g = TransDOFs(:,iel)
          vectorActiveDOFs(g) = 1
        END DO
        !Pressure B.C.
        !!DO i=1,SIZE(RHS)
        !!  IF( RestrictedDOFs(i) ) vectorActiveDOFs(i) = 0
        !!END DO
        counter = SUM(vectorActiveDOFs)
        
        ALLOCATE( activeDOFs(counter) )
        ALLOCATE( node(neq-counter), store(neq-counter) )
        node = 0
        counter = 0
        k = 0
        DO i=1,SIZE(RHS)
          IF(vectorActiveDOFs(i) > 0  )THEN
            counter = counter + 1
            activeDOFs(counter) = i
          ELSE
            k = k + 1
            node(k) = i
          END IF
        END DO

        !MEMORY ALLOCATION.
        ALLOCATE(d(0:neq),p(0:neq),u(0:neq),x(0:neq),xnew(0:neq),auxRealDb(0:neq))
        d=0.0
        p=0.0
        u=0.0
        x=0.0
        xnew=0.0
        auxRealDb=0.0
        
        ALLOCATE(elTransm(nod,nod))

        !--------------------------------------------------------------------------------------
        !PRE-CONDITIONER.

        nels = SIZE(Trans(1,1,:))
        
        ALLOCATE(diag_precon(0:neq))
        diag_precon = 0.0

        !Porous.
        DO i=1,numActiveElements
          iel = ActiveToFlowElements(i)
          g = TransDOFs(:,iel)
          DO k=1,nod
            diag_precon(g(k)) = diag_precon(g(k)) + Trans(k,k,iel)
          END DO
        END DO
        diag_precon(1:neq) = diag_precon(1:neq) + Trans_vector
        !Fracture.
        DO i=1,nFracElements
          iel = ListOfActiveFracElements(i)
          g(1:4) = FracTransDOFs(:,iel)
          elTransm(1:4,1:4) = FracTrans(:,:,iel)
          DO k=1,4
            diag_precon(g(k)) = diag_precon(g(k)) + elTransm(k,k)
          END DO
        END DO
        
        !Boundary conditions.
        DO i=1,neq
          IF(RestrictedDOFs(i))THEN
            diag_precon(i) = diag_precon(i) + penalty
            RHS(i) = diag_precon(i) * RHS(i)
          END IF
        END DO
        store = diag_precon(node)
        
        !Inverting.
        DO i=1,neq
          IF(diag_precon(i).NE.0.0)THEN
            diag_precon(i) = 1.0 / diag_precon(i)
          ELSE
            diag_precon(i) = 0.0
          END IF
        END DO

        !--------------------------------------------------------------------------------------
        d(1:neq) = diag_precon(1:neq) * RHS
        p = d
        x = 0.0

        !LOOPING ITERATIVO++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        cg_iters=0
        cg_converged=.TRUE.

        conjug_grad: DO
          cg_iters = cg_iters+1

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !u = k * p
          u = 0.0
          !CALL mkl_dcsrgemv('N', neq, LHSa, LHSia, LHSja, p, u)  !MKL - rotina blas2 para multiplicao matriz esparsa por vetor.
          !Porous.
          DO i=1,numActiveElements
            iel = ActiveToFlowElements(i)
            g = TransDOFs(:,iel)
            elTransm = Trans(:,:,iel)
            u(g) = u(g) + MATMUL( elTransm, p(g) )
          END DO
          u(1:neq) = u(1:neq) + Trans_vector * p(1:neq)
          !Fracture.
          DO i=1,nFracElements
            iel = ListOfActiveFracElements(i)
            g(1:4) = FracTransDOFs(:,iel)
            elTransm(1:4,1:4) = FracTrans(:,:,iel)
            u(g(1:4) ) = u(g(1:4) ) + MATMUL( elTransm(1:4,1:4), p(g(1:4) ) )
          END DO
          !BC.
          u(node) = p(node) * store
            
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          up=DOT_PRODUCT(RHS,d(1:neq))
          alpha=up/DOT_PRODUCT(p,u)
          xnew=x+p*alpha
          RHS=RHS-u(1:neq)*alpha
          d(1:neq)=diag_precon(1:neq)*RHS
          beta=DOT_PRODUCT(RHS,d(1:neq))/up
          p=d+p*beta

          temp1=MAXVAL(ABS(xnew-x))!MAXVAL(ABS(xnew(activeDOFs)-x(activeDOFs)))
          temp2=MAXVAL(ABS(xnew))
          convergence=temp1/temp2
          cg_converged=(convergence<=cg_tol)
          x=xnew

          IF(cg_converged.OR.cg_iters==cg_limit)EXIT

        END DO conjug_grad
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !GUARDANDO RESULTADO
        RHS = xnew(1:neq)

        !FINALIZANDO E ESCREVENDO RESULTADOS DA TELEMETRIA DO SOLVER
        IF (cg_iters==cg_limit)THEN
          WRITE(*,*)" "
          WRITE(*,*)"...linear solver did not converge... finalizing"
          READ(*,*)
          STOP
        END IF

        !!TELEMETRIA DE ITERACOES
        !WRITE(text,'(A40,I15)')'number of iterations: ',cg_iters
        !CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
        !convergence = 1.0_iwp - convergence
        !WRITE(text,'(A40,F15.8)')'convergence: ',convergence
        !CALL PrintMessage(.FALSE.,.FALSE.,TRIM(text))
        !temp1 = MAXVAL(ABS(xnew))
        !WRITE(text,'(A40,E15.4,A)')'maximum pressure variation: ',temp1/1000.0,' kPa'
        !CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )

        RETURN
    END SUBROUTINE

END MODULE

