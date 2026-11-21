!*****************************************************************************
! Substitutos das rotinas do Intel MKL usadas pelo iGeo (somente para gfortran).
!
! - mkl_dcsrgemv: produto matriz esparsa (CSR, base 1) por vetor, y = A*x.
!   Implementação equivalente à do MKL para transa = 'N', paralelizada
!   com OpenMP.
! - ddoti: produto escalar de um vetor esparso por um vetor cheio
!   (Sparse BLAS do MKL).
! - vdmul: multiplicação elemento a elemento, r = a*b (VML do MKL).
! - pardiso: o solver direto do MKL não tem substituto aqui. A única rotina
!   que o chama (LinearSolverPAR, em SolidMechanicsLib.f90) não é usada
!   pelo programa; se for chamada, a execução é interrompida com um erro.
!
! As rotinas BLAS/LAPACK padrão (dswap, dgeev) vêm da biblioteca ligada
! pelo Makefile (Accelerate no macOS, BLAS/LAPACK de referência ou
! OpenBLAS nos demais sistemas).
!
! As rotinas são externas (sem módulo), como no MKL, para que o código
! original funcione sem alterações.
!*****************************************************************************

REAL(8) FUNCTION ddoti(nz, x, indx, y)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: nz
  REAL(8), INTENT(IN):: x(*), y(*)
  INTEGER, INTENT(IN):: indx(*)

  INTEGER:: i

  ddoti = 0.0D0
  DO i = 1, nz
    ddoti = ddoti + x(i) * y(indx(i))
  END DO

END FUNCTION ddoti


SUBROUTINE vdmul(n, a, b, r)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: n
  REAL(8), INTENT(IN):: a(*), b(*)
  REAL(8), INTENT(OUT):: r(*)

  INTEGER:: i

  !$OMP PARALLEL DO SCHEDULE(STATIC)
  DO i = 1, n
    r(i) = a(i) * b(i)
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE vdmul


SUBROUTINE mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN):: transa
  INTEGER, INTENT(IN):: m
  REAL(8), INTENT(IN):: a(*), x(*)
  INTEGER, INTENT(IN):: ia(*), ja(*)
  REAL(8), INTENT(OUT):: y(*)

  INTEGER:: i, k
  REAL(8):: soma

  IF (transa /= 'N' .AND. transa /= 'n') THEN
    WRITE(*,*) 'mkl_dcsrgemv (compat): apenas transa = ''N'' é suportado.'
    ERROR STOP 1
  END IF

  !$OMP PARALLEL DO PRIVATE(k, soma) SCHEDULE(STATIC)
  DO i = 1, m
    soma = 0.0D0
    DO k = ia(i), ia(i+1) - 1
      soma = soma + a(k) * x(ja(k))
    END DO
    y(i) = soma
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE mkl_dcsrgemv


SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                   perm, nrhs, iparm, msglvl, b, x, error)
  IMPLICIT NONE
  INTEGER:: pt(*), maxfct, mnum, mtype, phase, n, ia(*), ja(*), &
            perm(*), nrhs, iparm(*), msglvl, error
  REAL(8):: a(*), b(*), x(*)

  WRITE(*,*) 'pardiso (compat): o solver PARDISO do Intel MKL não está'
  WRITE(*,*) 'disponível no build com gfortran.'
  error = -1
  ERROR STOP 1

END SUBROUTINE pardiso
