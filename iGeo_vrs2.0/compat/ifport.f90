!*****************************************************************************
! Compatibilidade com o módulo IFPORT da Intel (somente para gfortran).
!
! Rotinas do IFPORT usadas pelo iGeo:
! - SLEEP, IARGC, GETARG: o gfortran já as fornece como extensões GNU.
! - SORTQQ(LOC(v), n, SRT$INTEGER4): implementada aqui (apenas para
!   INTEGER(4), ordem crescente, que é o único uso no código).
!
! Requer -fdollar-ok, por causa do nome SRT$INTEGER4.
! Não use este arquivo com ifort/ifx: eles têm o IFPORT verdadeiro.
!*****************************************************************************
MODULE IFPORT
  USE, INTRINSIC:: ISO_C_BINDING, ONLY: C_PTR, C_INTPTR_T, C_F_POINTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: SORTQQ, SRT$INTEGER4

  ! Identificador do tipo INTEGER(4) para SORTQQ (o valor só importa aqui).
  INTEGER, PARAMETER:: SRT$INTEGER4 = 4

CONTAINS

  ! Ordena em ordem crescente os n primeiros elementos INTEGER(4) que
  ! começam no endereço adrarray (obtido com LOC).
  SUBROUTINE SORTQQ(adrarray, count, size)
    INTEGER(C_INTPTR_T), INTENT(IN):: adrarray
    INTEGER, INTENT(IN):: count, size

    TYPE(C_PTR):: ptr
    INTEGER(4), POINTER:: v(:)

    IF (size /= SRT$INTEGER4) THEN
      WRITE(*,*) 'SORTQQ (compat): apenas SRT$INTEGER4 é suportado.'
      ERROR STOP 1
    END IF
    IF (count < 2) RETURN

    ptr = TRANSFER(adrarray, ptr)
    CALL C_F_POINTER(ptr, v, [count])
    CALL HeapSortInt4(v)
  END SUBROUTINE SORTQQ

  SUBROUTINE HeapSortInt4(v)
    INTEGER(4), INTENT(IN OUT):: v(:)
    INTEGER:: n, i
    INTEGER(4):: tmp

    n = SIZE(v)
    DO i = n/2, 1, -1
      CALL SiftDown(v, i, n)
    END DO
    DO i = n, 2, -1
      tmp = v(1); v(1) = v(i); v(i) = tmp
      CALL SiftDown(v, 1, i-1)
    END DO
  END SUBROUTINE HeapSortInt4

  SUBROUTINE SiftDown(v, first, last)
    INTEGER(4), INTENT(IN OUT):: v(:)
    INTEGER, INTENT(IN):: first, last
    INTEGER:: root, child
    INTEGER(4):: tmp

    root = first
    DO WHILE (2*root <= last)
      child = 2*root
      IF (child < last) THEN
        IF (v(child) < v(child+1)) child = child + 1
      END IF
      IF (v(root) >= v(child)) RETURN
      tmp = v(root); v(root) = v(child); v(child) = tmp
      root = child
    END DO
  END SUBROUTINE SiftDown

END MODULE IFPORT
