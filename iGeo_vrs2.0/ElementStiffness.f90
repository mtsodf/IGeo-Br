
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************
  
  MODULE ElementStiffness
  
  
  CONTAINS
  
  
    SUBROUTINE ElementStiffnessCalculator(ndim,nip,nod,ndof,nst,dee,points,weights,coord,km)

      !USE IFPORT
      USE NumericalUtilities
      USE ConstitutiveGeoModels

      IMPLICIT NONE

      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER::i  !contadores
      INTEGER,INTENT(IN)::ndim,&  !n�mero de dimens�es do problema
          nip,&  !n�mero de pontos de integra��o
          nod,&  !n�mero de n�s por elemento
          ndof,&  !n�mero de graus de liberdade de deslocamento por elemento
          nst  !n�mero de componentes no tensor de tens�es
      REAL(iwp)::det  !determinantee do jacobiano
      REAL(iwp),INTENT(IN):: dee(6,6)

      !VARI�VEIS DIN�MICAS.
      REAL(iwp),ALLOCATABLE:: bee(:,:),&  !matriz de derivadas das fun��es de forma.
          jac(:,:),&  !matriz jacobiana das fun��es de forma.
          fun(:),&  !fun��es de forma.
          der(:,:),&  !derivadas das fun��es de forma dos deslocamentos em coordenadas naturais.
          deriv(:,:)  !derivadas das fun��es de forma dos deslocamentos em coordenadas globais.
      REAL(iwp),INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integra��o.
          weights(:),&  !peso dos pontos de integra��o.
          coord(:,:)  !coordenada dos n�s de um elemento.
      REAL(iwp),INTENT(IN OUT)::km(:,:)  !matriz de rigidez do elemento.

      !ALOCA��O DE MEM�RIA.
      ALLOCATE(bee(nst,ndof),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim))
      fun=0.0_iwp

      !INTEGRA��O NUM�RICA DA MATRIZ DE RIGIDEZ.
      km = 0.0_iwp
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
    
  
  END MODULE