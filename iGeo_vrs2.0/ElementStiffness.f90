
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
      INTEGER,INTENT(IN)::ndim,&  !número de dimensões do problema
          nip,&  !número de pontos de integração
          nod,&  !número de nós por elemento
          ndof,&  !número de graus de liberdade de deslocamento por elemento
          nst  !número de componentes no tensor de tensões
      REAL(iwp)::det  !determinantee do jacobiano
      REAL(iwp),INTENT(IN):: dee(6,6)

      !VARIÁVEIS DINÂMICAS.
      REAL(iwp),ALLOCATABLE:: bee(:,:),&  !matriz de derivadas das funções de forma.
          jac(:,:),&  !matriz jacobiana das funções de forma.
          fun(:),&  !funções de forma.
          der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais.
          deriv(:,:)  !derivadas das funções de forma dos deslocamentos em coordenadas globais.
      REAL(iwp),INTENT(IN)::points(:,:),&  !coordenadas dos pontos de integração.
          weights(:),&  !peso dos pontos de integração.
          coord(:,:)  !coordenada dos nós de um elemento.
      REAL(iwp),INTENT(IN OUT)::km(:,:)  !matriz de rigidez do elemento.

      !ALOCAÇÃO DE MEMÓRIA.
      ALLOCATE(bee(nst,ndof),der(ndim,nod),deriv(ndim,nod),fun(nod),jac(ndim,ndim))
      fun=0.0_iwp

      !INTEGRAÇÃO NUMÉRICA DA MATRIZ DE RIGIDEZ.
      km = 0.0_iwp
      DO i = 1, nip

        CALL shape_fun(fun,points,i)  !funções de forma no ponto de integração i
        CALL shape_der(der,points,i)  !derivada das funções de forma
        jac=MATMUL(der,coord)  !jacobiano do elemento no ponto de integração i
        CALL determinante(jac,det)  !determinante do jacobiano
        CALL invert(jac)  !inverso do jacobiano
        deriv=MATMUL(jac,der)  !derivadas das funções de forma em coordenadas reais
        CALL beemat(bee,deriv)  !Matriz de derivadas B

        km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)  !km=km+BT*D*B*det(J)*peso

      END DO

      RETURN
  END SUBROUTINE
    
  
  END MODULE