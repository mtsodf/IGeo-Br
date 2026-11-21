  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia
  !*****************************************************************************


  MODULE NumericalLib


  CONTAINS

    !Rotina para cálculo de determinante de matrizes pequenas.
    !Smith & Griffiths 2004.
    SUBROUTINE MatrixDeterminant(jac,det)

      ! This function returns the determinant of a 1x1, 2x2 or 3x3
      ! Jacobian matrix.

      IMPLICIT NONE
      REAL(8),INTENT(IN)::jac(:,:)
      REAL(8),INTENT(IN OUT)::det
      INTEGER::it
      it=UBOUND(jac,1)
      det=0.0
      SELECT CASE(it)
      CASE(1)
        det=1.0
      CASE(2)
        det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
      CASE(3)
        det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
        det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
        det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
        CASE DEFAULT
        WRITE(*,*)' wrong dimension for Jacobian matrix'
      END SELECT

      RETURN
    END SUBROUTINE MatrixDeterminant

    !Rotina para inverter pequenas matrizes.
    !Smith & Griffiths 2004.
    SUBROUTINE InvertMatrix(matrix)

      ! This subroutine inverts a small square matrix onto itself.

      IMPLICIT NONE
      REAL(8),INTENT(IN OUT)::matrix(:,:)
      REAL(8)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
      INTEGER::ndim,i,k
      ndim=UBOUND(matrix,1)
      IF(ndim==2)THEN
        det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
        j11=matrix(1,1)
        matrix(1,1)=matrix(2,2)
        matrix(2,2)=j11
        matrix(1,2)=-matrix(1,2)
        matrix(2,1)=-matrix(2,1)
        matrix=matrix/det
      ELSE IF(ndim==3)THEN
        det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
        det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
        det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
        j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
        j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
        j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
        j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
        j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
        j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
        j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
        j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
        j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
        matrix(1,1)=j11
        matrix(1,2)=j12
        matrix(1,3)=j13
        matrix(2,1)=j21
        matrix(2,2)=j22
        matrix(2,3)=j23
        matrix(3,1)=j31
        matrix(3,2)=j32
        matrix(3,3)=j33
        matrix=matrix/det
      ELSE
        DO k=1,ndim
          con=matrix(k,k)
          matrix(k,k)=1.0
          matrix(k,:)=matrix(k,:)/con
          DO i=1,ndim
            IF(i/=k)THEN
              con=matrix(i,k)
              matrix(i,k)=0.0
              matrix(i,:)=matrix(i,:)-matrix(k,:)*con
            END IF
          END DO
        END DO
      END IF

      RETURN
    END SUBROUTINE InvertMatrix

    !Rotina para cálculo de produtor escalar.
    FUNCTION DotProduct(VectorA,VectorB)

      IMPLICIT NONE
      REAL(8),INTENT(IN):: VectorA(3), VectorB(3)
      REAL(8):: DotProduct

      DotProduct = vectorA(1)*vectorB(1) + vectorA(2)*vectorB(2) + vectorA(3)*vectorB(3)

      END FUNCTION DotProduct

      !Rotina para cálculo de distância entre dois pontos.
      FUNCTION Distance(Coords1,Coords2)

      IMPLICIT NONE
      REAL(8):: Distance
      REAL(8),INTENT(IN):: Coords1(3),Coords2(3)

      Distance = SQRT(  ( Coords1(1)-Coords2(1) )**2.0 + &
                        ( Coords1(2)-Coords2(2) )**2.0 + &
                        ( Coords1(3)-Coords2(3) )**2.0 )

      END FUNCTION Distance

      !Rotina para cálculo de vetor entre 2 pontos.
      FUNCTION Vector2Pts(Coords1,Coords2)

      IMPLICIT NONE
      REAL(8):: Vector2Pts(3)
      REAL(8),INTENT(IN):: Coords1(3),Coords2(3)

      Vector2Pts(1) = Coords2(1) - Coords1(1)
      Vector2Pts(2) = Coords2(2) - Coords1(2)
      Vector2Pts(3) = Coords2(3) - Coords1(3)

    END FUNCTION Vector2Pts

    !Rotina para cálculo de vetor entre 2 pontos.
    FUNCTION Vector3Pts(Coords1,Coords2,Coords3)

      IMPLICIT NONE
      REAL(8):: Vector3Pts(3)
      REAL(8),INTENT(IN):: Coords1(3),Coords2(3),Coords3(3)

      REAL(8):: vectorA(3), VectorB(3)

      vectorA = Vector2Pts(Coords1,Coords2)
      vectorB = Vector2Pts(Coords1,Coords3)

      Vector3Pts(1) = vectorA(2)*vectorB(3) - vectorA(3)*vectorB(2)
      Vector3Pts(2) = vectorA(3)*vectorB(1) - vectorA(1)*vectorB(3)
      Vector3Pts(3) = vectorA(1)*vectorB(2) - vectorA(2)*vectorB(1)

    END FUNCTION Vector3Pts

    SUBROUTINE Triangle( coord, Area, normal )
      !Subrotina para cálculo de área e vetor normal a um triângulo a partir das
      !coordenadas nodais.
      !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

      !INPUTS
      INTEGER:: i
      REAL(iwp),INTENT(IN)::  coord(3,3)  !3 node coordenates (node,coordenate)
      REAL(iwp),INTENT(OUT):: Area,&      !triagle area
                              normal(3)   !normal vector

      REAL(iwp)::a(3), b(3),absProdutoVetorial, produtoVetorial(3)
      REAL(iwp):: one = 1.0_iwp, two = 2.0_iwp

      !Formação dos vetores
      DO i=1,3
        a(i) = coord(2,i) - coord(1,i)
        b(i) = coord(3,i) - coord(1,i)
      END DO

      !Produto Vetorial.
      produtoVetorial(1) = a(2)*b(3) - a(3)*b(2)
      produtoVetorial(2) = a(3)*b(1) - a(1)*b(3)
      produtoVetorial(3) = a(1)*b(2) - a(2)*b(1)

      !Vetor normal.
      absProdutoVetorial = SQRT( produtoVetorial(1)**2 + produtoVetorial(2)**2 + produtoVetorial(3)**2 )
      normal = one / absProdutoVetorial * produtoVetorial

      !Área triângulo = metada da área do paralelogramo formado pelos 2 vetores.
      Area = absProdutoVetorial / two

      RETURN
    END SUBROUTINE Triangle

    SUBROUTINE AngleVectos( vectorA, VectorB, theta )
  
      IMPLICIT NONE
      REAL(8),INTENT(IN)::  vectorA(3), VectorB(3)
      REAL(8),INTENT(IN OUT)::  theta
    
      REAL(8)::a,b
    
      a = SQRT( vectorA(1)**2.0 + vectorA(2)**2.0 + vectorA(3)**2.0 )
      b = SQRT( vectorB(1)**2.0 + vectorB(2)**2.0 + vectorB(3)**2.0 )
    
      theta = ( vectorA(1)*vectorB(1) + vectorA(2)*vectorB(2) + vectorA(3)*vectorB(3) ) / (a*b)
      theta = ACOS(theta)
      
      RETURN
    END SUBROUTINE

  END MODULE NumericalLib
