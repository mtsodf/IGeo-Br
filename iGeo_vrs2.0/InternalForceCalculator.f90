  
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE InternalForceCalculator

  CONTAINS

    !global -> integral(Bt,sigma)
    SUBROUTINE CalculateInternalForceVector(GeomechanicalModel)

      USE NumericalUtilities
      USE DiscretizedModel
      USE FEMlib
      USE MainConstants
      USE IO
      USE StressLoad

      IMPLICIT NONE

      !INPUTS / OUTPUTS
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel

      !VARIABLE.
      INTEGER:: iel !contadores.

      INTEGER,ALLOCATABLE:: num(:) ,&  !conectividade de um elemento
                            g(:) !vetor com os graus de liberdade de cada elemento.
      REAL(8),ALLOCATABLE:: points(:,:),&  !coordenadas dos pontos de integraÁ„o
                            weights(:),&  !peso dos pontos de integraÁ„o
                            coord(:,:),&  !coordenada dos nÛs de um elemento
                            g_coord(:,:)  !coordenada dos nÛs da malha.
      REAL(8),ALLOCATABLE:: sigma(:,:),elemNodalForce(:)

      INTEGER:: nels,nn,nip

      !CALL PrintMessage(.FALSE.,.FALSE.,'...ineternal force calculation.')

      GeomechanicalModel%InternalForce = 0.0

      element = GeomechanicalModel%element
      nels = GeomechanicalModel%NumElements
      nn = GeomechanicalModel%NumNodes
      ndim = GeomechanicalModel%NumDimensions
      nst = GeomechanicalModel%NumTensorComponents
      nip = GeomechanicalModel%NumIntegrationPoints
      nod = GeomechanicalModel%NumNodesPerElement
      ndof = GeomechanicalModel%NumDOFUperElement

      !MEMORY ALLOCATION.
      ALLOCATE( num(nod),g(ndof),points(nip,ndim),weights(nip),coord(nod,ndim),g_coord(ndim,nn),&
                sigma(nst,nip),elemNodalForce(ndof))

      !NODAL COORDENATES.
      g_coord=0.0
      g_coord(1,:) = GeomechanicalModel%X
      g_coord(2,:) = GeomechanicalModel%Y
      g_coord(3,:) = GeomechanicalModel%Z

      !INTEGRATION POINTS.
      CALL sample(element,points,weights)

      !COMPUTING NODAL FORCES.
      DO iel=1,nels

        sigma = GeomechanicalModel%EffectiveStressTensor(:,:,iel)

        num = GeomechanicalModel%ElementNodes(:,iel)
        coord = TRANSPOSE(g_coord(:,num))
        CALL ElementStressForce(nip,points,weights,coord,sigma,elemNodalForce)

        g = GeomechanicalModel%ElementDOFs(:,iel) !g_g(:,iel)
        GeomechanicalModel%InternalForce(g) = GeomechanicalModel%InternalForce(g) + elemNodalForce

      END DO

    END SUBROUTINE CalculateInternalForceVector

  END MODULE
