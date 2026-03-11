
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE InterfaceMechanics

    TYPE InterfaceStiffnessVector

      INTEGER:: NumEquations
        
      !DISPLACEMENT RESTRICTED DOFS=================
      INTEGER,ALLOCATABLE:: ListOfRestrictedDOFs(:)
      !=============================================

      !SDOF spring elements=========================
      INTEGER,ALLOCATABLE:: ListOfSupportedDOFs(:)    
      REAL(8),ALLOCATABLE:: Kvalues(:),Area(:)
      !=============================================

    END TYPE InterfaceStiffnessVector
    
    TYPE DisplacementBC
      
      INTEGER,ALLOCATABLE:: DOFs(:)           !DOFs with displacement function.
      REAL(8),ALLOCATABLE:: Displacement(:)   !displacement at dofs listed at DOFs(:)
      REAL(8),ALLOCATABLE:: RampFunction(:)   !ramp function for each time step.
      
      INTEGER,ALLOCATABLE:: ListYnodalRestriction(:)  !list of nodes to be restricted
                                                      !for negative Y displacement.
      
    END TYPE

    CONTAINS


    !Builds a list of restricted DOFs "ListOfRestrictedDOFs" base on a logical array "ListOfRestrictedDOFs".
    SUBROUTINE SupportBoundaryConditions(RestrictedDOFs,SupportedNodeDOFs,GeomechanicalModel,&
                                          ListOfRestrictedDOFs,ListOfSupportedDOFs,Kvalues)
    
      USE NearWellboreModel
      USE DiscretizedModel

      IMPLICIT NONE

      LOGICAL,INTENT(IN):: RestrictedDOFs(:,:),SupportedNodeDOFs(:,:)
      TYPE(Model),INTENT(IN):: GeomechanicalModel
      INTEGER,ALLOCATABLE,INTENT(IN OUT):: ListOfRestrictedDOFs(:),&
                                            ListOfSupportedDOFs(:)
      REAL(8),ALLOCATABLE,INTENT(IN OUT):: Kvalues(:)
        
      INTEGER::node,dof,numRestrictedDOFs,numSupportedDOFs
        
      dof = GeomechanicalModel%NumNodes * 3
      ALLOCATE(ListOfRestrictedDOFs(0:dof))
      ListOfRestrictedDOFs = 0
        
      !DISPLACEMENT BOUNDARY CONDITIONS-----------------------------------------
      numRestrictedDOFs = 0
      DO node=1,GeomechanicalModel%NumNodes
        IF(RestrictedDOFs(1,node)) THEN
          numRestrictedDOFs = numRestrictedDOFs + 1
          dof = GeomechanicalModel%NodeDOFs(1,node)
          ListOfRestrictedDOFs(numRestrictedDOFs) = dof
        END IF
        IF(RestrictedDOFs(2,node)) THEN
          numRestrictedDOFs = numRestrictedDOFs + 1
          dof = GeomechanicalModel%NodeDOFs(2,node)
          ListOfRestrictedDOFs(numRestrictedDOFs) = dof
        END IF
        IF(RestrictedDOFs(3,node))THEN
          numRestrictedDOFs = numRestrictedDOFs + 1
          dof = GeomechanicalModel%NodeDOFs(3,node)
          ListOfRestrictedDOFs(numRestrictedDOFs) = dof
        END IF
      END DO
      ListOfRestrictedDOFs(0) = numRestrictedDOFs
        
        
      !SDOF SPRING-MASS SYSTEM--------------------------------------------------
      ALLOCATE(ListOfSupportedDOFs(0:GeomechanicalModel%NumEquations),&
                Kvalues(GeomechanicalModel%NumEquations))
      Kvalues = 0.0
      numSupportedDOFs = 0
      DO node=1,GeomechanicalModel%NumNodes
        IF(SupportedNodeDOFs(1,node)) THEN
          numSupportedDOFs = numSupportedDOFs + 1
          dof = GeomechanicalModel%NodeDOFs(1,node)
          ListOfSupportedDOFs(numSupportedDOFs) = dof
        END IF
        IF(SupportedNodeDOFs(2,node)) THEN
          numSupportedDOFs = numSupportedDOFs + 1
          dof = GeomechanicalModel%NodeDOFs(2,node)
          ListOfSupportedDOFs(numSupportedDOFs) = dof
        END IF
        IF(SupportedNodeDOFs(3,node))THEN
          numSupportedDOFs = numSupportedDOFs + 1
          dof = GeomechanicalModel%NodeDOFs(3,node)
          ListOfSupportedDOFs(numSupportedDOFs) = dof
        END IF
      END DO
      ListOfSupportedDOFs(0) = numSupportedDOFs
    
    END SUBROUTINE SupportBoundaryConditions
    
    !Builds the support stiffness matrix
    
    
    
    !Casing boundary conditions.
    !SUBROUTINE CasingBC(GeomechanicalModel,ListOfRestrictedDOFs)
    !
    !    USE NearWellboreModel
    !    USE DiscretizedModel
    !
    !    IMPLICIT NONE
    !
    !    TYPE(Model),INTENT(IN):: GeomechanicalModel
    !    INTEGER,INTENT(IN OUT):: ListOfRestrictedDOFs(0:)
    !
    !    INTEGER:: node,dof,numRestrictedDOFs
    !    INTEGER:: iz,it,i
    !    
    !    node = 0
    !    DO iz=1,GeomechanicalModel%GridDivision(2)+1
    !      DO it=1,GeomechanicalModel%GridDivision(1)+1      !<<<<<<<<<<<<<<<
    !        node = node + 1
    !        DO i=1,3
    !          dof = GeomechanicalModel%NodeDOFs(i,node)
    !          IF(.NOT.ANY( ListOfRestrictedDOFs==dof ))THEN
    !            ListOfRestrictedDOFs(0) = ListOfRestrictedDOFs(0) + 1
    !            numRestrictedDOFs = ListOfRestrictedDOFs(0)
    !            ListOfRestrictedDOFs(numRestrictedDOFs) = dof
    !          END IF
    !        END DO
    !      END DO
    !    END DO
    !
    !END SUBROUTINE

END MODULE
