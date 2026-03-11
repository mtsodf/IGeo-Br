
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE StructuralModel

    TYPE StructuralSolver

        !Solver linear.
        REAL(8):: LS_tolerance
        INTEGER:: LS_maxIters

        !Solver nao-linear.
        REAL(8):: NLS_tolerance
        INTEGER:: NLS_maxIters
                
        CHARACTER(LEN=18):: StructuralSolver              !'Elastoplastic'
                                                          !'DynamicFEM'

    END TYPE

CONTAINS

    !DEFINING STRUCTURAL SOLVER
    SUBROUTINE ReadSolverParameters(Path, BaseFile, FreeChanelNum, StrucSolver )

        USE IO

        !VARIABLES.
        IMPLICIT NONE
        CHARACTER(LEN=320),INTENT(IN):: BaseFile, &   !file with basic information about the case.
                                        Path          !work directory where the case file is located and where temporary files will
                                                      !be created.
        INTEGER,INTENT(IN)::  FreeChanelNum
        TYPE (StructuralSolver), INTENT(IN OUT):: StrucSolver

        LOGICAL:: success

        !OPEN( FreeChanelNum, FILE = TRIM(TRIM(Path)//'/'//TRIM(BaseFile)) )

        !Default structural solver
        StrucSolver%StructuralSolver = 'Viscoplastic'

        CALL SearchInFile( FreeChanelNum, '#STRUCT_LSOLVER', success )
        IF(success)THEN
            READ(FreeChanelNum,*)StrucSolver%LS_tolerance
            READ(FreeChanelNum,*)StrucSolver%LS_maxIters
        END IF

        CALL SearchInFile( FreeChanelNum, '#STRUCT_SOLVER', success )
        IF(success)THEN
          READ(FreeChanelNum,*)StrucSolver%StructuralSolver   !'Viscoplastic'
                                                              !'DynamicFEM'
        END IF

        CALL SearchInFile( FreeChanelNum, '#MEDSolver', success )
        IF(success)THEN
          StrucSolver%StructuralSolver = 'DynamicFEM'
        END IF

        !CLOSE( FreeChanelNum )

    END SUBROUTINE ReadSolverParameters


END MODULE StructuralModel
