
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE MeshIO


  CONTAINS

    !Exporta grid em fornato padrão.
    SUBROUTINE MeshToFile( Path, FileName, FileNum, Coords, ElementIncidence )

        USE IO

        IMPLICIT NONE

        CHARACTER(LEN=*):: Path, FileName
        INTEGER,INTENT(IN):: FileNum

        REAL(8),INTENT(IN):: Coords(:,:)
        INTEGER,INTENT(IN):: ElementIncidence(:,:)

        INTEGER:: numNodes, numElements
        INTEGER:: i

        CALL PrintMessage( .FALSE., .FALSE., '...writing mesh file.' )

        OPEN(FileNum,FILE = TRIM( TRIM(Path)//'\'//TRIM(FileName)//'.out.msh' ))!,FORM="UNFORMATED" )

        numNodes = SIZE( Coords(1,:) )
        numElements = SIZE( ElementIncidence(1,:) )

        !Exportando coordenadas nodais.
        WRITE(FileNum,*)numNodes
        DO i=1, numNodes
            WRITE(FileNum,*)Coords(1,i),Coords(2,i),Coords(3,i)
        END DO

        !Exportando incidência dos elementos finitos.
        WRITE(FileNum,*)numElements
        DO i=1, numElements
            WRITE(FileNum,*)ElementIncidence(:,i)
        END DO

        RETURN
    END SUBROUTINE MeshToFile

END MODULE
