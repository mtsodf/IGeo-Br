
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  MODULE IO

  !*************************************************************************
  !DATA STRUCTURE.
  !*************************************************************************
  
  TYPE StepReport
    REAL(4):: time
    INTEGER:: globalStep,flowStep,structuralStep,fractureStep
    INTEGER:: linIters,nonLinIters
    INTEGER:: cycles                    !for dynamic solver
    REAL(4):: maxVAL
    CHARACTER(LEN=3):: stepResult       !P - proceed
                                        !R - repead step
                                        !F - failed
  END TYPE

  !*************************************************************************
  !FUNCTIONS AND SUBROUTINES.
  !*************************************************************************
  CONTAINS

  !Imprime mensagens na tela.
  SUBROUTINE PrintMessage( SpaceBefore, LineBefore, Message )

    LOGICAL,INTENT(IN):: SpaceBefore, LineBefore
    CHARACTER(LEN=*),INTENT(IN):: Message
    
    INTEGER:: Chanel_log=9

    IF(SpaceBefore) PRINT*,""
    IF(LineBefore)  PRINT*,"CGM: ----------------------------------------------------"
    WRITE(*,*)      "CGM: ",TRIM(Message)
    
    IF(SpaceBefore)     WRITE(Chanel_log,*)""
    IF(LineBefore)      WRITE(Chanel_log,*)"CGM: ----------------------------------------------------"
    WRITE(Chanel_log,*) "CGM: ",TRIM(Message)

  END SUBROUTINE

  !Imprime cabeÁalho do programa.
  SUBROUTINE PRINTHEAD( ChanelFile )

    INTEGER,INTENT(IN):: ChanelFile   !usualmente .log.

    PRINT*, " "
    PRINT*, "+++++++++++++++++++++++++++++++++++++++++++++"
    PRINT*, "                   CGM"
    PRINT*, "        Computational Geomechanics"
    PRINT*, "          UnB - Geotechnical Eng"
    PRINT*, "       Eng. Luis Carlos de Sousa Jr"
    PRINT*, "+++++++++++++++++++++++++++++++++++++++++++++"
    PRINT*, " "

    IF(ChanelFile.NE.0)THEN
      WRITE(ChanelFile,*) " "
      WRITE(ChanelFile,*) "++++++++++++++++++++++++++++++++++++++++++++++"
      WRITE(ChanelFile,*) "                   CGM"
      WRITE(ChanelFile,*) "        Computational Geomechanics"
      WRITE(ChanelFile,*) "          UnB - Geotechnical Eng"
      WRITE(ChanelFile,*) "        Eng. Luis Carlos de Sousa Jr"
      WRITE(ChanelFile,*) "++++++++++++++++++++++++++++++++++++++++++++++"
      WRITE(ChanelFile,*) " "
    END IF

  END SUBROUTINE

  !Procura por uma palavra chave no arquivo.
  SUBROUTINE SearchInFile( ChanelNumber, KeyWord, success )

    IMPLICIT NONE
    INTEGER,INTENT(IN):: ChanelNumber
    CHARACTER(LEN=*),INTENT(IN):: KeyWord
    LOGICAL,INTENT(IN OUT):: success

    INTEGER:: Reason
    CHARACTER(LEN=320):: text

    success = .FALSE.

    REWIND( ChanelNumber )
    Reason=1
    DO WHILE(Reason>=0)
      READ( ChanelNumber, *, IOSTAT=Reason ) text
      IF( TRIM(text)==KeyWord )THEN
        success = .TRUE.
        WRITE(text,*)'...',TRIM(KeyWord),' found.'
        CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
        EXIT
      END IF
    END DO

   !IF(.NOT.success) THEN
   ! WRITE(text,*)'...',KeyWord,' not found.'
   ! CALL PrintMessage( .FALSE., .FALSE., TRIM(text) )
   !END IF

  END SUBROUTINE SearchInFile

  !Routine to read a table separeted between '>' and '<'.
  !The file FileChanel reader should be in the postion of '>'
  SUBROUTINE ReadTable( FileChanel, numberOfParameters, table )

    IMPLICIT NONE
    INTEGER,INTENT(IN):: FileChanel, NumberOfParameters
    CHARACTER(LEN=50),ALLOCATABLE,INTENT(OUT):: table(:,:)

    CHARACTER(LEN=320):: text
    INTEGER:: i,count

    !Counting the positions.
    READ(FileChanel,*) text  !>
    count = -1
    DO WHILE( TRIM(text) .NE. '<' )
      count = count  + 1
      READ(FileChanel,*) text
    END DO
    ALLOCATE(table(NumberOfParameters,count))

    !READING DATA.
    DO i = 1, (count+1)
      BACKSPACE(FileChanel)
    END DO
    text=''
    DO i=1,count
      READ(FileChanel,*) table(:,i)
    END DO

  END SUBROUTINE


  SUBROUTINE PrintStepReport(time,globalStep,flowStep,structuralStep,fractureStep,&
                            linIters,nonLinIters,maxVAL,stepResult)
  
    IMPLICIT NONE
    REAL(4),INTENT(IN):: time
    INTEGER,INTENT(IN):: globalStep,flowStep,structuralStep,fractureStep
    INTEGER,INTENT(IN):: linIters,nonLinIters
    REAL(4),INTENT(IN):: maxVAL
    CHARACTER(LEN=3),INTENT(IN):: stepResult        !FW - move step forward
                                                    !R - repead step
                                                    !F - failed
    CHARACTER(LEN=100):: text
    
    WRITE(text,'(E11.4,I7,I7,I7,I7,I7,I7,E11.3,A4)')&
        time/60,&
        globalStep,flowStep,structuralStep,fractureStep,&
        linIters,nonLinIters,maxVAL,TRIM(stepResult)
    
    CALL PrintMessage( .FALSE., .FALSE., trim(text) )
    
  END SUBROUTINE
  
  SUBROUTINE PrintStepReportHEAD()
    
    IMPLICIT NONE
    
    CHARACTER(LEN=250):: text
    
    WRITE(text,*) '-----------------------------------------------------------------------'
    CALL PrintMessage( .FALSE., .FALSE., trim(text) )
    WRITE(text,*) '   time.min  Gstep   Fstep  Sstep FRstep Liter Nliter  Max.Val Result'
    CALL PrintMessage( .FALSE., .FALSE., trim(text) )
    WRITE(text,*) '-----------------------------------------------------------------------'
    CALL PrintMessage( .FALSE., .FALSE., trim(text) )
    
  END SUBROUTINE

  
  END MODULE IO
