
    !*****************************************************************************
    !Geomechanics
    !Eng. Luis Carlos de Sousa Junior
    !Geotechnical Engineering
    !University of Brasilia
    !*****************************************************************************

    PROGRAM IGeo

    USE IO
    USE OMP_LIB                 !bibliotecas OpenMP(R)
    USE IFPORT

    USE ProblemDefinition
    USE GeomechWell_main_these
    !USE Geotechnic
    USE WellStability

    !***************************************************************************
    !VARIABLES
    !***************************************************************************
    IMPLICIT NONE

    !Variables.
    INTEGER:: Chanel_log=9, Chanel_base=10
    CHARACTER(LEN=320)::  BaseFile, &     !file with basic information about the case.
                          Path            !work directory where the case file is located
    !and where temporary files will be created.
    CHARACTER(LEN=320)::buffer            !vari�vel para armazenamento de argumentos passados ao programa.
    INTEGER:: numARG
    LOGICAL:: test,&
              found
    CHARACTER(LEN=25):: ModelType     !Well
    !Soil

    test=.TRUE.
    
    !***************************************************************************
    !STARTING
    !***************************************************************************
    WRITE(*,*)'GEOTECHNICAL MODELING'
    WRITE(*,*)'STARTING'

    !***************************************************************************
    !FILES
    !***************************************************************************

    !Base file and path.
    numARG=0
    numARG=IARGC()

    IF(numARG>1)THEN

        CALL GETARG (1, buffer)
        Path = buffer
        CALL GETARG (2, buffer)
        BaseFile = buffer

    ELSEIF(test.EQV..FALSE.)THEN

        WRITE(*,*)'WTHCM: Entre com a pasta:'
        READ(*,*)Path
        WRITE(*,*)'WTHCM: Entre com o arquivo base:'
        READ(*,*)BaseFile

    ELSE
        !MAIN TESTES----------------------------------------------------------------
        !Path = 'C:\Users\cxnx\OneDrive - PETROBRAS\Documents\Pesquisas\Development\Geomechanics\Cases\BasicTests'
        !BaseFile = 'iGeo0.geo'                      !fluid compressibiliity system
        !BaseFile = 'iGeo1.geo'                      !transient linear flow test
        !BaseFile = 'iGeo2.geo'                      !transient radial flow test
        !BaseFile = 'iGeo3.geo'                      !transient radial fracture flow test
        !BaseFile = 'iGeo4.geo'                      !fracture propagation test
        !BaseFile = 'iGeo5a.geo'                     !poroelastic frature opening
        !BaseFile = 'iGeo5.geo'                      !poroelastic frature propagation
        !BaseFile = 'iGeo6.geo'                      !fracture containment
        !BaseFile = 'iGeo7.geo'                      !wellbore stability
        !BaseFile = 'iGeo8.geo'                      !radial flow for far field outflow option (equivalente to iGeo2.geo)
        !BaseFile = 'iGeo9.geo'                      !linear flow for far field outflow option (equivalente to iGeo1.geo)
        !BaseFile = 'iGeo10.geo'                     !fracture volume
        !BaseFile = 'iGeo11.geo'                     !imposed temperature variation (same data of iGeo6)
        !BaseFile = 'iGeo12.geo'                     !imposed porepressure variation (same data of iGeo6)
        !BaseFile = 'iGeo13.geo'                     !permeability damage due injected fluid (based on iGeo5.geo)
        !BaseFile = 'iGeo14.geo'                     !shear failure fracture - plasticity (based on iGeo5.geo)
        !Path = 'C:\Users\cxnx\OneDrive - PETROBRAS\Documents\Development\Geomechanics\Cases\HydraulicFracture'
        !BaseFile = 'MLS3B_injection.geo'
        
        !Path = 'C:\Users\cxnx\OneDrive - PETROBRAS\Documents\GMR\PWRI\RO152'
        !BaseFile = 'RO152.geo'
        
        !DRILLING TEST--------------------------------------------------------------
        Path = 'C:\Users\cxnx\OneDrive - PETROBRAS\Documents\Pesquisas\Development\Geomechanics\Cases\Drilling'
        !BaseFile = 'DrillingTest.geo'
        !BaseFile = 'MLS3B.geo'
        !BaseFile = 'B28.geo'
        BaseFile = 'SaltWell.geo'
        
        !ARTICLES-------------------------------------------------------------------
        !Path = 'C:\Users\cxnx\Desktop\Artigos\Sousa Jr artigo_2Dfrac\Model\Backup\Validation'
        !BaseFile = '2DvalidationCase.geo'
        !Path = 'C:\Users\cxnx\Desktop\Artigos\Sousa Jr artigo_2Dfrac\Model\Backup'
        !BaseFile = 'FlowTest_1.geo'
        !BaseFile = 'FlowTest_2.geo'
        !BaseFile = '2Dmodel_fracContainment.geo'
        
        !PETROBRAS------------------------------------------------------------------
        !Path = 'C:\Users\cxnx\Desktop\SEAP\HydraulicFraccture'
        !BaseFile = '2Da_Gc50_Ny50_k341_Ly160_R1.25_P48793.geo'
        
        !GEOTECHNICAL PROBLEM-------------------------------------------------------
        !Path = 'C:\Users\cxnx\OneDrive - PETROBRAS\Documents\Development\Geomechanics\Cases\FlowTest'
        !BaseFile = 'FlowTest_1.geo'
        !BaseFile = 'ThermalTest_Rad.geo'
        !BaseFile = 'SimpleFlowTest.geo'
        !BaseFile = '2D_BBG3_Analit.geo'
        
        !Path = 'C:\Users\cxnx\OneDrive - PETROBRAS\Documents\Development\Geomechanics\Cases\HydraulicFracture'
        !BaseFile = 'PKN_test.geo'
        !BaseFile = 'SES188_3D_TI.geo'
        !BaseFile = 'SES188_testRad.geo'
        !BaseFile = 'FracturePropagationTest.geo'
        !BaseFile = 'FracturePropagationTest_1.geo'
        !BaseFile = 'Validation2layersModel.geo'
        !BaseFile = 'test.geo'
        !BaseFile = 'SES188_test_TI.geo'
        !BaseFile = '2Dmodel_TI_SES188_KrelAbs_noPerm_Pressao_3D_TI_step2'
        !BaseFile = 'SoftSand.geo'
        
        !Path = 'C:\Users\cxnx\Desktop\FRATURA'
        !BaseFile = '3DRad_ATP3_PermQw.geo'
        !BaseFile = '3DRad_ATP3_PermQg6m.geo'
        
        !D.Sc. TESTS----------------------------------------------------------------
        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\BasicTests'
        !Path = 'C:\Users\cxnx\Documents\Development\Geomechanics\Cases'

        !BaseFile = 'TractionTest.geo'
        !BaseFile = 'UniaxialModel.geo'                       !OK 01/12/2017  30/09/2017   22/12/2017
        !BaseFile = 'UniaxialModel_Traction.geo'
        !BaseFile = 'UniaxialModel_Undrained.geo'             !OK 01/12/2017  01/09/2017
        !BaseFile = 'UniaxialModel_MC.geo'                    !OK 01/12/2017  30/09/2017   22/12/2017
        !BaseFile = 'UniaxialModel_HYP.geo'
        !BaseFile = 'UniaxialConstruction.geo'                !OK 01/12/2017
        !BaseFile = 'UniaxialConstructionOffshore.geo'        !OK 01/12/2017
        !BaseFile = 'UnxConsOffshoreMLS.geo'
        !BaseFile = 'Uniaxial_traction_1A.geo'                !OK 01/12/2017  30/09/2017   22/12/2017
        !BaseFile = 'Uniaxial_traction_1A_spring.geo'
        !BaseFile = '2DbearingCapacity.geo'                   !OK 01/12/2017
        !BaseFile = '2DbearingCapacity_finer.geo'
        !BaseFile = 'KIC_test.geo'
        !BaseFile = 'Beam_TractionTest.geo'
        !BaseFile = 'HydFracContainment_1.geo'
        !BaseFile = 'OffshoreSand_1.geo'                      !Flow only
        !BaseFile = 'HydFracContainment_2.geo'
        !BaseFile = 'Therzagui.geo'
        !BaseFile = 'FlowTest_1.geo'
        !BaseFile = 'FracContainment.geo'
        !BaseFile = 'Well_RAD.geo'
        !BaseFile = 'StressConcentration.geo'
        !BaseFile = 'B28_0.geo'
        !BaseFile = '2Dmodel_poromech_fracContainment.geo'
        !BaseFile = 'Embankment_test.geo'
        !BaseFile = 'Embankment.geo'
        !BaseFile = 'Brumadinho.geo'

        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\BasicTests\Thermal'
        !BaseFile = 'ThermalModel_1.geo'
        !BaseFile = 'ThermalModel_1_RAD.geo'
        !BaseFile = 'ThermalModel.geo'
        !BaseFile = 'test.geo'

        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\BasicTests\Uniaxial'
        !BaseFile = 'UniaxialModel.geo'

        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\BasicTests\Kic_tests'
        !BaseFile = 'KIC_test.geo'
        !BaseFile = 'KIC_containment.geo'
        !BaseFile = 'KIC_containment_perm.geo'

        !CONSTITTIVE MODELING-------------------------------------------------------
        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\BasicTests\ConstitutiveModel'
        !BaseFile = 'CAMCLAY_hydrostatic.geo'
        !BaseFile = 'CAMCLAY_UNIAXIAL_traj1.geo'
        !BaseFile = 'CAMCLAY_UNIAXIAL_traj2.geo'
        !BaseFile = 'CAMCLAY_UNIAXIAL_traj3.geo'
        !BaseFile = 'CAMCLAY_UNIAXIAL_traj4.geo'

        !D.Sc. TESTS----------------------------------------------------------------
        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\OnshoreSand'
        !BaseFile = 'OnshoreShallowSand.geo'
        !BaseFile = 'OnshoreShallowSand_1.geo'

        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\OffshoreSand'
        !BaseFile = '2D.geo'
        !BaseFile = 'test.geo'
        !BaseFile = 'OffshoreConstruction_1.geo'
        !BaseFile = 'OffshoreSand_1.geo'
        !BaseFile = 'OffshoreSand_1_BOX.geo'
        !BaseFile = 'OffshoreSand_2.geo'
        !BaseFile = 'OffshoreSand_3.geo'
        !BaseFile = 'OffshoreSand_horWell_1.geo'
        !BaseFile = 'OffshoreSand_horWell_3.geo'
        !BaseFile = 'OffshoreSand_VertWell_finner_flow.geo'

        !Path = 'C:\Users\cxnx\Documents\Doutorado\Codes\GeomecWELL\Cases\ShallowWell'
        !BaseFile = 'ShallowWell_1layer.geo'
        !BaseFile = 'ShallowWell.geo'

        !---------------------------------------------------------------------------
    END IF

    OPEN(Chanel_base,FILE=TRIM(TRIM(Path)//'/'//TRIM(BaseFile)))
    OPEN(Chanel_log,FILE = TRIM(TRIM(BaseFile))//".log")

    CALL PRINTHEAD(Chanel_log)

    CALL ReadModelType(Chanel_base,ModelType)
    IF( TRIM(ModelType) == 'WellModel0')THEN
        CALL GeomecW(Path,BaseFile)
    ELSEIF( TRIM(ModelType) == 'WellModel1')THEN
        CALL WellStabilityModule(Path,BaseFile)
    !ELSEIF( TRIM(ModelType) == 'Geotechnic')THEN
    !    CALL GeotechnicalSimulation(Path,BaseFile)
    ELSE
        CALL PrintMessage( .FALSE., .FALSE. , 'Wrong Model Type' )
        STOP
    END IF

    CALL PrintMessage( .FALSE., .FALSE., 'SIMULATION ENDED.' )
    CALL SLEEP(3)

    END PROGRAM IGeo

