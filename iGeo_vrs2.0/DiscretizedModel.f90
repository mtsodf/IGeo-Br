

  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  !******************************************
  !ESTRUTURA DE DADOS PARA REPRESENTAR UM
  !MODELO DISCRETO.
  !******************************************

  MODULE DiscretizedModel

  USE ConstitutiveGeoModels
  
  !*************************************************************************
  !DATA STRUCTURE.
  !*************************************************************************
  
  TYPE VariablesContainer

    !MAIN VARIABLES------------------------------------------------------
    REAL(8),ALLOCATABLE:: Displacement(:)               !dofs displacement.
    REAL(8),ALLOCATABLE:: EffectiveStressTensor(:,:,:)  !gauss points effective stress tensor.
    REAL(8),ALLOCATABLE:: StrainTensor(:,:,:)           !gauss points total strain tensor.
    REAL(8),ALLOCATABLE:: PlasticStrainTensor(:,:,:)    !gauss points plastic strain tensor.
    REAL(8),ALLOCATABLE:: DamageShear(:)                !Shear Strength damage.
    REAL(8),ALLOCATABLE:: DamageTraction(:)             !Shear Strength damage.
    REAL(8),ALLOCATABLE:: DamageStiffness(:)            !Young Modulus damage.
    LOGICAL,ALLOCATABLE:: IsFractured(:)                !Flag for fractured element.
    REAL(8),ALLOCATABLE:: Porepressure(:)               !nodal porepressure.
    REAL(8),ALLOCATABLE:: ElementAVGporepressure(:)     !element average porepressure.
    REAL(8),ALLOCATABLE:: Fracpressure(:)               !element fracture pressure.
    REAL(8),ALLOCATABLE:: Temperature(:)                !nodal temperature.
    REAL(8),ALLOCATABLE:: WaterSaturation(:)            !nodal water saturation.
    REAL(8),ALLOCATABLE:: PoreFluidBulkMod(:)           !pore fluid Bulk Modulus.
    REAL(8),ALLOCATABLE:: PoreStrain(:)                 !pore strain rate.
    REAL(8),ALLOCATABLE:: PoreCompressibility(:)
    REAL(8),ALLOCATABLE:: YieldFunc(:,:)                !Yield function at Gauss Points.
    REAL(8),ALLOCATABLE:: YieldEnvelopeTOL(:)
    REAL(8),ALLOCATABLE:: PlasticStrainWork(:)
    REAL(8),ALLOCATABLE:: StrainWorkDensity(:,:)
    REAL(8),ALLOCATABLE:: YYstrainDesityEnergy(:)
    REAL(8),ALLOCATABLE:: Jintegral(:)
    REAL(8),ALLOCATABLE:: JintegralCum(:,:)            !(:,iel) - (1,iel) - Zup; (2,iel) - Zdown; (3,iel) - X
    REAL(8),ALLOCATABLE:: Gamma(:)
    
    ! FOR DYNAMIC SOLVER-------------------------------------------------------
    ! Structural
    REAL(8),ALLOCATABLE:: Displacement0(:,:),Displacement1(:,:)   !nodes (3)
    ! Porous fluid
    CHARACTER(LEN=150),ALLOCATABLE:: ComponentName(:)             !name of fluid components
    REAL(8),ALLOCATABLE:: Porepressure0(:),Porepressure1(:)       !for each cell
    REAL(8),ALLOCATABLE:: Temperature0(:),Temperature1(:)         !for each cell
    REAL(8),ALLOCATABLE:: Sw0(:,:),Sw1(:,:)                       !cells - for each component
    REAL(8),ALLOCATABLE:: PorousVelocity(:,:)                     !for each interface (:,x) and component (x,:)
    ! Fracture
    REAL(8),ALLOCATABLE:: FractureNodes(:,:)                      !4 nodes element
    REAL(8),ALLOCATABLE:: FracturePressure(:,:)                   !fracture element's pressure
    REAL(8),ALLOCATABLE:: FractureSat(:,:)                        !fracture element's saturation for each componnent
    REAL(8),ALLOCATABLE:: FractureVelocity(:,:)                   !fracture element's (:,x) and component (x,:)
    
  END TYPE
  
  TYPE FracturedRock
    
    INTEGER:: ElementNumber
    INTEGER:: GaussPoint
    
    REAL(4):: EigenVectors(3,3)
    INTEGER:: EigenS3postion
    
    REAL(4):: Damage
    
  END TYPE
     
  TYPE ShearDamage
      
    INTEGER,ALLOCATABLE:: ElementList(:)
    REAL(8):: LimitPlasticDeviatoricStrain
    REAL(8):: ResidualDamage
    REAL(8):: ResidualStiffnessDamage
      
  END TYPE
    
  TYPE StiffnessDamage
      
    INTEGER,ALLOCATABLE:: ElementList(:)
    REAL(8):: MaxMeanEffectiveStress              !Liquefaction limit.
    REAL(8):: ShearDamage,ElasticDamage
      
  END TYPE
  
  !FRAC DATA.
  TYPE PermModif
      
    INTEGER,ALLOCATABLE:: ElementList(:)
    REAL(8):: PermIncrementByPlasticDilation      !mD/0.001
    REAL(8):: PermIncrementByPlasticShear         !mD/0.001
      
  END TYPE
  
  TYPE ContinuousDamageMechancis
    
    REAL(8):: ST
    REAL(8):: EyyMax
    
  END TYPE
  
  TYPE SoftSandModel
    
    INTEGER:: elementNumber
    REAL(8):: p_ref
    REAL(8):: m
    REAL(8):: eq0
    
  END TYPE

  TYPE BoundaryNode
    
    LOGICAL:: Xboundary = .FALSE.,Yboundary = .FALSE.,Zboundary = .FALSE.
    
  END TYPE
  
  TYPE CreepModel
      
      REAL(8):: mu0
      REAL(8):: Qc
      REAL(8):: N2
      
  END TYPE
  
  TYPE, PUBLIC:: Model

    CHARACTER(LEN=18):: SolverType              !'Elastoplastic'
                                                !'DynamicFEM'
                                                          
    !FOLDER AND FILES-----------------------------------------------------
    CHARACTER(LEN=320):: Path,&          !folder for work files.
                         BaseFile        !main ou base file.

    !8 NODE HEX MESH------------------------------------------------------
    INTEGER:: NumPartitions                 !mesh partitions for parallel processing.
    INTEGER:: NumElements
    INTEGER:: NumNodes
    INTEGER:: NumIntegrationPoints

    CHARACTER(LEN=320):: element            !('hexahedron' ou 'tetrahedron').

    REAL(8),ALLOCATABLE:: ElementVolume(:)
    INTEGER:: GridDivision(3)               !NX,NZ,NY for structured meshes.
    REAL(8),ALLOCATABLE:: X(:),Y(:),Z(:)    !nodal coordenates.
    INTEGER,ALLOCATABLE:: NodeDofs(:,:)     !(nf) node displacement DOFs.
    INTEGER,ALLOCATABLE:: ElementNodes(:,:) !(g_num) nodal incidence.
    INTEGER,ALLOCATABLE:: ElementDOFs(:,:)  !(g_g) element structural dofs.
    INTEGER:: NumEquations                  !(neq) number of structural dofs.
    INTEGER,ALLOCATABLE:: ElementIJK(:,:)   !I,J,K indexes for each element. (3,NumElements)
    
    REAL(4),ALLOCATABLE:: ElementMinLength(:)   !Min dimension for dynamic explicity solver.

    !STRUCTURAL SPECIFICATION--------------------------------------------
    INTEGER::   NumDimensions,&             !ndim = 3 for 3D problems.
                NumNodesPerElement,&        !nod - number of corner nodes per element (8).
                NumDOFUperElement,&         !ndof - number of displacement DOFs per element (8*3).
                NumDOFPperElement,&         !nodf - number of porepressure dofs per element (constant in all nodes - just for coupling).
                NumDOFUperNode,&            !nodof = 3 for 3D problems (displacement DOFs per node).
                NumTensorComponents         !nst = 6 for 3D problems (number of stress and strain components).

    REAL(8),ALLOCATABLE:: Volume(:)         !elements volumes.
        
    INTEGER,ALLOCATABLE:: ConstitutiveGeoModel(:)   !0 - Elastic
                                                    !1 - Traction damage (frac.) [Explicity Solver]
                                                    !2 - Mohr-Coulomb
                                                    !3 - CAM-CLAY modified
                                                    !4 - Double Mohr-Coulomb     [Explicity Solver]
                                                    !5 - Soil Model              [Explicity Solver]
    TYPE(SoftSandModel),ALLOCATABLE:: SandModel(:)
        
    LOGICAL,ALLOCATABLE:: Undrained(:)

    !Internal total stress load.
    REAL(8),ALLOCATABLE:: InternalForce(:)
        
    !Active Elements.
    LOGICAL,ALLOCATABLE:: ActiveElements(:)
    LOGICAL,ALLOCATABLE:: ActiveToFlowElements(:)   !active to porous flow.
  
    !Creep flag
    LOGICAL:: Creep = .FALSE.

    !ELEMENT PROPERTIES--------------------------------------------------
    !7 digits precision- REAL(8).
    !REAL(8),ALLOCATABLE:: YoungMod(:)
    !REAL(8),ALLOCATABLE:: Poisson(:)
    REAL(8),ALLOCATABLE:: BiotCoef(:)
    REAL(8),ALLOCATABLE:: ThermalExpCoeff(:)
    REAL(8),ALLOCATABLE:: ThermalConductivity(:)
    REAL(8),ALLOCATABLE:: RockSpecificHeatCapacity(:)
    REAL(8),ALLOCATABLE:: FluidSpecificHeatCapacity(:)
    REAL(8),ALLOCATABLE:: RockFluidHeatTransferCoef(:)
    REAL(8),ALLOCATABLE:: Cohesion(:)
    REAL(8),ALLOCATABLE:: FrictionAngle(:)
    REAL(8),ALLOCATABLE:: DilationAngle(:)
    INTEGER,ALLOCATABLE:: ElementTractionStrengthModel(:)
    REAL(8),ALLOCATABLE:: NTG(:)
    REAL(8),ALLOCATABLE:: POR(:)
    REAL(8),ALLOCATABLE:: KH(:)
    REAL(8),ALLOCATABLE:: KV(:)
    REAL(8),ALLOCATABLE:: Density(:)              !kg/m³
    REAL(4),ALLOCATABLE:: MCCparameters(:,:)      !Modified CAM-CLAY parameters (Pc,M,kappa,lambda,e0,P0,ec,evol_vp)
    REAL(4),ALLOCATABLE:: DAMAGEparameters(:,:)   !Damage parameters (traction strength St, strain at St, max strain)
    REAL(8),ALLOCATABLE:: Gc(:)                   !Critical release energy.
    LOGICAL,ALLOCATABLE:: CDMflag(:)
    INTEGER,ALLOCATABLE:: CDMmodelNum(:)
    TYPE(TransverseIsotroplicElasticity),ALLOCATABLE:: TransverIsotropicProps(:)
    INTEGER,ALLOCATABLE:: ElementCreep(:)         !pointer with the creep elements (position 0 - number of elements with creep
    TYPE(CreepModel),ALLOCATABLE:: CreepData(:)   !data for creep elements in the sequence ElementCreep(:)

    !MAIN VARIABLES------------------------------------------------------
    REAL(8),ALLOCATABLE:: Displacement(:)               !displacement.
    REAL(8),ALLOCATABLE:: EffectiveStressTensor(:,:,:)  !gauss points effective stress tensor.
    REAL(8),ALLOCATABLE:: StrainTensor(:,:,:)           !gauss points total strain tensor.
    REAL(8),ALLOCATABLE:: PlasticStrainTensor(:,:,:)    !gauss points plastic strain tensor.
    REAL(8),ALLOCATABLE:: CreepStrainIncrementTensor(:,:,:)    !gauss points incremental creep strain tensor.
    REAL(8),ALLOCATABLE:: CreepStrainTensor(:,:,:)      !gauss points creep strain tensor.
    REAL(8),ALLOCATABLE:: DamageShear(:)                !Shear Strength damage.
    REAL(8),ALLOCATABLE:: DamageTraction(:)             !Shear Strength damage.
    REAL(8),ALLOCATABLE:: DamageStiffness(:)            !Young Modulus damage.
    REAL(8),ALLOCATABLE:: DamagePerm(:)                 !Perm damage.
    LOGICAL,ALLOCATABLE:: IsFractured(:)                !Flag for fractured element.
    REAL(8),ALLOCATABLE:: Porepressure(:)               !nodal porepressure.
    REAL(8),ALLOCATABLE:: ElementAVGporepressure(:)     !element average porepressure.
    REAL(8),ALLOCATABLE:: Fracpressure(:)               !element fracture pressure.
    REAL(8),ALLOCATABLE:: FracAperture(:)               !element fracture aperture.
    REAL(8),ALLOCATABLE:: Temperature(:)                !nodal temperature.
    REAL(8),ALLOCATABLE:: WaterSaturation(:)            !nodal water saturation.
    REAL(8),ALLOCATABLE:: PoreFluidBulkMod(:)           !pore fluid Bulk Modulus.
    REAL(8),ALLOCATABLE:: PoreStrain(:)                 !pore strain rate.
    REAL(8),ALLOCATABLE:: PoreCompressibility(:)
    REAL(8),ALLOCATABLE:: YieldFunc(:,:)                !Yield function at Gauss Points.
    REAL(8),ALLOCATABLE:: YieldEnvelopeTOL(:)
    REAL(8),ALLOCATABLE:: PlasticStrainWork(:)
    REAL(8),ALLOCATABLE:: StrainWorkDensity(:,:)
    REAL(8),ALLOCATABLE:: YYstrainDesityEnergy(:)
    REAL(8),ALLOCATABLE:: Jintegral(:)
    REAL(8),ALLOCATABLE:: JintegralCum(:,:)             !(:,iel) - (1,iel) - Zup; (2,iel) - Zdown; (3,iel) - X
    REAL(8),ALLOCATABLE:: Gamma(:)
    REAL(8):: CreepWork(3)
        
    !DAMAGE MODEL-------------------------------------------------------
    !Model:
    !       prop = prop0 - decay * var
    !       var: plastic shear strain (Von Misses).
    !Shear strength decay.
    REAL(8),ALLOCATABLE:: ShearStrengthDecay(:,:)       !(1,element):
        
    !SECONDARY VARIABLES (POST-PROC)------------------------------------
    REAL(8),ALLOCATABLE:: CilindricalStress(:,:,:)      !gauss points stress tensor.: 1-rad, 2-circunferencial, 3-shear.
        
    !SCALING RELATIONS--------------------------------------------------
    INTEGER,ALLOCATABLE:: CoarseElement(:)              !Relation to coarse model elements.
    INTEGER,ALLOCATABLE:: NodesDependence(:,:)          !Border nodes dependence.
                                                        !(1:4,node at fine mesh) <- node 1.
                                                        !(2,node) <- node 2.
    REAL(8),ALLOCATABLE:: NodesDependenceWeight(:,:)    !dependeces weight for each node in NodesDependence.
        
    !DISCRETE FRACTURES-------------------------------------------------
    LOGICAL,ALLOCATABLE:: FracturedElementModel(:)      !(ELEMENT)
    !LOGICAL,ALLOCATABLE:: IsFractured(:,:)             !(GAUSS POINT, ELEMENT)
    INTEGER,ALLOCATABLE:: FracturedElementNumber(:,:)   !pointer for postion at FracturedRockElements.
    TYPE(FracturedRock),ALLOCATABLE:: FracturedRockElements(:)
    
    !DAMAGE MODEL-------------------------------------------------------
    TYPE (ContinuousDamageMechancis),ALLOCATABLE:: CDM(:)
    TYPE (ShearDamage),ALLOCATABLE:: ShearStrainDamage(:)
    TYPE (StiffnessDamage),ALLOCATABLE:: Liquefaction(:)
    
    !PERM MODIF DATA----------------------------------------------------
    TYPE (PermModif),ALLOCATABLE:: PermModel(:)
    
    !Soft sand model parameters-----------------------------------------
    REAL(8):: ShearHardeningPar(7)
    REAL(8):: Cp(7)
    REAL(8):: Cr(7)
    REAL(8):: Phi(7)
    REAL(8):: Phir(7)
    REAL(8):: DilAngle(7)
    
    !Material Group Model-----------------------------------------------
    !Aproach used for material behaviour.
    CHARACTER(LEN=20),ALLOCATABLE:: MaterialName(:)                 !material name for each element.
    INTEGER,ALLOCATABLE:: MaterialModelIdex(:)                      !material number in SoilMaterials(:)
    TYPE(SoilModel),ALLOCATABLE:: SoilModels(:)
    !-------------------------------------------------------------------
    
    !Boundary Condition-------------------------------------------------
    TYPE(BoundaryNode),ALLOCATABLE:: BoundaryNodes(:)
    INTEGER:: NumBoundaryNodes
    !-------------------------------------------------------------------
    
    !Duplicated Flow Variables------------------------------------------
    REAL(8),ALLOCATABLE:: InjRate(:)
    REAL(8),ALLOCATABLE:: Flow(:)
    
    !-------------------------------------------------------------------
    
    LOGICAL:: UpdateGeometry = .FALSE.
    
  END TYPE Model

  !*************************************************************************
  !FUNCOES E SUBROTINAS.
  !*************************************************************************
  CONTAINS

    !CONSTRUCAO DO MODELO GEOMECANICO - wellbore model.
    !#NWB_STRUCTUREDMESH card for mesh description.
    SUBROUTINE BuildWellGeomechanicalModel( DiscreteModel, Path, BaseFile, WellModel, WellData, Problem )

      USE NearWellboreModel
      USE ProblemDefinition
      USE IO
      USE NumericalUtilities
      USE FEMlib

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      CHARACTER(LEN=320),INTENT(IN):: Path,&                !diretOrio dos arquivos e de trabalho.
                                      BaseFile              !arquivo base de dados da simulacao.
      TYPE(NWB_Model),INTENT(IN OUT):: WellModel
      TYPE(WellboreData),INTENT(IN):: WellData
      TYPE(PD),INTENT(IN):: Problem

      INTEGER:: elementNum
      REAL(8):: Kg,Kb,biot
      REAL(8):: rho
        
      REAL(8):: stress(6)
      INTEGER:: i,j,k,ip,n,node,counter
      INTEGER:: irad,iz,itheta
        
      INTEGER,ALLOCATABLE:: sumVector(:),num(:),g(:)
        
      REAL(8),ALLOCATABLE:: pp(:)             !element nodal porepressure.
      REAL(8):: p
        
      INTEGER:: ndim,&                        !dimensões do problema
                nip,&                         !número de pontos de integração
                nodf                          !número de nós com grau de liberdade de poro-pressão.
      CHARACTER(LEN=320)::element
      REAL(8),ALLOCATABLE:: points(:,:),&     !coordenadas dos pontos de integração
                            weights(:),&      !peso dos pontos de integraçã
                            funf(:)           !funções de forma de poro-pressão
        
      REAL(8):: por,e,poisson
        
      REAL(8),ALLOCATABLE:: cw(:)
      
      LOGICAL,ALLOCATABLE:: allocatedNodes(:)
      
      REAL(8):: ShearHardeningPar(7)
      REAL(8):: Cp(7)
      REAL(8):: Cr(7)
      REAL(8):: Phi(7)
      REAL(8):: Phir(7)
      REAL(8):: DilAngle(7)
      
      LOGICAL:: fullStressTensor
        
      INTEGER:: I1,I2,J1,J2,K1,K2,index
      INTEGER:: index1_1,index1_2,index2_1,index2_2
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !MODEL BASIC DATA.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      CALL SetFiniteElementType( DiscreteModel )

      !Dados de malha.
      DiscreteModel%NumElements = WellModel%NumElements
      DiscreteModel%NumIntegrationPoints = 1  !8
      DiscreteModel%NumNodes = WellModel%NumNodes
      DiscreteModel%GridDivision(1) = WellModel%NRAD!WellModel%NTHETA
      DiscreteModel%GridDivision(2) = WellModel%NZ
      DiscreteModel%GridDivision(3) = WellModel%NTHETA!WellModel%NRAD
      !A condicao de contorno de deslocamento e imposta em matriz de C.C.
      DiscreteModel%NumEquations = MAXVAL(WellModel%NodeDOFs) !3*DiscreteModel%NumNodes 
      
      !Arquivos.
      DiscreteModel%Path = Path
      DiscreteModel%BaseFile = BaseFile

      !Coordenadas nodais.
      ALLOCATE( DiscreteModel%X(WellModel%NumNodes), DiscreteModel%Y(WellModel%NumNodes), &
                  DiscreteModel%Z(WellModel%NumNodes) )
      DiscreteModel%X = WellModel%Coords(1,:)
      DiscreteModel%Y = WellModel%Coords(2,:)
      DiscreteModel%Z = WellModel%Coords(3,:)
        
      ALLOCATE(num(8),g(24))

      !Elements Volume.
      ALLOCATE( DiscreteModel%ElementVolume(WellModel%NumElements) )
      DiscreteModel%ElementVolume = 0.0
        
      !Element Index.
      ALLOCATE( DiscreteModel%ElementIJK(3,WellModel%NumElements) )
      DiscreteModel%ElementIJK = 0
        
      !Active Elements.
      ALLOCATE( DiscreteModel%ActiveElements(WellModel%NumElements) )
      DiscreteModel%ActiveElements = .TRUE.
      ALLOCATE( DiscreteModel%ActiveToFlowElements(WellModel%NumElements) )
      DiscreteModel%ActiveToFlowElements = .TRUE.
      
      ALLOCATE( DiscreteModel%ElementMinLength(WellModel%NumElements) )
      DiscreteModel%ElementMinLength = 0.0

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !PROPERTIES.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
      !Lendo algumas propriedades que são muito utilizadas.
      !ALLOCATE( & !DiscreteModel%YoungMod(WellModel%NumElements), &
      !    DiscreteModel%Poisson(WellModel%NumElements), 
      ALLOCATE( DiscreteModel%TransverIsotropicProps(WellModel%NumElements),&
          DiscreteModel%BiotCoef(WellModel%NumElements), &
          DiscreteModel%ThermalExpCoeff(WellModel%NumElements), &
          DiscreteModel%Cohesion(WellModel%NumElements), DiscreteModel%FrictionAngle(WellModel%NumElements),&
          DiscreteModel%DilationAngle(WellModel%NumElements), & !DiscreteModel%TractionStrength(WellModel%NumElements),&
          DiscreteModel%ElementTractionStrengthModel(WellModel%NumElements),&
          DiscreteModel%NTG(WellModel%NumElements),&
          DiscreteModel%POR(WellModel%NumElements),DiscreteModel%KH(WellModel%NumElements),&
          DiscreteModel%ThermalConductivity(WellModel%NumElements),&
          DiscreteModel%RockSpecificHeatCapacity(WellModel%NumElements),&
          DiscreteModel%FluidSpecificHeatCapacity(WellModel%NumElements),&
          DiscreteModel%RockFluidHeatTransferCoef(WellModel%NumElements),&
          DiscreteModel%KV(WellModel%NumElements),DiscreteModel%MCCparameters(8,WellModel%NumElements),&
          DiscreteModel%DAMAGEparameters(3,WellModel%NumElements), DiscreteModel%Gc(WellModel%NumElements),&
          DiscreteModel%CDMflag(WellModel%NumElements), DiscreteModel%CDMmodelNum(WellModel%NumElements)) !,&
          !DiscreteModel%TRACflag(WellModel%NumElements),DiscreteModel%TRACmodelNum(WellModel%NumElements) )
      DiscreteModel%NTG = 0.0
      DiscreteModel%POR = 0.0
      DiscreteModel%KH = 0.0
      DiscreteModel%KV = 0.0
      DiscreteModel%ThermalConductivity = 0.0
      DiscreteModel%ThermalExpCoeff = 0.0
      DiscreteModel%RockSpecificHeatCapacity = 0.0
      DiscreteModel%FluidSpecificHeatCapacity = 0.0
      DiscreteModel%RockFluidHeatTransferCoef = 0.0
      DiscreteModel%Cohesion = 0.0
      DiscreteModel%FrictionAngle = 0.0
      DiscreteModel%DilationAngle = 0.0
      !DiscreteModel%TractionStrength = 0.0
      DiscreteModel%ElementTractionStrengthModel = 0
      DiscreteModel%MCCparameters = 0.0
      DiscreteModel%DAMAGEparameters = 0.0
      DiscreteModel%Gc = 0.0
      DiscreteModel%CDMflag = .FALSE.
      DiscreteModel%CDMmodelNum = 0
      !DiscreteModel%TRACflag = .FALSE.
      !DiscreteModel%TRACmodelNum = 0
      
      REWIND(13)
      elementNum = 0
      DO itheta=1, WellModel%NTHETA
        DO iz=1, WellModel%NZ
          DO irad=1, WellModel%NRAD
            elementNum = elementNum + 1

            DiscreteModel%POR(elementNum) = WellData%Porosity(iz)
            DiscreteModel%NTG(elementNum) = 1.0
            DiscreteModel%KH(elementNum) = WellData%EffectivePerms(1,iz) * 9.869D-16    !Kh
            DiscreteModel%KV(elementNum) = WellData%EffectivePerms(2,iz) * 9.869D-16    !Kv
            !DiscreteModel%YoungMod(elementNum) = WellData%YoungModulus(1,iz)*1.0E3      !kPa -> Pa
            DiscreteModel%TransverIsotropicProps(elementNum)%E1 = WellData%YoungModulus(1,iz)*1.0E3      !kPa -> Pa
            DiscreteModel%TransverIsotropicProps(elementNum)%E2 = WellData%YoungModulus(2,iz)*1.0E3      !kPa -> Pa
            !DiscreteModel%Poisson(elementNum) = WellData%Poisson(1,iz)
            DiscreteModel%TransverIsotropicProps(elementNum)%n1 = WellData%Poisson(1,iz)
            DiscreteModel%TransverIsotropicProps(elementNum)%n2 = WellData%Poisson(2,iz)
            IF(  WellData%EffectivePerms(1,iz)<1.0E-15 ) DiscreteModel%ActiveToFlowElements(elementNum) = .FALSE.  ! 1md cutoff

            IF( DiscreteModel%TransverIsotropicProps(elementNum)%E1 < DiscreteModel%TransverIsotropicProps(elementNum)%E2 ) THEN
                DiscreteModel%TransverIsotropicProps(elementNum)%minE = DiscreteModel%TransverIsotropicProps(elementNum)%E1
            ELSE
                DiscreteModel%TransverIsotropicProps(elementNum)%minE = DiscreteModel%TransverIsotropicProps(elementNum)%E2
            END IF
            IF( DiscreteModel%TransverIsotropicProps(elementNum)%n1 < DiscreteModel%TransverIsotropicProps(elementNum)%n2 ) THEN
                DiscreteModel%TransverIsotropicProps(elementNum)%minNi = DiscreteModel%TransverIsotropicProps(elementNum)%n1
            ELSE
                DiscreteModel%TransverIsotropicProps(elementNum)%minNi = DiscreteModel%TransverIsotropicProps(elementNum)%n2
            END IF
            
            Kg = WellData%Kg(iz)*1.0E9
            e = 0.5 * (DiscreteModel%TransverIsotropicProps(elementNum)%E1 + DiscreteModel%TransverIsotropicProps(elementNum)%E2)
            poisson = 0.5 * (DiscreteModel%TransverIsotropicProps(elementNum)%n1 + DiscreteModel%TransverIsotropicProps(elementNum)%n2)
            Kb = e / 3.0 / (1.0-2.0*poisson)
            DiscreteModel%BiotCoef(elementNum) = 1.0 - Kb / Kg
            IF(DiscreteModel%POR(elementNum)<1.0E-5) DiscreteModel%BiotCoef(elementNum) = 0.0
                    
            DiscreteModel%ElementIJK(1,elementNum) = irad
            DiscreteModel%ElementIJK(2,elementNum) = iz
            DiscreteModel%ElementIJK(3,elementNum) = itheta

          END DO
        END DO
      END DO

      !CONECTIVIDADE NODAL.
      ALLOCATE( DiscreteModel%ElementNodes(DiscreteModel%NumNodesPerElement,&
                                                      WellModel%NumElements) )
      DiscreteModel%ElementNodes = WellModel%ElementIncidence
        
      !NEAR WELLBORE PERMEABILITY.
      DO elementNum = 1, DiscreteModel%NumElements
        num = DiscreteModel%ElementNodes(:,elementNum)
        IF(ANY(WellData%WellBoreWallFlag(num)) .EQV. .TRUE.)THEN
          DiscreteModel%KH(elementNum) = DiscreteModel%KH(elementNum) * WellData%WellPermMult
          DiscreteModel%KV(elementNum) = DiscreteModel%KV(elementNum) * WellData%WellPermMult
        END IF
      END DO

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !DEGREES OF FREEDON.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
      ALLOCATE( DiscreteModel%NodeDofs(DiscreteModel%NumDOFUperNode,SIZE(WellModel%NodeDOFs(1,:))),&
                DiscreteModel%ElementDOFs(DiscreteModel%NumDOFUperElement,DiscreteModel%NumElements) )
      DiscreteModel%NodeDofs = WellModel%NodeDOFs   != 1
                                                    !CALL formnf(DiscreteModel%NodeDofs)
      DiscreteModel%ElementDOFs = WellModel%ElementDOFs
      !CALL BuildDOFmatrix(DiscreteModel%NumElements,DiscreteModel%NumNodesPerElement,&
      !                    DiscreteModel%NumDOFUperElement,DiscreteModel%ElementNodes,&
      !                    DiscreteModel%NodeDofs,DiscreteModel%ElementDOFs)

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !DYNAMIC PROPERTIES.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      ALLOCATE( DiscreteModel%Displacement(DiscreteModel%NumEquations), &
                DiscreteModel%EffectiveStressTensor(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                DiscreteModel%PlasticStrainTensor(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                DiscreteModel%CreepStrainIncrementTensor(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                DiscreteModel%CreepStrainTensor(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                DiscreteModel%YieldFunc(DiscreteModel%NumIntegrationPoints,&
                                        DiscreteModel%NumElements),&
                DiscreteModel%YieldEnvelopeTOL(DiscreteModel%NumElements),&
                DiscreteModel%StrainTensor(6,DiscreteModel%NumIntegrationPoints,&
                                            DiscreteModel%NumElements),&
                DiscreteModel%Porepressure(WellModel%NumNodes),&
                DiscreteModel%ElementAVGporepressure(DiscreteModel%NumElements),&
                DiscreteModel%Temperature(WellModel%NumNodes),&
                DiscreteModel%WaterSaturation(WellModel%NumNodes),&
                DiscreteModel%CilindricalStress(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                DiscreteModel%DamageShear(DiscreteModel%NumElements),&
                DiscreteModel%DamageTraction(DiscreteModel%NumElements),&
                DiscreteModel%DamageStiffness(DiscreteModel%NumElements),&
                DiscreteModel%DamagePerm(DiscreteModel%NumElements),&
                DiscreteModel%IsFractured(DiscreteModel%NumElements),&
                DiscreteModel%Fracpressure(DiscreteModel%NumElements),&
                DiscreteModel%FracAperture(DiscreteModel%NumElements),&
                DiscreteModel%PoreFluidBulkMod(DiscreteModel%NumElements),&
                DiscreteModel%PoreStrain(WellModel%NumNodes),&
                DiscreteModel%PoreCompressibility(DiscreteModel%NumElements),&
                DiscreteModel%PlasticStrainWork(DiscreteModel%NumElements),&
                DiscreteModel%StrainWorkDensity(DiscreteModel%NumIntegrationPoints,DiscreteModel%NumElements),&
                DiscreteModel%YYstrainDesityEnergy(DiscreteModel%NumElements),&
                DiscreteModel%Jintegral(DiscreteModel%NumElements),&
                DiscreteModel%JintegralCum(3,DiscreteModel%NumElements),&
                DiscreteModel%Gamma(DiscreteModel%NumElements),&
                DiscreteModel%InjRate(WellModel%NumNodes),&
                DiscreteModel%Flow(WellModel%NumNodes) )
      DiscreteModel%Displacement = 0.0
      DiscreteModel%EffectiveStressTensor = 0.0
      DiscreteModel%StrainTensor = 0.0
      DiscreteModel%PlasticStrainTensor = 0.0
      DiscreteModel%CreepStrainIncrementTensor = 0.0
      DiscreteModel%CreepStrainTensor = 0.0
      DiscreteModel%YieldFunc = 0.0
      DiscreteModel%YieldEnvelopeTOL = 0.0
      DiscreteModel%Porepressure = 0.0
      DiscreteModel%ElementAVGporepressure = 0.0
      DiscreteModel%Temperature = 0.0
      DiscreteModel%WaterSaturation = 0.0
      DiscreteModel%CilindricalStress = 0.0
      DiscreteModel%DamageShear = 0.0
      DiscreteModel%DamageTraction = 0.0
      DiscreteModel%DamageStiffness = 0.0
      DiscreteModel%DamagePerm = 0.0
      DiscreteModel%IsFractured = .FALSE.
      DiscreteModel%Fracpressure = 0.0
      DiscreteModel%FracAperture = 0.0
      DiscreteModel%PoreFluidBulkMod = 0.0
      DiscreteModel%PoreStrain = 0.0
      DiscreteModel%PoreCompressibility = 0.0
      DiscreteModel%PlasticStrainWork = 0.0
      DiscreteModel%StrainWorkDensity = 0.0
      DiscreteModel%YYstrainDesityEnergy = 0.0
      DiscreteModel%Jintegral = 0.0
      DiscreteModel%JintegralCum = 0.0
      DiscreteModel%Gamma = 1.0
      DiscreteModel%InjRate = 0.0
      DiscreteModel%Flow = 0.0
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !STRESS AND POREPRESSURE.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      element = DiscreteModel%element
      ndim = DiscreteModel%NumDimensions
      nip = DiscreteModel%NumIntegrationPoints
      nodf = 8  !porepressure dofs per element for integration by f.e..
      ALLOCATE(pp(nodf),funf(nodf),points(nip,ndim),weights(nip))  
      CALL sample(element,points,weights)

      IF( SIZE(WellData%Porepressure)<1 )THEN
        DO i=1,DiscreteModel%NumNodes
          DiscreteModel%Porepressure(i) = Problem%OWpressure - &
                                               (DiscreteModel%Z(i)-Problem%OWdepth) * Problem%nwGrad
        END DO
      ELSE
        
        elementNum = 0
        DO i=1, WellModel%NTHETA
          DO k=1, WellModel%NZ
            DO j=1, WellModel%NRAD
              elementNum = elementNum + 1
              num = DiscreteModel%ElementNodes(:,elementNum)
              DiscreteModel%Porepressure(num) = WellData%Porepressure(k)
            END DO
          END DO
        END DO
        
      END IF
      DiscreteModel%Porepressure = DiscreteModel%Porepressure * 1000.0
      
      fullStressTensor = .FALSE.
      IF( SIZE(WellData%MaxHorizontalStress)<1) fullStressTensor = .TRUE.
      
      elementNum = 0
      DO i=1, WellModel%NTHETA
        DO k=1, WellModel%NZ
              
          stress = 0.0
          
          IF(fullStressTensor)THEN
            stress = WellData%StressTensor(:,k)*1000.0
          ELSE
            stress(1) = WellData%MaxHorizontalStress(k)*1000.0
            stress(2) = WellData%MinHorizontalStress(k)*1000.0
            stress(3) = WellData%VerticalStress(k)*1000.0
          END IF
                
          DO j=1, WellModel%NRAD
            elementNum = elementNum + 1
                                
            num = DiscreteModel%ElementNodes(:,elementNum)
            biot = DiscreteModel%BiotCoef(elementNum)
            pp = DiscreteModel%Porepressure(num)
            
            DO ip=1,nip
              CALL shape_function(funf,points,ip)
              p = DOT_PRODUCT(funf,pp)
              DiscreteModel%EffectiveStressTensor(1,ip,elementNum) = -stress(1) + biot*p
              DiscreteModel%EffectiveStressTensor(2,ip,elementNum) = -stress(2) + biot*p
              DiscreteModel%EffectiveStressTensor(3,ip,elementNum) = -stress(3) + biot*p
              DiscreteModel%EffectiveStressTensor(4,ip,elementNum) =  stress(4)
              DiscreteModel%EffectiveStressTensor(5,ip,elementNum) =  stress(5)
              DiscreteModel%EffectiveStressTensor(6,ip,elementNum) =  stress(6)
            END DO
                    
            !Element average pressure.
            DiscreteModel%ElementAVGporepressure(elementNum) = SUM(DiscreteModel%Porepressure(num))/nodf
                    
          END DO
        END DO
      END DO

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
      ALLOCATE(DiscreteModel%InternalForce(DiscreteModel%NumEquations))
      DiscreteModel%InternalForce = 0.0
        
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !ELEMENT VOLUME.
      ALLOCATE(DiscreteModel%Volume(DiscreteModel%NumElements))
      DiscreteModel%Volume = 0.0
        
      !DENSITY.
      ALLOCATE(DiscreteModel%Density(DiscreteModel%NumElements))
      DiscreteModel%Density = 2500.0
!DENSITY BY LAYER'S VERTICAL STRESS
      IF( SIZE(WellData%VerticalStress)>0 )THEN
        elementNum = 0
        DO itheta=1, WellModel%NTHETA
          DO iz=1, WellModel%NZ
            DO irad=1, WellModel%NRAD
              elementNum = elementNum + 1
              IF(iz<WellModel%NZ)THEN
                rho = WellData%VerticalStress(iz+1) - WellData%VerticalStress(iz)
                rho = rho / (-WellData%LayersTops(iz+1) + WellData%LayersTops(iz))
                rho = rho * 1000.0
              ELSE
                !keep rho calculated previouslly.
              END IF
              IF(rho<1.0E-5) rho = 2500.0
              DiscreteModel%Density(elementNum) = rho
            END DO
          END DO
        END DO
      END IF
!END DENSITY

      !INITIAL PORE COMPRESSBILITY.
      !Element pore compressibility.
      DO i=1,WellModel%NumElements
        IF(DiscreteModel%POR(i)>0.0)THEN
          biot = DiscreteModel%BiotCoef(i)
          por = DiscreteModel%POR(i)
          e = 0.5 * (DiscreteModel%TransverIsotropicProps(i)%E1 + DiscreteModel%TransverIsotropicProps(i)%E2)
          !poisson = DiscreteModel%Poisson(i)
          poisson = 0.5 * (DiscreteModel%TransverIsotropicProps(i)%n1 + DiscreteModel%TransverIsotropicProps(i)%n2)
          IF(e>0) DiscreteModel%PoreCompressibility(i) = biot**2/por/e*(1+poisson)*(1-2*poisson)/(1-poisson)
        END IF
      END DO
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !CONSTITUTIVE MODEL.
      !0 - Elastic
      !1 - Traction damage (frac.)
      !2 - Mohr-Coulomb
      !3 - CAM-CLAY (modified)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      ALLOCATE(DiscreteModel%ConstitutiveGeoModel(DiscreteModel%NumElements))
      DiscreteModel%ConstitutiveGeoModel = 0  !ELASTIC.
      
      !DAMAGE.
      !TYPE(ContinuousDamageMechancis),ALLOCATABLE:: CDM(:)
      IF(SIZE(Problem%CDM(:))>0) THEN
        n = SIZE(Problem%CDM)
        ALLOCATE(DiscreteModel%CDM(n))
        
        DO i=1,n
          
          DiscreteModel%CDM(i)%ST = Problem%CDM(i)%ST
          DiscreteModel%CDM(i)%EyyMax = Problem%CDM(i)%EyyMax
          
          elementNum = 0
          DO itheta=1, WellModel%NTHETA
            DO iz=1, WellModel%NZ
              DO irad=1, WellModel%NRAD
                elementNum = elementNum + 1
                IF( iz>=Problem%CDM(i)%K1 .AND. iz<=Problem%CDM(i)%K2 .AND. itheta==1) THEN
                  DiscreteModel%CDMflag(elementNum) = .TRUE.
                  DiscreteModel%CDMmodelNum(elementNum) = i
                  DiscreteModel%Gc(elementNum) = Problem%CDM(i)%Gc
                END IF
              END DO
            END DO
          END DO
            
        END DO
        
      END IF

      DO i=1,SIZE(WellData%ConsModel(1,:))
        
        elementNum = 0
        n = 0
        DO itheta=1, WellModel%NTHETA
          DO iz=1, WellModel%NZ
            DO irad=1, WellModel%NRAD
              elementNum = elementNum + 1
                
              IF((iz>=WellData%ConsModel(2,i)).AND.(iz<=WellData%ConsModel(3,i)) )THEN
                
                IF(WellData%ConsModel(1,i)==1 .AND. &     !DAMAGE (no shear strength limit)   <- Explicit Solver.
                   DiscreteModel%CDMflag(elementNum) .EQV. .TRUE.)THEN
                  DiscreteModel%ConstitutiveGeoModel(elementNum) = 1
                
                ELSEIF(WellData%ConsModel(1,i)==2)THEN    !Mohr-Coulomb.
                  DiscreteModel%ConstitutiveGeoModel(elementNum) = 2
                  DiscreteModel%Cohesion(elementNum) = WellData%ConsModelParameters(1,i) !*1000.0
                  DiscreteModel%FrictionAngle(elementNum) = WellData%ConsModelParameters(2,i)
                  DiscreteModel%DilationAngle(elementNum) = WellData%ConsModelParameters(3,i)
                  
                ELSEIF(WellData%ConsModel(1,i)==3)THEN    !Modified CAM-CLAY.
                  DiscreteModel%ConstitutiveGeoModel(elementNum) = 3
                  DiscreteModel%MCCparameters(1:6,elementNum) = WellData%ConsModelParameters(1:6,i)
                  !Pc(kPa)  M  kappa  lambda   e0   P0(kPa)
                  DiscreteModel%MCCparameters(1,elementNum) = - DiscreteModel%MCCparameters(1,elementNum)
                  DiscreteModel%MCCparameters(6,elementNum) = - DiscreteModel%MCCparameters(6,elementNum)
                  
                ELSEIF(WellData%ConsModel(1,i)==4)THEN    !Hyperbolic Model.
                  DiscreteModel%ConstitutiveGeoModel(elementNum) = 4
                  DiscreteModel%Cohesion(elementNum) = WellData%ConsModelParameters(1,i) !*1000.0
                  DiscreteModel%FrictionAngle(elementNum) = WellData%ConsModelParameters(2,i)
                  
                ELSEIF(WellData%ConsModel(1,i)==5)THEN    !Double Mohr-Coulomb.
                  DiscreteModel%ConstitutiveGeoModel(elementNum) = 5
                  n = n + 1
                  DiscreteModel%Cohesion(elementNum) = WellData%ConsModelParameters(1,i) !*1000.0
                  DiscreteModel%FrictionAngle(elementNum) = WellData%ConsModelParameters(2,i)

                ELSE
                END IF
                
              END IF
                
            END DO
          END DO
        END DO
        
      END DO
    
      !Soft Sand Model.
      ALLOCATE( DiscreteModel%SandModel(n) )
      
      elementNum = 0
      n = 0
      DO itheta=1, WellModel%NTHETA
        DO iz=1, WellModel%NZ
          DO irad=1, WellModel%NRAD
            elementNum = elementNum + 1
                      
            DO i=1,SIZE(WellData%ConsModel(1,:))
              IF( (WellData%ConsModel(1,i)==5) .AND.&
                  (iz>=WellData%ConsModel(2,i)).AND.(iz<=WellData%ConsModel(3,i)) )THEN
                n = n + 1
                DiscreteModel%SandModel(n)%elementNumber = elementNum
                DiscreteModel%SandModel(n)%p_ref = WellData%ConsModelParameters(3,i)
                DiscreteModel%SandModel(n)%m = WellData%ConsModelParameters(4,i)
                DiscreteModel%SandModel(n)%eq0 = WellData%ConsModelParameters(5,i)
                
              END IF
            END DO

          END DO
        END DO
      END DO
      
      IF(n>0)THEN
      
        DATA ShearHardeningPar /   0,	0.003	 , 0.03	 ,   0.05,	  0.07,  0.1,	0.2  /
        DATA Cp                /2600,	 2365	 , 1650	 , 870	 ,  333	  , 43.3,	10   /
        DATA Cr                /2600,	  814  ,  568  ,  300	 ,  115	  , 10	, 10   /
        DATA Phi               /   0,	   11	 ,   28.1,	 31.8,   33	  , 34	, 34   /
        DATA Phir              /   0,	   14.2,   37.3,	 42.7,	 45	  , 46	, 46   /
        DATA DilAngle          /   0,      0 ,    0	 ,    3.2,	  6.3	, 13.8,	13.8 /
  
        DiscreteModel%ShearHardeningPar = ShearHardeningPar
        DiscreteModel%Cp = Cp
        DiscreteModel%Cr = Cr
        DiscreteModel%Phi = Phi
        DiscreteModel%Phir = Phir
        DiscreteModel%DilAngle = DilAngle
      END IF

      
      !UNDRAINED BEHAVIOUR.
      ALLOCATE(DiscreteModel%Undrained(DiscreteModel%NumElements),cw(DiscreteModel%NumElements))
      cw = 0.0
      DiscreteModel%Undrained = .FALSE.                
      elementNum = 0
      DO itheta=1, WellModel%NTHETA
        DO iz=1, WellModel%NZ
          DO irad=1, WellModel%NRAD
            elementNum = elementNum + 1
            IF( ANY(WellData%ListOfUndrainedLayers(1:)==iz) )THEN
              DiscreteModel%Undrained(elementNum) = .TRUE.
              cw(elementNum) = WellData%UndrainedLayerParameters(1,1)
            END IF
          END DO
        END DO
      END DO
      CALL UpdatePoreFluidCompressibility(DiscreteModel,cw)
      
      !CREEP
      !  TYPE CreepModel
      !     LOGICAL:: Creep
      !     REAL(8):: mu0
      !     REAL(8):: Qc
      !     REAL(8):: N2
      IF( SIZE(Problem%CreepTable(1,:))>0 )THEN
          DiscreteModel%Creep = .TRUE.
          
          ALLOCATE( DiscreteModel%CreepData( SIZE(Problem%CreepTable(1,:))),&
              DiscreteModel%ElementCreep(DiscreteModel%NumElements) )
          DiscreteModel%ElementCreep = 0
          
          counter = 0
          DO i=1, SIZE(Problem%CreepTable(1,:))
              READ(Problem%CreepTable(1,i),*) I1
              READ(Problem%CreepTable(2,i),*) I2
              READ(Problem%CreepTable(3,i),*) J1
              READ(Problem%CreepTable(4,i),*) J2
              READ(Problem%CreepTable(5,i),*) K1
              READ(Problem%CreepTable(6,i),*) K2
              READ(Problem%CreepTable(7,i),*) DiscreteModel%CreepData(i)%mu0
              READ(Problem%CreepTable(8,i),*) DiscreteModel%CreepData(i)%Qc
              READ(Problem%CreepTable(9,i),*) DiscreteModel%CreepData(i)%N2
              
              elementNum = 0
              DO itheta=1, WellModel%NTHETA
                  DO iz=1, WellModel%NZ
                      DO irad=1, WellModel%NRAD
                          elementNum = elementNum + 1
                              
                          IF( (itheta>=J1) .AND. (itheta<=J2) .AND. &
                              (iz>=K1) .AND. (iz<=K2) .AND. &
                              (irad>=I1) .AND. (irad<=I2) )THEN
                              counter = counter + 1
                              DiscreteModel%ElementCreep(elementNum) = i
                          END IF
                              
                      END DO
                  END DO
              END DO
          END DO
          
      END IF
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !PERMEABILITY ENHANCEMENT MODEL.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ALLOCATE( DiscreteModel%FracturedElementModel(DiscreteModel%NumElements) )
      !ALLOCATE( DiscreteModel%IsFractured(nip,DiscreteModel%NumElements) )
      ALLOCATE( DiscreteModel%FracturedElementNumber(nip,DiscreteModel%NumElements) )
      DiscreteModel%FracturedElementModel = .FALSE.
      !DiscreteModel%IsFractured = .FALSE.
      DiscreteModel%FracturedElementNumber = 0
      
      !DO i=1,SIZE(WellData%FractureInterval(1,:))
      !  
      !  elementNum = 0
      !  DO itheta=1, WellModel%NTHETA
      !    DO iz=1, WellModel%NZ
      !      DO irad=1, WellModel%NRAD
      !        elementNum = elementNum + 1
      !        
      !        IF( (irad>=WellData%FractureInterval(1,i)).AND.(irad<=WellData%FractureInterval(2,i)).AND. &
      !            (iz>=WellData%FractureInterval(3,i)).AND.(iz<=WellData%FractureInterval(4,i)).AND. &
      !            (itheta>=WellData%FractureInterval(5,i)).AND.(itheta<=WellData%FractureInterval(6,i)) )THEN
      !          DiscreteModel%FracturedElementModel(elementNum)=.TRUE.
      !          DiscreteModel%TractionStrength(elementNum) = WellData%TractionStrengthRatio(i)*DiscreteModel%Cohesion(elementNum)
      !        END IF
      !            
      !      END DO
      !    END DO
      !  END DO
      !  
      !END DO
      
      !***************************************************************************
      !DAMAGE.
      !***************************************************************************
      !TYPE (ShearDamage):: ShearStrainDamage
      !Sempre alocado (tamanho 0 quando a opção não é usada), pois é
      !percorrido com SIZE() em outros módulos.
      ALLOCATE(DiscreteModel%ShearStrainDamage(SIZE(Problem%ShearStrainDamage)))
      IF(SIZE(Problem%ShearStrainDamage(:))>0) THEN
        n = SIZE(Problem%ShearStrainDamage)
        DO i=1,n
          
          DiscreteModel%ShearStrainDamage(i)%LimitPlasticDeviatoricStrain = &
                          Problem%ShearStrainDamage(i)%LimitPlasticDeviatoricStrain
          DiscreteModel%ShearStrainDamage(i)%ResidualDamage =  &
                          Problem%ShearStrainDamage(i)%ResidualDamage
          DiscreteModel%ShearStrainDamage(i)%ResidualStiffnessDamage =  &
                          Problem%ShearStrainDamage(i)%ResidualDamage
          
          ALLOCATE( DiscreteModel%ShearStrainDamage(i)%ElementList(0:DiscreteModel%NumElements))
          DiscreteModel%ShearStrainDamage(i)%ElementList = 0        
          j = 0
          elementNum = 0
          DO itheta=1, WellModel%NTHETA
            DO iz=1, WellModel%NZ
              DO irad=1, WellModel%NRAD
                elementNum = elementNum + 1
                IF(iz>=Problem%ShearStrainDamage(i)%K1.AND.IZ<=Problem%ShearStrainDamage(i)%K2)THEN
                  j = j + 1
                  DiscreteModel%ShearStrainDamage(i)%ElementList(j) = elementNum
                END IF
              END DO
            END DO
          END DO
          DiscreteModel%ShearStrainDamage(i)%ElementList(0) = j
        
        END DO
        
      END IF
      
      !TYPE(StiffnessDamage):: Liquefaction
      !Sempre alocado (tamanho 0 quando a opção não é usada), pois é
      !percorrido com SIZE() em outros módulos.
      ALLOCATE(DiscreteModel%Liquefaction(SIZE(Problem%Liquefaction)))
      IF(SIZE(Problem%Liquefaction(:))>0) THEN
        n = SIZE(Problem%Liquefaction)
        
        DO i=1,n
          
          DiscreteModel%Liquefaction(i)%MaxMeanEffectiveStress = Problem%Liquefaction(i)%MaxMeanEffectiveStress
          DiscreteModel%Liquefaction(i)%ShearDamage = Problem%Liquefaction(i)%ShearDamage
          DiscreteModel%Liquefaction(i)%ElasticDamage = Problem%Liquefaction(i)%ElasticDamage
          
          ALLOCATE( DiscreteModel%Liquefaction(i)%ElementList(0:DiscreteModel%NumElements))
          DiscreteModel%Liquefaction(i)%ElementList = 0
          j = 0
          elementNum = 0
          DO itheta=1, WellModel%NTHETA
            DO iz=1, WellModel%NZ
              DO irad=1, WellModel%NRAD
                elementNum = elementNum + 1
                IF(iz>=Problem%Liquefaction(i)%K1.AND.IZ<=Problem%Liquefaction(i)%K2)THEN
                  j = j + 1
                  DiscreteModel%Liquefaction(i)%ElementList(j) = elementNum
                END IF
              END DO
            END DO
          END DO
          DiscreteModel%Liquefaction(i)%ElementList(0) = j
          
        END DO
      
      END IF

      !TYPE (PermModif):: PermMultModel
      !Sempre alocado (tamanho 0 quando a opção não é usada), pois é
      !percorrido com SIZE() em outros módulos.
      ALLOCATE(DiscreteModel%PermModel(SIZE(Problem%Frac)))
      IF(SIZE(Problem%Frac(:))>0) THEN
        n = SIZE(Problem%Frac)
        
        DO i=1,n
          
          DiscreteModel%PermModel(i)%PermIncrementByPlasticDilation = Problem%Frac(i)%PermIncrementByPlasticDilation
          DiscreteModel%PermModel(i)%PermIncrementByPlasticShear = Problem%Frac(i)%PermIncrementByPlasticShear
          
          ALLOCATE( DiscreteModel%PermModel(i)%ElementList(0:DiscreteModel%NumElements))
          DiscreteModel%PermModel(i)%ElementList = 0
          j = 0
          elementNum = 0
          DO itheta=1, WellModel%NTHETA
            DO iz=1, WellModel%NZ
              DO irad=1, WellModel%NRAD
                elementNum = elementNum + 1
                IF(iz>=Problem%Frac(i)%K1.AND.IZ<=Problem%Frac(i)%K2)THEN
                  j = j + 1
                  DiscreteModel%PermModel(i)%ElementList(j) = elementNum
                END IF
              END DO
            END DO
          END DO
          DiscreteModel%PermModel(i)%ElementList(0) = j
        END DO
        
      END IF
      
      !TYPE (TensionDamage):: TractionStrengthModel
      !IF(SIZE(Problem%TracDamage(:))>0) THEN
      !  n = SIZE(Problem%TracDamage)
      !  ALLOCATE(DiscreteModel%TractionStrengthModel(n))
      !  
      !  DO i=1,n
      !    
      !    DiscreteModel%TractionStrengthModel(i)%T0 = Problem%TracDamage(i)%T0
      !    DiscreteModel%TractionStrengthModel(i)%MaxPlasticDilation = Problem%TracDamage(i)%MaxPlasticDilation
      !    !DiscreteModel%TractionStrengthModel(i)%FracDirection = Problem%TracDamage(i)%FracDirection
      !    
      !    elementNum = 0
      !    DO itheta=1, WellModel%NTHETA
      !      DO iz=1, WellModel%NZ
      !        DO irad=1, WellModel%NRAD
      !          elementNum = elementNum + 1
      !          IF( iz>=Problem%TracDamage(i)%K1 .AND. IZ<=Problem%TracDamage(i)%K2 .AND. itheta==1) THEN
      !            DiscreteModel%ElementTractionStrengthModel(elementNum) = i
      !            DiscreteModel%FracturedElementModel(elementNum) = .TRUE.
      !            DiscreteModel%Gc(elementNum) = Problem%TracDamage(i)%Gc
      !            DiscreteModel%TRACflag(elementNum) = .TRUE.
      !            DiscreteModel%TRACmodelNum(elementNum) = I
      !          END IF
      !        END DO
      !      END DO
      !    END DO
      !  END DO
      !  
      !END IF
      
      !***************************************************************************
      !THERMAL MODEL.
      !***************************************************************************
      !REAL(8),ALLOCATABLE:: ThermalExpCoeff(:)
      !REAL(8),ALLOCATABLE:: ThermalConductivity(:)
      IF(SIZE(Problem%ThermalProps(1,:))>0) THEN
        
        n = SIZE(Problem%ThermalProps(1,:))
        
        DO i=1,n
          
          elementNum = 0
          DO itheta=1, WellModel%NTHETA
            DO iz=1, WellModel%NZ
              DO irad=1, WellModel%NRAD
                elementNum = elementNum + 1
                
                IF(iz>=Problem%ThermalPropsKinterval(1,i).AND.iz<=Problem%ThermalPropsKinterval(2,i))THEN
                  DiscreteModel%ThermalExpCoeff(elementNum) = Problem%ThermalProps(1,i)
                  DiscreteModel%ThermalConductivity(elementNum) = Problem%ThermalProps(2,i)
                  DiscreteModel%RockSpecificHeatCapacity(elementNum) = Problem%ThermalProps(3,i)
                  DiscreteModel%FluidSpecificHeatCapacity(elementNum) = Problem%ThermalProps(4,i)
                  DiscreteModel%RockFluidHeatTransferCoef(elementNum) = Problem%ThermalProps(5,i)
                END IF
                
              END DO
            END DO
          END DO
          
        END DO
        
      END IF
      
      !***************************************************************************
      !OTHER FUNCTIONALITIES.
      !***************************************************************************
      DiscreteModel%UpdateGeometry = Problem%UpdateGeometry
      
      
    RETURN
    END SUBROUTINE


    !DEFINITE O TIPO DE ELEMENTO FINITO DA MALHA.
    !DiscreteModel%element="hexahedron"            !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
    !DiscreteModel%NumDimensions = 3               !ndim - número de dimensões do problema.
    !DiscreteModel%NumNodesPerElement = 8          !nod - número de nós de deslocamento por elemento.
    !DiscreteModel%NumDOFUperElement = 24          !ndof - número de graus de liberdade de deslocamento por elemento.
    !DiscreteModel%NumDOFPperElement = 8           !nodf - número de nós com grau de liberdade de poro-pressão (mas será usado como constante).
    !DiscreteModel%NumDOFUperNode = 3              !nodof - número de graus de liberdade de deslocamento por nó.
    !DiscreteModel%NumTensorComponents = 6
    SUBROUTINE SetFiniteElementType( DiscreteModel )

        IMPLICIT NONE
        TYPE(Model),INTENT(IN OUT) :: DiscreteModel

        DiscreteModel%element="hexahedron"            !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
        DiscreteModel%NumDimensions = 3               !ndim - número de dimensões do problema.
        DiscreteModel%NumNodesPerElement = 8          !nod - número de nós de deslocamento por elemento.
        DiscreteModel%NumDOFUperElement = 24          !ndof - número de graus de liberdade de deslocamento por elemento.
        DiscreteModel%NumDOFPperElement = 8           !nodf - número de nós com grau de liberdade de poro-pressão (mas será usado como constante).
        DiscreteModel%NumDOFUperNode = 3              !nodof - número de graus de liberdade de deslocamento por nó.
        DiscreteModel%NumTensorComponents = 6         !nst - número de componentes no tensor de tensões

        RETURN
    END SUBROUTINE

    
    SUBROUTINE BuildDOFmatrix(nels,nod,ndof,g_num,nf,g_g)
      !Monta matriz de graus de liberdade dos nós dos elementos.

      USE NumericalUtilities

      !variaveis.
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

      INTEGER,ALLOCATABLE,INTENT(IN)::g_num(:,:),&  !conectividade dos elementos da malha
                                      nf(:,:)       !matriz de graus de liberdade
      INTEGER,ALLOCATABLE,INTENT(IN OUT)::g_g(:,:)  !matriz com os graus de liberdade por elemento
      INTEGER,INTENT(IN)::nels,&                    !número de elementos
                          nod,&                     !número de nós de deslocamento por elemento
                          ndof                      !número de graus de liberdade de deslocamento por elemento

      INTEGER::iel                                  !contador
      INTEGER,ALLOCATABLE:: num(:),&                !conectividade de um elemento
                            g(:)                    !graus de liberdade de 1 elemento

      ALLOCATE(g(ndof),num(nod))

      DO iel=1,nels
        num = g_num(:,iel)
        CALL num_to_g(num,nf,g)
        g_g(:,iel) = g
      END DO

      RETURN
    END SUBROUTINE

    
    SUBROUTINE UpdateGeomechanicalModel(DiscreteModel,PressureBC_flag,&
                                        Porepressure,Temperature,WaterSaturation,RAD)
      USE FEMlib
      USE NumericalLib
      USE SolidMechanicsLib
      USE MainConstants
      USE StressLoad

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      LOGICAL,INTENT(IN):: PressureBC_flag(:)
      REAL(8),INTENT(IN):: Porepressure(:)               !element porepressure.
      REAL(8),INTENT(IN):: Temperature(:)                !element temperature.
      REAL(8),INTENT(IN):: WaterSaturation(:)            !element water saturation.
      LOGICAL,INTENT(IN):: RAD

      INTEGER:: i,ip,nip,nels,nn

      INTEGER,ALLOCATABLE:: num(:),g(:)
      REAL(8),ALLOCATABLE:: coord(:,:),g_coord(:,:)
      REAL(8):: n(3),coordsTRI(3,3),area1,area2,normal1(3),normal2(3),theta
      REAL(8):: stress(6),sxx,syy,sxy,normal,shear
        
      REAL(8),ALLOCATABLE:: points(:,:),&                 !coordenadas dos pontos de integração
                            weights(:)                    !peso dos pontos de integração
      REAL(8),ALLOCATABLE:: sigma(:,:)

      element = DiscreteModel%element
      nels = DiscreteModel%NumElements
      nn = DiscreteModel%NumNodes
      ndim = DiscreteModel%NumDimensions
      nst = DiscreteModel%NumTensorComponents
      nip = DiscreteModel%NumIntegrationPoints
      nod = DiscreteModel%NumNodesPerElement
      ndof = DiscreteModel%NumDOFUperElement
        
      !MEMORY ALLOCATION.
      ALLOCATE( num(nod),g(ndof),points(nip,ndim),weights(nip),coord(nod,ndim),g_coord(ndim,nn),&
                sigma(nst,nip))
        
      !NODAL COORDINATES.
      g_coord=0.0
      g_coord(1,:) = DiscreteModel%X
      g_coord(2,:) = DiscreteModel%Y
      g_coord(3,:) = DiscreteModel%Z

      ! NODAL TEMPERATURE, POREPRESSURE AND WATER SATURATION.
      DiscreteModel%Temperature = Temperature
      DO i=1,DiscreteModel%NumNodes
        DiscreteModel%Porepressure(i) = Porepressure(i)
        DiscreteModel%WaterSaturation(i) = WaterSaturation(i)
      END DO

      ! AVERAGE ELEMENT POREPRESSURE.
      DO i=1,DiscreteModel%NumElements                   
          IF(DiscreteModel%ActiveToFlowElements(i))THEN
              num = DiscreteModel%ElementNodes(:,i)
              DiscreteModel%ElementAVGporepressure(i) = SUM(Porepressure(num))/8.0
          END IF
      END DO
        
      ! NODAL PORE STRAIN.
      ! ...
      

    RETURN
    END SUBROUTINE

    !ENERGY CALCULATIONS.    
    SUBROUTINE UpdateStrainWork(DiscreteModel,Stress,DisplacementIncrement)
    
      USE NumericalUtilities
      USE FEMlib
      USE ConstitutiveGeoModels
    
      IMPLICIT NONE
    
      ! INPUT.
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      REAL(8):: Stress(:,:,:)
      REAL(8),INTENT(IN):: DisplacementIncrement(0:)
      
      INTEGER::i,iel  !contadores
      REAL(8)::det  !determinante do jacobiano
    
      INTEGER:: nels,&  !número de elementos.
                nn,&    !número de nós.
                ndim,&  !número de dimensões do problema.
                nst,&   !número de componentes no tensor de tensões.
                nip,&   !número de pontos de integração.
                nod,&   !número de nós por elemento.
                ndof,&  !número de graus de liberdade por elemento.
                neq
      REAL(8),ALLOCATABLE:: bee(:,:),&  !matriz de derivadas das funções de forma
                            jac(:,:),&  !matriz jacobiana das funções de forma
                            fun(:),&    !funções de forma
                            der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                            deriv(:,:)  !derivadas das funções de forma dos deslocamentos em coordenadas globais
      REAL(8),ALLOCATABLE:: eld(:),eps(:)
      REAL(8),ALLOCATABLE:: g_coord(:,:),&  !coordenadas dos nós da malha.
                            coord(:,:),&    !coordenada dos nós de um elemento.
                            points(:,:),&   !coordenadas dos pontos de integração.
                            weights(:)      !peso dos pontos de integração.
      REAL(8),ALLOCATABLE:: st(:),&         !delta de tensão 
                            st_elastic(:),& !delta de tensão elástica
                            dee(:,:)        !stress-strain matrix
      INTEGER,ALLOCATABLE:: g(:),&  !vetor com graus de liberdade de um elemento.
                            num(:)  !conectividade de um elemento.
      
      REAL(8):: biot,pp
      
      REAL(8):: ws_pl
      REAL(8):: e,v

      CHARACTER(LEN=320)::element  !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
      
      element = DiscreteModel%element
      nels = DiscreteModel%NumElements
      nn = DiscreteModel%NumNodes
      ndim = DiscreteModel%NumDimensions
      nst = DiscreteModel%NumTensorComponents
      nip = DiscreteModel%NumIntegrationPoints
      nod = DiscreteModel%NumNodesPerElement
      ndof = DiscreteModel%NumDOFUperElement
      
      neq = DiscreteModel%NumEquations
    
      ALLOCATE(bee(nst,ndof),jac(ndim,ndim),fun(nod),der(ndim,nod),deriv(ndim,nod),&
               g_coord(ndim,nn),st(nst),st_elastic(nst),points(nip,ndim),weights(nip),g(ndof),&
               num(nod),coord(nod,ndim),eld(ndof),eps(nst),dee(nst,nst))
      
      g_coord=0.0
      g_coord(1,:) = DiscreteModel%X
      g_coord(2,:) = DiscreteModel%Y
      g_coord(3,:) = DiscreteModel%Z
      
      CALL sample(element,points,weights)
      
      DO iel=1,nels
 
        num = DiscreteModel%ElementNodes(:,iel)
        coord = TRANSPOSE(g_coord(:,num))
        g = DiscreteModel%ElementDOFs(:,iel)
        eld = DisplacementIncrement(g)
        
        biot = DiscreteModel%BiotCoef(iel)
        pp = DiscreteModel%ElementAVGporepressure(iel)
          
        !e = DiscreteModel%YoungMod(iel)
        !v = DiscreteModel%Poisson(iel)
        !CALL D_elastic_homogenious(dee,e,v)
        dee = D_TransverseIsotropic(DiscreteModel%TransverIsotropicProps(iel)%E1,&
            DiscreteModel%TransverIsotropicProps(iel)%E2,&
            DiscreteModel%TransverIsotropicProps(iel)%n1,&
            DiscreteModel%TransverIsotropicProps(iel)%n2,&
            DiscreteModel%TransverIsotropicProps(iel)%strike,&
            DiscreteModel%TransverIsotropicProps(iel)%dip)
          
        ws_pl = 0.0
        
        DO i=1,nip
          
          !FINITE ELEMENT STUFF------------------------------------------------
          CALL shape_der(der,points,i)  !derivada das funções de forma
          jac = MATMUL(der,coord)       !jacobiano do elemento no ponto de integração i
          CALL determinante(jac,det)    !determinante do jacobiano
          CALL invert(jac)              !inverso do jacobiano
          deriv = MATMUL(jac,der)       !derivadas das funções de forma em coordenadas reais
          CALL beemat(bee,deriv)        !matriz de derivadas B
          eps = MATMUL(bee,eld)         !strain
          st = Stress(:,i,iel)          !- biot * pp  for total stress.
          
          ! STRAIN WORK (TOTAL: ELASTIC + PLASTIC)-----------------------------
          DiscreteModel%StrainWorkDensity(i,iel) = DiscreteModel%StrainWorkDensity(i,iel) + &
                                                                      0.5 * DOT_PRODUCT(st,eps)
          
          ! PLASTIC WORK-------------------------------------------------------
          st_elastic = MATMUL(dee,eps)
          ws_pl = ws_pl + 0.5 * ( DOT_PRODUCT(st_elastic,eps) - DOT_PRODUCT(st,eps) ) * det * weights(i)
                      
          ! SYY*EYY J-INTEGRAL-------------------------------------------------
          !IF( DiscreteModel%ConstitutiveGeoModel(iel).NE.1 )THEN      !1- Damage (EXPLICITY) - calculated internally.
            !IF( st(2)>0.0 )&
              !st(2) = st(2) - DiscreteModel%BiotCoef(iel) * DiscreteModel%ElementAVGporepressure(iel)
              DiscreteModel%YYstrainDesityEnergy(iel) = DiscreteModel%YYstrainDesityEnergy(iel) + st(2) * eps(2)
                                                                                    !0.5 * st(2) * eps(2)
          !END IF
                    
        END DO
          
        DiscreteModel%PlasticStrainWork(iel) = DiscreteModel%PlasticStrainWork(iel) + ws_pl
        
      END DO

    RETURN
    END SUBROUTINE
                                        
    
    SUBROUTINE UpdatePoreFluidCompressibility(DiscreteModel,Cwater)

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      REAL(8),INTENT(IN):: Cwater(:)                 !water compressibility.
    
      DiscreteModel%PoreFluidBulkMod(:) = 1.0 / Cwater(:)
    
    RETURN
    END SUBROUTINE
     
    
    SUBROUTINE InternalEffectiveStressForce(DiscreteModel,Force)
      !       fumbal(# nós) = int( transpose(B)* sigma' * dv
    
      USE NumericalUtilities
      USE FEMlib
    
      IMPLICIT NONE
    
      ! INPUT.
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      REAL(8),INTENT(OUT):: Force(0:)  !matriz de rigidez do elemento
      
      INTEGER::i,iel  !contadores
      REAL(8)::det  !determinante do jacobiano
    
      INTEGER:: nels,&  !número de elementos.
                nn,&    !número de nós.
                ndim,&  !número de dimensões do problema.
                nst,&   !número de componentes no tensor de tensões.
                nip,&   !número de pontos de integração.
                nod,&   !número de nós por elemento.
                ndof,&  !número de graus de liberdade por elemento.
                neq
      REAL(8),ALLOCATABLE:: bee(:,:),&  !matriz de derivadas das funções de forma
                            jac(:,:),&  !matriz jacobiana das funções de forma
                            fun(:),&    !funções de forma
                            der(:,:),&  !derivadas das funções de forma dos deslocamentos em coordenadas naturais
                            deriv(:,:)  !derivadas das funções de forma dos deslocamentos em coordenadas globais
      REAL(8),ALLOCATABLE:: g_coord(:,:),&  !coordenadas dos nós da malha.
                            coord(:,:),&    !coordenada dos nós de um elemento.
                            points(:,:),&   !coordenadas dos pontos de integração.
                            weights(:)      !peso dos pontos de integração.
      REAL(8),ALLOCATABLE:: stress(:)       !delta de tensão 
      INTEGER,ALLOCATABLE:: g(:),&  !vetor com graus de liberdade de um elemento.
                            num(:)  !conectividade de um elemento.

      CHARACTER(LEN=320)::element  !tipo de elemento finito ('hexahedron' ou 'tetrahedron').
      
      element = DiscreteModel%element
      nels = DiscreteModel%NumElements
      nn = DiscreteModel%NumNodes
      ndim = DiscreteModel%NumDimensions
      nst = DiscreteModel%NumTensorComponents
      nip = DiscreteModel%NumIntegrationPoints
      nod = DiscreteModel%NumNodesPerElement
      ndof = DiscreteModel%NumDOFUperElement
      
      neq = DiscreteModel%NumEquations
    
      ALLOCATE(bee(nst,ndof),jac(ndim,ndim),fun(nod),der(ndim,nod),deriv(ndim,nod),&
               g_coord(ndim,nn),stress(nst),points(nip,ndim),weights(nip),g(ndof),&
               num(nod),coord(nod,ndim))
      fun=0.0
      
      g_coord=0.0
      g_coord(1,:) = DiscreteModel%X
      g_coord(2,:) = DiscreteModel%Y
      g_coord(3,:) = DiscreteModel%Z
    
      !INTEGRAÇÃO NUMÉRICA DA TENSÃO DESBALANCEADA
      !     int( transpose(B)* (sigma*-sigma) * dv ) = 
      !                           = SOMA(pontos de Gauss, transpose(B) * (sigma*-sigma) * detJacobiano * pesoPtoGauss
      !Force = 0.0
      
      CALL sample(element,points,weights)
      
      DO iel=1,nels
        IF( DiscreteModel%ActiveElements(iel) .EQV. .TRUE. ) THEN       
          num = DiscreteModel%ElementNodes(:,iel)
          coord = TRANSPOSE(g_coord(:,num))
          g = DiscreteModel%ElementDOFs(:,iel)  !=g_g(:,iel)
          
          DO i=1,nip
            CALL shape_fun(fun,points,i)  !funções de forma no ponto de integração i
            CALL shape_der(der,points,i)  !derivada das funções de forma
            jac = MATMUL(der,coord)       !jacobiano do elemento no ponto de integração i
            CALL determinante(jac,det)    !determinante do jacobiano
            CALL invert(jac)              !inverso do jacobiano
            deriv = MATMUL(jac,der)       !derivadas das funções de forma em coordenadas reais
            CALL beemat(bee,deriv)        !matriz de derivadas B

            stress = DiscreteModel%EffectiveStressTensor(:,i,iel)
          
            Force(g) = Force(g) + MATMUL( TRANSPOSE( bee ) , stress ) * det * weights(i)
          
          END DO
          
        END IF
      END DO
      
      Force = -Force

    RETURN
    END SUBROUTINE

    
    SUBROUTINE UpdateCompressibility(DiscreteModel,ElementVolume,&
                                      DeltaVolumetricStrain,DeltaPorepressure,&
                                      NodalPoreCompressbility)
    
      IMPLICIT NONE
    
      ! INPUT.
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      REAL(8),INTENT(IN):: ElementVolume(:),DeltaVolumetricStrain(:),DeltaPorepressure(:)
      REAL(8),INTENT(IN OUT):: NodalPoreCompressbility(:)
      
      INTEGER:: i,num(8),nod
      REAL(8):: biot,por,evol,vol,e,poisson,deltaP,minDeltaP,tol=1.0E-2
      
      REAL(8),ALLOCATABLE:: sumVector1(:),sumVector2(:)
      
      CHARACTER(LEN=30):: method='StressCoupling'  !StrainCoupling

      ALLOCATE(sumVector1(DiscreteModel%NumNodes),sumVector2(DiscreteModel%NumNodes))
      sumVector1 = 0.0
      sumVector2 = 0.0
      
      nod = DiscreteModel%NumNodesPerElement
    
      !Element compressiblity.
      DO i=1,DiscreteModel%NumElements
        num = DiscreteModel%ElementNodes(:,i)
        vol = 0.0
        IF( (DiscreteModel%POR(i)>0.0).AND.&
            (DiscreteModel%KH(i)>0.0).AND.(DiscreteModel%KV(i)>0.0))THEN

          IF(TRIM(method)=='StressCoupling')THEN
            biot = DiscreteModel%BiotCoef(i)
            por = DiscreteModel%POR(i)
            e = 0.5 * (DiscreteModel%TransverIsotropicProps(i)%E1 + DiscreteModel%TransverIsotropicProps(i)%E2)
            poisson = 0.5 * (DiscreteModel%TransverIsotropicProps(i)%n1 + DiscreteModel%TransverIsotropicProps(i)%n2)
 
            IF(e>0.0)THEN
              DiscreteModel%PoreCompressibility(i) = biot**2/por/e*(1+poisson)*(1-2*poisson)/(1-poisson)
            ELSE
              DiscreteModel%PoreCompressibility(i) = 0.0
            END IF

          ELSEIF(TRIM(method)=='StrainCoupling')THEN
            deltaP = SUM(DeltaPorepressure(num)) / nod
            minDeltaP = tol * MAXVAL(ABS(DeltaPorepressure))
            IF(ABS(deltaP)>minDeltaP)THEN
              biot = DiscreteModel%BiotCoef(i)
              por = DiscreteModel%POR(i)
              evol = DeltaVolumetricStrain(i)
              DiscreteModel%PoreCompressibility(i) = biot/por*evol/deltaP
            END IF

          ELSE
          END IF
          vol = ElementVolume(i)
        ELSE
          DiscreteModel%PoreCompressibility(i) = 0.0
        END IF

        sumVector1(num) = sumVector1(num) + DiscreteModel%PoreCompressibility(i)*vol/8
        sumVector2(num) = sumVector2(num) + vol/8
        
      END DO
      
      !Average nodal pore compressibility.
      NodalPoreCompressbility = 0.0
      DO i=1,DiscreteModel%NumNodes
        IF(sumVector2(i)>0.0) NodalPoreCompressbility(i) = sumVector1(i) / sumVector2(i)
      END DO
    
    RETURN
    END SUBROUTINE
    
                                      
    SUBROUTINE CalculateVolumetricStrain(DiscreteModel,VolumetricStrain)
    
      IMPLICIT NONE
    
      ! INPUT.
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      REAL(8),INTENT(IN OUT):: VolumetricStrain(:)  !Element volumetric strain.
      
      INTEGER:: i,ip,nip
      REAL(8):: evol
      
      nip = DiscreteModel%NumIntegrationPoints
    
      VolumetricStrain = 0.0
      
      DO i=1,DiscreteModel%NumElements
        DO ip=1,nip
          evol = SUM(DiscreteModel%StrainTensor(1:3,ip,i))
          VolumetricStrain(i) = VolumetricStrain(i) + evol
          DiscreteModel%PoreStrain = VolumetricStrain(i) * DiscreteModel%BiotCoef(i)  !+
        END DO
        VolumetricStrain(i) = VolumetricStrain(i) / nip
      END DO
            
    RETURN
    END SUBROUTINE


    SUBROUTINE AddNewFracturedRock(DiscreteModel,FracturedRockElements,ElementNumber,GaussPoint,&
                                    EigenVectors,EigenS3postion,Damage)
    
      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      TYPE(FracturedRock),ALLOCATABLE,INTENT(IN OUT):: FracturedRockElements(:)
      INTEGER,INTENT(IN):: ElementNumber
      INTEGER,INTENT(IN):: GaussPoint
      REAL(4),INTENT(IN):: EigenVectors(3,3)
      INTEGER:: EigenS3postion
      REAL(4),INTENT(IN):: Damage
    
      TYPE(FracturedRock),ALLOCATABLE:: temp(:)
      
      INTEGER:: OLDsize
      
      ALLOCATE(temp(OLDsize))
      temp = FracturedRockElements
      DEALLOCATE(FracturedRockElements)
      ALLOCATE(FracturedRockElements(OLDsize+1))
      FracturedRockElements(1:OLDsize) = temp
      
      !New position.
      FracturedRockElements(OLDsize+1)%ElementNumber = ElementNumber
      FracturedRockElements(OLDsize+1)%GaussPoint = GaussPoint
      FracturedRockElements(OLDsize+1)%EigenVectors = EigenVectors
      FracturedRockElements(OLDsize+1)%EigenS3postion = EigenS3postion
      FracturedRockElements(OLDsize+1)%Damage = Damage
      
      DiscreteModel%FracturedElementNumber(GaussPoint,ElementNumber) = OLDsize+1
      
    RETURN
  END SUBROUTINE
      
  
    SUBROUTINE InterpolateSandModelTable(DiscreteModel,eps,Cp,Cr,Phi,Phir,DilAngle)
  
    
      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      REAL(8),INTENT(IN):: eps
      REAL(8),INTENT(OUT):: Cp
      REAL(8),INTENT(OUT):: Cr
      REAL(8),INTENT(OUT):: Phi
      REAL(8),INTENT(OUT):: Phir
      REAL(8),INTENT(OUT):: DilAngle
    
      INTEGER:: i,n
      REAL(8):: factor
    
      n = 7
    
      DO i=2,n
      
        IF( ( eps >= DiscreteModel%ShearHardeningPar(i-1) ).AND.(eps < DiscreteModel%ShearHardeningPar(i)) )THEN
        
          factor = (eps-DiscreteModel%ShearHardeningPar(i-1)) /  &
                          ( DiscreteModel%ShearHardeningPar(i)-DiscreteModel%ShearHardeningPar(i-1) )
          Cp = DiscreteModel%Cp(i-1) + factor * (DiscreteModel%Cp(i)-DiscreteModel%Cp(i-1))
          Cr = DiscreteModel%Cr(i-1) + factor * (DiscreteModel%Cr(i)-DiscreteModel%Cr(i-1))
          Phi = DiscreteModel%Phi(i-1) + factor * (DiscreteModel%Phi(i)-DiscreteModel%Phi(i-1))
          Phir = DiscreteModel%Phir(i-1) + factor * (DiscreteModel%Phir(i)-DiscreteModel%Phir(i-1))
          DilAngle = DiscreteModel%DilAngle(i-1) + factor * (DiscreteModel%DilAngle(i)-DiscreteModel%DilAngle(i-1))
        
          EXIT
        
        END IF
      
      END DO
    
  
    RETURN
    END SUBROUTINE
  
    
    SUBROUTINE SaveVariables(DiscreteModel,Variables)
    
      IMPLICIT NONE
      TYPE(Model),INTENT(IN):: DiscreteModel
      TYPE(VariablesContainer),INTENT(IN OUT):: Variables
    
      Variables%Displacement = DiscreteModel%Displacement
      Variables%EffectiveStressTensor = DiscreteModel%EffectiveStressTensor
      Variables%StrainTensor = DiscreteModel%StrainTensor
      Variables%PlasticStrainTensor = DiscreteModel%PlasticStrainTensor
      Variables%DamageShear = DiscreteModel%DamageShear
      Variables%DamageTraction = DiscreteModel%DamageTraction
      Variables%DamageStiffness = DiscreteModel%DamageStiffness
      !Variables%IsFractured = DiscreteModel%IsFractured
      Variables%Porepressure = DiscreteModel%Porepressure
      Variables%ElementAVGporepressure = DiscreteModel%ElementAVGporepressure
      Variables%Fracpressure = DiscreteModel%Fracpressure
      Variables%Temperature = DiscreteModel%Temperature
      Variables%WaterSaturation = DiscreteModel%WaterSaturation
      Variables%PoreFluidBulkMod = DiscreteModel%PoreFluidBulkMod
      Variables%PoreStrain = DiscreteModel%PoreStrain
      Variables%PoreCompressibility = DiscreteModel%PoreCompressibility
      Variables%YieldFunc = DiscreteModel%YieldFunc
      Variables%YieldEnvelopeTOL = DiscreteModel%YieldEnvelopeTOL
      Variables%PlasticStrainWork = DiscreteModel%PlasticStrainWork
      Variables%StrainWorkDensity = DiscreteModel%StrainWorkDensity
      Variables%YYstrainDesityEnergy = DiscreteModel%YYstrainDesityEnergy
      Variables%Jintegral = DiscreteModel%Jintegral
      Variables%JintegralCum = DiscreteModel%JintegralCum
      Variables%Gamma = DiscreteModel%Gamma
      
      RETURN
    END SUBROUTINE
    
    
    SUBROUTINE RecoverVariables(DiscreteModel,Variables)
    
      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      TYPE(VariablesContainer),INTENT(IN):: Variables
    
      DiscreteModel%Displacement = Variables%Displacement
      DiscreteModel%EffectiveStressTensor = Variables%EffectiveStressTensor
      DiscreteModel%StrainTensor = Variables%StrainTensor
      DiscreteModel%PlasticStrainTensor = Variables%PlasticStrainTensor
      DiscreteModel%DamageShear = Variables%DamageShear
      DiscreteModel%DamageTraction = Variables%DamageTraction
      DiscreteModel%DamageStiffness = Variables%DamageStiffness
      !DiscreteModel%IsFractured = Variables%IsFractured
      DiscreteModel%Porepressure = Variables%Porepressure
      DiscreteModel%ElementAVGporepressure = Variables%ElementAVGporepressure
      DiscreteModel%Fracpressure = Variables%Fracpressure
      DiscreteModel%Temperature = Variables%Temperature
      DiscreteModel%WaterSaturation = Variables%WaterSaturation
      DiscreteModel%PoreFluidBulkMod = Variables%PoreFluidBulkMod
      DiscreteModel%PoreStrain = Variables%PoreStrain
      DiscreteModel%PoreCompressibility = Variables%PoreCompressibility
      DiscreteModel%YieldFunc = Variables%YieldFunc
      DiscreteModel%YieldEnvelopeTOL = Variables%YieldEnvelopeTOL
      DiscreteModel%PlasticStrainWork = Variables%PlasticStrainWork
      DiscreteModel%StrainWorkDensity = Variables%StrainWorkDensity
      DiscreteModel%YYstrainDesityEnergy = Variables%YYstrainDesityEnergy
      DiscreteModel%Jintegral = Variables%Jintegral
      DiscreteModel%JintegralCum = Variables%JintegralCum
      DiscreteModel%Gamma = Variables%Gamma
      
      RETURN
    END SUBROUTINE
    
    
    SUBROUTINE SetupVariables(DiscreteModel,Variables)
    
      IMPLICIT NONE
      TYPE(Model),INTENT(IN):: DiscreteModel
      TYPE(VariablesContainer),INTENT(IN OUT):: Variables
    
      ALLOCATE( Variables%Displacement(DiscreteModel%NumEquations), &
                Variables%EffectiveStressTensor(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                Variables%PlasticStrainTensor(6,DiscreteModel%NumIntegrationPoints,&
                                                  DiscreteModel%NumElements),&
                Variables%YieldFunc(DiscreteModel%NumIntegrationPoints,&
                                        DiscreteModel%NumElements),&
                Variables%YieldEnvelopeTOL(DiscreteModel%NumElements),&
                Variables%StrainTensor(6,DiscreteModel%NumIntegrationPoints,&
                                            DiscreteModel%NumElements),&
                Variables%Porepressure(DiscreteModel%NumNodes),&
                Variables%ElementAVGporepressure(DiscreteModel%NumElements),&
                Variables%Temperature(DiscreteModel%NumNodes),&
                Variables%WaterSaturation(DiscreteModel%NumNodes),&
                Variables%DamageShear(DiscreteModel%NumElements),&
                Variables%DamageTraction(DiscreteModel%NumElements),&
                Variables%DamageStiffness(DiscreteModel%NumElements),&
                Variables%IsFractured(DiscreteModel%NumElements),&
                Variables%Fracpressure(DiscreteModel%NumElements),&
                Variables%PoreFluidBulkMod(DiscreteModel%NumElements),&
                Variables%PoreStrain(DiscreteModel%NumNodes),&
                Variables%PoreCompressibility(DiscreteModel%NumElements),&
                Variables%PlasticStrainWork(DiscreteModel%NumElements),&
                Variables%StrainWorkDensity(DiscreteModel%NumIntegrationPoints,DiscreteModel%NumElements),&
                Variables%YYstrainDesityEnergy(DiscreteModel%NumElements),&
                Variables%Jintegral(DiscreteModel%NumElements),&
                Variables%JintegralCum(3,DiscreteModel%NumElements),&
                Variables%Gamma(DiscreteModel%NumElements))
      Variables%Displacement = 0.0
      Variables%EffectiveStressTensor = 0.0
      Variables%StrainTensor = 0.0
      Variables%PlasticStrainTensor = 0.0
      Variables%YieldFunc = 0.0
      Variables%YieldEnvelopeTOL = 0.0
      Variables%Porepressure = 0.0
      Variables%ElementAVGporepressure = 0.0
      Variables%Temperature = 0.0
      Variables%WaterSaturation = 0.0
      Variables%DamageShear = 0.0
      Variables%DamageTraction = 0.0
      Variables%DamageStiffness = 0.0
      Variables%IsFractured = .FALSE.
      Variables%Fracpressure = 0.0
      Variables%PoreFluidBulkMod = 0.0
      Variables%PoreStrain = 0.0
      Variables%PoreCompressibility = 0.0
      Variables%PlasticStrainWork = 0.0
      Variables%StrainWorkDensity = 0.0
      Variables%YYstrainDesityEnergy = 0.0
      Variables%Jintegral = 0.0
      Variables%JintegralCum = 0.0
      Variables%Gamma = 1.0
      
      RETURN
    END SUBROUTINE
       
    
    SUBROUTINE EditProperties(DiscreteModel,Problem)
    
      USE ProblemDefinition

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: DiscreteModel
      TYPE (PD), INTENT(IN):: Problem
      
      INTEGER:: n,i,itheta,iz,irad,elementNum
    
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !PROPERTIES EDIT
      !REAL(8),ALLOCATABLE:: EditProps_values(:)
      !REAL(8),ALLOCATABLE:: EditProps_IJKinterval(:,:)
      !CHARACTER(LEN=40),ALLOCATABLE:: EditProps_PropsNames(:)
      n= SIZE(Problem%EditProps_PropsNames(:))
      
      DO i=1,n
          
        IF(TRIM(Problem%EditProps_PropsNames(i))=='KH')THEN
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  DiscreteModel%KH(elementNum) = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO
            
        END IF
        
        IF(TRIM(Problem%EditProps_PropsNames(i))=='KV')THEN
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  DiscreteModel%KV(elementNum) = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO

        END IF
                
        IF(TRIM(Problem%EditProps_PropsNames(i))=='Young')THEN
          
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  !DiscreteModel%YoungMod(elementNum) = Problem%EditProps_values(i)
                  DiscreteModel%TransverIsotropicProps(elementNum)%E1 = Problem%EditProps_values(i)
                  DiscreteModel%TransverIsotropicProps(elementNum)%E2 = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO
            
        END IF

        IF(TRIM(Problem%EditProps_PropsNames(i))=='Poisson')THEN
          
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  !DiscreteModel%Poisson(elementNum) = Problem%EditProps_values(i)
                    DiscreteModel%TransverIsotropicProps(elementNum)%n1 = Problem%EditProps_values(i)
                    DiscreteModel%TransverIsotropicProps(elementNum)%n2 = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO
            
        END IF
        
        IF(TRIM(Problem%EditProps_PropsNames(i))=='Cohesion')THEN
          
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  DiscreteModel%Cohesion(elementNum) = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO
            
        END IF
                
        IF(TRIM(Problem%EditProps_PropsNames(i))=='FrictionAngle')THEN
          
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  DiscreteModel%FrictionAngle(elementNum) = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO
            
        END IF
                        
        IF(TRIM(Problem%EditProps_PropsNames(i))=='DilationAngle')THEN
          
          elementNum = 0
          DO itheta=1,DiscreteModel%GridDivision(3)    !WellModel%NTHETA
            DO iz=1,DiscreteModel%GridDivision(2)      !WellModel%NZ
              DO irad=1,DiscreteModel%GridDivision(1)  !WellModel%NRAD
                elementNum = elementNum + 1
                  
                IF( irad>=Problem%EditProps_IJKinterval(1,i).AND.irad<=Problem%EditProps_IJKinterval(2,i).AND.&
                    iz>=Problem%EditProps_IJKinterval(3,i).AND.iz<=Problem%EditProps_IJKinterval(4,i).AND.&
                    itheta>=Problem%EditProps_IJKinterval(5,i).AND.itheta<=Problem%EditProps_IJKinterval(6,i)  )THEN
                  DiscreteModel%DilationAngle(elementNum) = Problem%EditProps_values(i)
                END IF
                  
              END DO
            END DO
          END DO
            
        END IF
        
      END DO
      
    RETURN
    END SUBROUTINE
        
    
  END MODULE DiscretizedModel
