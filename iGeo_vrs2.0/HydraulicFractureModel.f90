  
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************


  MODULE HydraulicFractureModel

      TYPE, PUBLIC:: HydraulicFracture

        INTEGER:: NumElements
          
        LOGICAL:: ContainsFracModel = .FALSE.
        !INTEGER:: StartDOF  !after fluid pressure.
        
        !INTEGER,ALLOCATABLE:: FracElsNodes(:,:)
        !REAL(8),ALLOCATABLE:: CellCenter(:,:)
        REAL(4),ALLOCATABLE:: AreaJ(:)
        REAL(4),ALLOCATABLE:: NodeArea(:)      !number of nodes of Geomechanical Model.
        REAL(4),ALLOCATABLE:: NormalJ(:,:)
        REAL(4),ALLOCATABLE:: DJ(:)
        INTEGER,ALLOCATABLE:: fracIndex(:,:)  !(ix,iz) fracture element for a given ix,iz index
        INTEGER,ALLOCATABLE:: FracToSolidElement(:),SolidToFracElement(:)
        
        !TRANSMISSIBILITY FRACTURE TO FRACTURE---------------------------------------------------------------
        REAL(8),ALLOCATABLE:: Trans(:)                !(4x4 Transmissibilities, frac Cell Number)
        INTEGER,ALLOCATABLE:: TransDOFs(:,:)          !(2 dofs,position at Trans)
        LOGICAL,ALLOCATABLE:: OpenFracElements(:)
        REAL(8),ALLOCATABLE:: ConectionFlowRates(:)   !stores internodes flow rates.
          
        INTEGER:: NumActiveElements
        INTEGER,ALLOCATABLE:: ListOfActiveElements(:)
        INTEGER,ALLOCATABLE:: ListOfyyNodes(:,:)
        
        !SPRINT STIFFNESS------------------------------------------------------------------------------------
        REAL(8),ALLOCATABLE:: Kn(:)
        INTEGER,ALLOCATABLE:: Kn_dofs(:)
                                
        !DYNAMIC PROPERTIES----------------------------------------------------------------------------------
        REAL(8),ALLOCATABLE:: Aperture(:),Aperture0(:)      !elemento based
        REAL(8),ALLOCATABLE:: NodeAperture0(:)              !nodal initial aperture
        REAL(8),ALLOCATABLE:: FluidPressure(:)
        INTEGER,ALLOCATABLE:: JfractureFront(:),KfractureFront(:)
        REAL(8),ALLOCATABLE:: Jintegral(:),J_h(:)
        REAL(8),ALLOCATABLE:: cumJintegral(:,:,:)           !(1,...NZ,...NX) - Z up J-Integral
                                                            !(2,...NZ,...NX) - Z down J-Integral
                                                            !(3,...NZ,...NX) - X J-Integral
        LOGICAL,ALLOCATABLE:: ActiveToFlow(:)               !define if the fracture element has flow.
          
        !FRACTURE DAMAGE MODEL -Contitinous damage mechanics-------------------------------------------------
        REAL(8),ALLOCATABLE:: Damage(:)                     !Damage variable for each element.
        
        LOGICAL:: FracClosure= .FALSE.                     !.TRUE. for fracture closing
          
      END TYPE

  CONTAINS
  
  
    ! UPDATES THE HYDRAULIC FRACTURE MODEL - CONTINUOS DAMAGE MODEL.
    SUBROUTINE UpdateHydraulicFracture(HFmodel,GeomechanicalModel,Load)
      
      USE NumericalLib
      USE NumericalUtilities
      USE IO
      USE DiscretizedModel

      IMPLICIT NONE
      TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),INTENT(IN OUT):: Load(0:)
      
      INTEGER,ALLOCATABLE:: g(:),num(:)
      REAL(8):: normal(3)
      REAL(8):: damageCutoff = 0.50,&
                minAperture = 1.0E-8, &
                maxAperture = 1.0E-2
      
      INTEGER:: i,j,k,iel,nip,dof
      
      REAL(8):: area, pp,  p, fy, force(3)
      REAL(8):: syy,sxy,szy,sn,an
      
      INTEGER:: ix,iz,NX,NZ
      
      REAL(8):: youngNormal
      !LOGICAL:: fractureClosing
      
      ALLOCATE( g(24),num(8) )
      
      nip = GeomechanicalModel%NumIntegrationPoints
      
      
      ! FRACTURE STRESS ANALYSIS-----------------------------------------
      ! Update element constitutive model.
      ! Update fracture load.
      DO i=1,HFmodel%NumElements
        
        iel = HFmodel%FracToSolidElement(i)
        
        !J-INTEGRAL------------------------------------------------------
        !HFmodel%Jintegral(i) = HFmodel%J_h(i) * GeomechanicalModel%YYplasticStrainDesityenergy(iel)
        !IF( HFmodel%Jintegral(i) > GeomechanicalModel%Gc(iel) )THEN
        !  HFmodel%OpenFracElements(i) = .TRUE.
        !  GeomechanicalModel%IsFractured(iel) = .TRUE.
        !END IF
        
        !STRUTURAL LOAD AND B.C.-----------------------------------------
        !IF( HFmodel%OpenFracElements(i) .AND. GeomechanicalModel%ElementIJK(3,iel)==1 )THEN
        IF( HFmodel%OpenFracElements(i) .AND. HFmodel%ActiveToFlow(i) .OR. HFmodel%Aperture0(i)>0.0 )THEN
            
          num = GeomechanicalModel%ElementNodes(:,iel)    
          area = HFmodel%AreaJ(i)                             !element YY stress.
          p = HFmodel%FluidPressure(i)                        !fracture pressure.
          pp = GeomechanicalModel%ElementAVGporepressure(iel) !!porous media average pressure.
          
          !normal stress (YY).
          sn = 0.0
          DO j=1,GeomechanicalModel%NumIntegrationPoints
            sn = sn + GeomechanicalModel%EffectiveStressTensor(2,j,iel)
          END DO
          sn = sn/GeomechanicalModel%NumIntegrationPoints - GeomechanicalModel%BiotCoef(iel) * pp
          
          !Y DOFs
          g(1) = GeomechanicalModel%NodeDOFs(2,num(1))
          g(2) = GeomechanicalModel%NodeDOFs(2,num(2))
          g(3) = GeomechanicalModel%NodeDOFs(2,num(3))
          g(4) = GeomechanicalModel%NodeDOFs(2,num(4))
          
          !Aperture 
          an = SUM(GeomechanicalModel%Displacement(g(1:4))) / 4.
          HFmodel%Aperture(i) = HFmodel%Aperture0(i) + an
          !IF(HFmodel%Aperture(i)<1.E-5) HFmodel%Aperture(i) = 1.E-5
          GeomechanicalModel%FracAperture(iel) = HFmodel%Aperture(i)
          
          !Umbalanced Force on fracture wall.
          !>>> negative fy force could cause non preciesly fracture propagation.
          fy = (p + sn) * area * HFmodel%Damage(i)                  !residual force to be applyed (over the total stress).
          IF ( .NOT.HFmodel%FracClosure .AND. fy<0.0 ) fy = 0.0     !to avoid oscilating solution in prescribed presseure constrained simulations.
          force = 0
          force(2) = fy
          Load(g(1:4)) = Load(g(1:4)) + force(2)/4.0
          
        END IF
        
      END DO
      
    RETURN
    END SUBROUTINE

    
    ! UPDATE HYDRAULIC TRANMISSIBILITY BASEAD ON DAMAGE AND APERTURE.
    SUBROUTINE UpdateHFtransmissibilityMatrix(HFmodel,GeomechanicalModel)

      USE NumericalLib
      USE NumericalUtilities
      USE FEMlib
      USE IO
      USE DiscretizedModel

      IMPLICIT NONE
      TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
          
      INTEGER:: fracElement,i,ip,iel,nFracsEls,j,counter

      REAL(8):: Kabs, aperture, &
                minAperture = 1.0E-8, &
                maxAperture = 2.0E-2
          
      CHARACTER(LEN=320)::elementType         !'hexahedron' ou 'tetrahedron'.
      INTEGER,ALLOCATABLE:: num(:)            !element incidence.
      REAL(8):: det
      REAL(8),ALLOCATABLE:: coord(:,:)        !element nodal coordenates.
      REAL(8),ALLOCATABLE:: points(:,:),&     !coordenadas dos pontos de integração.
                            weights(:)        !peso dos pontos de integração.
      REAL(8),ALLOCATABLE:: jac(:,:),der(:,:),deriv(:,:)
          
      REAL(8):: PermTensor(2,2)
      REAL(8),ALLOCATABLE:: elTransm(:,:)
      
      REAL(8):: damageCutoff = 0.95
      
      INTEGER:: ndim = 3                      !number of dimensions: 2D -> 2.
      INTEGER:: nod = 8                       !number of nodes per element.
      INTEGER:: nip = 8                       !number of integration points.
      
      REAL(8):: LX,LY,LZ,AX,AZ
      
      nFracsEls = HFmodel%NumElements
      
      !FEM CALCULATION.
      ALLOCATE(points(nip,ndim),coord(nod,ndim),                        &
                jac(ndim,ndim),weights(nip),num(nod),der(ndim,nod),     &
                deriv(ndim,nod),elTransm(nod,nod) )
      
      counter = 0.0
      DO i=1,nFracsEls
        fracElement = i
        iel = HFmodel%FracToSolidElement(fracElement)
              
        IF( HFmodel%OpenFracElements(fracElement) .OR. HFmodel%Aperture0(fracElement)>0.0 )THEN
        
          !   2--------3   
          !   |        |
          !   1--------4   --->Y
              
          num = GeomechanicalModel%ElementNodes(:,iel) !HFmodel%FracElsNodes(:,fracElement)
              
          coord(:,1) = GeomechanicalModel%X(num)
          coord(:,2) = GeomechanicalModel%Z(num)
          
          aperture = HFmodel%Aperture(fracElement) !- SUM(GeomechanicalModel%Displacement(g))/nod
          IF(aperture < minAperture .AND. aperture>0.0)THEN
            aperture = minAperture
          ELSEIF(aperture > maxAperture)THEN
            aperture = maxAperture
          ELSE
          END IF
          
          PermTensor = 0.0
          PermTensor(1,1) = aperture**2.0 / 12
          PermTensor(2,2) = PermTensor(1,1)
          
          LX = ABS(coord(1,1)-coord(4,1))
          LZ = ABS(coord(1,2)-coord(2,2))
          AX = aperture * LZ / 2.0
          AZ = aperture * LX / 2.0
          
          !1-2
          counter=counter+1
          HFmodel%Trans(counter) = PermTensor(2,2)*AZ/LZ
          HFmodel%TransDOFs(1,counter) = num(1)             !num(1)
          HFmodel%TransDOFs(2,counter) = num(1+1)           !num(2)
          !2-3
          counter=counter+1
          HFmodel%Trans(counter) = PermTensor(1,1)*AX/LX
          HFmodel%TransDOFs(1,counter) = num(1+1)           !num(2)
          HFmodel%TransDOFs(2,counter) = num(1+2)           !num(3)
          !3-4
          counter=counter+1
          HFmodel%Trans(counter) = PermTensor(2,2)*AZ/LZ
          HFmodel%TransDOFs(1,counter) = num(1+2)           !num(3)
          HFmodel%TransDOFs(2,counter) = num(1+3)           !num(4)
          !1-4
          counter=counter+1
          HFmodel%Trans(counter) = PermTensor(1,1)*AX/LX
          HFmodel%TransDOFs(1,counter) = num(1)             !num(1)
          HFmodel%TransDOFs(2,counter) = num(1+3)           !num(4)
          
          HFmodel%ActiveToFlow(fracElement) = .TRUE.
                  
        END IF
            
      END DO
      
      j = 0
      DO i=1,SIZE(HFmodel%Trans)
        IF(HFmodel%Trans(i).NE.0.0)THEN
          j=j+1
          HFmodel%ListOfActiveElements(j) = i
        END IF
      END DO
      HFmodel%ListOfActiveElements(0) = j
      
    RETURN
    END SUBROUTINE
    
    
    ! UPDATE HYDRAULIC FRACTURE PRESSURE WITH NODAL POREPRESSURE FROM
    ! GEOMECHANICAL MODEL.
    SUBROUTINE UpdateHFfluidPressure(GeomechanicalModel,HFmodel,Porepressure)

      USE NumericalLib
      USE NumericalUtilities
      USE FEMlib
      USE IO
      USE DiscretizedModel

      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
      REAL(8),INTENT(IN):: Porepressure(:)
      
      INTEGER:: num(4)            !element incidence.
    
      REAL(8):: damageCutoff = 0.95
      
      INTEGER:: iel,fracElement,nFracsEls
      
      nFracsEls = HFmodel%NumElements
      
      DO fracElement=1,nFracsEls
        
        IF(HFmodel%OpenFracElements(fracElement))THEN
          
          iel = HFmodel%FracToSolidElement(fracElement)
          num = HFmodel%ListOfyyNodes(:,fracElement)
          
          HFmodel%FluidPressure(fracElement) = SUM(Porepressure(num(1:4)))/4            !MAXVAL(Porepressure(num(1:4)))
          GeomechanicalModel%Fracpressure(iel) = HFmodel%FluidPressure(fracElement)
        
        END IF
      END DO
    
    RETURN
    END SUBROUTINE

    
    ! PROPAGATE FRACTURE.
    ! CRITERIA: traction damage based on plastic strain.
    SUBROUTINE PropagateFracture_TractionDamage(HFmodel,GeomechanicalModel,RestrictedDOFs,Stress)!,&
                                                !UnbalancedLoad)
    
      USE FEMlib
      USE IO
      USE DiscretizedModel
      USE MainConstants
    
      IMPLICIT NONE
      TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      INTEGER,INTENT(IN OUT):: RestrictedDOFs(0:)
      REAL(8),INTENT(IN):: Stress(:,:,:)
      !REAL(8),INTENT(IN OUT):: UnbalancedLoad(0:)
      
      INTEGER,ALLOCATABLE:: g(:),num(:)
      REAL(8):: normal(3)
      REAL(8):: damageCutoff = 0.50
      
      INTEGER:: i,j,k,iel,nip,dof
      
      REAL(8):: area, pp,  p, fy, force(3)
      REAL(8):: syy,sxy,szy,sn,an
      
      REAL(8):: J_integral
      
      INTEGER:: ix,iz,NX,NZ,numNodes
      
      INTEGER,ALLOCATABLE:: freeNodeCounter(:),nodeCounter(:)
      
      ALLOCATE( g(24),num(8) )
      
      nip = GeomechanicalModel%NumIntegrationPoints
      
      numNodes = GeomechanicalModel%NumNodes
      ALLOCATE( freeNodeCounter(numNodes),nodeCounter(numNodes) )
      freeNodeCounter = 0
      nodeCounter = 0
      
      !UnbalancedLoad = 0.0
      
      ! J-INTEGRAL AND FRACTURE PROPAGATION------------------------------
      HFmodel%cumJintegral = 0.
      DO ix=1,GeomechanicalModel%GridDivision(1)
          DO iz=1,GeomechanicalModel%GridDivision(2)
              i = HFmodel%fracIndex(iz,ix)
              iel = HFmodel%FracToSolidElement(i)
              
              IF(GeomechanicalModel%CDMflag(iel))THEN
                sn = Stress(2,1,iel)
                IF( sn>0 ) J_integral = sn * GeomechanicalModel%StrainTensor(2,1,iel) * HFmodel%DJ(i)                                         !OLD-> GeomechanicalModel%YYstrainDesityEnergy(iel) * HFmodel%DJ(i)
                IF( J_integral < 0 ) J_integral = 0.    
                
                GeomechanicalModel%Jintegral(iel) = J_integral                                                                              !OLD-> GeomechanicalModel%Jintegral(iel) + J_integral
                HFmodel%cumJintegral(1,iz,ix) = GeomechanicalModel%Jintegral(iel)
              END IF
              
          END DO
      END DO
      HFmodel%cumJintegral(2,:,:) = HFmodel%cumJintegral(1,:,:)
      HFmodel%cumJintegral(3,:,:) = HFmodel%cumJintegral(1,:,:)
      
      !VERTICAL CONTOUR (UP)
      DO ix=1,GeomechanicalModel%GridDivision(1)
          DO iz=2,GeomechanicalModel%GridDivision(2)
              i = HFmodel%fracIndex(iz,ix)
              iel = HFmodel%FracToSolidElement(i)
              HFmodel%cumJintegral(1,iz,ix) = HFmodel%cumJintegral(1,iz,ix) + HFmodel%cumJintegral(1,iz-1,ix)  !cummulative
              GeomechanicalModel%JintegralCum(1,iel) = HFmodel%cumJintegral(1,iz,ix)
              
              J_integral = HFmodel%cumJintegral(1,iz-1,ix)
              !IF( HFmodel%cumJintegral(1,iz,ix)>0. ) J_integral = J_integral + HFmodel%cumJintegral(1,iz,ix)
              J_integral = 2.*J_integral !symmetry condition
              !IF( HFmodel%OpenFracElements(i) .AND. J_integral>GeomechanicalModel%Gc(iel) )THEN !&                     !(8/8/2021) 
              IF( HFmodel%OpenFracElements(i) )THEN !&
                  !2.0 * HFmodel%cumJintegral(1,iz-1,ix) > GeomechanicalModel%Gc(iel))THEN
                  !2.0 * HFmodel%cumJintegral(1,iz,ix) > GeomechanicalModel%Gc(iel))THEN
                  
                  HFmodel%cumJintegral(1,iz,ix) = 0.
                  GeomechanicalModel%JintegralCum(1,iel) = 0.
                  !Frac previous element.
                  IF(J_integral>GeomechanicalModel%Gc(iel))THEN                                                         !(8/8/2021) 
                      i = HFmodel%fracIndex(iz-1,ix)
                      iel = HFmodel%FracToSolidElement(i)
                      HFmodel%OpenFracElements(i) = .TRUE.
                      GeomechanicalModel%IsFractured(iel) = .TRUE.
                      HFmodel%Damage(i) = 1.00
                      num = GeomechanicalModel%ElementNodes(:,iel)
                      freeNodeCounter(num(1:4)) = freeNodeCounter(num(1:4)) + 1
                      nodeCounter(num(1:4)) = nodeCounter(num(1:4)) + 1
                      !EXIT
                  END IF
              END IF
          END DO
      END DO
      
      !VERTICAL CONTOUR (DOWN)
      DO ix=1,GeomechanicalModel%GridDivision(1)
          DO iz=GeomechanicalModel%GridDivision(2)-1,1,-1
              i = HFmodel%fracIndex(iz,ix)
              iel = HFmodel%FracToSolidElement(i)
              HFmodel%cumJintegral(2,iz,ix) = HFmodel%cumJintegral(2,iz,ix) + HFmodel%cumJintegral(2,iz+1,ix)
              GeomechanicalModel%JintegralCum(2,iel) = HFmodel%cumJintegral(2,iz,ix)
                            
              J_integral = HFmodel%cumJintegral(2,iz+1,ix)
              !IF( HFmodel%cumJintegral(1,iz,ix)>0. ) J_integral = J_integral + HFmodel%cumJintegral(1,iz,ix)
              J_integral = 2.*J_integral !symmetry condition
              !IF( HFmodel%OpenFracElements(i) .AND. J_integral>GeomechanicalModel%Gc(iel) )THEN !&                         !(8/8/2021)
              IF( HFmodel%OpenFracElements(i) )THEN !&
                  !2.0 * HFmodel%cumJintegral(1,iz+1,ix) > GeomechanicalModel%Gc(iel))THEN
                  
                  HFmodel%cumJintegral(2,iz,ix) = 0.
                  GeomechanicalModel%JintegralCum(2,iel) = 0.
                  !Frac previous element.
                  IF( J_integral>GeomechanicalModel%Gc(iel) )THEN                                                           !(8/8/2021) 
                      i = HFmodel%fracIndex(iz+1,ix)
                      iel = HFmodel%FracToSolidElement(i)
                      HFmodel%OpenFracElements(i) = .TRUE.
                      GeomechanicalModel%IsFractured(iel) = .TRUE.
                      HFmodel%Damage(i) = 1.00
                      num = GeomechanicalModel%ElementNodes(:,iel)
                      freeNodeCounter(num(1:4)) = freeNodeCounter(num(1:4)) + 1
                      nodeCounter(num(1:4)) = nodeCounter(num(1:4)) + 1
                      !EXIT
                  END IF
              END IF
          END DO       
      END DO
      !HORIZONTAL CONTOUR (RIGHT)
      DO iz=1,GeomechanicalModel%GridDivision(2)
          DO ix=GeomechanicalModel%GridDivision(1)-1,1,-1
              i = HFmodel%fracIndex(iz,ix)
              iel = HFmodel%FracToSolidElement(i)
              HFmodel%cumJintegral(3,iz,ix) = HFmodel%cumJintegral(3,iz,ix) + HFmodel%cumJintegral(3,iz,ix+1)
              GeomechanicalModel%JintegralCum(3,iel) = HFmodel%cumJintegral(3,iz,ix)
                                          
              J_integral = HFmodel%cumJintegral(3,iz,ix+1)
              !IF( HFmodel%cumJintegral(1,iz,ix)>0. ) J_integral = J_integral + HFmodel%cumJintegral(1,iz,ix)
              J_integral = 2.*J_integral !symmetry condition
              !IF( HFmodel%OpenFracElements(i) .AND. J_integral>GeomechanicalModel%Gc(iel) )THEN !&                         !(8/8/2021)
              IF( HFmodel%OpenFracElements(i) )THEN    
                  !2.0 * HFmodel%cumJintegral(1,iz,ix+1) > GeomechanicalModel%Gc(iel))THEN
                  
                  HFmodel%cumJintegral(3,iz,ix) = 0.
                  GeomechanicalModel%JintegralCum(3,iel) = 0.
                  !Frac previous element.
                  IF( J_integral>GeomechanicalModel%Gc(iel) )THEN                                                           !(8/8/2021)
                      i = HFmodel%fracIndex(iz,ix+1)
                      iel = HFmodel%FracToSolidElement(i)
                      HFmodel%OpenFracElements(i) = .TRUE.
                      GeomechanicalModel%IsFractured(iel) = .TRUE.
                      HFmodel%Damage(i) = 1.00
                      num = GeomechanicalModel%ElementNodes(:,iel)
                      freeNodeCounter(num(1:4)) = freeNodeCounter(num(1:4)) + 1
                      nodeCounter(num(1:4)) = nodeCounter(num(1:4)) + 1
                      !EXIT
                  END IF
              END IF
              
          END DO
      END DO
      
      !UPDATING DISPLACEMENTS CONTRAINTS VECTOR--------------------------
      DO i=1,numNodes
        IF(nodeCounter(i).NE.0 .AND. freeNodeCounter(i)==nodeCounter(i) )THEN
          dof = GeomechanicalModel%NodeDOFs(2,i)
          DO j=1,RestrictedDOFs(0)
            IF(RestrictedDOFs(j) == dof)THEN
              RestrictedDOFs(j) = 0
              !UnbalancedLoad(dof) = -GeomechanicalModel%InternalForce(dof)
              EXIT
            END IF
          END DO
        END IF
      END DO
      
      j = 0
      DO i=1,RestrictedDOFs(0)
        IF( RestrictedDOFs(i) .NE. 0 )THEN
          j = j + 1
          RestrictedDOFs(j) = RestrictedDOFs(i)
        END IF
      END DO
      RestrictedDOFs(0) = j
      
    RETURN
    END SUBROUTINE
    
   
    ! BUILDS HYDRAULIC FRACTURE MODEL.   
    SUBROUTINE BuildHydraulicFractureModel( HFmodel,Problem,&
                                            GeomechanicalModel,RestrictedDOFs)
      
      USE ProblemDefinition
      USE NumericalLib
      USE NumericalUtilities
      USE FEMlib
      USE IO
      USE DiscretizedModel

      IMPLICIT NONE
      TYPE(HydraulicFracture),INTENT(IN OUT):: HFmodel
      TYPE (PD), INTENT(IN):: Problem
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      INTEGER,INTENT(IN OUT):: RestrictedDOFs(0:)
          
      INTEGER:: Nrad,Nz,Ntheta,nNodes,node
          
      INTEGER:: element,fracElement,irad,iz,itheta,i,j,ip,n(4),g(4),nPorousCells,nModels
      INTEGER:: K1,K2
      REAL(8):: Xlength
      REAL(8):: area,L1,L2,cellCenter,x

      REAL(8):: Kabs
          
      CHARACTER(LEN=320)::elementType         !'hexahedron' ou 'tetrahedron'.
      INTEGER:: num(8)                        !element incidence.
      REAL(8):: det
      REAL(8):: coord(3,3)        !element nodal coordenates.
      REAL(8),ALLOCATABLE:: points(:,:),&     !coordenadas dos pontos de integração.
                            weights(:)        !peso dos pontos de integração.
      REAL(8),ALLOCATABLE:: jac(:,:),der(:,:),deriv(:,:)
          
      REAL(8):: PermTensor(2,2)
      REAL(8),ALLOCATABLE:: elTransm(:,:)
      
      INTEGER:: ndim = 2                      !number of dimensions: 3D -> 3.
      INTEGER:: nod = 4                       !number of nodes per element.
      INTEGER:: nip                           !number of integration points.
      
      REAL(8):: st,maxStrain,Jc
      
      REAL(8):: coords1(3),coords2(3),coords3(3),coords(3,3),dx,dz
      
      LOGICAL:: found
      
      REAL(8):: e
      
      LOGICAL,ALLOCATABLE:: dofRestriction(:)
      INTEGER,ALLOCATABLE:: copyRestrictedDOFs(:)
      INTEGER:: dof
      
      REAL(8),ALLOCATABLE:: NodalForce(:)
      
      !LOGICAL,ALLOCATABLE:: springNodes(:)
      INTEGER:: numNodes
      
      INTEGER,ALLOCATABLE:: freeNodeCounter(:),nodeCounter(:)

      Ntheta = GeomechanicalModel%GridDivision(3)
      Nz = GeomechanicalModel%GridDivision(2)
      Nrad = GeomechanicalModel%GridDivision(1)
          
      ALLOCATE( HFmodel%SolidToFracElement(GeomechanicalModel%NumElements),&
                HFmodel%fracIndex(Nz,Nrad) )
      HFmodel%SolidToFracElement = 0
      HFmodel%fracIndex = 0
      
      !HYDRAULIC FRACTURE MODEL--------------------------------------------------------------------------------
      element = 0
      fracElement = 0
      itheta=1
      DO iz=1,Nz
        DO irad=1,Nrad
          element = element + 1
          !IF( GeomechanicalModel%CDMflag(element) ) THEN   !.OR. GeomechanicalModel%TRACflag(element) ) THEN
            fracElement = fracElement + 1
            HFmodel%SolidToFracElement(element) = fracElement
            HFmodel%fracIndex(iz,irad) = fracElement
          !END IF
        END DO
      END DO
      
      HFmodel%NumElements = fracElement
      IF(fracElement>0)THEN
        HFmodel%ContainsFracModel = .TRUE.
      ELSE
        HFmodel%ContainsFracModel = .FALSE.
      END IF
            
      !TRANSMISSIBILITIES--------------------------------------------------------------------------------------
      ALLOCATE(HFmodel%Trans(4*HFmodel%NumElements),HFmodel%TransDOFs(2,4*HFmodel%NumElements),&
               HFmodel%ConectionFlowRates(4*HFmodel%NumElements) )
      HFmodel%Trans = 0.0
      HFmodel%TransDOFs = 0
      HFmodel%ConectionFlowRates = 0.0
      
      ALLOCATE(HFmodel%FracToSolidElement(fracElement))
      HFmodel%FracToSolidElement = 0
      DO i=1,GeomechanicalModel%NumElements
        IF( HFmodel%SolidToFracElement(i)>0 ) &
                  HFmodel%FracToSolidElement(HFmodel%SolidToFracElement(i)) = i
      END DO
      
      ALLOCATE( HFmodel%Aperture(HFmodel%NumElements),&
                HFmodel%Aperture0(HFmodel%NumElements),&
                HFmodel%NodeAperture0(GeomechanicalModel%NumNodes),&
                HFmodel%OpenFracElements(HFmodel%NumElements),&
                HFmodel%FluidPressure(HFmodel%NumElements),&
                HFmodel%ListOfActiveElements(0:SIZE(HFmodel%Trans)),&
                HFmodel%AreaJ(HFmodel%NumElements),&
                HFmodel%NodeArea(GeomechanicalModel%NumNodes),&
                HFmodel%NormalJ(3,HFmodel%NumElements),&
                HFmodel%DJ(HFmodel%NumElements),&
                HFmodel%ListOfyyNodes(4,HFmodel%NumElements),&
                HFmodel%ActiveToFlow(HFmodel%NumElements)    )
      HFmodel%Aperture = 0.0
      HFmodel%Aperture0 = 0.0
      HFmodel%NodeAperture0 = 0.0
      HFmodel%OpenFracElements = .FALSE.
      HFmodel%FluidPressure = 0.0
      HFmodel%ListOfActiveElements = 0
      HFmodel%AreaJ = 0.0
      HFmodel%NodeArea = 0.
      HFmodel%NormalJ = 0.0
      HFmodel%DJ = 0.0
      HFmodel%ListOfyyNodes = 0
      HFmodel%ActiveToFlow = .FALSE.
      
      numNodes = GeomechanicalModel%NumNodes
      ALLOCATE( freeNodeCounter(numNodes),nodeCounter(numNodes) )
      freeNodeCounter = 0
      nodeCounter = 0
      
      ALLOCATE(HFmodel%cumJintegral(3,GeomechanicalModel%GridDivision(2),GeomechanicalModel%GridDivision(1)))
      HFmodel%cumJintegral = 0.
      
      !BUILDING FRACTURE ELEMENTS AND GEOMETRICAL PROPERTIES--------------------------------------------------
      !ALLOCATE( springNodes((Ntheta+1)*(Nz+1)) )
      !springNodes = .FALSE.
      HFmodel%NumActiveElements = 0
      DO fracElement=1,HFmodel%NumElements
          
        element = HFmodel%FracToSolidElement(fracElement)
        
        !
        !   2--------3   
        !   |        |
        !   1--------4   --->X
        !
        n(1) = GeomechanicalModel%ElementNodes(1,element)
        n(2) = GeomechanicalModel%ElementNodes(2,element)
        n(3) = GeomechanicalModel%ElementNodes(3,element)
        n(4) = GeomechanicalModel%ElementNodes(4,element)
        nodeCounter(n) = nodeCounter(n) + 1
        
        !springNodes(n) = .TRUE.
                        
        HFmodel%ListOfyyNodes(1,fracElement) = GeomechanicalModel%ElementNodes(1,element)
        HFmodel%ListOfyyNodes(2,fracElement) = GeomechanicalModel%ElementNodes(2,element)
        HFmodel%ListOfyyNodes(3,fracElement) = GeomechanicalModel%ElementNodes(3,element)
        HFmodel%ListOfyyNodes(4,fracElement) = GeomechanicalModel%ElementNodes(4,element)
            
        !Normal Area.
        coords1(1) = GeomechanicalModel%X( HFmodel%ListOfyyNodes(1,fracElement) )
        coords1(2) = GeomechanicalModel%Y( HFmodel%ListOfyyNodes(1,fracElement) )
        coords1(3) = GeomechanicalModel%Z( HFmodel%ListOfyyNodes(1,fracElement) )
        coords2(1) = GeomechanicalModel%X( HFmodel%ListOfyyNodes(4,fracElement) )
        coords2(2) = GeomechanicalModel%Y( HFmodel%ListOfyyNodes(4,fracElement) )
        coords2(3) = GeomechanicalModel%Z( HFmodel%ListOfyyNodes(4,fracElement) )
        dx = Distance(coords1,coords2)
        coords3(1) = GeomechanicalModel%X( HFmodel%ListOfyyNodes(2,fracElement) )
        coords3(2) = GeomechanicalModel%Y( HFmodel%ListOfyyNodes(2,fracElement) )
        coords3(3) = GeomechanicalModel%Z( HFmodel%ListOfyyNodes(2,fracElement) )
        dz = Distance(coords1,coords3)
        HFmodel%AreaJ(fracElement) = dx * dz
        HFmodel%NodeArea(n(1:4)) = HFmodel%NodeArea(n(1:4)) + HFmodel%AreaJ(fracElement) / 4.
                            
        !Normal vector.
        coord(1,:) = coords1
        coord(2,:) = coords2
        coord(3,:) = coords3
        CALL Triangle( coord, dx, coords1 )
        HFmodel%NormalJ(:,fracElement) = coords1
            
        !DJ size
        n(1) = GeomechanicalModel%ElementNodes(1,element)
        n(2) = GeomechanicalModel%ElementNodes(5,element)
        n(3) = GeomechanicalModel%ElementNodes(4,element)
        n(4) = GeomechanicalModel%ElementNodes(8,element)
        coords1(1) = 0.5 * ( GeomechanicalModel%X( n(1) ) + GeomechanicalModel%X( n(3) ) )
        coords1(2) = 0.5 * ( GeomechanicalModel%Y( n(1) ) + GeomechanicalModel%Y( n(3) ) )
        coords1(3) = 0.5 * ( GeomechanicalModel%Z( n(1) ) + GeomechanicalModel%Z( n(3) ) )
        coords2(1) = 0.5 * ( GeomechanicalModel%X( n(2) ) + GeomechanicalModel%X( n(4) ) )
        coords2(2) = 0.5 * ( GeomechanicalModel%Y( n(2) ) + GeomechanicalModel%Y( n(4) ) )
        coords2(3) = 0.5 * ( GeomechanicalModel%Z( n(2) ) + GeomechanicalModel%Z( n(4) ) )
        HFmodel%DJ(fracElement) = Distance(coords1,coords2)
            
        HFmodel%FluidPressure(fracElement) = GeomechanicalModel%Porepressure(element)

      END DO
      
      !PRE-FRACTURE DATA---------------------------------------------------------------------------------------
      
      ALLOCATE(HFmodel%Damage(HFmodel%NumElements))
      HFmodel%Damage = 0.0
      
      !ALLOCATE( HFmodel%JfractureFront(Nz),HFmodel%KfractureFront(Ntheta))
      !HFmodel%JfractureFront = 0
      !HFmodel%KfractureFront = Nz
      
      element = 0
      IF(Problem%WellData%PreFractureExists) THEN
        
        DO itheta=1,Ntheta
          DO iz=1,Nz
            DO irad=1,Nrad
              element = element + 1
              
              IF(itheta==1)THEN
                num = GeomechanicalModel%ElementNodes(:,element)
                x = SUM( GeomechanicalModel%X(num) )/8.0
                
                IF( x <= Problem%WellData%PreFractureLengh(iz) )THEN
                  
                  fracElement = HFmodel%SolidToFracElement(element)
                  HFmodel%OpenFracElements(fracElement) = .TRUE.
                  HFmodel%Damage(fracElement) = 1.00
                  HFmodel%Aperture0(fracElement) = Problem%WellData%PreFractureAperture(iz)
                  HFmodel%Aperture(fracElement) = Problem%WellData%PreFractureAperture(iz)
                  HFmodel%NodeAperture0(num) = 0  !Problem%WellData%PreFractureAperture(iz)
                  GeomechanicalModel%FracAperture(element) = HFmodel%Aperture(fracElement)
                  !GeomechanicalModel%DamageTraction(element) = 1.0
                  GeomechanicalModel%IsFractured(element) = .TRUE.
                  
                  !Node restriction.
                  freeNodeCounter(num(1:4)) = freeNodeCounter(num(1:4)) + 1
                  
                  !Fracture Front track.
                  !IF( itheta==1 .AND. irad>HFmodel%JfractureFront(iz) ) HFmodel%JfractureFront(iz) = irad
                  !IF( itheta==1 .AND. iz<HFmodel%KfractureFront(irad) ) HFmodel%KfractureFront(irad) = iz

                END IF
                
              END IF
              
            END DO
          END DO
        END DO
        
        !UPDATING DISPLACEMENTS CONTRAINTS VECTOR.
        DO i=1,numNodes
          IF(nodeCounter(i).NE.0 .AND. freeNodeCounter(i)==nodeCounter(i) )THEN
            dof = GeomechanicalModel%NodeDOFs(2,i)
            DO j=1,RestrictedDOFs(0)
              IF(RestrictedDOFs(j) == dof)THEN
                RestrictedDOFs(j) = 0
                EXIT
              END IF
            END DO
          END IF
        END DO
      
        j = 0
        DO i=1,RestrictedDOFs(0)
          IF( RestrictedDOFs(i) .NE. 0 )THEN
            j = j + 1
            RestrictedDOFs(j) = RestrictedDOFs(i)
          END IF
        END DO
        RestrictedDOFs(0) = j
        
      END IF
      
      HFmodel%FracClosure = Problem%FracClosure
      
    RETURN
  END SUBROUTINE

  END MODULE
