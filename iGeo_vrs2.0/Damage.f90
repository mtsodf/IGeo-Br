
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  
  MODULE Damage
  
  CONTAINS
  
  
    !*************************************************************************
    !DAMAGE.
    !*************************************************************************
    
    ! Damage Models.
    SUBROUTINE DamageUpdate(GeomechanicalModel,EffectiveStress,PlasticStrain)
  
      USE DiscretizedModel
      USE SolidMechanicsLib
  
      IMPLICIT NONE
      TYPE(Model),INTENT(IN OUT):: GeomechanicalModel
      REAL(8),INTENT(IN):: EffectiveStress(:,:,:),PlasticStrain(:,:,:)
      
      INTEGER:: iel,i,n,nip,ip,J
      REAL(8):: newCohesion,newFrictionAngle,strain,factor
      REAL(8):: sigm,dsbar,theta
      LOGICAL:: failure
      
      nip = GeomechanicalModel%NumIntegrationPoints
      
      !Liquefaction.
      IF(ALLOCATED(GeomechanicalModel%Liquefaction))THEN
        
        DO i=1,SIZE(GeomechanicalModel%Liquefaction)
          
          DO j=1,GeomechanicalModel%Liquefaction(i)%ElementList(0)
          
            iel = GeomechanicalModel%Liquefaction(i)%ElementList(j)
            failure = .FALSE.
            DO ip = 1, nip
              CALL invar(GeomechanicalModel%EffectiveStressTensor(:,ip,iel)+EffectiveStress(:,ip,iel),&
                         sigm,dsbar,theta)
              IF( sigm>GeomechanicalModel%Liquefaction(i)%MaxMeanEffectiveStress ) failure = .TRUE.
            END DO
            
            !Failure.
            IF(failure)THEN
              GeomechanicalModel%DamageShear(iel) = 1.0
              GeomechanicalModel%DamageStiffness(iel) = 0.99
            END IF
            
          END DO
        END DO
      
      END IF
        
      ! Plastic Shear Damage.
      n= SIZE(GeomechanicalModel%ShearStrainDamage)
      IF(n>0)THEN
          
        DO i=1,n
            
          DO j=1,GeomechanicalModel%ShearStrainDamage(i)%ElementList(0)
            iel = GeomechanicalModel%ShearStrainDamage(i)%ElementList(j)
          
            !plastic strain (deviatory)
            strain = 0.0
            DO ip = 1, nip
              CALL invar(GeomechanicalModel%PlasticStrainTensor(:,ip,iel)+PlasticStrain(:,ip,iel),&
                         sigm,dsbar,theta)
              strain = strain + dsbar
            END DO
            strain = strain / nip
              
            !Damage calculation.
            factor = 1. - strain / GeomechanicalModel%ShearStrainDamage(i)%LimitPlasticDeviatoricStrain * &
                                    (1.-GeomechanicalModel%ShearStrainDamage(i)%ResidualDamage)
            IF( factor < GeomechanicalModel%ShearStrainDamage(i)%ResidualDamage )&
                            factor = GeomechanicalModel%ShearStrainDamage(i)%ResidualDamage
            
            GeomechanicalModel%DamageShear(iel) = 1. - factor
            !GeomechanicalModel%DamageStiffness(iel) = factor
              
          END DO
            
        END DO
        
      END IF
      
      ! Reform Stiffness.
      
      
      RETURN
  END SUBROUTINE
    
  END MODULE
  