
  !*****************************************************************************
  !COMPUTATIONAL GEOMECHANICS
  !Eng. Luis Carlos de Sousa Junior
  !Programa de Pos-graduação em Geotecnia.
  !University of Brasilia


  !*****************************************************************************

  
  MODULE ConstitutiveGeoModels
    
    TYPE SoilModel
      
      CHARACTER(LEN=50):: Name
      
      REAL(8):: &
        density,&
        porosity,&
        poisson
      
      !Shear Modulus.
      !G=Ag+(s'v/Pa)^Bg  (kPa)
      REAL(8):: Ag,Bg
      
      !Shear Strength (Mohr-Coulomb)
      REAL(8):: cohesion, frictionAngle, dilationAngle
      
    END TYPE
    
    TYPE TransverseIsotroplicElasticity
        REAL(8):: E1,E2                     !young
        REAL(8):: n1,n2                     !poisson
        REAL(8):: strike,dip                !degrees
        REAL(8):: minE,minNi
    END TYPE
  
  CONTAINS
  
    !*************************************************************************
    !SOIL MODEL
    !*************************************************************************
    
    SUBROUTINE ViscoSoilStress( Strain, CurrentStress, Soil, &
                                StressIncrement, PlasticStrainIncrement,&
                                UmbalanceStressIncrement,&
                                Young, Yield, DeltaTime )
    
      USE SolidMechanicsLib
      
      IMPLICIT NONE
      REAL(8),INTENT(IN):: Strain(6),CurrentStress(6),Young
      TYPE(SoilModel),INTENT(IN):: Soil
      REAL(8),INTENT(IN OUT):: StressIncrement(6),PlasticStrainIncrement(6),&
                               UmbalanceStressIncrement(6),Yield,DeltaTime
        
      INTEGER:: nst = 6
      REAL(8):: dee(6,6),v,g
      REAL(8):: c,phi,psi
      REAL(8):: sigm,dsbar,lode_theta,f,dq1,dq2,dq3
      REAL(8),ALLOCATABLE:: m1(:,:),m2(:,:),m3(:,:),flow(:,:)
      
      ALLOCATE( m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst) )
      
      !Elaastic Stress---------------------------------------------------------
      v = Soil%Poisson
      CALL D_elastic_homogenious(dee,Young,v)
      StressIncrement = MATMUL(dee,Strain)
      
      !Mohr-Coulomb Envelop----------------------------------------------------
      c = Soil%cohesion    !coesão (kPa)
      phi = Soil%frictionAngle  !ângulo de atrito (graus)
      psi = Soil%dilationAngle  !ângulo de dilatação (graus)
      
      !Invariantes and yield function
      CALL invar(StressIncrement+CurrentStress,sigm,dsbar,lode_theta)
      CALL mocouf(phi,c,sigm,dsbar,lode_theta,f)
  
      !Strain rate.
      !dQ/dsigma = dQ/dsigm*dsigm/ds + dQ/dJ2*dJ2/ds + dQ/dJ3*dJ3/ds
      !devp/td = f * dQ/dsigma
      IF(f>0.0)THEN
        
        CALL mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)                       !dq1=dQ/dsigm, dq2=dQ/dJ2, dq3=dQ/dJ3
        CALL formm(StressIncrement+CurrentStress,m1,m2,m3)                  !m1=dsigm/ds, m2=dJ2/ds, m3=dJ3/ds
        flow = f * ( m1 * dq1 + m2 * dq2 + m3 * dq3 )
        
        PlasticStrainIncrement = MATMUL(flow,StressIncrement+CurrentStress) !f*dQ/dsigma
        PlasticStrainIncrement = PlasticStrainIncrement * DeltaTime
        UmbalanceStressIncrement =  MATMUL(dee,PlasticStrainIncrement)
        
      ELSE
        PlasticStrainIncrement = 0.0
      END IF
      
      Yield = f
      
    RETURN
    END SUBROUTINE
                                
                           
    SUBROUTINE SoilStress( Strain, CurrentStress, Soil, &
                           StressIncrement, Erate, Young )
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: Strain(6),CurrentStress(6)
      TYPE(SoilModel),INTENT(IN):: Soil
      REAL(8),INTENT(IN OUT):: StressIncrement(6),Erate(6),Young
      
      REAL(8):: dee(6,6),v,g
      
      Erate = 0.0
      
      v = Soil%Poisson
      
      CALL D_elastic_homogenious(dee,Young,v)
      
      StressIncrement = MATMUL(dee,Strain)
      
    RETURN
    END SUBROUTINE                       
    
    !G=Ag+(s'v/Pa)^Bg  (originally in kPa)
    FUNCTION ShearModulus(Model,EffetiveStress)
        
      USE SolidMechanicsLib
    
      IMPLICIT NONE
      REAL(8):: ShearModulus
      TYPE(SoilModel),INTENT(IN):: Model
      REAL(8),INTENT(IN):: EffetiveStress(:)
      
      REAL(8):: eigenValues(3),S1,sigm,dsbar,theta,pi = 3.1415926535897932384626433832795
      
      !CALL TensorEigenvalues(EffetiveStress,eigenValues)
      CALL invar(EffetiveStress,sigm,dsbar,theta)
      S1 = sigm + 2.0/3.0 * dsbar * SIN( theta - 2.0 * pi / 3.0 )
      IF(S1>0.0) S1=0                                             !traction test.
      
      ShearModulus = ( Model%Ag + (ABS(S1)/1.0) ** (Model%Bg) )*1.0E3  !(MPa)
    
    END FUNCTION

    !*************************************************************************
    !TRACTION STRENGTH MODEL
    !*************************************************************************
    
    !STRESS STRENGTH MODEL.
    !Fracture plane:: X
    SUBROUTINE StressDamage(SA,EA,StrainIncrement,Young,Poisson,Strength,Damage,SD)
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: SA(6)              !Current stress tensor.
      REAL(8),INTENT(IN):: EA(6)              !Current strain tensor (total).
      REAL(8),INTENT(IN):: StrainIncrement(6) !Strain increment.
      REAL(8),INTENT(IN):: Young, Poisson     !Homogenious elastic properties.
      REAL(8),INTENT(IN OUT):: Strength(3)    !1 - traction stregth.
                                              !2 - traction strain at 1.
                                              !3 - maximum traction strain.
      REAL(8),INTENT(IN OUT):: Damage         !Damage variable.
      REAL(8),INTENT(IN OUT):: SD(6)          !New Stress Tensor (total: current + increment)
    
      REAL(8):: SB(6),EB(6),&                 !stress and strain states at point B (trial).
                S0(6),E0(6),&                 !stress and strain states at point 0 (unload).
                SC(6),&                       !stress and strain states at point C (reaload).
                syy_eff,&
                e_unload(6),&                 !unload strain.
                deltaStress(6),&              !delta stress.
                YoungNormal,&                 !damaged young modulus.
                eyy_max,&                     !maximum YY strain.
                syy_max,&                     !traction strength.
                Fy,&                          !yield function.
                newFy,&
                factor,d
      
      REAL(8):: one = 1.0, zero = 0.0
      
      INTEGER:: i
      
      ! CURRENT STATE+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      eyy_max = Strength(3)
      syy_max = Strength(1)
              
      IF(Damage>0.98) THEN
        Damage = 0.98
        YoungNormal = (1.0-Damage) * Young
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,StrainIncrement)
        SD = SA + deltaStress
        RETURN
      END IF
      
      YoungNormal = (1.0-Damage) * Young

      SB = 0.0
      SC = 0.0
      SD = 0.0
      E0 = 0.0
      e_unload = 0.0
      Fy = 0.0
      
      EB = EA + StrainIncrement
      
      ! TRIAL STRESS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      IF(SA(2)<=0.0)THEN
        deltaStress = StressElastic(Young,Poisson,StrainIncrement)                      !stress increment.
      ELSE
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,StrainIncrement)
      END IF
      
      SB = SA + deltaStress
      
      syy_eff = SB(2)
      !IF(Damage<0.99)syy_eff = syy_eff / (1.0-Damage)
      Fy = TractionLinearYieldFunction(syy_eff,EB(2),syy_max,eyy_max) 
      
      ! TRIVIAL CASE (COMPRESSION TO COMPRESSION) (1.3 and 2.3)+++++++++++++++++++++++++
      IF(SA(2)<=0.0.AND.SB(2)<=0.0)THEN
        SD = SB
        RETURN
      END IF
      
      ! CASE 1.1 (TRACTION TO COMPRESSION)++++++++++++++++++++++++++++++++++++++++++++++
      IF(SA(2)>0.0.AND.SB(2)<0.0)THEN
        
        !RETURN TO (A->0).
        !Return normal do YY.
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,StrainIncrement)
        factor = RootFactor(zero,SA(2),one,SA(2)+deltaStress(2))
        e_unload = factor * StrainIncrement
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,e_unload)
        S0 = SA + deltaStress
        
        !0 TO NEW B.
        e_unload = (1-factor) * StrainIncrement
        deltaStress = StressElastic(Young,Poisson,e_unload) 
        SB = S0 + deltaStress
        SD = SB
        
        RETURN
      END IF
      
      ! CASE 2.1 (COMPRESSION TO TRACTION)++++++++++++++++++++++++++++++++++++++++++++++
      IF( SA(2)<0.0 .AND. SB(2)>0.0 .AND. Fy>0.0 )THEN
        
        !A TO 0.
        deltaStress = StressElastic(Young,Poisson,1.0*StrainIncrement)
        factor = RootFactor( zero,SA(2), one,SA(2)+deltaStress(2) )
        deltaStress = StressElastic(Young,Poisson,factor*StrainIncrement)
        S0 = SA + deltaStress
        
        !0 TO B (new damage).
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,(1-factor)*StrainIncrement)
        SB = S0 + deltaStress
        syy_eff = SB(2)
        !IF(Damage<0.99)syy_eff = SB(2) / (1.0-Damage)
        Fy = TractionLinearYieldFunction(syy_eff,EB(2),syy_max,eyy_max)
        !-> continue for next IF (TRACTION TO FAILURE).
      END IF  
      
      ! CASE 2.2 (TRACTION TO FAILURE)++++++++++++++++++++++++++++++++++++++++++++++++++
      IF( SA(2)>=0.0 .AND. SB(2)>0.0 .AND. Fy>0.0 )THEN
        
        !B TO 0.
        !Return normal do YY.
        e_unload = -StrainIncrement
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,e_unload)
        factor = RootFactor(zero,SB(2),one,SB(2)+deltaStress(2))
        e_unload = - (1.0-factor)*e_unload
        deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,e_unload)
        S0 = SB + deltaStress
        
        !0 TO C (new damage).
        DO i=1,2
          d = Damage + 0.001
          YoungNormal = (1.0-d) * Young
          deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,-e_unload)
          syy_eff = deltaStress(2)
          !IF(Damage<0.99)syy_eff = deltaStress(2) / (1.0-Damage)
          newFy = TractionLinearYieldFunction(syy_eff,EB(2),syy_max,eyy_max)
          factor = RootFactor(  Damage, Fy,&
                                d,      newFy )
          Damage = Damage + factor * (d-Damage)
          IF(Damage>0.98) Damage = 0.98
          YoungNormal = (1.0-Damage) * Young
          deltaStress = StressElasticOrtoYY(Young,Poisson,YoungNormal,-e_unload)
        END DO
        SC = S0 + deltaStress
        SD = SC
!write(15,*)SD(2)/1000.0/98.1,EB(2)
        RETURN

      ELSE
        SD = SB
      END IF 
      

      RETURN
    END SUBROUTINE
    
    
    ! Find root factor for 2 points line.
    ! RootFrator = (Xroot-X1)/(X2-X1)
    FUNCTION RootFactor(X1,F1,X2,F2)
    
      IMPLICIT NONE
      REAL(8):: RootFactor
      REAL(8),INTENT(IN):: X1,F1,X2,F2
    
      REAL(8):: xroot,tang
      
      tang = (F2-F1)/(X2-X1)
      
      xroot = X1 - F1 / tang
      
      RootFactor = (xroot-X1)/(X2-X1)
    
    END FUNCTION
    
    FUNCTION TractionLinearYieldFunction(StressYY,strainYY,TractionStrength,MaxStrainYY)
    
        IMPLICIT NONE
        REAL(8):: TractionLinearYieldFunction
        REAL(8),INTENT(IN):: StressYY,strainYY,TractionStrength,MaxStrainYY
        
        REAL(8):: s0,t
        
        ! Fy = Syy - ( s0 - t * StrainYY )
        
        ! Yield limit.
        s0 = TractionStrength
        t = -s0 / MaxStrainYY
        TractionLinearYieldFunction = StressYY - (s0 + t * strainYY)
    
      RETURN
    END FUNCTION
    
    FUNCTION CalcTractionStrength(T0,MaxPlasticDilation,PlasticDilation)
    
      IMPLICIT NONE
      REAL(8):: CalcTractionStrength
      REAL(8),INTENT(IN):: T0,MaxPlasticDilation,PlasticDilation
      
      IF(PlasticDilation>MaxPlasticDilation)THEN
        CalcTractionStrength = 0.0
      ELSE
        CalcTractionStrength = T0 * (1.0-PlasticDilation/MaxPlasticDilation)
      END IF
      
    RETURN
    END FUNCTION
    
    !Potencial de escoamento para o critério de Mohr-Coulomb.
    SUBROUTINE tractionDQDS(tensor_s,ST,DQDS)
      ! This subroutine forms the derivatives of a Traction potential
      ! dq1=dQ/dsigm, dq2=dQ/dJ2, dq3=dQ/dJ3
      
      USE SolidMechanicsLib
  
      IMPLICIT NONE
      REAL(8),INTENT(IN)::ST,tensor_s(:)
      REAL(8),INTENT(OUT)::DQDS(6)
      
      REAL(8):: sigm, dsbar, lode_theta, S3, S3_normal(3)  !invariantes de tensão.
      REAL(8):: eigenValues(3),eigenVectors(3,3),tensor(6),tensorRotated(6)
      INTEGER:: position(1)
      
      CALL TensorEig(tensor_s,eigenValues,eigenVectors)
      position = MAXLOC(eigenValues(1:3))

      tensor = 0.0
      tensor(position(1)) = 1.0
      
      CALL RotateStress(tensor,eigenVectors,tensorRotated)
      
      DQDS = tensorRotated ! * tensor_s
      
      RETURN
    END SUBROUTINE
  
    
    !DAMAGE SINGLE SPRING MODEL. NOT TESTED.
    !Fracture plane: X
    SUBROUTINE SpringStress(DeltaY,Area,Kn0,Stmax,UYmax,Damage,Fy,DeltaFy)
    
      IMPLICIT NONE
      REAL(8),INTENT(IN)::  DeltaY
      REAL(8),INTENT(IN OUT)::  Fy,&      !Current force.
                                DeltaFy,&
                                Damage    !current damage.
      REAL(8),INTENT(IN)::  Area,&        !area of influence (stress<->force).
                            Kn0,&         !original stiffness.
                            Stmax,&       !maximum traction stress.
                            UYmax          !maximum displacement.
      
      REAL(8):: kn,SyyC,SyyD,u0,St,deltaForce
      REAL(8):: currentDamage
      
      kn = (1.0-Damage) * Kn0
      DeltaFy = Fy + Kn * DeltaY
      
      IF( DeltaFy > 0.0 )THEN
        
        SyyC = DeltaFy / Area
        u0 = DeltaY - DeltaFy / kn
        St = Stmax - Stmax/UYmax * (deltaY-u0)
        
        IF( SyyC>St)THEN
          SyyD = St
          DeltaFy = SyyD * Area
          Damage = 1.0 - Fy / (deltaY-u0) / kn0
        END IF
        
      END IF
      
      DeltaFy = DeltaFy - Fy
      
      RETURN
    END SUBROUTINE 
    
    !*************************************************************************
    !MOHR-COULOMB.
    !*************************************************************************
    
    !MOHR-COULOMB STRESS CALCULATOR.
    FUNCTION StressElastoplasticMC( Young,Poisson,CurrentStress,&
                                    Cohesion,FrictionAngle,DilationAngle,&
                                    Strain,StrainEP,YieldFunc)
    
      USE NumericalUtilities
      USE SolidMechanicsLib
    
      IMPLICIT NONE
      REAL(8):: StressElastoplasticMC(6)
      REAL(8),INTENT(IN):: Young,Poisson,CurrentStress(6)
      REAL(8),INTENT(IN):: Cohesion,FrictionAngle,DilationAngle
      REAL(8),INTENT(IN):: Strain(6)
      REAL(8),INTENT(IN OUT):: StrainEP(6),YieldFunc
      
      REAL(8):: dee(6,6),Depl(6,6),stressIncrement(6),&
                SC(6) !stress at failure envelope.
      REAL(8):: sigm, dsbar, lode_theta       !stress invariants.
      REAL(8):: phi,c,dil,f                   !MC paramters.
      REAL(8):: f_old, factor
      
      StrainEP = 0.0
      
      !Elastic Stress---------------------------------------------------------------
      CALL D_elastic_homogenious(dee,Young,Poisson)
      stressIncrement = MATMUL(dee,Strain)
      
      !Check failure envelope and calculate viscoplastic strain---------------------
      CALL invar(CurrentStress+stressIncrement, sigm, dsbar, lode_theta)
      c = Cohesion
      phi = FrictionAngle
      CALL mocouf(phi,c,sigm,dsbar,lode_theta,f)
      
      !Elastoplastic calculation-----------------------------------------------------
      
      IF(f>0.0)THEN
        
        !previous stress state.
        CALL invar(CurrentStress, sigm, dsbar, lode_theta)
        CALL mocouf(phi,c,sigm,dsbar,lode_theta,f_old)
        
        !factor to return to yield envelope.
        factor = f / (f - f_old)
        
        !stress at failure envelop.
        SC = CurrentStress + (1.0-factor) * stressIncrement
        
        !elastoplastic D matrix.
        dil = DilationAngle
        CALL mcdpl(phi,dil,dee,SC,Depl)
        Depl = factor * Depl
        StrainEP = factor * Strain
      
        !new stress.
        stressIncrement = stressIncrement - MATMUL(Depl,Strain)
        
        !final yield function.
        CALL invar(CurrentStress+stressIncrement, sigm, dsbar, lode_theta)
        CALL mocouf(phi,c,sigm,dsbar,lode_theta,f)
        
      END IF
      
      !Final stress increment-------------------------------------------------------
      StressElastoplasticMC = stressIncrement
      
      !Final Yield------------------------------------------------------------------
      YieldFunc = f
    
    END FUNCTION
    

    !Yield function for Mohr-Coulomb
    SUBROUTINE mocouf(phi,c,sigm,dsbar,theta,f)
      ! This subroutine calculates the value of the yield function
      ! for a Mohr-Coulomb material (phi in degrees).

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::phi,c,sigm,dsbar,theta   
      REAL(iwp),INTENT(OUT)::f
      REAL(iwp):: phir,snph,csph,csth,snth,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,   &
                  d180=180.0_iwp
                  phir=phi*d4*ATAN(one)/d180
                  snph=SIN(phir) 
                  csph=COS(phir) 
                  csth=COS(theta)
                  snth=SIN(theta)
      f=snph*sigm+dsbar*(csth/SQRT(d3)-snth*snph/d3)-c*csph
      RETURN
    END SUBROUTINE mocouf

    !Rotina de Smith e Griffiths para cálculo da função potencial de escoamento para o critério de Mohr-Coulomb.
    SUBROUTINE mocouq(psi,dsbar,theta,dq1,dq2,dq3)
      ! This subroutine forms the derivatives of a Mohr-Coulomb potential
      ! function with respect to the three invariants (psi in degrees).
  
      IMPLICIT NONE
      REAL(8),INTENT(IN)::psi,dsbar,theta
      REAL(8),INTENT(OUT)::dq1,dq2,dq3
      REAL(8):: psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth,zero=0.0,     &
                pt49=0.49,pt5=0.5,one=1.0,d3=3.0,d4=4.0,  &
                d180=180.0
                psir=psi*d4*ATAN(one)/d180 
                snth=SIN(theta) 
                snps=SIN(psir)
                sq3=SQRT(d3)  
                dq1=snps
      IF(ABS(snth).GT.pt49)THEN
        c1=one
        IF(snth.LT.zero)c1=-one
        dq2=(sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar 
        dq3=zero
      ELSE
        csth=COS(theta)
        cs3th=COS(d3*theta)
        tn3th=TAN(d3*theta)
        tnth=snth/csth
        dq2=sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5
        dq3=pt5*d3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
      END IF
      RETURN
    END SUBROUTINE mocouq
  
    SUBROUTINE mcdpl(phi,psi,dee,stress,pl)
      ! This subroutine forms the plastic stress/strain matrix
      ! for a Mohr-Coulomb material (phi,psi in degrees).
      
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(8),INTENT(IN)::stress(:),dee(:,:),phi,psi  
      REAL(8),INTENT(OUT)::pl(:,:)
      
      REAL(iwp),ALLOCATABLE::dfds(:),dqds(:),ddqds(:),dfdsd(:)
      REAL(iwp)::t1,t2,t3,t4,t5,t6,t8,t10,t12,t13,t14,t15,t16,t17,t18,t19,t20, &
        t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,   &
        t38,t39,t40,t41,t42,t43,t44,t45,t46,t48,t49,t50,t51,t53,t54,t55,t56,   &
        t60,t61,t63,t64,t68,t69,t70,t71,t73,t74,t77,t79,t80,t82,t83,t84,t85,   &
        t86,t89,t92,t93,t94,t97,t98,t101,t103,t106,t110,t111,t113,t122,t129,   &
        t133,t140,t145,t152,t166,t186,t206,pm,pi,phir,snph,snth,sq3,sx,sy,sz,  &
        txy,tyz,tzx,zero=0.0_iwp,pt49=0.49_iwp,one=1.0_iwp,two=2.0_iwp,        &
        d3=3.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d180=180.0_iwp,snps,psir
      
      REAL(iwp)::denom
      INTEGER::i,j,ih
      ih=SIZE(stress)
      
      ALLOCATE(dfds(ih),dqds(ih),ddqds(ih),dfdsd(ih))
      
      pi=ACOS(-one) 
      phir=phi*pi/d180
      snph=SIN(phir)
      psir=psi*pi/d180
      snps=SIN(psir)
      sq3=SQRT(d3)
      
      sx=stress(1)
      sy=stress(2)
      sz=stress(3)   
      txy=stress(4) 
      tyz=stress(5) 
      tzx=stress(6) 
      t3=one/d3
      t4=(sx+sy+sz)*t3
      t5=sx-t4
      t6=sy-t4
      t8=sz-t4
      t10=tyz**2
      t12=tzx**2
      t14=txy**2
      t23=(sx-sy)**2
      t25=(sy-sz)**2
      t27=(sz-sx)**2
      t32=SQRT((t23+t25+t27)/d6+t14+t10+t12)
      t33=t32**2
      t37=d3*sq3*(t5*t6*t8-t5*t10-t6*t12-t8*t14+two*txy*tyz*tzx)/two/t33/t32
      IF(t37>one)t37=one
      IF(t37<-one)t37=-one
      t38=ASIN(t37)
      t40=SIN(t38*t3)
      snth=-t40
      IF(ABS(snth).GT.pt49)THEN
        pm=-one
        IF(snth.LT.zero)pm=one  
        t2=snph/d3
        t4=(sx-sy)**2
        t6=(sy-sz)**2
        t8=(sz-sx)**2
        t10=one/d6
        t12=txy**2
        t13=tyz**2
        t14=tzx**2
        t16=SQRT((t4+t6+t8)*t10+t12+t13+t14)
        t19=one/t16/sq3
        t20=one/two
        t21=t19*t20
        t23=d3+pm*snph
        t25=two*sy
        t26=two*sz
        t33=two*sx
        t48=t20*t23
        dfds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
        dfds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
        dfds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
        dfds(4)=t19*t48*txy
        dfds(5)=t19*t48*tyz
        dfds(6)=t19*t48*tzx
        t2=snps/d3
        t23=d3+pm*snps
        t48=t20*t23
        dqds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
        dqds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
        dqds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
        dqds(4)=t19*t48*txy
        dqds(5)=t19*t48*tyz
        dqds(6)=t19*t48*tzx
      ELSE
        t1=one/d3
        t2=snph*t1
        t4=(sx-sy)**2
        t6=(sy-sz)**2
        t8=(sz-sx)**2
        t10=one/d6
        t12=txy**2
        t13=tyz**2
        t14=tzx**2
        t15=(t4+t6+t8)*t10+t12+t13+t14
        t16=SQRT(t15)
        t18=d3*sq3
        t20=(sx+sy+sz)*t1
        t21=sx-t20
        t22=sy-t20
        t23=t21*t22
        t24=sz-t20
        t29=two*txy
        t32=t23*t24-t21*t13-t22*t14-t24*t12+t29*tyz*tzx
        t33=one/two
        t35=t16**2
        t37=one/t35/t16
        t39=t18*t32*t33*t37
        IF(t39>one)t39=one
        IF(t39<-one)t39=-one
        t40=ASIN(t39)
        t41=t40*t1
        t42=COS(t41)
        t43=SIN(t41)
        t45=one/sq3
        t48=one/t16*(t42+t43*snph*t45)
        t50=two*sy
        t51=two*sz
        t53=(d4*sx-t50-t51)*t10
        t56=one-t1
        t60=t21*t1*t24
        t61=t23*t1
        t63=t1*t14
        t64=t1*t12
        t69=t18*t32
        t70=t35**2
        t73=t33/t70/t16
        t77=t18*(t56*t22*t24-t60-t61-t56*t13+t63+t64)*t33*t37-               &
          d3/two*t69*t73*t53
        t79=d3**2
        t80=sq3**2
        t82=t32**2
        t83=two**2
        t86=t15**2
        t92=SQRT(one-t79*t80*t82/t83/t86/t15)
        t93=one/t92
        t94=t93*t1
        t98=t2*t45
        t103=two*sx
        t106=(-t103+d4*sy-t51)*t10
        t110=t1*t22*t24
        t113=t1*t13
        t122=t18*(-t110+t21*t56*t24-t61+t113-t56*t14+t64)*t33*t37-           &
          d3/two*t69*t73*t106
        t133=(-t50+d4*sz-t103)*t10
        t145=t18*(-t110-t60+t23*t56+t113+t63-t56*t12)*t33*t37-               &
          d3/two*t69*t73*t133
        t166=t18*(-two*t24*txy+two*tyz*tzx)*t33*t37-d3*t69*t73*txy
        t186=t18*(-two*t21*tyz+t29*tzx)*t33*t37-d3*t69*t73*tyz
        t206=t18*(-two*t22*tzx+t29*tyz)*t33*t37-d3*t69*t73*tzx
        dfds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
        dfds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
        dfds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
        dfds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
        dfds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
        dfds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
        t2=snps*t1
        t48=one/t16*(t42+t43*snps*t45)
        t98=t2*t45
        dqds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
        dqds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
        dqds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
        dqds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
        dqds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
        dqds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
      END IF

      ddqds=MATMUL(dee,dqds)
      dfdsd=MATMUL(dee,dfds) 
      denom=DOT_PRODUCT(dfdsd,dqds)
      
      DO i=1,ih
        DO j=1,ih
          pl(i,j)=ddqds(i)*dfdsd(j)/denom
        END DO
      END DO
      DEALLOCATE(dfds,dqds,ddqds,dfdsd)
      
    RETURN
    END SUBROUTINE mcdpl

    !*************************************************************************
    !HYPERBOLIC MODEL.
    !*************************************************************************
    
    SUBROUTINE StressHYP( Poisson,E_ref,P_ref,M,Cp,Cr,Phi,Phir,DilAngle,&
                          CurrentStress,Strain,StressIncrement,StrainEP,YieldFunc)
    
      USE SolidMechanicsLib
      
      IMPLICIT NONE
      REAL(8),INTENT(IN):: Strain(6),CurrentStress(6)
      REAL(8),INTENT(IN):: Poisson,E_ref,P_ref,M,Cp,Cr,Phi,Phir,DilAngle
      REAL(8),INTENT(IN OUT):: StressIncrement(6),StrainEP(6),YieldFunc
      
      REAL(8):: dee(6,6),Depl(6,6),SC(6)
      REAL(8):: sigm, dsbar, lode_theta       !stress inveriants.
      REAL(8):: e,v,s3,Cohesion,FrictionAngle,f_old,factor
      REAL(8):: pi = 3.1415926535897932384626433832795
      REAL(8):: pa = 101325.

      !Update Elastic properties----------------------------------------------------
      
      CALL invar( CurrentStress, sigm, dsbar, lode_theta )
      s3 = sigm + 2.0/3.0 * dsbar * SIN( lode_theta + 2.0 * pi / 3.0)
      
      e = E_ref * ( (-s3)/pa )**M
      v = Poisson
      
      !Elastic Stress---------------------------------------------------------------
      
      CALL D_elastic_homogenious(dee,e,v)
      stressIncrement = MATMUL(dee,Strain)
      CALL invar( CurrentStress+stressIncrement, sigm, dsbar, lode_theta )
      
      !Envelope parameters----------------------------------------------------------
      
      CALL mocouf(Phir,Cr,sigm,dsbar,lode_theta,YieldFunc)
      
      IF( YieldFunc>0.0 )THEN   !priotiry: low stress envelope.
        !Low stress.
        Cohesion = Cr
        FrictionAngle = Phir
      ELSE
        !High stress.
        Cohesion = Cp
        FrictionAngle = Phi
        CALL mocouf(Phi,Cp,sigm,dsbar,lode_theta,YieldFunc)
      END IF
      
      !Elastoplastic calculation-----------------------------------------------------
      
      IF(YieldFunc>0.0)THEN
        
        !previous stress state.
        CALL invar( CurrentStress, sigm, dsbar, lode_theta )
        CALL mocouf(FrictionAngle,Cohesion,sigm,dsbar,lode_theta,f_old)
        
        !factor to return to yield envelope.
        factor = YieldFunc / (YieldFunc - f_old)
        
        !stress at failure envelop.
        SC = CurrentStress + (1.0-factor) * stressIncrement
        
        !elastoplastic D matrix.
        CALL mcdpl(FrictionAngle,DilAngle,dee,SC,Depl)
        Depl = factor * Depl
        StrainEP = factor * Strain
      
        !new stress.
        stressIncrement = stressIncrement - MATMUL(Depl,Strain)
        
        !final yield function.
        CALL invar(CurrentStress+stressIncrement,sigm,dsbar,lode_theta)
        CALL mocouf(FrictionAngle,Cohesion,sigm,dsbar,lode_theta,YieldFunc)
        
      END IF
      
    END SUBROUTINE
      
    !*************************************************************************
    !CAM-CLAY.
    !*************************************************************************
    
    FUNCTION StressMCC(Strain,Young,Poisson,CurrentStress,&
                       Pc,M,kappa,lambda,e0,P0,StrainEP)
    
      USE SolidMechanicsLib
      
      IMPLICIT NONE
      REAL(8):: StressMCC(6)
      REAL(8),INTENT(IN):: Strain(6)
      REAL(8),INTENT(IN):: M,kappa,lambda,e0,P0,CurrentStress(6)
      REAL(8),INTENT(IN OUT):: Pc,Young,Poisson
      REAL(8),INTENT(IN OUT):: StrainEP(6)
      
      REAL(8):: dee(6,6),Depl(6,6),stressIncrement(6),&
                SC(6),& !stress at failure envelope.
                evol,evol_vp
      REAL(8):: sigm, dsbar, lode_theta       !stress inveriants.
      REAL(8):: f, f_old, factor
      REAL(8):: ni
      
      evol = 0.0
      
      !Current elastic properties---------------------------------------------------
      evol_vp = SUM(StrainEP(1:3))
      evol =  SUM(Strain(1:3))
      CALL mccdProps(M,e0,P0,lambda,kappa,SUM(CurrentStress(1:3))/3.0,&
                     evol,evol_vp,Pc,Young,Poisson)
    
      !Elastic Stress---------------------------------------------------------------
      CALL D_elastic_homogenious(dee,Young,Poisson)
      stressIncrement = MATMUL(dee,Strain)
      
      !Check failure envelope and calculate viscoplastic strain---------------------
      CALL invar(CurrentStress+stressIncrement, sigm, dsbar, lode_theta)
      f = mccf(M,Pc,sigm,dsbar)

      !Elastoplastic calculation-----------------------------------------------------
      IF(f>1.0E-6*ABS(sigm))THEN
        
        !previous stress state.
        CALL invar(CurrentStress, sigm, dsbar, lode_theta)
        f_old = mccf(M,Pc,sigm,dsbar)

        !factor to return to yield envelope.
        IF(f.NE.f_old)THEN
          factor = f / (f - f_old)
        ELSE
          factor=0.0
        END IF
        
        !stress at failure envelop.
        SC = CurrentStress + (1.0-factor) * stressIncrement
        
        !elastoplastic D matrix.
        CALL mccdpl(M,Pc,e0,P0,lambda,kappa,dee,SC,Depl)
        Depl = factor * Depl
        
        !new stress.
        stressIncrement = stressIncrement - MATMUL(Depl,Strain)
        
        StrainEP = StrainEP + factor * Strain

!SC = CurrentStress + stressIncrement
!call invar(SC,sigm,dsbar,lode_theta)
!f = mccf(M,Pc,sigm,dsbar)
!write(16,'(4F10.2)')Pc/1e5,sigm/1e5,dsbar/1e5,f/1e5

      END IF
      
      !Final stress increment-------------------------------------------------------
      
      StressMCC = stressIncrement
    
!write(9,'(5E15.4)')e0,P0,sigm,dsbar,evol

    RETURN
    END FUNCTION
    
    !MODIFIED CAM-CLAY YIELD FUNCTION.                   
    FUNCTION mccf(M,Pc,P,Q)
    
      IMPLICIT NONE
      REAL(8):: mccf
      REAL(8),INTENT(IN)::  M,Pc,P,Q
      
      !mccf = (Q/M)**2.0 - P * (Pc-P)
      !mccf = -M**2.0 * P + M**2.0 * Pc - Q**2.0 / P
      mccf = M**2 * P**2.0 - M**2.0 * Pc * P + Q**2.0
      
      RETURN
    END FUNCTION mccf
                   
    !MODIFIED CAM-CLAY ELASTOPLASTIC STRESS-STRAIN MATRIX.
    SUBROUTINE mccdpl(M,Pc,e0,P0,lambda,kappa,dee,stress,pl)
      
      USE SolidMechanicsLib
      USE NumericalUtilities
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: M,e0,P0,lambda,kappa
      REAL(8),INTENT(IN OUT):: Pc
      REAL(8),INTENT(IN):: stress(:),dee(:,:) 
      REAL(8),INTENT(IN OUT):: pl(:,:)
      
      REAL(8),ALLOCATABLE::dfds(:),dqds(:),dpds(:),dJ2ds(:),d_dqds(:),dfds_d(:)
      
      REAL(8):: p, q, lode_theta       !stress inveriants.
      REAL(8):: x
      REAL(8):: dfdp, dfdq, dfdj3
      
      REAL(8)::denom
      INTEGER::i,j,ih
      
      REAL(8):: ni,h(6),hp       !hardening modulus.
      REAL(8):: e
      
      ih = SIZE(stress)
      ALLOCATE(dfds(ih),dqds(ih),dpds(ih),dJ2ds(ih),d_dqds(ih),dfds_d(ih))
      
      CALL invar(stress, p, q, lode_theta)
      
      ! DERIVATIVES.
      !F = M^2 * P^2 - M^2 * P * Pc + Q^2.0
      
      !dF/dP, dF/dQ, dF/dJ3
      dfdp = 2.0*M**2.0*P - M**2.0 * Pc
      dfdq = 2.0 * Q
      dfdj3 = 0.0
      
      !dP/dsigma
      dpds = 0.0
      dpds(1:3) = 1.0/3.0
      
      !dQ/dsigma
      IF(q>0.0)THEN
        x = 2.0 * q**2.0
        dJ2ds(1) = ( 2.0*(stress(1)-stress(2))-2.0*(stress(3)-stress(1)) )
        dJ2ds(2) = ( 2.0*(stress(2)-stress(3))-2.0*(stress(1)-stress(2)) )
        dJ2ds(3) = ( 2.0*(stress(3)-stress(1))-2.0*(stress(2)-stress(3)) )
        dJ2ds(4) = ( 12.0*stress(4) )
        dJ2ds(5) = ( 12.0*stress(5) )
        dJ2ds(6) = ( 12.0*stress(6) )
        dJ2ds = dJ2ds * 1.0/SQRT(2.0) * 1.0/2.0 / SQRT(x)
      ELSE
        dJ2ds = 0.0
      END IF
      
      !dF/ds
      dfds = dfdp*dpds + dfdq*dJ2ds     !<- dFunction/ds
      dqds = dfds                       !<- dPotencial/ds
      
      !HARDENING.
      h = 0.0
      hp = 0.0
                                                
      !F = M^2 * P^2 - M^2 * Pc * P + Q^2                                          
      !hp = - dF/dK * {dK/devol_vp} . {dF/dp}
      ni = 1+e0                                 !specific volume
      ni = ni - lambda*LOG(Pc/p0)               !Pc updated
      !h(1:3) = -p * (-P0*ni/(lambda-kappa))
      !hp = DOT_PRODUCT(h,dqds)
      IF(lambda.NE.kappa) hp = -(-M**2.0*P) * (ni*Pc/(lambda-kappa)) * ( 2*M**2*P-M**2*Pc )
      
      !FINAL EQUATION.
      !Dee = [D]*{dQds} * {dFds}*[D] / ( {dFds}.[D].{dQds} + h )
      d_dqds = MATMUL(dee,dqds)
      dfds_d = MATMUL(dee,dfds) 
      denom = DOT_PRODUCT(d_dqds,dqds) - hp
      pl = 0.0
      DO i=1,ih
        DO j=1,ih
          pl(i,j)=d_dqds(i)*dfds_d(j)/denom
        END DO
      END DO
      
    RETURN
    END SUBROUTINE
    
    !CURRENT ELASTIC PROPERTIES.
    SUBROUTINE mccdProps(M,e0,P0,lambda,kappa,P,evol,evol_vp,Pc,Young,Poisson)
      
      USE SolidMechanicsLib
      USE NumericalUtilities
    
      IMPLICIT NONE
      REAL(8),INTENT(IN):: M,e0,P0,lambda,kappa
      REAL(8),INTENT(IN):: P,&                      !mean effective stress
                           evol,evol_vp              !volumetric strain times biot
      REAL(8),INTENT(IN OUT):: Young,Poisson,Pc     !updated properties.
      
      REAL(8):: e,ni,ni0,phi,phi0,dPhi,K,G
      
      !Constant parameter.
      !G = Young / 2.0 / (1.0+Poisson)              !constant.
      
      IF(p.NE.0.0)THEN
        
        phi0 = e0/(1.0+e0)
        dPhi = evol
        phi = phi0 + dPhi
        e = phi / (1.0-phi)                         !void ratio
        ni = 1 + e                                  !specific volume
        
        !ni0 = 1 + e0
        
        !IF(lambda>kappa) Pc = P0 + P0 * ni * (-evol_vp) / (lambda-kappa)
        IF(lambda>kappa) Pc = P0 * EXP( ni*(-evol_vp) / (lambda-kappa) )

        K = - ni / kappa * p
        
        !Poisson = (3.0*K-2.0*G)/2.0/(3.0*K+G)
        Young = K/3.0/(1.0-2.0*Poisson)  !9.0*K*G/(3*K+G)
        
      END IF

!if(poisson>0.45)THEN
!  READ(*,*)
!end if
      
    RETURN
    END SUBROUTINE
        
    !Função potencial de escoamento para o critério Cam-Clay Modificado.
    !dq1=dQ/dsigm, dq2=dQ/dJ2, dq3=dQ/dJ3
    SUBROUTINE mcc_q(Pc,M,P,Q,dq1,dq2,dq3)
      ! This subroutine forms the derivatives of a Mohr-Coulomb potential
      ! function with respect to the three invariants (psi in degrees).
  
      IMPLICIT NONE
      REAL(8),INTENT(IN)::Pc,M,P,Q
      REAL(8),INTENT(OUT)::dq1,dq2,dq3
      
      !mccf = M**2 * P**2.0 - M**2.0 * Pc * P + Q**2.0
      dq1 = 2.0 * M**2.0 * P - M**2.0 * Pc
      dq2 = 2.0 * Q
      dq3 = 0.0
      
      RETURN
    END SUBROUTINE mcc_q
  
    !Hardening
    SUBROUTINE mccHardening(Pc0,Evol0,Pc,M,Evol_vp)
      ! This subroutine forms the derivatives of a Mohr-Coulomb potential
      ! function with respect to the three invariants (psi in degrees).
  
      IMPLICIT NONE
      REAL(8),INTENT(IN)::Pc0,Evol0,M,Evol_vp
      REAL(8),INTENT(OUT)::Pc
      
      Pc = Pc0
      
      RETURN
    END SUBROUTINE mccHardening
    
    !Pc intercept
    FUNCTION Pc_intercept(M,P,Q)
      ! This subroutine forms the derivatives of a Mohr-Coulomb potential
      ! function with respect to the three invariants (psi in degrees).
  
      IMPLICIT NONE
      REAL(8)::Pc_intercept
      REAL(8),INTENT(IN)::M,P,Q
      
      REAL(8):: Pc,f1,f2,p1,p2,dp,df,t
      
      !ITERATION 1.
      p1 = P
      f1 = mccf(M,p1,P,Q)
      p2 = 0.99 * P
      f2 = mccf(M,p2,P,Q)
      dp = p2-p1
      df = f2-f1
      t = -df/dp
      
      !INTERATION 2.
      p1 = p2
      f1 = f2
      p2 = p1 + f1/t
      f2 = mccf(M,p2,P,Q)
      dp = p2 - p1
      p2 = p1 + f1/t
      t = -df/dp
      
      !INTERATION 2.
      p1 = p2
      f1 = f2
      p2 = p1 + f1/t
      f2 = mccf(M,p2,P,Q)
      dp = p2 - p1
      p2 = p1 + f1/t
      t = -df/dp
      
      p1 = p2
      f1 = f2
      p2 = p1 + f1/t
      Pc_intercept = p2
      
      RETURN
    END FUNCTION Pc_intercept
    
    !CAM-CLAY: He + Hp
    FUNCTION mcc_HeHp(M,Pc,e0,P0,lambda,kappa,dee,stress)
      
      USE SolidMechanicsLib
      USE NumericalUtilities
    
      IMPLICIT NONE
      REAL(8):: mcc_HeHp
      REAL(8),INTENT(IN):: M,e0,P0,lambda,kappa
      REAL(8),INTENT(IN OUT):: Pc
      REAL(8),INTENT(IN):: stress(:),dee(:,:)
      
      REAL(8),ALLOCATABLE::dfds(:),dqds(:),dpds(:),dJ2ds(:),d_dqds(:),dfds_d(:)
      
      REAL(8):: p, q, lode_theta       !stress inveriants.
      REAL(8):: x
      REAL(8):: dfdp, dfdq, dfdj3
      
      REAL(8)::denom
      INTEGER::i,j,ih
      
      REAL(8):: ni,h(6),hp       !hardening modulus.
      REAL(8):: e
      
      ih = SIZE(stress)
      ALLOCATE(dfds(ih),dqds(ih),dpds(ih),dJ2ds(ih),d_dqds(ih),dfds_d(ih))
      
      CALL invar(stress, p, q, lode_theta)
      
      ! DERIVATIVES.
      !F = M^2 * P^2 - M^2 * P * Pc + Q^2.0
      
      !dF/dP, dF/dQ, dF/dJ3
      dfdp = 2.0*M**2.0*P - M**2.0 * Pc
      dfdq = 2.0 * Q
      dfdj3 = 0.0
      
      !dP/dsigma
      dpds = 0.0
      dpds(1:3) = 1.0/3.0
      
      !dQ/dsigma
      IF(q>0.0)THEN
        x = 2.0 * q**2.0
        dJ2ds(1) = ( 2.0*(stress(1)-stress(2))-2.0*(stress(3)-stress(1)) )
        dJ2ds(2) = ( 2.0*(stress(2)-stress(3))-2.0*(stress(1)-stress(2)) )
        dJ2ds(3) = ( 2.0*(stress(3)-stress(1))-2.0*(stress(2)-stress(3)) )
        dJ2ds(4) = ( 12.0*stress(4) )
        dJ2ds(5) = ( 12.0*stress(5) )
        dJ2ds(6) = ( 12.0*stress(6) )
        dJ2ds = dJ2ds * 1.0/SQRT(2.0) * 1.0/2.0 / SQRT(x)
      ELSE
        dJ2ds = 0.0
      END IF
      
      !dF/ds
      dfds = dfdp*dpds + dfdq*dJ2ds     !<- dFunction/ds
      dqds = dfds                       !<- dPotencial/ds
      
      !HARDENING.
      h = 0.0
      hp = 0.0
      
      !F = M^2 * P^2 - M^2 * Pc * P + Q^2                                          
      !hp = - dF/dK * {dK/devol_vp} . {dF/dp}
      ni = 1+e0                                 !specific volume
      ni = ni - lambda*LOG(Pc/p0)               !Pc updated
      !h(1:3) = -p * (-P0*ni/(lambda-kappa))
      !hp = DOT_PRODUCT(h,dqds)
      IF(lambda.NE.kappa) hp = (-M**2.0*P) * (ni*Pc/(lambda-kappa)) * ( 2*M**2*P-M**2*Pc )
      
      !FINAL EQUATION.
      !{dFds}.[D].{dQds} + hp )
      d_dqds = MATMUL(dee,dqds)
      dfds_d = MATMUL(dee,dfds) 
      mcc_HeHp = DOT_PRODUCT(d_dqds,dqds) + hp
      
    RETURN
    END FUNCTION
    
    !*************************************************************************
    !ELASTIC RELATIONS.
    !*************************************************************************

    ! ELASTIC STRESS CALCULATION FROM STRAIN TENSOR - HOMOGENIOUS MATERIALS.
    FUNCTION StressElastic(Young,Poisson,Strain)
    
      IMPLICIT NONE
      REAL(8):: StressElastic(6)
      REAL(8),INTENT(IN):: Young,Poisson
      REAL(8),INTENT(IN):: Strain(6)
      
      REAL(8):: dee(6,6)
      
      CALL D_elastic_homogenious(dee,Young,Poisson)
      
      StressElastic = MATMUL(dee,Strain)
    
    RETURN
    END FUNCTION

    ! [D] MATRIX FOR ELASTIC HOMOGENIOUS MATERIALS.
    SUBROUTINE D_elastic_homogenious(dee,e,v)
      
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(8),INTENT(IN)::e,v
      REAL(8),INTENT(IN OUT)::dee(:,:)
      REAL(8)::v1,v2,c,vv,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp
      INTEGER::i
      
      !QUALITY.
      IF((v.GT.0.5_iwp).OR.(v.LT.zero))THEN
        WRITE(*,*)"ERRO NA MONTAGEM DA MATRIZ ELASTICA [D]"
        WRITE(*,*)"POISSON FORA DO LIMITE ACEITAVEL {0.0 , 0.5}"
        WRITE(*,*)"subrotina: D_HomogeneousElasticMatrix (D_HomogeneousElasticMatrix.dll)"
        WRITE(*,*)"FINALIZANDO..."
        STOP
      END IF
      !IF(e.LT.zero)THEN
      !  WRITE(*,*)"ERRO NA MONTAGEM DA MATRIZ ELASTICA [D]"
      !  WRITE(*,*)"MODULO DE YOUNG NEGATIVO"
      !  WRITE(*,*)"subrotina: D_HomogeneousElasticMatrix (D_HomogeneousElasticMatrix.dll)"
      !  WRITE(*,*)"FINALIZANDO..."
      !  STOP
      !END IF
      
      ![D] MATRIX.
      dee=zero
      v1=one-v
      c=e/((one+v)*(one-two*v))
      v2=v/(one-v)
      vv=(one-two*v)/(one-v)*pt5
      DO i=1,3
        dee(i,i)=one
      END DO
      DO i=4,6
        dee(i,i)=vv
      END DO
      dee(1,2)=v2
      dee(2,1)=v2
      dee(1,3)=v2
      dee(3,1)=v2
      dee(2,3)=v2
      dee(3,2)=v2
      dee=dee*e/(two*(one+v)*vv)
      
    RETURN
    END SUBROUTINE
   
    ! ELASTIC STRESS CALCULATION FROM STRAIN TENSOR - ANYSOTROPIC MATERIALS. 
    FUNCTION StressElasticOrtoYY(Young,Poisson,YoungNormal,Strain)
    
      IMPLICIT NONE
      REAL(8):: StressElasticOrtoYY(6)
      REAL(8),INTENT(IN):: Young,YoungNormal,Poisson
      REAL(8),INTENT(IN):: Strain(6)
      
      REAL(8):: dee(6,6)
      REAL(8):: angle = 90.0
      
      !***********************************************************************      
      ! Módulo de Young, no plano de fraqueza.                               !
      ! Módulo de Young, normal ao plano de fraqueza                         !
      ! Poisson para deformação no plano de fraqueza                         !
      !     devido às tensões atuantes na normal ao plano de fraqueza.       !
      ! Poisson para deformação no plano de deposição                        !
      !     devido às tensões atuantes no mesmo plano.                       !
      ! Ângulo azimute da estrutura geológica planar (strike)                !
      !     em relação ao eixo y (norte) e no sentido horário                !
      ! Ângulo de inclinação da estrutura geológica planar (dip)             !           
      !     em relação ao plano horizontal xy e no sentido horário           !
      ! Matriz constitutiva elástica global                                  !
      !***********************************************************************
      CALL D_HomogeneousElasticMatrixTransvIso_PUCtecgraf(Young,YoungNormal,&
                                                Poisson,Poisson,angle,angle,dee)

      StressElasticOrtoYY = MATMUL(dee,Strain)
    
    RETURN
    END FUNCTION
    
    ! INTERFACE FOR MATRIX FOR ELASTIC TRANSVERSE ISOTROPIC MATERIALS (D_HomogeneousElasticMatrixTransvIso_PUCtecgraf).
    FUNCTION D_TransverseIsotropic(E1,E2,N1,N2,STRIKE,DIP)
        
        IMPLICIT NONE
        REAL(8):: D_TransverseIsotropic(6,6)
        REAL(8),INTENT(IN)::E1,E2,N1,N2,STRIKE,DIP
        
        REAL(8) :: Ea, Ec, niac, niab, sf, df
        
        Ea = E2     !no plano
        Ec = E1     !normal ao plano
        niac = N2   !no plano
        niab = N1   !normal ao plano
        sf = STRIKE !azimute
        df = DIP    !dip
        
        CALL D_HomogeneousElasticMatrixTransvIso_PUCtecgraf( Ea, Ec, niac, niab, sf, df, D_TransverseIsotropic )
        
    END FUNCTION
    
    ! [D] MATRIX FOR ELASTIC TRANSVERSE ISOTROPIC MATERIALS.
    SUBROUTINE D_HomogeneousElasticMatrixTransvIso_PUCtecgraf( Ea, Ec, niac, niab, sf, df, MG )
  
      !***********************************************************************
      !																	     !
      ! Tecnologia em Computação Gráfica - TECGRAF                           !
      ! Pontifícia Universidade Católica do Rio de Janeiro - PUC-Rio         !
      !                                                                      !
      ! Grupo de Geomecânica Computacional                                   !
      ! Coordenadora: Deane Roehl                                            !
      ! Contato: geomec@tecgraf.puc-rio.br                                   !
      !                                                                      !
      !***********************************************************************
      !                                                                      !
      ! Rotina:  rot_mat_const                                               !
      ! Arquivo: rot_mat_const.f90 										     !
      ! Data:    26/09/2011												     !
      ! Autor:   Roberto J. Q. Quispe 									     !
      !																	     !
      ! Modificações													     !
      !   Data:															     !
      !																	     !
      !***********************************************************************
      ! 																	 !
      ! Descrição: Rotacão da matriz constitutiva elástica para um           !
      !            um meio geologico transversalmente isotrópico com  	     !
      !            vetor de tensões: (Sxx, Syy, Szz, Sxy, Sxz, Syz)          !
      !                                                                      !
      ! Ref: Kim, J-M. Vectorial formulation of direction cosines for        !
      !       anisotropic geologic structures from their geologic angle      !
      !       measurements. Mathematical Geology, 37 (8), 2005.              !
      ! 																	 !
      ! Argumentos:														     !
      !   Ea       - Módulo de Young, no plano de deposição ab               !
      !   Ec       - Módulo de Young, na direção de deposição c              !
      !   niac     - Módulo de Poisson para deformação no plano de deposição !
      !              devido às tensões atuantes na direção de deposição      !
      !   niab     - Módulo de Poisson para deformação no plano de deposição !
      !              devido às tensões atuantes no mesmo plano               !
      !   sf[DEG]  - Ângulo azimute da estrutura geológica planar (strike)   !
      !              em relação ao eixo y (norte) e no sentido horário       !
      !	df[DEG]	 - Ângulo de inclinação da estrutura geológica planar (dip)  ! 
      !              em relação ao plano horizontal xy e no sentido horário  !
      !   ML       - Matriz constitutiva elástica local                      !
      !   R        - Matriz de rotação do sistema local para o global        !
      !   MG       - Matriz constitutiva elástica global                     !
      !***********************************************************************
      IMPLICIT NONE

      ! Variáveis
      REAL(8), intent(in) :: Ea, Ec, niac, niab, sf, df
      REAL(8) :: sf_rad, df_rad
      REAL(8) :: Eb, nibc, Gab, Gbc, Gca
      REAL(8) :: PI, ZER, ONE, TWO, CL
      ! Matriz constitutiva elástica local
      REAL(8), dimension(6,6) :: ML, R
      ! Vetores unitários e vetores auxiliares
      REAL(8), dimension(3) :: a, b, c
      ! Matriz constitutiva elástica global 
      REAL(8), dimension(6,6), intent(out) :: MG
  
      ! Definição de constantes numéricas
      PI  = DACOS(-1.00D0)
      ZER = 0.00D0
      ONE = 1.00D0
      TWO = 2.00D0
 
      ! Transformação dos ângulos de graus para radianos
      sf_rad = sf*PI/180.00D0 !strike
      df_rad = df*PI/180.00D0 !dip
 
      ! Constantes para material transversalmente isotrópico
      Eb   = Ea
      nibc = niac
      Gca  = 0.25D0*(Ea+Ec)/(ONE+niac) !Aproximação
      Gbc  = Gca
      Gab  = Ea/(TWO*(ONE+niab))

      ! Montagem da matriz constitutiva elástica local
      CL = ONE/(Ea*Eb*Ec) - niab**TWO/(Ea**TWO*Eb) - nibc**TWO/(Ea*Eb**TWO) &
                          - niac**TWO/(Ea**TWO*Eb) - TWO*niab*nibc*niac/(Ea**TWO*Eb)

      ML(1,1) = (ONE/(Eb*Ec)-(nibc/Eb)**TWO)*CL**(-ONE)
      ML(1,2) = (niab/(Ea*Ec)+niac*nibc/(Ea*Eb))*CL**(-ONE)
      ML(1,3) = (niac/(Ea*Eb)+niab*nibc/(Ea*Eb))*CL**(-ONE)
      ML(1,4) = ZER
      ML(1,5) = ZER
      ML(1,6) = ZER
      ML(2,1) = ML(1,2)
      ML(2,2) = (ONE/(Ea*Ec)-(niac/Ea)**TWO)*CL**(-ONE)
      ML(2,3) = (niab*niac/Ea**TWO + nibc/(Ea*Eb))*CL**(-ONE)
      ML(2,4) = ZER
      ML(2,5) = ZER
      ML(2,6) = ZER
      ML(3,1) = ML(1,3)
      ML(3,2) = ML(2,3)
      ML(3,3) = (ONE/(Ea*Eb)-(niab/Ea)**TWO)*CL**(-ONE)
      ML(3,4) = ZER
      ML(3,5) = ZER
      ML(3,6) = ZER
      ML(4,1) = ML(1,4)
      ML(4,2) = ML(2,4)
      ML(4,3) = ML(3,4)
      ML(4,4) = Gab
      ML(4,5) = ZER
      ML(4,6) = ZER
      ML(5,1) = ML(1,5)
      ML(5,2) = ML(2,5)
      ML(5,3) = ML(3,5)
      ML(5,4) = ML(4,5)
      ML(5,5) = Gca
      ML(5,6) = ZER
      ML(6,1) = ML(1,6)
      ML(6,2) = ML(2,6)
      ML(6,3) = ML(3,6)
      ML(6,4) = ML(4,6)
      ML(6,5) = ML(5,6)
      ML(6,6) = Gbc

      ! Vetores dos cossenos diretores entre o sistema local e global
      a = (/DCOS(sf_rad)*DCOS(df_rad), -DSIN(sf_rad)*DCOS(df_rad), -DSIN(df_rad)/) 
      b = (/DSIN(sf_rad), DCOS(sf_rad), 0.00D0/)
      c = (/DCOS(sf_rad)*DSIN(df_rad), -DSIN(sf_rad)*DSIN(df_rad), DCOS(df_rad)/)

      !Matriz de rotação 
      R(1,:) = (/a(1)**TWO,     a(2)**TWO,     a(3)**TWO,     a(1)*a(2),           a(3)*a(1),           a(2)*a(3)          /)
      R(2,:) = (/b(1)**TWO,     b(2)**TWO,     b(3)**TWO,     b(1)*b(2),           b(3)*b(1),           b(2)*b(3)          /)
      R(3,:) = (/c(1)**TWO,     c(2)**TWO,     c(3)**TWO,     c(1)*c(2),           c(3)*c(1),           c(2)*c(3)          /)
      R(4,:) = (/TWO*a(1)*b(1), TWO*a(2)*b(2), TWO*a(3)*b(3), a(1)*b(2)+b(1)*a(2), a(3)*b(1)+b(3)*a(1), a(2)*b(3)+b(2)*a(3)/)
      R(5,:) = (/TWO*c(1)*a(1), TWO*c(2)*a(2), TWO*c(3)*a(3), c(1)*a(2)+a(1)*c(2), c(3)*a(1)+a(3)*c(1), c(2)*a(3)+a(2)*c(3)/)
      R(6,:) = (/TWO*b(1)*c(1), TWO*b(2)*c(2), TWO*b(3)*c(3), b(1)*c(2)+c(1)*b(2), b(3)*c(1)+c(3)*b(1), b(2)*c(3)+c(2)*b(3)/)

      !Matriz constitutiva elástica global 
      MG =  MATMUL(TRANSPOSE(R),MATMUL(ML,R))

      RETURN

    END SUBROUTINE
    
    ! FLUID HYDROSTATIC STRESS FROM STRAIN TENSOR.
    FUNCTION FluidHydrostaticStress(Kf,BiotCoef,Strain)
    
      IMPLICIT NONE
      REAL(8):: FluidHydrostaticStress(6)
      REAL(8),INTENT(IN):: Kf, BiotCoef
      REAL(8),INTENT(IN):: Strain(6)
      
      REAL(8):: dee(6,6)
      
      dee = 0.0
      dee(1:3,1:3) = Kf

      FluidHydrostaticStress = BiotCoef * MATMUL(dee,Strain)
    
    RETURN
    END FUNCTION


    !*************************************************************************
    !AUXILIARES FUNCTIONS.
    !*************************************************************************
    
    !Smith e Griffiths' routine.
    SUBROUTINE formm(stress,m1,m2,m3)
       IMPLICIT NONE
       INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
       REAL(iwp),INTENT(IN)::stress(:)
       REAL(iwp),INTENT(OUT)::m1(:,:),m2(:,:),m3(:,:)
       
       REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
         two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
       INTEGER::i,j
       
      sx=stress(1)
      sy=stress(2)    
      sz=stress(3)
      txy=stress(4)  
      tyz=stress(5) 
      tzx=stress(6)
      sigm=(sx+sy+sz)/three
      dx=sx-sigm  
      dy=sy-sigm 
      dz=sz-sigm
      m1=zero
      m2=zero
      m1(1:3,1:3)=one/(three*sigm)
      DO i=1,3 
        m2(i,i)=two 
        m2(i+3,i+3)=six 
      END DO
      m2(1,2)=-one
      m2(1,3)=-one 
      m2(2,3)=-one
      m3(1,1)=dx
      m3(1,2)=dz 
      m3(1,3)=dy 
      m3(1,4)=txy  
      m3(1,5)=-two*tyz
      m3(1,6)=tzx 
      m3(2,2)=dy 
      m3(2,3)=dx 
      m3(2,4)=txy
      m3(2,5)=tyz 
      m3(2,6)=-two*tzx 
      m3(3,3)=dz
      m3(3,4)=-two*txy
      m3(3,5)=tyz 
      m3(3,6)=tzx
      m3(4,4)=-three*dz 
      m3(4,5)=three*tzx
      m3(4,6)=three*tyz
      m3(5,5)=-three*dx
      m3(5,6)=three*txy 
      m3(6,6)=-three*dy
      DO i=1,6 
        DO j=i+1,6
          m1(j,i)=m1(i,j) 
          m2(j,i)=m2(i,j)   
          m3(j,i)=m3(i,j)
        END DO
      END DO
      m1=m1/three
      m2=m2/three
      m3=m3/three
    RETURN   
    END SUBROUTINE formm
    
    FUNCTION GetSoilProp(SoilModels,SoilName,PropName)
    
      IMPLICIT NONE
      REAL(8):: GetSoilProp
      TYPE(SoilModel),INTENT(IN):: SoilModels(:)
      CHARACTER(LEN=20),INTENT(IN):: SoilName
      CHARACTER(LEN=*),INTENT(IN):: PropName
      
      INTEGER:: i
      
      GetSoilProp = 0.0
      
      DO i=1,SIZE(SoilModels)
        IF(SoilModels(i)%Name==SoilName)THEN
          IF(PropName=='density') GetSoilProp = SoilModels(i)%density
          IF(PropName=='poisson') GetSoilProp = SoilModels(i)%poisson
          EXIT
        END IF
      END DO
      
    END FUNCTION
    
    FUNCTION GetSoilPropNumber(SoilModels,SoilName)
    
      IMPLICIT NONE
      REAL(8):: GetSoilPropNumber
      TYPE(SoilModel),INTENT(IN):: SoilModels(:)
      CHARACTER(LEN=20),INTENT(IN):: SoilName
      
      INTEGER:: i
      
      GetSoilPropNumber = 0.0
      
      DO i=1,SIZE(SoilModels)
        IF(SoilModels(i)%Name==SoilName)THEN
          GetSoilPropNumber=i
          EXIT
        END IF
      END DO
      
    END FUNCTION
    
  END MODULE