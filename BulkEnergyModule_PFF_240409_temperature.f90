!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! bulk energy module
!
! Martin Olbricht, TU Bergakademie Freiberg, 12.01.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BulkEnergyModule

USE ABQINTERFACE_PF
USE FLOATNUMBERS
USE FreeEnergyModule
USE AliasModulePF

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS


!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_eps_elast_d_eps()
    ! elastic strain w.r.t. strain
    
	  USE TensorModule
      
      IMPLICIT NONE
      !
      !
      DIMENSION d_eps_elast_d_eps(3,3,3,3)
      !
      !
      d_eps_elast_d_eps = SymIdent4th()
      
      
    END FUNCTION d_eps_elast_d_eps
    
!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_eps_elast_d_temperature(eps,dTemperature,nThermalpar,parThermalMatrixPhase)
    ! elastic strain w.r.t. strain
      
      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      !
      REAL(kind=AbqRK) :: therm_exp
      INTEGER(kind=AbqIK) :: i1
      !
      DIMENSION :: d_eps_elast_d_temperature(3,3)
      !
      !
      therm_exp = parThermalMatrixPhase(5)
      !
      FORALL (i1=1:3)
        d_eps_elast_d_temperature(i1,i1) = -therm_exp
      END FORALL
      
      
    END FUNCTION d_eps_elast_d_temperature

!------------------------------------------------------------------------------------!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION HFEDtens_H(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! bulk energy density
      
      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      !
      HFEDtens_H = HFEDpos(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      
      
    END FUNCTION HFEDtens_H
    
!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION bulkED(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! bulk energy density

      
      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      REAL(kind=AbqRK) :: HFEDtens, HFEDcomp
      REAL(kind=AbqRK) :: degD
      !
	  !
      HFEDtens = HFEDpos(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      HFEDcomp = HFEDneg(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      bulkED = degD*HFEDtens + HFEDcomp
      
      
    END FUNCTION bulkED
    
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               first derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. strain
    
      USE TensorModule

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3)
      REAL(kind=AbqRK) :: d_HFEDcomp_d_eps_e(3,3)
      REAL(kind=AbqRK) :: degD
      !
      DIMENSION d_bulkED_d_eps(3,3)
      !
      d_eps_e_d_eps = d_eps_elast_d_eps()
      d_eps_e_d_temperature = d_eps_elast_d_temperature(eps,dTemperature,nThermalpar,parThermalMatrixPhase)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_eps = degD*DOUBLECONTRACTIONTwoFour(d_HFEDtens_d_eps_e,d_eps_e_d_eps) + DOUBLECONTRACTIONTwoFour(d_HFEDcomp_d_eps_e,d_eps_e_d_eps)
      
      
    END FUNCTION d_bulkED_d_eps

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_damage_H(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. damage

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      REAL(kind=AbqRK) :: d_degD_d_damage_H, HFEDtens
      !
      HFEDtens = H
      d_degD_d_damage_H = d_Degradation_d_damage(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_damage_H = d_degD_d_damage_H*HFEDtens
      
      
    END FUNCTION d_bulkED_d_damage_H
    

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_temperature(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. temperature
	  
	  USE TensorModule

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3)
      REAL(kind=AbqRK) :: d_HFEDcomp_d_eps_e(3,3)
      REAL(kind=AbqRK) :: degD
      !
      DIMENSION d_bulkED_d_temperature(3,3)
      !
      !
      d_eps_e_d_temperature = d_eps_elast_d_temperature(eps,dTemperature,nThermalpar,parThermalMatrixPhase)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_temperature = degD*DOUBLECONTRACTIONTwoTwo(d_HFEDtens_d_eps_e,d_eps_e_d_temperature) + DOUBLECONTRACTIONTwoTwo(d_HFEDcomp_d_eps_e,d_eps_e_d_temperature)
      
      
    END FUNCTION d_bulkED_d_temperature

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------!
!                               second derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!


    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps_d_eps(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! doubled partial derivative of bulk energy density w.r.t. strain strain

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      !
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e_d_eps_e(3,3,3,3)
      REAL(kind=AbqRK) :: d_HFEDcomp_d_eps_e_d_eps_e(3,3,3,3)
      !
      REAL(kind=AbqRK) :: degD
      !
      DIMENSION d_bulkED_d_eps_d_eps(3,3,3,3)
      !
      !
      d_HFEDtens_d_eps_e_d_eps_e = d_HFEDpos_d_eps_e_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      d_HFEDcomp_d_eps_e_d_eps_e = d_HFEDneg_d_eps_e_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      !
      d_bulkED_d_eps_d_eps = degD*d_HFEDtens_d_eps_e_d_eps_e + d_HFEDcomp_d_eps_e_d_eps_e
      

    END FUNCTION d_bulkED_d_eps_d_eps

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps_d_damage(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. strain, damage
	  
	  USE TensorModule

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3)
      REAL(kind=AbqRK) :: d_HFEDcomp_d_eps_e(3,3)
      REAL(kind=AbqRK) :: d_degD_d_damage
      !
      DIMENSION d_bulkED_d_eps_d_damage(3,3)
      !
      !
      d_eps_e_d_eps = d_eps_elast_d_eps()
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      !
      d_degD_d_damage = d_Degradation_d_damage(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_eps_d_damage = d_degD_d_damage*DOUBLECONTRACTIONTwoFour(d_HFEDtens_d_eps_e,d_eps_e_d_eps)
      
      
    END FUNCTION d_bulkED_d_eps_d_damage

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps_d_temperature(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. strain, temperature
	  
	  USE TensorModule

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      !
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e_d_eps_e(3,3,3,3)
      REAL(kind=AbqRK) :: d_HFEDcomp_d_eps_e_d_eps_e(3,3,3,3)
      !
      REAL(kind=AbqRK) :: degD
      !
      DIMENSION d_bulkED_d_eps_d_temperature(3,3)
      !
      !
      d_eps_e_d_eps = d_eps_elast_d_eps()
      d_eps_e_d_temperature = d_eps_elast_d_temperature(eps,dTemperature,nThermalpar,parThermalMatrixPhase)
      !
      d_HFEDtens_d_eps_e_d_eps_e = d_HFEDpos_d_eps_e_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      d_HFEDcomp_d_eps_e_d_eps_e = d_HFEDneg_d_eps_e_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_eps_d_temperature = degD*DOUBLECONTRACTIONTwoTwo(d_HFEDtens_d_eps_e_d_eps_e,d_eps_e_d_temperature) + DOUBLECONTRACTIONTwoTwo(d_HFEDtens_d_eps_e_d_eps_e,d_eps_e_d_temperature)
      
      
    END FUNCTION d_bulkED_d_eps_d_temperature

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_damage_d_damage_H(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! partial derivative of the bulk energy w.r.t. damage and damage
    !

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      REAL(kind=AbqRK) :: d_degD_d_damage_d_damage_H, HFEDtens
      !
      HFEDtens = H
      d_degD_d_damage_d_damage_H = d_Degradation_d_damage_d_damage(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_damage_d_damage_H = d_degD_d_damage_d_damage_H * HFEDtens

    END FUNCTION d_bulkED_d_damage_d_damage_H

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_damage_d_temperature(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. strain

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      !
      REAL(kind=AbqRK) :: degD
      !
      d_eps_e_d_eps = d_eps_elast_d_eps()
      d_eps_e_d_temperature = d_eps_elast_d_temperature(eps,dTemperature,nThermalpar,parThermalMatrixPhase)
      !
      !
      d_bulkED_d_damage_d_temperature = zero
      
      
    END FUNCTION d_bulkED_d_damage_d_temperature

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_temperature_d_temperature(eps,damage,dTemperature,nThermalpar,parThermalMatrixPhase,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. strain
	  
	  USE TensorModule

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      INTEGER(kind=AbqIK), INTENT(IN) :: nThermalpar
      REAL(kind=AbqRK), INTENT(IN) :: parThermalMatrixPhase(nThermalpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: dTemperature
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_betaPF
      REAL(kind=AbqRK), INTENT(IN) :: H
      !
      REAL(kind=AbqRK) :: d_eps_e_d_eps(3,3,3,3)
      REAL(kind=AbqRK) :: d_eps_e_d_temperature(3,3)
      !
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e_d_eps_e(3,3,3,3)
      REAL(kind=AbqRK) :: d_HFEDcomp_d_eps_e_d_eps_e(3,3,3,3)
      !
      REAL(kind=AbqRK) :: degD, temp1, temp2
      !
      d_eps_e_d_eps = d_eps_elast_d_eps()
      d_eps_e_d_temperature = d_eps_elast_d_temperature(eps,dTemperature,nThermalpar,parThermalMatrixPhase)
      !
      !
      d_HFEDtens_d_eps_e_d_eps_e = d_HFEDpos_d_eps_e_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      d_HFEDcomp_d_eps_e_d_eps_e = d_HFEDneg_d_eps_e_d_eps_e(eps,dTemperature,nHFEDpar,parHFEDMatrixPhase,nThermalpar,parThermalMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      ! Term 1
      temp1 = DOUBLECONTRACTIONSTwoFourTwo(d_eps_e_d_temperature,d_HFEDtens_d_eps_e_d_eps_e,d_eps_e_d_temperature)
      ! Term 2
      temp2 = DOUBLECONTRACTIONSTwoFourTwo(d_eps_e_d_temperature,d_HFEDcomp_d_eps_e_d_eps_e,d_eps_e_d_temperature)
      
      d_bulkED_d_temperature_d_temperature = degD*temp1 + temp2
      
      
    END FUNCTION d_bulkED_d_temperature_d_temperature

!------------------------------------------------------------------------------------!

END MODULE BulkEnergyModule
