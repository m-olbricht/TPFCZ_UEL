!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! viscous dissipation energy module
!
! Martin Olbricht, TU Bergakademie Freiberg, 14.03.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



MODULE ViscousDissipationModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS
  
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VDED_GL(phase_rate,dtime,nVDEDpar,parVDEDMatrixPhase)
    ! viscous dissipation energy
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase_rate
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: eta
      !
      !
      eta              = parVDEDMatrixPhase(1)
      !
      VDED_GL = zero
      !
      VDED_GL = half*eta*phase_rate**two
      
      END FUNCTION VDED_GL

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_GL_d_phase(phase_rate,dtime,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase_rate
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: eta
      !
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_VDED_GL_d_phase = zero
      !
      d_VDED_GL_d_phase = eta * phase_rate/dtime
      
      
      END FUNCTION d_VDED_GL_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_GL_d_phase_d_phase(phase_rate,dtime,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase_rate
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: eta
      !
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_VDED_GL_d_phase_d_phase = zero
      !
      d_VDED_GL_d_phase_d_phase = eta / dtime**two
      
      
      END FUNCTION d_VDED_GL_d_phase_d_phase

!------------------------------------------------------------------------------------  

! ---------------------------- Modified Viscous Energy ------------------------------

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VDED_mod(phase_rate,dtime,nVDEDpar,parVDEDMatrixPhase)
    ! viscous dissipation energy
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase_rate
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: eta
      !
      !
      eta              = parVDEDMatrixPhase(1)
      !
      VDED_mod = zero
      !
      VDED_mod = half*eta*dtime*phase_rate**two
      
      END FUNCTION VDED_mod

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_mod_d_phase(phase_rate,dtime,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase_rate
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: eta
      !
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_VDED_mod_d_phase = zero
      !
      d_VDED_mod_d_phase = eta * phase_rate
      
      
      END FUNCTION d_VDED_mod_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_mod_d_phase_d_phase(phase_rate,dtime,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase_rate
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: eta
      !
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_VDED_mod_d_phase_d_phase = zero
      !
      d_VDED_mod_d_phase_d_phase = eta / dtime
      
      
      END FUNCTION d_VDED_mod_d_phase_d_phase

!------------------------------------------------------------------------------------

END MODULE ViscousDissipationModule
