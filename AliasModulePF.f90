MODULE AliasModulePF

!  Alias Modul
!  
!  Speichert alle Umbenennungen, damit sie von verschiedenen Subroutinen aufgerufen werden können

   ! ==========================================================
   ! Viscous Dissipation PF
   ! ==========================================================
!~    USE ViscousDissipationModule, ONLY: VDED => VDED_GL, &
!~                                                d_VDED_d_damage => d_VDED_GL_d_phase, &
!~                                                d_VDED_d_damage_d_damage => d_VDED_GL_d_phase_d_phase

   USE ViscousDissipationModule, ONLY: VDED => VDED_mod, &
                                               d_VDED_d_damage => d_VDED_mod_d_phase, &
                                               d_VDED_d_damage_d_damage => d_VDED_mod_d_phase_d_phase

   ! ==========================================================
   ! Interface Energy (AT2 aktiv)
   ! Zum Wechseln auf AT1 einfach hier umstellen
   ! ==========================================================
   USE InterfaceEnergyModule, ONLY: &
   
        IED  => IED_AT2, &
        d_IED_d_damage => d_IED_AT2_d_damage, &
        d_IED_d_damage_d_damage => d_IED_AT2_d_damage_d_damage, &
        d_IED_d_grad_damage => d_IED_AT2_d_grad_damage, &
        d_IED_d_grad_damage_d_grad_damage => &
        d_IED_AT2_d_grad_damage_d_grad_damage, &
        d_IED_d_damage_d_grad_damage => &
        d_IED_AT2_d_damage_d_grad_damage
        
        
   ! ==========================================================
   ! Degradation Function
   ! ==========================================================

!~ 		USE DegradationModule, ONLY: Degradation => DegradationFlex, &
!~                                               d_Degradation_d_damage => d_DegradationFlex_d_damage, &
!~                                               d_Degradation_d_damage_d_damage => d_DegradationFlex_d_damage_d_damage 
                                              
		USE DegradationModule, ONLY: Degradation => DegradationQuadratic, &
                                              d_Degradation_d_damage => d_DegradationQuadratic_d_damage, &
                                              d_Degradation_d_damage_d_damage => d_DegradationQuadratic_d_damage_d_damage 
                                               
   ! ==========================================================
   ! Energy Split
   ! ==========================================================
   
   
  USE SplitModulePF, ONLY: HFEDpos => HFEDposNoSplit, &
                                               d_HFEDpos_d_eps_e => d_HFEDposNoSplit_d_eps_e, &
                                               d_HFEDpos_d_eps_e_d_eps_e => d_HFEDposNoSplit_d_eps_e_d_eps_e, &
                                               HFEDneg => HFEDnegNoSplit, &
                                               d_HFEDneg_d_eps_e => d_HFEDnegNoSplit_d_eps_e, &
                                               d_HFEDneg_d_eps_e_d_eps_e => d_HFEDnegNoSplit_d_eps_e_d_eps_e
  USE SplitModulePF, ONLY: eps_elast 



        
   IMPLICIT NONE
   
   PUBLIC :: VDED, d_VDED_d_damage, d_VDED_d_damage_d_damage
   PUBLIC :: IED, d_IED_d_damage, d_IED_d_damage_d_damage
   PUBLIC :: Degradation, d_Degradation_d_damage, d_Degradation_d_damage_d_damage
   PUBLIC :: HFEDpos, d_HFEDpos_d_eps_e, d_HFEDpos_d_eps_e_d_eps_e, HFEDneg, d_HFEDneg_d_eps_e, d_HFEDneg_d_eps_e_d_eps_e
   PUBLIC :: eps_elast
   

END MODULE AliasModulePF
