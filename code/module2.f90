	module moduletwo
	
	implicit none

	real(KIND=8), parameter :: Distance_nearestneighbour_MoS2=0.182442685, Distance_nearestneighbour_WS2=0.182038540	! nm
	
	integer(KIND=4) :: SetTransitionMetalDichalcogenides=1         !  1 --- Mos2,  2----Ws2

	
	
!*****************************************************************************************************	
	
	!~ real(KIND=8), parameter :: Delta0_MoS2=-1.512, Delta1_MoS2=0., Delta2_MoS2=-3.025, DeltaP_MoS2=-1.276, DeltaZ_MoS2=-8.236
	!~ !Mo-S
	!~ real(KIND=8), parameter :: V_pd_sigma_MoS2=-2.619, V_pd_pi_MoS2=-1.396
	!~ !Mo-Mo
	!~ real(KIND=8), parameter :: V_dd_sigma_MoS2=-0.933, V_dd_pi_MoS2=-0.478, V_dd_delta_MoS2=-0.442
	!~ !S-S
	!~ real(KIND=8), parameter :: V_pp_sigma_MoS2=0.696, V_pp_pi_MoS2=0.278
	
	!~ real(KIND=8), parameter :: Delta0_MoS2=-1.096, Delta1_MoS2=0., Delta2_MoS2=-1.512, DeltaP_MoS2=-3.560, DeltaZ_MoS2=-6.886
	!~ !Mo-S
	!~ real(KIND=8), parameter :: V_pd_sigma_MoS2=3.688, V_pd_pi_MoS2=-1.241
	!~ !Mo-Mo
	!~ real(KIND=8), parameter :: V_dd_sigma_MoS2=-0.895, V_dd_pi_MoS2=0.252, V_dd_delta_MoS2=0.228
	!~ !S-S
	!~ real(KIND=8), parameter :: V_pp_sigma_MoS2=1.225, V_pp_pi_MoS2=-0.467
	!SOC
	!~ real(KIND=8), parameter :: SOC_dd_MoS2=0.075, SOC_pp_MoS2=0.052
	!~ !S-S Interlayer
	!~ real(KIND=8), parameter :: V_pp_sigma_interlayer_MoS2=-0.774, V_pp_pi_interlayer_MoS2=0.123
	
	real(KIND=8), parameter :: Delta0_MoS2=-1.09353d0, Delta1_MoS2=5.5d0, Delta2_MoS2=-1.51187d0,&
    &DeltaP_MoS2=-3.55909d0, DeltaZ_MoS2=-6.88559d0
	!Mo-S
	real(KIND=8), parameter :: V_pd_sigma_MoS2=3.68886d0, V_pd_pi_MoS2=-1.24057d0
	!Mo-Mo
	real(KIND=8), parameter :: V_dd_sigma_MoS2=-0.895078d0, V_dd_pi_MoS2=0.252318d0, V_dd_delta_MoS2=0.228446d0
	!S-S
	real(KIND=8), parameter :: V_pp_sigma_MoS2=1.22524d0, V_pp_pi_MoS2=-0.467313d0
	!SOC
	real(KIND=8), parameter :: SOC_dd_MoS2=0.075d0, SOC_pp_MoS2=0.05d0
	!S-S Interlayer
	real(KIND=8), parameter :: V_pp_sigma_interlayer_MoS2=-0.774d0, V_pp_pi_interlayer_MoS2=0.123d0
!------------------------------------------------------------------------------------------------------------------------------------------------------------------	

!	real(KIND=8), parameter :: Delta0_MoS2=0.0, Delta1_MoS2=0.0, Delta2_MoS2=0.0, DeltaP_MoS2=0.0, DeltaZ_MoS2=0.0
	!Mo-S
!	real(KIND=8), parameter :: V_pd_sigma_MoS2=0.0, V_pd_pi_MoS2=0.0
	!Mo-Mo
!	real(KIND=8), parameter :: V_dd_sigma_MoS2=0.0, V_dd_pi_MoS2=0.0, V_dd_delta_MoS2=0.0
	!S-S
!	real(KIND=8), parameter :: V_pp_sigma_MoS2=0.0, V_pp_pi_MoS2=0.0
	!SOC
!	real(KIND=8), parameter :: SOC_dd_MoS2=0.0, SOC_pp_MoS2=0.0
	!S-S Interlayer
!	real(KIND=8), parameter :: V_pp_sigma_interlayer_MoS2=0.0, V_pp_pi_interlayer_MoS2=0.0
	
	
!------------------------------------------------------------------------------------------------------------------------------------------------------------------	
	
	real(KIND=8), parameter :: Delta0_WS2=-0.872d0, Delta1_WS2=0.42d0, Delta2_WS2=-2.065d0, DeltaP_WS2=-3.468d0, DeltaZ_WS2=-3.913d0
	!W-S
	real(KIND=8), parameter :: V_pd_sigma_WS2=3.603d0, V_pd_pi_WS2=-0.942d0
	!W-W
	real(KIND=8), parameter :: V_dd_sigma_WS2=-1.216d0, V_dd_pi_WS2=0.177d0, V_dd_delta_WS2=0.243d0
	!S-S
	real(KIND=8), parameter :: V_pp_sigma_WS2=0.749d0, V_pp_pi_WS2=0.236d0
	!SOC
	real(KIND=8), parameter :: SOC_dd_WS2=0.215d0, SOC_pp_WS2=0.057d0
	!S-S Interlayer
	real(KIND=8), parameter :: V_pp_sigma_interlayer_WS2=-0.55d0, V_pp_pi_interlayer_WS2=-0.6d0
	
	!real(kind=8), Save :: V_pd_sigma,V_pd_pi,V_dd_sigma,V_dd_pi,V_dd_delta,V_pp_sigma,V_pp_pi,V_pp_sigma_interlayer,V_pp_pi_interlayer,SOC_dd,SOC_pp,V_pp_sigma_temp,V_pp_pi_temp
	
!*****************************************************************************************************	
	






	End module moduletwo
