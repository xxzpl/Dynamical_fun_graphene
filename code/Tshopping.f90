!!!!!!!!-----------------------------Subrountine: Initialize the Exchanging Interaction-------one layer---------------------------------------------------------------------------

subroutine Hoppingvalues(Tswitch, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ)
use moduleone
use moduletwo
	
Implicit real (kind=8) (a-h,o-z) 

integer(KIND=4)::Tswitch

complex(kind=8):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
real(KIND=8)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
real(KIND=8)::V_pd_sigma, V_pd_pi, V_dd_sigma, V_dd_pi, V_dd_delta, V_pp_sigma, V_pp_pi, V_pp_sigma_interlayer
real(KIND=8)::V_pp_pi_interlayer, SOC_dd,SOC_pp


	
!	MoS2
	if (Tswitch==1) then
	Delta0=Delta0_MoS2
	Delta1=Delta1_MoS2
	Delta2=Delta2_MoS2
	DeltaP=DeltaP_MoS2
	DeltaZ=DeltaZ_MoS2
	V_pd_sigma=V_pd_sigma_MoS2
	V_pd_pi=V_pd_pi_MoS2
	V_dd_sigma=V_dd_sigma_MoS2
	V_dd_pi=V_dd_pi_MoS2
	V_dd_delta=V_dd_delta_MoS2
	V_pp_sigma=V_pp_sigma_MoS2
	V_pp_pi=V_pp_pi_MoS2
	V_pp_sigma_interlayer=V_pp_sigma_interlayer_MoS2
	V_pp_pi_interlayer=V_pp_pi_interlayer_MoS2
	SOC_dd=SOC_dd_MoS2
	SOC_pp=SOC_pp_MoS2	
!	WS2
	else if (Tswitch==2) then
	Delta0=Delta0_WS2
	Delta1=Delta1_WS2
	Delta2=Delta2_WS2
	DeltaP=DeltaP_WS2
	DeltaZ=DeltaZ_WS2
	V_pd_sigma=V_pd_sigma_WS2
	V_pd_pi=V_pd_pi_WS2
	V_dd_sigma=V_dd_sigma_WS2
	V_dd_pi=V_dd_pi_WS2
	V_dd_delta=V_dd_delta_WS2
	V_pp_sigma=V_pp_sigma_WS2
	V_pp_pi=V_pp_pi_WS2
	V_pp_sigma_interlayer=V_pp_sigma_interlayer_WS2
	V_pp_pi_interlayer=V_pp_pi_interlayer_WS2
	SOC_dd=SOC_dd_WS2
	SOC_pp=SOC_pp_WS2		
	end if
	
!-------------------------------------------------------------------------------------------------------------------------------
	
!	Next Nearest Hopping
!	Mo-Mo
!	Z2-1, X2-2, Y2-3
!		6(alpha)	5(-gemma)
!	1(-beta)			4(beta)
!		2(gemma)	3(-alpha)
!	alpha
	tdd(1,1,6)=1./16*(9*V_dd_sigma+4*V_dd_pi+3*V_dd_delta)
	tdd(1,2,6)=0.
	tdd(1,3,6)=0.
	tdd(1,4,6)=sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(1,5,6)=3./8*(V_dd_sigma-V_dd_delta)
	
	tdd(2,1,6)=0.
	tdd(2,2,6)=1./4*(3*V_dd_pi+V_dd_delta)
	tdd(2,3,6)=-sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(2,4,6)=0.
	tdd(2,5,6)=0.
	
	tdd(3,1,6)=0.
	tdd(3,2,6)=-sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(3,3,6)=1./4*(V_dd_pi+ 3*V_dd_delta)
	tdd(3,4,6)=0.
	tdd(3,5,6)=0.	

	tdd(4,1,6)=sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(4,2,6)=0.
	tdd(4,3,6)=0.
	tdd(4,4,6)=1./16*(3*V_dd_sigma+12*V_dd_pi+V_dd_delta)
	tdd(4,5,6)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	
	tdd(5,1,6)=3./8*(V_dd_sigma-V_dd_delta)
	tdd(5,2,6)=0.
	tdd(5,3,6)=0.
	tdd(5,4,6)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	tdd(5,5,6)=1./4*(V_dd_sigma+3*V_dd_delta)

	
!	beta
	tdd(1,1,4)=V_dd_pi
	tdd(1,2,4)=0.
	tdd(1,3,4)=0.
	tdd(1,4,4)=0.
	tdd(1,5,4)=0.
	
	tdd(2,1,4)=0.
	tdd(2,2,4)=V_dd_delta
	tdd(2,3,4)=0.
	tdd(2,4,4)=0.
	tdd(2,5,4)=0.
	
	tdd(3,1,4)=0.
	tdd(3,2,4)=0.
	tdd(3,3,4)=V_dd_pi
	tdd(3,4,4)=0.
	tdd(3,5,4)=0.	

	tdd(4,1,4)=0.
	tdd(4,2,4)=0.
	tdd(4,3,4)=0.
	tdd(4,4,4)=1./4*(3*V_dd_sigma+V_dd_delta)
	tdd(4,5,4)=-sqrt(3.)/4*(V_dd_sigma-V_dd_delta)
	
	tdd(5,1,4)=0.
	tdd(5,2,4)=0.
	tdd(5,3,4)=0.
	tdd(5,4,4)=-sqrt(3.)/4*(V_dd_sigma-V_dd_delta)
	tdd(5,5,4)=1./4*(V_dd_sigma+3*V_dd_delta)
	
	
!	gamma
	tdd(1,1,2)=1./16*(9*V_dd_sigma+4*V_dd_pi+3*V_dd_delta)
	tdd(1,2,2)=0.
	tdd(1,3,2)=0.
	tdd(1,4,2)=-sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(1,5,2)=-3./8*(V_dd_sigma-V_dd_delta)
	
	tdd(2,1,2)=0.
	tdd(2,2,2)=1./4*(3*V_dd_pi+V_dd_delta)
	tdd(2,3,2)=sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(2,4,2)=0.
	tdd(2,5,2)=0.
	
	tdd(3,1,2)=0.
	tdd(3,2,2)=sqrt(3.)/4*(V_dd_pi-V_dd_delta)
	tdd(3,3,2)=1./4*(V_dd_pi+ 3*V_dd_delta)
	tdd(3,4,2)=0.
	tdd(3,5,2)=0.	

	tdd(4,1,2)=-sqrt(3.)/16*(3*V_dd_sigma-4*V_dd_pi+V_dd_delta)
	tdd(4,2,2)=0.
	tdd(4,3,2)=0.
	tdd(4,4,2)=1./16*(3*V_dd_sigma+12*V_dd_pi+V_dd_delta)
	tdd(4,5,2)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	
	tdd(5,1,2)=-3./8*(V_dd_sigma-V_dd_delta)
	tdd(5,2,2)=0.
	tdd(5,3,2)=0.
	tdd(5,4,2)=sqrt(3.)/8*(V_dd_sigma-V_dd_delta)
	tdd(5,5,2)=1./4*(V_dd_sigma+3*V_dd_delta)

	tdd(:,:,1)=tdd(:,:,4)
	tdd(:,:,3)=tdd(:,:,6)
	tdd(:,:,5)=tdd(:,:,2)
	
	!~ tdd(:,:,1)=dconjg(tdd(:,:,4))
	!~ tdd(:,:,3)=dconjg(tdd(:,:,6))
	!~ tdd(:,:,5)=dconjg(tdd(:,:,2))
	
!	S-S
!	x-1, y-2, z-3
!		6(alpha)	5(-gemma)
!	1(-beta)			4(beta)
!		2(gemma)	3(-alpha)
!	alpha
	tpp(1,1,6)=1./4*(3*V_pp_pi+V_pp_sigma)
	tpp(1,2,6)=sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(1,3,6)=0.
	tpp(2,1,6)=sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(2,2,6)=1./4*(V_pp_pi+3*V_pp_sigma)
	tpp(2,3,6)=0.
	tpp(3,1,6)=0.
	tpp(3,2,6)=0.
	tpp(3,3,6)=V_pp_pi

!	beta
	tpp(1,1,4)=V_pp_sigma
	tpp(1,2,4)=0.
	tpp(1,3,4)=0.
	tpp(2,1,4)=0.
	tpp(2,2,4)=V_pp_pi
	tpp(2,3,4)=0.
	tpp(3,1,4)=0.
	tpp(3,2,4)=0.
	tpp(3,3,4)=V_pp_pi
	
!	gemma
	tpp(1,1,2)=1./4*(3*V_pp_pi+V_pp_sigma)
	tpp(1,2,2)=-sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(1,3,2)=0.
	tpp(2,1,2)=-sqrt(3.)/4*(V_pp_pi-V_pp_sigma)
	tpp(2,2,2)=1./4*(V_pp_pi+3*V_pp_sigma)
	tpp(2,3,2)=0.
	tpp(3,1,2)=0.
	tpp(3,2,2)=0.
	tpp(3,3,2)=V_pp_pi 

	tpp(:,:,1)=tpp(:,:,4)
	tpp(:,:,3)=tpp(:,:,6)
	tpp(:,:,5)=tpp(:,:,2)
	
	
!	MoS2
	if (Tswitch==1) then	
!--------------------------------------------------------------------------------------------------------------------

!	Nearest Hoping
!	Mo-S	Bottom layer
!    2 (gemma)		3 (alpha)
!		Mo
!		1 (beta)
!	(Mo,S,Direction)  
	a=sqrt(1./7)/7
!	beta	
	tdp(1,1,1)=14*a*V_pd_pi
	tdp(1,2,1)=0.
	tdp(1,3,1)=7*a*sqrt(3.)*V_pd_pi
	tdp(1,4,1)=0.
	tdp(1,5,1)=0.
	
	tdp(2,1,1)=0.
	tdp(2,2,1)=sqrt(3.)*a*(4*sqrt(3.)*V_pd_sigma-V_pd_pi )
	tdp(2,3,1)=0.
	tdp(2,4,1)=-1*a*(4*sqrt(3.)*V_pd_sigma+6*V_pd_pi )
	tdp(2,5,1)=2*a*(V_pd_sigma- 3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,1)=0.
	tdp(3,2,1)=2*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,1)=0.
	tdp(3,4,1)=-2*a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,1)=sqrt(3.)*a*(V_pd_sigma+ 4*sqrt(3.)*V_pd_pi )
	
!	gemma	
	tdp(1,1,2)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,2)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,2)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,2)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,2)=sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,2)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,2)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,2)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,2)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,2)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,2)=-3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,2)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,2)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,2)=a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,2)=sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	
!	alpha
	tdp(1,1,3)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,3)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,3)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,3)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,3)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,3)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,3)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,3)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,3)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,3)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,3)=3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,3)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,3)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,3)=a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,3)=sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	

!	Nearest Hoping
!	Mo-S Top Layer
!    5 (5gemma)		6 (6alpha)
!		Mo
!		4 (4beta)
!	(Mo,S,Direction)  


	a=sqrt(1./7)/7
!	4beta	
	tdp(1,1,4)=14*a*V_pd_pi
	tdp(1,2,4)=0.
	tdp(1,3,4)=-7*a*sqrt(3.)*V_pd_pi
	tdp(1,4,4)=0.
	tdp(1,5,4)=0.
	
	tdp(2,1,4)=0.
	tdp(2,2,4)=-sqrt(3.)*a*(4*sqrt(3.)*V_pd_sigma-V_pd_pi )
	tdp(2,3,4)=0.
	tdp(2,4,4)=-1*a*(4*sqrt(3.)*V_pd_sigma+6*V_pd_pi )
	tdp(2,5,4)=2*a*(V_pd_sigma- 3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,4)=0.
	tdp(3,2,4)=2*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,4)=0.
	tdp(3,4,4)=2*a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,4)=-sqrt(3.)*a*(V_pd_sigma+ 4*sqrt(3.)*V_pd_pi )
	
!	gemma	
	tdp(1,1,5)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,5)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,5)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,5)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,5)=sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,5)=sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,5)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,5)=3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,5)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,5)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,5)=3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,5)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,5)=sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,5)=-a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,5)=-sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	
!	alpha
	tdp(1,1,6)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,2,6)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(1,3,6)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(1,4,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(1,5,6)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(2,1,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,2,6)=-sqrt(3.)*a*(sqrt(3.)*V_pd_sigma+5*V_pd_pi )
	tdp(2,3,6)=-3*a*(sqrt(3.)*V_pd_sigma -2*V_pd_pi )
	tdp(2,4,6)=-sqrt(3.)*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	tdp(2,5,6)=-1*a*(V_pd_sigma-3*sqrt(3.)*V_pd_pi )
	
	tdp(3,1,6)=-3*a*(sqrt(3.)*V_pd_sigma- 2*V_pd_pi )
	tdp(3,2,6)=-1*a*(3*sqrt(3.)*V_pd_sigma+V_pd_pi )
	tdp(3,3,6)=-sqrt(3.)*a*(3*sqrt(3.)*V_pd_sigma +V_pd_pi )
	tdp(3,4,6)=-a*(3*V_pd_sigma-2*sqrt(3.)*V_pd_pi )
	tdp(3,5,6)=-sqrt(3.)*a*(V_pd_sigma+4*sqrt(3.)*V_pd_pi )
	
!--------------------------------------------- interlayer hopping  S-S ------------------------------- 	
!	Interlayer Hoping
!	S-S	
!    8 (u)		9 (v)
!		  S
!		7 (lambda)
!	(Mo,S,Direction)  
	
!	lambda
	tpp(1,1,7)=V_pp_pi
	tpp(1,2,7)=0.
	tpp(1,3,7)=0.
	tpp(2,1,7)=0.
	tpp(2,2,7)=0.214*V_pp_sigma+0.786*V_pp_pi
	tpp(2,3,7)=0.410*V_pp_sigma-0.410*V_pp_pi
	tpp(3,1,7)=0.
	tpp(3,2,7)=0.410*V_pp_sigma-0.410*V_pp_pi
	tpp(3,3,7)=0.785*V_pp_sigma+0.215*V_pp_pi
	
!	u
	tpp(1,1,8)=0.161*V_pp_sigma+0.839*V_pp_pi
	tpp(1,2,8)=-0.093*V_pp_sigma+ 0.093*V_pp_pi
	tpp(1,3,8)=0.355*V_pp_sigma- 0.355*V_pp_pi
	tpp(2,1,8)=-0.093*V_pp_sigma+ 0.093*V_pp_pi
	tpp(2,2,8)=0.054*V_pp_sigma+0.946*V_pp_pi
	tpp(2,3,8)=-0.206*V_pp_sigma+ 0.206*V_pp_pi
	tpp(3,1,8)=0.355*V_pp_sigma- 0.355*V_pp_pi
	tpp(3,2,8)=-0.206*V_pp_sigma+ 0.206*V_pp_pi
	tpp(3,3,8)=0.785*V_pp_sigma+0.215*V_pp_pi
	
!	v
	tpp(1,1,9)=0.161*V_pp_sigma+0.839*V_pp_pi
	tpp(1,2,9)=0.093*V_pp_sigma-0.093*V_pp_pi
	tpp(1,3,9)=-0.355*V_pp_sigma+0.355*V_pp_pi
	tpp(2,1,9)=0.093*V_pp_sigma-0.093*V_pp_pi
	tpp(2,2,9)=0.054*V_pp_sigma+0.946*V_pp_pi
	tpp(2,3,9)=-0.206*V_pp_sigma+ 0.206*V_pp_pi
	tpp(3,1,9)=-0.355*V_pp_sigma+0.355*V_pp_pi
	tpp(3,2,9)=-0.206*V_pp_sigma+ 0.206*V_pp_pi
	tpp(3,3,9)=0.785*V_pp_sigma+0.215*V_pp_pi	

!	WS2
	else if (Tswitch==2) then

!--------------------------------------------------------------------------------------------------------------------

!	Nearest Hoping
!	Mo-S	
!    2 (gemma)		3 (alpha)
!		Mo
!		1 (beta)
!	(Mo,S,Direction)  
	
!	beta	
	tdp(1,1,1)=0.75706*V_pd_pi
	tdp(1,2,1)=0.
	tdp(1,3,1)=0.653346*V_pd_pi
	tdp(1,4,1)=0.
	tdp(1,5,1)=0.
	
	tdp(2,1,1)=0.
	tdp(2,2,1)=0.648581*V_pd_sigma-0.0955712*V_pd_pi 
	tdp(2,3,1)=0.
	tdp(2,4,1)=-0.37577*V_pd_sigma-0.323159*V_pd_pi 
	tdp(2,5,1)=0.106209*V_pd_sigma-0.559728*V_pd_pi 
	
	tdp(3,1,1)=0.
	tdp(3,2,1)=0.559728*V_pd_sigma+0.110741*V_pd_pi 
	tdp(3,3,1)=0.
	tdp(3,4,1)=-0.324291*V_pd_sigma+0.374459*V_pd_pi 
	tdp(3,5,1)=0.0916587*V_pd_sigma+0.648581*V_pd_pi 
	
!	gemma	
	tdp(1,1,2)=-0.281827*V_pd_sigma-0.0531043*V_pd_pi 
	tdp(1,2,2)=-0.280844*V_pd_sigma+0.324291*V_pd_pi 
	tdp(1,3,2)= 0.486436*V_pd_sigma +0.0916584*V_pd_pi 
	tdp(1,4,2)= 0.162713*V_pd_sigma +0.467748*V_pd_pi 
	tdp(1,5,2)= 0.0919795*V_pd_sigma -0.484739*V_pd_pi 
	
	tdp(2,1,2)=0.162713*V_pd_sigma +0.467748*V_pd_pi 
	tdp(2,2,2)=0.162145*V_pd_sigma +0.466117*V_pd_pi 
	tdp(2,3,2)=-0.280844*V_pd_sigma +0.324291*V_pd_pi 
	tdp(2,4,2)=-0.0939423*V_pd_sigma+0.487005*V_pd_pi 
	tdp(2,5,2)=-0.0531044*V_pd_sigma+0.279864*V_pd_pi 
	
	tdp(3,1,2)=-0.280844*V_pd_sigma +0.324291*V_pd_pi 
	tdp(3,2,2)=-0.279864*V_pd_sigma -0.0553706*V_pd_pi 
	tdp(3,3,2)=0.484739*V_pd_sigma +0.0959047*V_pd_pi 
	tdp(3,4,2)=0.162145*V_pd_sigma -0.187229*V_pd_pi 
	tdp(3,5,2)=0.0916587*V_pd_sigma +0.648581*V_pd_pi 
	
!	alpha
	tdp(1,1,3)=-0.281827*V_pd_sigma -0.0531043*V_pd_pi 
	tdp(1,2,3)=0.280844*V_pd_sigma -0.324291*V_pd_pi 
	tdp(1,3,3)=0.486436*V_pd_sigma +0.0916584*V_pd_pi 
	tdp(1,4,3)=-0.162713*V_pd_sigma -0.467748*V_pd_pi 
	tdp(1,5,3)=-0.0919795*V_pd_sigma +0.484739*V_pd_pi 
	
	tdp(2,1,3)=-0.162713*V_pd_sigma -0.467748*V_pd_pi 
	tdp(2,2,3)=0.162145*V_pd_sigma +0.466117*V_pd_pi 
	tdp(2,3,3)=0.280844*V_pd_sigma -0.324291*V_pd_pi 
	tdp(2,4,3)=-0.0939423*V_pd_sigma +0.487005*V_pd_pi 
	tdp(2,5,3)=-0.0531044*V_pd_sigma +0.279864*V_pd_pi 
	
	tdp(3,1,3)=0.280844*V_pd_sigma -0.324291*V_pd_pi 
	tdp(3,2,3)=-0.279864*V_pd_sigma -0.0553706*V_pd_pi 
	tdp(3,3,3)=-0.484739*V_pd_sigma -0.0959047*V_pd_pi 
	tdp(3,4,3)=0.162145*V_pd_sigma -0.187229*V_pd_pi 
	tdp(3,5,3)=0.0916587*V_pd_sigma +0.648581*V_pd_pi 
	

!	Nearest Hoping
!	Mo-S	
!    5 (5gemma)		6 (6alpha)
!		Mo
!		4 (4beta)
!	(Mo,S,Direction)  


	
!	4beta	
	tdp(1,1,4)=0.75706*V_pd_pi
	tdp(1,2,4)=0.
	tdp(1,3,4)=-0.653346*V_pd_pi
	tdp(1,4,4)=0.
	tdp(1,5,4)=0.
	
	tdp(2,1,4)=0.
	tdp(2,2,4)=-0.648581*V_pd_sigma +0.0955712*V_pd_pi 
	tdp(2,3,4)=0.
	tdp(2,4,4)=-0.37577*V_pd_sigma -0.323159*V_pd_pi 
	tdp(2,5,4)=0.106209*V_pd_sigma -0.559728*V_pd_pi 
	
	tdp(3,1,4)=0.
	tdp(3,2,4)=0.559728*V_pd_sigma +0.110741*V_pd_pi 
	tdp(3,3,4)=0.
	tdp(3,4,4)=0.324291*V_pd_sigma -0.374459*V_pd_pi 
	tdp(3,5,4)=-0.0916587*V_pd_sigma -0.648581*V_pd_pi 
	
!	gemma	
	tdp(1,1,5)=-0.281827*V_pd_sigma -0.0531043*V_pd_pi 
	tdp(1,2,5)=0.280844*V_pd_sigma -0.324291*V_pd_pi 
	tdp(1,3,5)=-0.486436*V_pd_sigma -0.0916584*V_pd_pi 
	tdp(1,4,5)=0.162713*V_pd_sigma +0.467748*V_pd_pi 
	tdp(1,5,5)=0.0919795*V_pd_sigma -0.484739*V_pd_pi 
	
	tdp(2,1,5)=0.162713*V_pd_sigma +0.467748*V_pd_pi 
	tdp(2,2,5)=-0.162145*V_pd_sigma -0.466117*V_pd_pi 
	tdp(2,3,5)=0.280844*V_pd_sigma -0.324291*V_pd_pi 
	tdp(2,4,5)=-0.0939423*V_pd_sigma +0.487005*V_pd_pi 
	tdp(2,5,5)=-0.0531044*V_pd_sigma +0.279864*V_pd_pi 
	
	tdp(3,1,5)=0.280844*V_pd_sigma -0.324291*V_pd_pi 
	tdp(3,2,5)=-0.279864*V_pd_sigma -0.0553706*V_pd_pi 
	tdp(3,3,5)=0.484739*V_pd_sigma +0.0959047*V_pd_pi 
	tdp(3,4,5)=-0.162145*V_pd_sigma +0.187229*V_pd_pi 
	tdp(3,5,5)=-0.0916587*V_pd_sigma -0.648581*V_pd_pi 
	
!	alpha
	tdp(1,1,6)=-0.281827*V_pd_sigma -0.0531043*V_pd_pi 
	tdp(1,2,6)=-0.280844*V_pd_sigma +0.324291*V_pd_pi 
	tdp(1,3,6)=-0.486436*V_pd_sigma -0.0916584*V_pd_pi 
	tdp(1,4,6)=-0.162713*V_pd_sigma -0.467748*V_pd_pi 
	tdp(1,5,6)=-0.0919795*V_pd_sigma +0.484739*V_pd_pi 
	
	tdp(2,1,6)=-0.162713*V_pd_sigma -0.467748*V_pd_pi 
	tdp(2,2,6)=-0.162145*V_pd_sigma -0.466117*V_pd_pi 
	tdp(2,3,6)=-0.280844*V_pd_sigma +0.324291*V_pd_pi 
	tdp(2,4,6)=-0.0939423*V_pd_sigma +0.487005*V_pd_pi 
	tdp(2,5,6)=-0.0531044*V_pd_sigma +0.279864*V_pd_pi 
	
	tdp(3,1,6)=-0.280844*V_pd_sigma +0.324291*V_pd_pi 
	tdp(3,2,6)=-0.279864*V_pd_sigma -0.0553706*V_pd_pi 
	tdp(3,3,6)=-0.484739*V_pd_sigma -0.0959047*V_pd_pi 
	tdp(3,4,6)=-0.162145*V_pd_sigma +0.187229*V_pd_pi 
	tdp(3,5,6)=-0.0916587*V_pd_sigma -0.648581*V_pd_pi 	
	

!--------------------------------------------- interlayer hopping  S-S ------------------------------- 	
!	Interlayer Hoping
!	S-S	
!    8 (u)		9 (v)
!		S
!		7 (lambda)
!	(Mo,S,Direction)  
	
!	lambda
	tpp(1,1,7)=V_pp_pi
	tpp(1,2,7)=0.
	tpp(1,3,7)=0.
	tpp(2,1,7)=0.
	tpp(2,2,7)=0.210467*V_pp_sigma+0.789533*V_pp_pi
	tpp(2,3,7)=0.407641*V_pp_sigma-0.407641*V_pp_pi
	tpp(3,1,7)=0.
	tpp(3,2,7)=0.407641*V_pp_sigma-0.407641*V_pp_pi
	tpp(3,3,7)=0.789534*V_pp_sigma+0.210466*V_pp_pi
	
!	u
	tpp(1,1,8)=0.15785*V_pp_sigma+0.84215*V_pp_pi
	tpp(1,2,8)=-0.0911352*V_pp_sigma+ 0.0911352*V_pp_pi
	tpp(1,3,8)=0.353027*V_pp_sigma- 0.353027*V_pp_pi
	tpp(2,1,8)=-0.0911352*V_pp_sigma+0.0911352*V_pp_pi
	tpp(2,2,8)=0.052617*V_pp_sigma+0.947383*V_pp_pi
	tpp(2,3,8)=-0.203821*V_pp_sigma+ 0.203821*V_pp_pi
	tpp(3,1,8)=0.353027*V_pp_sigma- 0.353027*V_pp_pi
	tpp(3,2,8)=-0.203821*V_pp_sigma+ 0.203821*V_pp_pi
	tpp(3,3,8)=0.789534*V_pp_sigma+0.210466*V_pp_pi
	
!	v
	tpp(1,1,9)=0.15785*V_pp_sigma+0.84215*V_pp_pi
	tpp(1,2,9)=0.0911352*V_pp_sigma-0.0911352*V_pp_pi
	tpp(1,3,9)=-0.353027*V_pp_sigma+0.353027*V_pp_pi
	tpp(2,1,9)=0.0911352*V_pp_sigma-0.0911352*V_pp_pi
	tpp(2,2,9)=0.052617*V_pp_sigma+0.947383*V_pp_pi
	tpp(2,3,9)=-0.203821*V_pp_sigma+ 0.203821*V_pp_pi
	tpp(3,1,9)=-0.353027*V_pp_sigma+0.353027*V_pp_pi
	tpp(3,2,9)=-0.203821*V_pp_sigma+ 0.203821*V_pp_pi
	tpp(3,3,9)=0.789534*V_pp_sigma+0.210466*V_pp_pi	
end if

!--------------------------------------------- hopping  S(PT)-S(PB) ------------------------------- 	
!	
!	S-S	
!                    S
!     10
!		S
!	
	tpp(1,1,10)=V_pp_pi
	tpp(1,2,10)=0.
	tpp(1,3,10)=0.
	tpp(2,1,10)=0.
	tpp(2,2,10)=V_pp_pi
	tpp(2,3,10)=0.
	tpp(3,1,10)=0.
	tpp(3,2,10)=0.
	tpp(3,3,10)=V_pp_sigma	
	
!===============================================================================

!~ open (unit =UNIT3, file = FILENAME3,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)
!~ write(UNIT3, '( 60("-") )')
!~ write (UNIT3, '("#TDP(1:3,1:5,1:6)", 3(" "), "tdp(1,k,i)", 3(" "), "tdp(2,k,i)", 3(" "), "tdp(3,k,i)" )')
!~ write(UNIT3, '( 60("-") )')
!~ write(UNIT3, '(a)') 
 
   !~ Do i=1,6
          !~ Do k=1,5
           !~ write(UNIT3, '(1x,3(3x2E18.11) )') tdp(1,k,i), tdp(2,k,i), tdp(3,k,i)
	  !~ end do
  !~ end do

!~ write(UNIT3, '(a)')   
!~ write(UNIT3, '( 60("-") )')
!~ write (UNIT3, '("#TPP(1:3,1:3,6|4|2)", 3(" "), "tpp(1,k,i)", 3(" "), "tpp(2,k,i)", 3(" "), "tpp(3,k,i)" )')
!~ write(UNIT3, '( 60("-") )')
!~ write(UNIT3, '(a)') 
 
   !~ Do i=6, 2, -2
          !~ Do k=1,3
           !~ write(UNIT3, '(1x,3(3x2E18.11) )') tpp(1,k,i), tpp(2,k,i), tpp(3,k,i)
	  !~ end do
  !~ end do 
  
!~ write(UNIT3, '(a)')   
!~ write(UNIT3, '( 60("-") )')
!~ write (UNIT3, '("#TDD(1:5,1:5,6|4|2)", 3(" "), "tdd(1,k,i)", 3(" "), "tdd(2,k,i)", 3(" "), "tdd(3,k,i)" &
                         !~ ,3(" "), "tdd(4,k,i)" , 3(" "), "tdd(5,k,i)" )')
!~ write(UNIT3, '( 60("-") )')
!~ write(UNIT3, '(a)') 
!~ ! 
   !~ Do i=6, 2, -2
          !~ Do k=1,5
           !~ write(UNIT3, '(1x,5(3x2E18.11) )') tdd(1,k,i), tdd(2,k,i), tdd(3,k,i), tdd(4,k,i), tdd(5,k,i)
	  !~ end do
  !~ end do 

!~ close(UNIT3)
end 




	
	
	