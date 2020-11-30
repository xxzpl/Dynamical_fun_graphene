program DynamicalpolarizationPI
use moduleone
use moduletwo
Implicit real (kind=8) (a-h,o-z)
real(kind=8)::  Qx, Qy, Omega, Deltasmall,   kx, ky,   Vdeta12, Vdeta23
complex(kind=8):: DynamicalPoPI
real(KIND=8), ALLOCATABLE, dimension(:) :: rankx, ranky, RCKS

real(kind=8),dimension(Hdimension):: eig
complex(kind=8),dimension(Hdimension,Hdimension):: vec

complex(kind=8), Allocatable, dimension(:):: Cks
complex(kind=8), Allocatable, dimension(:, :, :):: ckvec
real(kind=8), ALLOCATABLE, dimension( : , :  )::  ckeig
complex(kind=8)::  CKQ
integer:: Local(1)

Omega=0.d0
Deltasmall=0.5d0*10.0d0**(-3)
alpha=2.0/(4*pi**2)


!~ ===================================================
    open (unit =UNIT11, file = FILENAME11,  status="UNKNOWN" ,  iostat=ierror)   !! wavenumber.dat
    do i=1,2                 !//Skip the first 2 lines
        read(UNIT11,*) line
    end do
    Lens=0
	Do
	    read(UNIT11, *, IOSTAT=ierr) line
		if(ierr/=0) exit
		Lens=Lens+1
	end do
    rewind(UNIT11)
	
    !~ write(*,  '(  1xI8 )' ) Lens
    allocate( rankx(1:Lens) )
    allocate( ranky(1:Lens) )
	allocate(  RCks(1:Lens) )
	allocate(    Cks(1:Lens) )
	
	allocate(ckvec(1:Lens, Hdimension,Hdimension) )
	allocate(ckeig(1:Lens, Hdimension) )
	do i=1,2                 !//Skip the first 2 lines
        read(UNIT11,*) line
    end do
    do j=1, Lens           !//Read the data from wavenumbers.dat
        read(UNIT11, '(13x, 2(1xES23.15E3) )' ) rankx(j), ranky(j)
    end do
	close(UNIT11)
	!~ |deltak|
   !~ write(*, '(13x, 2(1xES23.15E3) )' ) rankx(12)-rankx(11), ranky(13)-ranky(12)
 


    write(*, '( a )' )
	write(*, '(13x, 2(1xES23.15E3) )' ) rankx(12), ranky(12)
    Cks=CMPLX( rankx, ranky)
    write(*, '( a )' )
	write(*, * ) Cks(12)

Qx= 4.54902617;      Qy=0.00000000
!~ Qx= 4.52389342;      Qy=0.41887902

   write(*, '( a )' )
   write(*, '( a, I6 )'  ) "Lens is::", lens
   
kx=rankx(12)+Qx
ky=ranky(12)+Qy
write(*,*)
write(*,*) kx, ky

call wavenumberlocal( Lens, rankx(1:Lens), ranky(1:Lens),  rankx(12),ranky(12), Qx, Qy,  M ) 
   write(*, '( a )' )
   write(*, '( a, I6 )'  ) "The value of M is::", M


Bvector1x=(2.0_DP*Pi)/(sq3*carbondistance)
Bvector1y=(2.0_DP*Pi)/(3.0_DP*carbondistance)
Xbound=(4.0_DP*Pi)/(3.0_DP*sq3*carbondistance)
Ybound=(2.0_DP*Pi)/(3.0_DP*carbondistance)

	          if(    ky>= sq3*kx          .and. &
					   ky> Ybound              .and. &
					   ky>= -sq3*kx   )   then
		    ky=ky-2.0_DP*Bvector1y			
		else if (  ky>=0.0d0                 .and. &
					   ky<= sq3*kx         .and. &
					   ky> 2.0_DP*Ybound-sq3*kx   ) then		
             kx=kx-Bvector1x
		     ky=ky-Bvector1y 
		else if (  ky<=0.0d0                 .and. &
					   ky>=-sq3*kx         .and. &
					   ky< -2.0_DP*Ybound+sq3*kx   ) then
             kx=kx-Bvector1x
		     ky=ky+Bvector1y
		else if (  ky<= sq3*kx          .and. &
					   ky< -Ybound              .and. &
					   ky<= -sq3*kx   )   then					   
		     ky=ky+2.0_DP*Bvector1y			 
		else if (  ky<=0.0d0                 .and. &
					   ky>= sq3*kx         .and. &
					   ky< -2.0_DP*Ybound-sq3*kx   ) then		
             kx=kx+Bvector1x
		     ky=ky+Bvector1y	
		else if (  ky>=0.0d0                 .and. &
					   ky<=-sq3*kx         .and. &
					   ky> 2.0_DP*Ybound+sq3*kx   ) then
             kx=kx+Bvector1x
		     ky=ky-Bvector1y			 
		end if	 


CKQ=CMPLX(kx, ky)

   write(*, '( a )' )
   write(*, *  ) CKQ

Cks=Cks-CKQ



 
do i=1, Lens
        RCks(i)=abs(Cks(i) )
end do

Local=MINLOC( RCks )
     
   write(*, '( a )' )
   write(*, '( I6 )'  ) Local(1)   
    write(*, '( a )' )
   write(*, '( a )' )	

   
Vdeta12=0.0d0
Vdeta23=0.0d0
   
call Hoppingvalues(SetTransitionMetalDichalcogenides, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ)   !! get the hopping values ts

 	   call Hamiltonianofxs2( kx, ky, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)
               
write(*, *  ) eig
write(*,'(a)' )
write(*,'(a)' )


i=Local(1)
 	   call Hamiltonianofxs2(rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)
 
write(*,'(a)' ) "Polarization eig" 
write(*, *  ) eig
write(*,'(a)' )
write(*,'(a)' )

!~ i=M
!~ write(*,*) M  ! 214796
 	   !~ call Hamiltonianofxs2(rankx(214796), ranky(214796), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)
 
!~ write(*,'(a)' ) "Subrouting  eig(M)" 
!~ write(*, *  ) eig


!~ 546   245  2.155759751314448E+000 -4.188790204786397E-002   215088

!~ CQ=CMPLX(Qx, Qy)

!~ complex(kind=8), Allocatable, dimension(:, :, :):: ckvec
!~ real(kind=8), ALLOCATABLE, dimension( : , :  )::  ckeig


!~ ===================================================


!~ end function
end




