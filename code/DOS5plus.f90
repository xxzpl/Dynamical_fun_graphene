!===========================================================================
subroutine Counter(A, probability, Densityofstate)
use moduleone
Implicit real (kind=8) (a-h,o-z) 

real(KIND=8) :: A,  probability
real(KIND=8), dimension(0:Ndeltax) :: Densityofstate

Densityofstate=0.d0
kc=Nint((A- Al)/deltaX) 

Densityofstate(kc)=Densityofstate(kc)+probability

end	subroutine Counter

!===========================================================================

  subroutine Rnorm(nstates,psi_R,r0)
        implicit real(kind=8) (a-h,o-z)
        real(kind=8) psi_R(0:nstates-1)

        r0=0.0
        do j=0,nstates-1
        r0=r0+psi_R(j)
        enddo
end  subroutine Rnorm


!===========================================================================

subroutine Vfun(nprime, carriernumber, distanceofinterlayer, Sunitcell,  Vdeta)
 implicit none
 
 integer(KIND=4), intent(in):: carriernumber
 
 real(KIND=8), intent(in):: nprime, Sunitcell
 
 real(KIND=8):: Vdeta, electron, epsilong, kappa, distanceofinterlayer
 
 
  electron=1.60217653*10.0**(-19) !!   unit ( C )

  kappa = (3.90d0+1.0d0)/2.0d0 ;     

  epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V
    
Vdeta=electron*distanceofinterlayer*carriernumber*nprime/(epsilong*kappa*Sunitcell)    ! UNIT ~ eV

!~ AA    Vdeta=electron*distanceofinterlayer*nprime*carriernumber*10.d0**(-17)/(epsilong*kappa*Sunitcell) ! UNIT ~ eV  


end subroutine Vfun



subroutine Efun(nprime, Efiled)    !! nprime  UNIT   (     10^{13} cm^{-2}   )
 implicit none
 real(KIND=8), intent(in):: nprime     !! UNIT   (     10^{13} cm^{-2}   )
 
 real(KIND=8), intent(out):: Efiled
 
 real(KIND=8) :: electron, epsilong, kappa
 
  electron=1.60217653*10.0**(-19) !!   unit ( C )

  kappa = (3.90d0+1.0d0)/2.0d0 ;     

  epsilong = 8.854187817*10.0**(-12)   !! unit (  F/m  )    F=C/V
  
  
Efiled=electron*nprime*10.d0**(8)/(2.0d0*epsilong*kappa)    !! UNIT  (    V/nm   )


end subroutine Efun







subroutine   area( Densityofstate, N0, Nend, S )
use moduleone
Implicit real (kind=8) (a-h,o-z) 

real (kind=8) , dimension(0:Ndeltax) :: Densityofstate

integer(kind=4):: N0, Nend

real(KIND=8) :: S

   
         S=0.0_DP 
              do i=N0+1, Nend-1
                  S=S+Densityofstate(i)
              end do
	 S=S+ (Densityofstate(N0)+Densityofstate(Nend))/2.0_DP

end  subroutine area


!~ ========================================
                       !~ wavenumberlocal
!~ ========================================


subroutine wavenumberlocal( Lens, vectorx, vectory,  wkx, wky, Qx, Qy,  M ) 

use moduleone
Implicit real (kind=8) (a-h,o-z)

integer(kind=4), intent(IN):: Lens
real(kind=8),dimension(1:Lens), intent(in):: vectorx, vectory
real(kind=8)::  kx, ky, Qx, Qy, wkx, wky
integer(kind=4),intent(out):: M 

integer:: Local(1)
complex(kind=8), dimension(Lens):: Cks
real(KIND=8), dimension(Lens) :: RCKS
complex(kind=8)::  CKQ


!~ write(*,*)
!~ write(*, '(13x, 2(1xES23.15E3) )' ) vectorx(12),  vectory(12)
!~ write(*,*) 


Cks=CMPLX( vectorx, vectory)

!~ write(*,*)
!~ write(*,*) Cks(12)

!~ write(*, '( a )' )
!~ write(*, '( a, I6 )'  ) "Subroutine Lens is::", lens



Bvector1x=(2.0_DP*Pi)/(sq3*carbondistance)
Bvector1y=(2.0_DP*Pi)/(3.0_DP*carbondistance)
Xbound=(4.0_DP*Pi)/(3.0_DP*sq3*carbondistance)
Ybound=(2.0_DP*Pi)/(3.0_DP*carbondistance)

!~ write(*,*)
!~ write(*, '( 2(1xES23.15E3) )') kx, ky

kx=wkx+Qx
ky=wky+Qy

!~ write(*,*)
!~ write(*,*) kx, ky

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

!~ write(*,*)
!~ write(*,*) kx, ky

CKQ=CMPLX(kx, ky)
   !~ write(*, '( a )' )
   !~ write(*, *  ) CKQ
Cks=Cks-CKQ

do i=1, Lens
        RCKS(i)=abs(Cks(i) )
end do

Local=MINLOC( RCks )
M=Local(1)

   !~ write(*, '( a )' )
   !~ write(*, '( I6 )'  ) Local
   !~ write(*, '( a )' )
   !~ write(*, '( I6 )'  ) Local(1)   
   
end subroutine wavenumberlocal


function FermiDistribution (E,Ef,Beta)
	implicit none
	real*8 FermiDistribution,E,Ef,Beta
	integer imax
	imax=100
	FermiDistribution=Beta*(E-Ef)
	if (FermiDistribution>imax) then
	FermiDistribution=0.
	else if (FermiDistribution<-imax) then
	FermiDistribution=1.
	else
	FermiDistribution=1.+dexp(FermiDistribution)
	FermiDistribution=1./FermiDistribution
	end if
	end



