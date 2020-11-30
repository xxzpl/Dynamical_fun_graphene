subroutine wavenumbers(M, N, deltak, Lens)
use moduleone
!==================================================
! wavenumbers in BZ
!
!  Kprime (2Pi/3sqrt(3)a ,  2pi/3a  ) ;  K  (4pi/3asqrt3, 0 )
!                      -----
!                   /          \
!                   \          / 
!                       ----
!
!!!          BZ for graphene
!
!  peiliang write at Feb. 2018
!==================================================
Implicit real (kind=DP) (a-h,o-z)  

integer(kind=4), intent(in) ::M, N
real(KIND=DP), dimension(0:M,0:N) :: rankx, ranky
real(Kind=DP), intent(out):: deltaK
integer(kind=4), intent(OUT) :: Lens


Xbound=(4.0_DP*Pi)/(3.0_DP*sq3*carbondistance)
Ybound=(2.0_DP*Pi)/(3.0_DP*carbondistance)

deltak=2.0_DP*Ybound/real(N,DP)   !!  dependent on the value of N


!## open file  'wavenumber.dat'
open (unit =UNIT11, file = FILENAME11,  status="UNKNOWN" ,  iostat=ierror)   
write(UNIT11,'(69("=") )')
write(UNIT11, '( 6(" "), a, 5(" "), a, 12(" "), a, 20(" "), a )' ) "M", "N", "Kx", "Ky"

Lens=0
do i=0,M
    do j=0,N
	   ranky(i,j)=real(j,DP)*deltak-Ybound
	   
	   rankx(i,j)=real(i,DP)*deltak-Xbound
	   
        if (    (ranky(i,j)>= -2.0_DP*Ybound+sq3*rankx(i,j) ) .and. &
				 (ranky(i,j)<= 2.0_DP*Ybound-sq3*rankx(i,j))  .and. &
		         ( ranky(i,j)>=-2.0_DP*Ybound-sq3*rankx(i,j)) .and. &
		         ( ranky(i,j)<= 2.0_DP*Ybound+sq3*rankx(i,j) ) .and. &
		         (rankx(i,j)<=Xbound)         ) then
		  
		  write(UNIT11, '(1x, I6, I6, 2(1xES23.15E3) )') i, j, rankx(i,j), ranky(i,j)
		  Lens=Lens+1
		  !~ write(UNIT11, '(1x, I10, I6, I6, 2(1xES23.15E3) )') i, j, rankx(i,j), ranky(i,j)		  
		 end if 
	end do
end do
close(UNIT11)

end subroutine wavenumbers


!-------------------------------------------------------------------------

!-------------------------------------------------------------------------


subroutine BandKs(Blens, Bkx, Bky)
use moduleone
Implicit real (kind=DP) (a-h,o-z)  
integer(kind=4):: Blens   !  839
real   (KIND=DP), dimension(0:Blens-1) :: Bkx, Bky


!~ sq3=dsqrt(3.0_DP)
sin30=sin(pi/6.0d0)
cos30=cos(pi/6.0d0)

Bkxbound=(4.0_DP*Pi)/(3.0_DP*sq3*carbondistance)
deltaBkx=Bkxbound*(1.50d0+0.5d0*sq3)/real(514-2)

BMX=pi/(sq3*carbondistance)
BMY=Pi/(3.0_DP*carbondistance)

  
open (unit =UNIT6, file = FILENAME6,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)
!write (UNIT6, '(60("="))' ) 
!write(UNIT6,'(a)')
!write(UNIT6,'( a)') 'Bkxs       and           Bkys' 
!write(UNIT6,'( 60("-") )')

do i=0, 839

If (i <= 186 ) then

Bkx(i)=real(i)*deltaBkx*cos30
Bky(i)=real(i)*deltaBkx*sin30
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)
  
else if ( i==187 ) then

Bkx(i)=BMX
Bky(i)=BMY
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

else if ( (i.GE. 188)  .and.  (i .LE. 295)  ) then

Bkx(i)=BMX+ real(i-187)*deltaBkx*sin30
Bky(i)=BMY- real(i-187)*deltaBkx*cos30
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

else if ( i==296 ) then

Bkx(i)=Bkxbound
Bky(i)=0.0d0
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

else if ( (i.GE. 297)  .and.  (i .LE. 512)  ) then

Bkx(i)=Bkxbound- real(i-296)*deltaBkx
Bky(i)=0.0d0
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)


else if ( (i.GE. 513)  .and.  (i .LE. 729)  ) then

Bkx(i)=real(i-513)*deltaBkx*sin30
Bky(i)=real(i-513)*deltaBkx*cos30
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

else if ( i==730 ) then

Bkx(i)=Bkxbound*0.5d0
Bky(i)=BMY*2.0d0
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

else if ( (i.GE. 731)  .and.  (i .LE. 838)  ) then

Bkx(i)=Bkxbound*0.5d0+ real(i-730)*deltaBkx*sin30
Bky(i)=BMY*2.0d0- real(i-730)*deltaBkx*cos30
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

else if ( i==839 ) then

Bkx(i)=BMX
Bky(i)=BMY
write(UNIT6, '(1x,2(1xE18.11) )') Bkx(i),  Bky(i)

end if
end do 
close(UNIT6)
end subroutine

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------




	
     