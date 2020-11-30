program main
use moduleone
use moduletwo

Implicit real (kind=DP) (a-h,o-z) 

real(KIND=8), ALLOCATABLE, dimension(:) :: rankx, ranky, rankxQ, rankyQ

integer(kind=4):: numberofcarrier, judge, KNumberatfermi
real(KIND=8), ALLOCATABLE, dimension(:) :: rankFSx, rankFSy, rankFSxnew, rankFSynew
real(KIND=DP), dimension(0:Blengths-1) :: Bkx, Bky
real(KIND=8)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
complex(kind=8):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
real(kind=8),dimension(Hdimension):: eig,  prob1, prob2, prob3
complex(kind=8),dimension(Hdimension,Hdimension):: vec

real(kind=8) :: RdistanceX,  RdistanceY  !    Rdistance(  1 )=Rdistance_X     Rdistance( 2 )=Rdistance_y

real(kind=8),dimension(3*Nlayer,Hdimension)::problayerU, problayerD
real(kind=8),dimension(3*Nlayer, 0:Ndeltax)::DensityofstatelayerU, DensityofstatelayerD, Densityofstatelayer, DoslayerU, DoslayerD

real(KIND=8), dimension(0:Ndeltax) :: Densityofstate, Densityofstate1, Densityofstate2, Densityofstate3, Densityofstatecheck
real(KIND=8), dimension(0:Ndeltax) :: Dos0, Dos1, Dos2, Dos3
real(kind=8),dimension(3*Nlayer)::  RNprimelayer, Potentialonlayer, checkprimelayer
real(kind=8),dimension(1:3*Nlayer-1):: carriydensity

real(KIND=8):: prob0, RNprime, distanceoflayer,  distanceoflayerW
real(kind=8):: EMiu, EMiu0, Eerror

character (len=50) :: test_name
character (len=4) :: material_name
character (len=*), parameter:: material_name1="MoS2", material_name2="WS2"

real(Kind=8), allocatable, dimension(:,:)::   FermiKS, QKS
complex(kind=8), Allocatable, dimension(:, :, :):: ckvec, ckvecQ
real(kind=8), ALLOCATABLE, dimension( : , :  )::  ckeig, ckeigQ

real(kind=8)::  Omega, Deltasmall
complex(KIND=8):: DynamicalpolarizationPI 

real(kind=8), external:: FermiDistribution    !  Function~FermiDistribution (E,Ef,Beta) 
real(kind=8), allocatable, dimension(:,:,:)::  Fenzi
integer(kind=4), allocatable, dimension(:,:,:):: MBianhua

complex(kind=8):: complexcoeff, phase
real(kind=8):: corelationeff


!-------------------------------------------------------------------------------------------------------------------

open (unit =UNIT9, file = FILENAME9,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror) !!  open informatjion.out
print*, "=================================================="
call timestamp (6)
call timestamp (UNIT9)
write(UNIT9, '( 79("=")  )') 
!write(UNIT9,'( "#", 79("-"))')
write(6,  '(a)' ) "The program is running"
write(6,  '(a)' ) "Waiting ......."
!------------------------------------------------------------------------------------------------------------------
IF (SetTransitionMetalDichalcogenides==1 ) then
    numberofcarrier=neNumberofMos2
    Sunitcell=Sunitcellmos2
	distanceoflayer=distanceoflayermos2
	material_name=material_name1
    L=4    ! lens of the material name
	write(UNIT9, '(a)')
	write(UNIT9, '(a)')"The object of material is MoS2. "
else if (SetTransitionMetalDichalcogenides==2 ) then
    numberofcarrier=neNumberofWs2
    Sunitcell=SunitcellWs2
	distanceoflayer=distanceoflayerWs2
	material_name=material_name2
    L=3    ! lens of the material name
	write(UNIT9, '(a)')
	write(UNIT9, '(a)')"The object of material is WS2 "
end if
material_name=adjustl(material_name)
material_name=trim(material_name)

!------------------------------------------------------------------------------------------------------------------
ConstantMoverN=1.2d0              !!!  (2.0d0/sq3=1.1547)   W/H

N=500*4           !   Y axis
M=NINT(N*ConstantMoverN)
write(UNIT9,'(a)')
write(UNIT9,'(a,2I6)')"M and N are:: ", M, N 
write(UNIT9,'(a)')

!----------------------------------------------------  Getting Ks  -----------------------------------------------------------------------
!! lens is the number of Ks in BZ
!!! real(Kind=DP), intent(out):: deltaK

call wavenumbers(M, N, deltak, Lens) 

!~ integer(kind=4), intent(in) ::M, N
!~ real(KIND=DP), dimension(0:M,0:N) :: rankx, ranky
!~ real(Kind=DP), intent(out):: deltaK
!~ integer(kind=4), intent(OUT) :: Lens

write(UNIT9,'(a)')
write(UNIT9,'(a,1xES23.15E3, a, I6)')"deltak value is 2.0*Ybound/real(N,DP)= ", deltak, "   and Lens=", Lens 
write(UNIT9,'(a)')
allocate(rankx(1:Lens))
allocate(ranky(1:Lens))
allocate(rankxQ(1:Lens))
allocate(rankyQ(1:Lens))
allocate(ckvec(Hdimension,Hdimension, 1:Lens) )
allocate(ckeig( Hdimension, 1:Lens) )
allocate(ckvecQ(Hdimension,Hdimension, 1:Lens) )
allocate(ckeigQ( Hdimension, 1:Lens) )


allocate(Fenzi(Hdimension,Hdimension, 1:Lens) )
allocate(MBianhua(Hdimension,Hdimension, 1:Lens) )


open (unit =UNIT11, file = FILENAME11,  status="UNKNOWN" ,  iostat=ierror)   !! wavenumber.dat

 !//Skip the first 2 lines 
    do i=1,2                
        read(UNIT11,*) line
    end do

!//Read the data from wavenumbers.dat
    do j=1, Lens           
        read(UNIT11, '(13x, 2(1xES23.15E3) )' ) rankx(j), ranky(j)
    end do
	close(UNIT11)
write(UNIT9,'(a)') "please have a check rankx and ranky in wavenumber.out"
write(UNIT9,'(a)')
write(UNIT9,'(a, 2(1xES23.15E3) )') "rankx(Lens), ranky(Lens) are:: ",  rankx(Lens), ranky(Lens)
write(UNIT9,'(a)')

write(*,'(a)')
write(*,'(a, I6)') "Lens=", Lens 
write(*,'(a)')
write(*,'(a)')
write(*,'(a, 2(1xES23.15E3) )') "rankx(10), ranky(10) are:: ",  rankx(10), ranky(10)

!-----------------------------------------------------------------------------------------------------------------
!~                                                               End of Getting Ks 
!-----------------------------------------------------------------------------------------------------------------

open (unit =UNIT7, file = FILENAME7,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)         !!===   Bandstructure0.dat 
open (unit =UNIT200, file = FILENAME200,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  !!=== dos0.dat

!~ call Hoppingvalues(SetTransitionMetalDichalcogenides, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ)   !! get the hopping values ts

call BandKs(Blengths, Bkx, Bky)   !! Blenghts defined in moduleone
 
 !~ Vdeta12=0.0d0
 !~ Vdeta23=0.0d0

!~ do i=0, Blengths-1
        !~ call Hamiltonianofxs2(Bkx(i), Bky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)    
	     !~ do j=1,11	     
                !~ write(UNIT7, '(1x,2(1xE18.11) )') real(i),  eig(j)      !!!===   Bandstructure0.dat 
	     !~ end do	     
!~ end do


!~ do i=0, Blengths-1
        !~ call Hamiltonianof1layergraphene(Bkx(i), Bky(i), eig, vec) 
	     !~ do j=1,2	     
                !~ write(UNIT7, '(1x,2(1xE18.11) )') real(i),  eig(j)      !!!===   Bandstructure0.dat 
	     !~ end do	     
!~ end do

!~ write(8, '( 2a )'  ) " ckeig(1,LK) " ,   " ckeig(2,LK) "	
          !~ write(8, '( 2ES23.15E3 )'  ) ckeig(1,i)  ,   ckeig(2,i)

 !-------------------------------------------------     Density of state    -------------------------------------------------------------------     
  
 Densityofstate=0.d0
 
 !++!$OMP  PARALLEL DO PRIVATE(I) REDUCTION(+:   Densityofstate) 
 	do i=1, Lens
	   call Hamiltonianof1layergraphene(rankx(i), ranky(i), eig, vec) 
 		   do j=1,2       
 		                 prob0=1.0d0
 		                  call counter( eig(j), prob0, Dos0 )			
			         Densityofstate=Dos0+Densityofstate

						 ckeig( j, i )=eig(j)                          !  ckeig(1, i)= - abs(E)  ckeig(2, i)=+E
						 ckvec( :, j, i )=vec(:, j)
                        
						!  ckeig(1, i)= - abs(E)  ckeig(2, i)=+E					 
		   end do	   
	end do	
!~ ============= ==============
!~ eigenfunction should be got after time phase

RdistanceX=0.0d0;     RdistanceY=  -1.0d0 
!~ phase=exp(AI*(Qx*RdistanceX+Qy*RdistanceY ))
!~ ckvecQ( 2 , :, :)=ckvecQ( 2 , :, :)*phase
!~ ================================================



!~ ===============================
!~ Dynamical polarization function pi(q, w)
!~ 2018 June 5 peiliang zhao nijmegen
!~ ===============================

theta=0.d0*pi/6.0d0
!~ theta=pi/6.0d0

QQ=1.0d0
Qx=QQ*cos(theta)
Qy=QQ*sin(theta)
rankxQ=rankx+Qx
rankyQ=ranky+Qy

 	do i=1, Lens
	   call Hamiltonianof1layergraphene(rankxQ(i), rankyQ(i), eig, vec) 
 		   do j=1,2       
						 ckeigQ( j, i )=eig(j)
						 ckvecQ( :, j, i )=vec(:, j)					 
		   end do	     
	end do	
!~ ============= ====================
!~ eigenfunction should be got after time phase
!~ =================================
open (unit =UNIT18, file = FILENAME18,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  !    IMageDynaPol.dat

Omega=0.d0
Wmax=24.0d0
Wstep=Wmax/1024.0d0
Fmiu=0.d0
Deltasmall=0.005d0
alpha=2.0/(4*pi**2)
alpha=(8.0d0*pi**2/(3.0d0*sq3*Lens))*alpha
!~ print*,""
!~ print*, "alpha=", alpha

DynamicalpolarizationPI=cmplx(0._DP, 0._DP) 

DO Lk=1, Lens     
        do j=1, 2						
                Do i=1, 2
				 fermiIJ=FermiDistribution ( ckeig(i,Lk), Fmiu, BoltzmannBETA )-FermiDistribution( ckeigQ(j,Lk), Fmiu, BoltzmannBETA  )
				 complexcoeff= dot_product(ckvec( : , i, Lk) , ckvecQ( : , j, LK))
                 corelationeff=abs(complexcoeff)**2
                 Fenzi(i, j, Lk)=corelationeff*fermiIJ
		        end do
        end do
end do   ! LK

Do LW=0, 1023
Omega=LW*Wstep
    DO Lk=1, Lens                  
        do j=1, 2
                Do i=1, 2
                 DynamicalpolarizationPI=DynamicalpolarizationPI+ Fenzi(i, j, Lk)/(AI*Deltasmall+Omega+ckeig(i,Lk)-ckeigQ(j,Lk))		 
		        end do
        end do
    end do   ! LK

!~ write(UNIT18, '(  a, 3ES26.15E3  )' )   " Omega,   DynamicalpolarizationPI= " ,  Omega,   DynamicalpolarizationPI
!~ write(*, '(a)' )
!~ write(*, '(  a, 3ES26.15E3  )' )   " Omega,   DynamicalpolarizationPI= " ,  Omega,   DynamicalpolarizationPI
	
DynamicalpolarizationPI= alpha*DynamicalpolarizationPI
write(UNIT18, '( 1x, 3(3xe20.10)  )' )  Omega/2.7d0,   -DynamicalpolarizationPI*2.7d0

DynamicalpolarizationPI=cmplx(0._DP, 0._DP) 
end do   ! LW
close(UNIT18)
!~ ===============================
 !~ oooooooooooooooooooooooooo 
!~ Dynamical polarization function pi(q, w)
!~ ===============================
 	
call Rnorm(Ndeltax+1, Densityofstate,r0)  !! DOS5plus.f90
Densityofstate=Densityofstate*2.0_DP/r0

write(  UNIT9,    '( a )'   ) 
write(  UNIT9,    '( a, 1xES23.15E3 )'   ) "The separation of Density of state [deltax] is ", deltaX  !! valued in module1.f90

rdos=0.d0		 
 do i=0, NE0
		  r=Densityofstate(i)
 	      rdos=rdos+r
 end do	     
write(  UNIT9,    '( a )'   ) ''
write(  UNIT9,    '( a, F14.10 )'   )  'rdos is ', rdos

rdos2=0.d0	 
 do i=NE0, Ndeltax
               r=Densityofstate(i)
 	      rdos2=rdos2+r
end do	
write(  UNIT9,    '( a )'   ) ''
write(  UNIT9,    '( a, F14.10 )'   )  'rdos2 is ' ,rdos2
write(  UNIT9,    '( a )'   ) ''
write(  UNIT9,    '( a, F14.10 )'   )  'rdos+rods2 should be 2 here = ',  rdos+rdos2

Densityofstate=Densityofstate/deltaX
!!write( UNIT200,   '(120("="))' )  !!=== dos0.out                     !!write( UNIT200,    '( 120("-") )')

write( UNIT200,    '(5x,"Energy(E)" ,6x, "Densityofstate" )' )   !!=== dos0.out
do i=0,  Ndeltax
    write(  UNIT200,   '(1x,E16.8,3x, E16.8 )' )   AL+real(i)*deltax ,  Densityofstate(i)       !!=== dos0.out     
end do   
 
 !-------------------------------------------------     End of Density of state    -----------------------------------------------   
 
 
write(  *,  '(  a  )' )''
call timestamp (6)	
write(  UNIT9,  '(  a  )' )''
call timestamp (UNIT9)
!stop

end program
