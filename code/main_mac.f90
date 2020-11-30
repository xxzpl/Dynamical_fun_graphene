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
stop 		

 
 
 
 
 
 !-------------------------------------------------   Reading Carrier Density   -----------------------------------------------
 
open (unit =UNIT10, file = FILENAME10,  status="UNKNOWN" , ACTION="READ", iostat=ios)   ! input.txt
   if(ios.ne.0) then
        stop ' I/O error, where is input file for reading'
   endif
	! input-data-file interpreter
	ios=0
	read(UNIT10,'(a)',iostat=ios) line
	read(UNIT10,'(a)',iostat=ios) line	
	
!####	
	read(UNIT10,'(a)',iostat=ios) line		
	if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='RNprime (unit/~E13.cm^{-2} )='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (RNprime (unit/~E13.cm^{-2} )=)'
        endif
	word=adjustl(word)
	read(word,*,iostat=ios) RNprimereal
        endif
        endif
        endif

call  Efun(RNprimereal, Efiled)    !! nprime  UNIT   (     10^{13} cm^{-2}   )

RNprime=RNprimereal*Sunitcell*10.0d0**(17)/real(Hdimension)  

write(UNIT9, '( a)') ""
write(UNIT9, '( a,  f14.8,  a )') "The given carrier density=",  RNprimereal , "*E13/cm^2  "
write(*, '( a)') ""
write(*, '( a,  f14.8,  a )') "The given carrier density=",  RNprimereal , "*E13/cm^2  "
write(UNIT9, '( a)') ""
write(UNIT9, '( a, 1xES23.15E3, a )'  ) "The given carrier density=",  RNprime,&
  &", which is  [ n in electrons per unit cell ] / [ the number of orbits ]."
write(UNIT9, '( a)') ""
write(UNIT9, '( a,  1xES23.15E3,  " V/nm"   )' )  "The outside electric filed that is&
  & caused by a positively charged gate E=", Efiled
write(UNIT9, '( a)') ""

read(UNIT10,'(a)',iostat=ios) line
	if(ios==0) then
        line=adjustl(line)
        if(index(line(1:1),"!").eq.0 .and. index(line(1:1),"#").eq.0) then ! skip comment
        word='Eerror of Precision='
        if(index(line,trim(word)).ne.0) then
        read(line(index(line,trim(word))+len_trim(word):),'(a)',iostat=ios) word
        if(ios.ne.0) then
        stop ' I/O error during input: (Eerror of Precision=)'
        endif
	word=adjustl(word)
	read(word,*,iostat=ios) Eerror
        endif
        endif
        endif
		
close(UNIT10)		
		
EMiu=0.0_DP     
Eerror=Eerror*deltaX  
write( UNIT9,  '( a )'   )
write(UNIT9, '( a, 1xES23.15E3)') 'Eerror=', Eerror

!! mkdir a folder 
NYU=(RNprimereal-int(RNprimereal))*100 
If ( int(RNprimereal) <10.0  .and.   NYU<10   ) then
write(test_name,  '( a , "Nprime", I1,".",I1, a )' ) material_name(1:L),  int(RNprimereal), NYU, "0"
else if( int(RNprimereal) >=10.0 .and. NYU<10     ) then
write(test_name,  '( a , "Nprime", I2,".",I1, a )' ) material_name(1:L),  int(RNprimereal), NYU, "0"
else if( int(RNprimereal) >=10.0  .and. NYU>=10  ) then
write(test_name,  '( a , "Nprime", I2,".",I2 )' ) material_name(1:L),  int(RNprimereal), NYU
else if( int(RNprimereal) <10.0  .and. NYU>=10  ) then
write(test_name,  '( a , "Nprime", I1,".",I2 )' ) material_name(1:L),  int(RNprimereal), NYU
end if

call system ("checkexits.bat "//test_name)
!~ call system ( "mkdir " //test_name )
!~ write(*,'(a)') test_name

!~ ------------------------------------------------------------------------------

!~ ------------------------------------------------------------------------------


RNprime1=RNprime*0.30d0
RNprime2=RNprime*0.35d0
RNprime3=RNprime*0.35d0

!~ RNprime1=RNprime*0.3d0
!~ RNprime2=RNprime*0.35d0
!~ RNprime3=RNprime*0.35d0

 !-------------------------------------------------    End of reading Carrier Density     -----------------------------------------------






 !-------------------------------------------------    neutrality point  and  carry density    -----------------------------------------------

do  M=1,50     !! DO LOOP III
 !======================   Loops   ============ 
write(UNIT9, '(a, I2, a)') "---------------- LOOP ", M,   " ---------------------------"

!~ call Vfun(nprime, carriernumber, distanceoflayer, Sunitcell,  Vdeta)
call Vfun(RNprime2+RNprime3, Hdimension, distanceoflayer,  Sunitcell, Vdeta12)  ! same to PRB79,035421 model
call Vfun(RNprime3, Hdimension, distanceoflayer, Sunitcell, Vdeta23)

write( UNIT9,  '( a )'   )
write( UNIT9,  '( a,  2x,E16.9 )'   ) "Vdeta12= ", Vdeta12
write( UNIT9,  '( a )'   )
write( UNIT9,  '( a,  2x,E16.9 )'   ) "Vdeta23= ", Vdeta23     !!  v2=0 eV   V3=V2+Vdeta23  V1=-Vdeta12


rdosprime0=rdos

    Densityofstate=0.d0
    Densityofstate1=0.d0
    Densityofstate2=0.d0
    Densityofstate3=0.d0 
	
	do i=1, Lens
	   call Hamiltonianofxs2(rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)
		   do j=1,11
		       prob1(j)= abs(vec(1,j))**2+abs(vec(2,j))**2+abs(vec(3,j))**2
		       prob2(j)= abs(vec(4,j))**2+abs(vec(5,j))**2+abs(vec(6,j))**2+abs(vec(7,j))**2+abs(vec(8,j))**2		       
		       prob3(j)= abs(vec(9,j))**2+abs(vec(10,j))**2+abs(vec(11,j))**2		       
		       prob0=1.0d0
		       call counter( eig(j), prob1(j), Dos1 )			
			Densityofstate1=Dos1+Densityofstate1
		       call counter( eig(j), prob2(j), Dos2 )			
			Densityofstate2=Dos2+Densityofstate2 
		       call counter( eig(j), prob3(j), Dos3 )			
			Densityofstate3=Dos3+Densityofstate3 
		       call counter( eig(j), prob0, Dos0 )			
			Densityofstate=Dos0+Densityofstate
		   end do	     
	end do	
   call Rnorm(Ndeltax+1, Densityofstate,r0)
   
     Densityofstate=Densityofstate*2.0_DP/(deltaX*r0*rdosprime0)
   Densityofstate1=Densityofstate1*2.0_DP/(deltaX*r0*rdosprime0)
   Densityofstate2=Densityofstate2*2.0_DP/(deltaX*r0*rdosprime0)  
   Densityofstate3=Densityofstate3*2.0_DP/(deltaX*r0*rdosprime0)  

   Biggest_inter_Value=0.0_DP
   do i=1, Ndeltax
	  B=  deltaX*abs( Densityofstate(i)-Densityofstate(i-1) )
	  if( B>= Biggest_inter_Value ) then
	     Biggest_inter_Value=B
	   else 			        
		 cycle
	  end if	  
   end do
Serror=1.1_DP*Biggest_inter_Value
write( UNIT9,  '(a)'  )
write( UNIT9,  '(a, 1xES23.15E3)'  ) "The Serror is equal to ",  Serror

rdosprime0= 1.0d0
   
Densityofstatecheck=Densityofstate-Densityofstate1- Densityofstate2- Densityofstate3
   N0=0
   	 call  area( Densityofstatecheck, N0, Ndeltax, S )
	    S=S*deltaX	
write( UNIT9,  '(a)'  )
write( UNIT9,  '(a, 1xES23.15E3)'  )    'Zero should be here = ', S 

Do kndelta=10, Ndeltax  
	 N0=0
	 S=0.0_DP 
	 call  area( Densityofstate, N0, kndelta, S )
	    S=S*deltaX	
!write(  *,    '( a, F14.10 )'   )  'S= ' ,S
	 if ( S-rdosprime0 > Serror) then
	         write( UNIT9, '( a, F14.10 )'   )  'S= ' ,S
	         write( UNIT9, '(  a  )'   )
             write( UNIT9, '(  a  )'   ) 'The program was stopped, please have a check Serror. '	 
 	         write( UNIT9, '(  a  )'   )	    	 
	         stop  
	  else if ( abs(S-rdosprime0) .LE. Serror ) then
                 s2=0.0_DP
		 call  area( Densityofstate, N0, kndelta+1, s2 )
		 s2=s2*deltaX		 
			     if( abs(s2-rdosprime0)  .GT. abs(S-rdosprime0)   ) then
					 write( UNIT9, '(  a  )'   )
					 write( UNIT9, '(  a, 1xES23.15E3  )'   )	 ' S should be close to one = ' ,S				 
				 				neutralitypoint=kndelta
					 write( UNIT9, '(  a  )'   )
					 write( UNIT9, '(  a, I6  )'   )	'The neutralitypoint is   ', neutralitypoint			
				 exit
			    else 			        
          		      cycle
          		 end if 
		end if			 
end do		        

!-------------------------------------------------------------	
!write(  *,  '(  a  )' )''
!call timestamp (6)
!stop	
!-------------------------------------------------------------		

    Do kndelta=neutralitypoint+3, Ndeltax   
	    S=0.0_DP 
	 call  area( Densityofstate, neutralitypoint, kndelta, S )
	    S=S*deltaX	 
		
	 if ( S-RNprime > Serror) then
	         write( UNIT9, '( a, F14.10 )'   )  'S= ' ,S
	         write( UNIT9, '(  a  )'   )
			 write( UNIT9, '(  a  )'   )'The program was stopped, please have a check Serror. '	 
 	         write( UNIT9, '(  a  )'   )    	 
	         stop  
	  else if ( abs(S-RNprime) .LE. Serror ) then
                 s2=0.0_DP
		 call  area( Densityofstate, neutralitypoint, kndelta+1, s2 )
		 s2=s2*deltaX		 
			     if( abs(s2-RNprime)  .GT. abs(S-RNprime)   ) then
					 write( UNIT9, '(  a  )'   )
					 write( UNIT9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu)-Sprime(neutralitypoint)= ", S	
				 				Kmiu=kndelta
					 write( UNIT9, '(  a  )'   )
					 write( UNIT9, '(  a, I6  )'   )	 'The Kmiu is   ', Kmiu			
				 exit
			    else 			        
          		      cycle
          		 end if 
		end if			 
end do		 	



	 call  area( Densityofstate, neutralitypoint, Kmiu+1, S )
	    S=S*deltaX			
	write( Unit9, '(  a  )'   )		
	write( Unit9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu+1)-Sprime(neutralitypoint)= ", S	
	
	call  area( Densityofstate, neutralitypoint, Kmiu-1, S )
	    S=S*deltaX		
	write( Unit9, '(  a  )'   )	
	write( Unit9, '(  a, 1xES23.15E3  )'   ) "Sprime(Kmiu-1)-Sprime(neutralitypoint)= ", S	
	
	
	EMiu0=AL+real(neutralitypoint)*deltaX
	EMiu=AL+real(Kmiu)*deltaX   
	
	write (UNIT9, '(a)') ''
    write (UNIT9, '(a, F16.8)') 'EMiu0=  ', EMiu0	    
    write (UNIT9, '(a)') ''      
    write (UNIT9, '(a, F16.8)') 'EMiu=  ', EMiu	    

	    A1=0.0_DP 
	 call  area( Densityofstate1, neutralitypoint, KMiu, A1)
	    A1=A1*deltaX	 	 

	    A2=0.0_DP 
	 call  area( Densityofstate2, neutralitypoint, KMiu, A2)
	    A2=A2*deltaX	 

	    A3=0.0_DP 
	 call  area( Densityofstate3, neutralitypoint, KMiu, A3)
	    A3=A3*deltaX

Errordelta=Biggest_inter_Value*0.9_DP

If (   ( abs(A1-RNprime1)<= Errordelta  .and. abs(A2-RNprime2)<= Errordelta) .and.  (abs(A3-RNprime3)<= Errordelta)  ) then	 

     write(UNIT9, '(a)')
	 write(UNIT9, '(a, 1xES23.15E3)')  'Nprime1=  ', A1
	 write(UNIT9, '(a, 1xES23.15E3)')  'Nprime2=  ', A2	 
	 write(UNIT9, '(a, 1xES23.15E3)')  'Nprime3=  ', A3
	 write(UNIT9, '(a)')
	 write(UNIT9, '(a, 1xES23.15E3)')  'Nprimenew=  ', A3+A2+A1		

         exit 
else  
	     RNprime1=(RNprime1+A1)/2.0d0
         RNprime2=(RNprime2+A2)/2.0d0
         RNprime3=(RNprime3+A3)/2.0d0 
	 
	 write(UNIT9, '(a)')""
	 write(UNIT9, '(a,  1xES23.15E3)')  'New_Nprime1=  ', RNprime1
	 write(UNIT9, '(a,  1xES23.15E3)')  'New_Nprime2=  ', RNprime2
	 write(UNIT9, '(a,  1xES23.15E3)')  'New_Nprime3=  ', RNprime3
	 write(UNIT9, '(a)')""	 
	 
end if 	
end do   !! DO LOOP III

 !-------------------------  End of  neutrality point  and  carry density    ----------------------
 
 !-----------------------------------------------------------------------------------------------------


precisionofnprime=  abs( A3+A2+A1-RNprime )*100.0d0/RNprime                                  ! percent

	 write(UNIT9, '(    a     )' )""
	 write(UNIT9, '(    a, f10.6, a   )' ) 'The precisionofnprime = ', precisionofnprime, '%.'
	 write(UNIT9, '(    40("-") , a, 40("-")    )' ) "  Loop End  "
	 !!write( UNIT200,   '(120("="))' )  

open (unit =UNIT8, file = FILENAME8,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)    ! gap1.dat

    do i=0, Blengths-1
        call Hamiltonianofxs2(Bkx(i), Bky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)    
	     do j=1,11	     
                write(UNIT8, '(1x,2(1xE18.11) )') real(i),  eig(j)
	     end do	     
    end do

Y0=0.2 ;  y1=1.2
call showgap1plot (y0, y1, RNprimereal, EMiu0, EMiu  )   ! gap1.dat.pdf

Y0=-3 ;  y1=3
call showBandstructure0plot (Y0, Y1)  ! Bandstructure0.pdf

Densityofstate=Densityofstate*rdos
Densityofstate1=Densityofstate1*rdos
Densityofstate2=Densityofstate2*rdos
Densityofstate3=Densityofstate3*rdos
  	 N0=0
	 S=0.0_DP 
	 call  area( Densityofstate, N0, Ndeltax, S )
	    S=S*deltaX	
write(  UNIT9,    '( a )'   ) ''	    
write(  UNIT9,    '( a, F14.10 )'   )  'The area of Densityofstate should be 2 at the place:: ' ,S


open (unit =UNIT5, file = FILENAME5,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)  ! dos.dat
!~ write (UNIT5, '(120("="))' ) 
write(UNIT5,'(5x,"Energy(E)" ,6x, "Density of state",2x," Sum of DOS(123)", 4x,&
  &" Densityofstate1 ", 2x," Densityofstate2 ", 2x," Densityofstate3 "   )' )
!~ write(UNIT5,'( 120("-") )')

   do i=0,Ndeltax
write(UNIT5, '(1x,E16.8,3x, E16.8, 4(3x,E16.8) )' ) AL+real(i)*deltax , Densityofstate(i), Densityofstate1(i)+&
  &Densityofstate2(i)+Densityofstate3(i),  Densityofstate1(i), Densityofstate2(i), Densityofstate3(i)
   end do
close(UNIT5)


!~ ===============================
!~ K points at fermi level
!~ 2018 June 5 peiliang zhao nijmegen
!~ ===============================

!~ allocate(ckvec(Hdimension,Hdimension, 1:Lens) )
!~ allocate(ckeig( Hdimension, 1:Lens) )

open (unit =UNIT201, file = FILENAME201,  status="UNKNOWN" , ACTION="READWRITE", iostat=ierror)  !K_ferimilevel.dat      READWRITE?

open (unit =UNIT12, file = FILENAME12,  status="UNKNOWN" , ACTION="READWRITE", iostat=ierror)  !  'qvectors.dat '  

KNumberoffermi=0
!===================================================================
do i=1, Lens
 	   call Hamiltonianofxs2(rankx(i), ranky(i), tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)	   
	          judge=0
		   do k=1,11
                         If ( abs(eig(k)-EMiu) <= Fermibound ) Then    ! Fermibound   defined in module1 file
				             judge=judge+1
			             end if	
						 
						 ckeig( k, i )=eig(k)
						 ckvec( :, k, i )=vec(:, k)
						 
		   end do
	   If (Judge.GE. 1) then 
	   KNumberoffermi=KNumberoffermi+1
	   write (UNIT201, '(1x,2f16.8)')  rankx(i), ranky(i)      !! FILENAME201='K_ferimilevel.dat'
	   end if	   
end do   
     		
!~ checking  vkvec

!~ ckeig(  k, Lens)=eig(k)
!~ ckvec( :, k, Lens)=vec(:, k)

r=0.d0
do i=1,11
r=r+abs(ckvec( i, 3, 9))**2
end do

	 write( UNIT9,  '(a)' )
	 write( UNIT9, '(a, 3f16.8 )' ) "The Vec_norm 1 should be here::",  r,  dot_product( ckvec( :, 3, 9),  ckvec( :, 3, 9) )
	 
	 	 
!~ dot_product(conjg(ckvec( :, 3, 9))*ckvec( :, 3, 9))	 
			
	 write( UNIT9,  '(a)' )
	 write( UNIT9,  '(a, I6)' ) "The number of the K points at Fermi level First time is ",  KNumberoffermi

	 
     rewind(UNIT201)  ! Put the reading position at the first letter in the first line	
	 
allocate(FermiKS( KNumberoffermi , 2))
allocate(  QKS(KNumberoffermi*(KNumberoffermi-1)/2,   2)  )

do j=1, KNumberoffermi    
	 Read (UNIT201, '(1x,2f16.8)') FermiKS(j , :)
end do
close(UNIT201)      !K_ferimilevel.dat  
call Kferimilevelplot (  )




!~ ===============      QK    =================

write (UNIT12, '(9x,a, 15x,a )') "qx", "qy"
write (UNIT12, '( a )') 

Lqks=0
do j=KNumberoffermi, 2, -1   
    do i=1,j-1
	 Lqks=Lqks+1
	 QKS(Lqks, :)=FermiKS(j , :)-FermiKS(i , :)
	 write (UNIT12, '(1x,2f16.8)') FermiKS(j , :)-FermiKS(i , :)   !! 'qvectors.dat ' 	  
	end do 
end do
write(UNIT9,  '(a)' )
write(UNIT9, '( a,  1x,2f16.8  )' ) "Random check QKS(13) in qvector.dat file at 15th line::",  QKS(13, :)

!~ ===============   END   QK    =================


!~ ===============================
!~ Dynamical polarization function pi(q, w)
!~ 2018 June 5 peiliang zhao nijmegen
!~ ===============================
 
!~ Rdistance( 1:3 , 1:3, :)= 0.0d0                                                               ! ST -->  ST
!~ Rdistance( 1:3 , 4:8, 1)=0.0d0;       Rdistance( 1:3 , 4:8, 2)= -1.0d0        ! Mo -->  ST        !      Rdistance(  : , : ,  1 )=Rdistance_X     Rdistance(  : , : ,  2 )=Rdistance_y
!~ Rdistance( 1:3 , 9:11, :)= 0.0d0                                                             ! SB -->  ST 

!~ Rdistance( 4:8 , 1:3, 1)= 0.0d0;     Rdistance( 4:8 , 1:3, 2)= 1.0d0;          !ST -->  Mo
!~ Rdistance( 4:8 , 4:8, :)= 0.0d0                                                                !Mo -->  Mo
!~ Rdistance( 4:8 , 9:11, 1)= 0.0d0;   Rdistance( 4:8 , 9:11, 2)= 1.0d0;        !SB -->  Mo

!~ Rdistance( 9:11 , 1:3, 1)=0.0d0;       Rdistance( 9:11 , 1:3, 2)=  0.0d0     !ST -->  SB
!~ Rdistance( 9:11 , 4:8, 1)=0.0d0;       Rdistance( 9:11 , 4:8, 2)= -1.0d0     !Mo -->  SB
!~ Rdistance( 9:11 , 9:11, 1)=0.0d0;     Rdistance( 9:11 , 9:11, 2)= 0.0d0     !SB -->  SB


!~ Qx=QKS(1, 1)
!~ Qy=QKS(1, 2)

!~ DynamicalpolarizationPI=(0._DP, 0._DP) 
!~ Fmiu=EMiu-EMiu0

!~ Omega=0.d0
!~ Deltasmall=0.5d0*10.0d0**(-3)
!~ alpha=2.0/(4*pi**2)

!~ DO Lk=1, Lens                  
        !~ wkx=rankx(Lk)
        !~ wky=ranky(Lk)
        !~ call wavenumberlocal( Lens, rankx(1:Lens), ranky(1:Lens),  wkx, wky, Qx, Qy,  M ) 
        write(UNIT9, '( a )' )
        write(UNIT9, '( a, I6 )'  ) "The value of relocal M is::", M
        !~ do j=1, 11
                !~ Do i=1, 11
				 !~ fermidistribution=1.0d0/(1.0+exp(BoltzmannBETA*(ckeig(i,Lk)-Fmiu)))&
								   !~ -1.0d0/(1.0+exp(BoltzmannBETA*(ckeig(j, M)-Fmiu)))	
								   
				 !~ phase=exp(-(Qx*Rdistance( i, j, 1)+Qy*Rdistance( i, j, 2)))	

                 !~ corelationeff= abs(dot_product(ckvec( : , i, Lk) , ckvec( : , j, M))*phase)**2
				 
                 !~ DynamicalpolarizationPI=corelationeff*fermidistribution/(AI*Deltasmall+Omega+ckeig(i,Lk)-ckeig(j,M))+DynamicalpolarizationPI
		        !~ end do
        !~ end do
!~ end do   

!~ DynamicalpolarizationPI= alpha*DynamicalpolarizationPI


!~ write(UNIT9,  '(a)' )
!~ write(UNIT9, '( a, 2f16.8  )' ) "DynamicalpolarizationPI( Q,W )=", DynamicalpolarizationPI

!~ dot_product(a, b)

!~ write(UNIT9, '(a)' )
!~ write(UNIT9, *) ckeig(  :, M )
!~ call Hamiltonianofxs2(rankx(13)+Qx, ranky(13)+Qy, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,Vdeta12, Vdeta23, eig, vec)	
!~ write(UNIT9, '(a)' )		
!~ write(UNIT9, *  ) eig

!~ ===============================
                 !~ END
!~ Dynamical polarization function pi(q, w)
!~ 2018 June 5 peiliang zhao nijmegen
!~ ===============================

write( UNIT9,  '(a)' )
write(UNIT9, '( 79("=")  )') 
call timestamp (UNIT9)
PRINT*
call timestamp (6)
print*, "=================================================="

close(UNIT12)   !  'qvectors.dat '  
close(UNIT8)  ! gap1.dat
close(UNIT7) ! Bandstructure0.dat
close(UNIT200 )  !dos0.dat
close(UNIT9)

call system ("move.bat " // test_name)
!~ call system ("./move.batch " // test_name)

end program
