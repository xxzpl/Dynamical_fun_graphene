!~ ==========================================================

                                                       !~ Hamiltonian of XS2

!~ ==========================================================
subroutine Hamiltonianofxs2(XK, Yk, tdp, tdd, tpp,Delta0,Delta1,Delta2,DeltaP,DeltaZ,v12,v23, eig, vec) !  ,v1,v2,v3
use moduleone

Implicit real (kind=8) (a-h,o-z) 
!----------------------------------------------------------------------------------------------------------------------
        integer,parameter:: neigenvalues=11
        real(kind=8),dimension(0:neigenvalues-1):: eig
        complex(kind=8):: vec(0:neigenvalues-1,0:neigenvalues-1)
        complex(kind=8):: H(0:(neigenvalues*(neigenvalues+1))/2-1)
        complex(kind=8):: cwork(0:2*neigenvalues-1)
        real(kind=8):: aux(0:7*neigenvalues-1)
        integer:: iwork(0:5*neigenvalues-1),ifail(0:neigenvalues-1)
        real(kind=8),parameter:: abstol=2*tiny(abstol)  ! see ZHPEVX
	
complex(kind=8):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
real(KIND=8)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
real(KIND=8):: v12, v23
!----------------------------------------------------------------------------------------------------------------------
real(KIND=8):: XK, Yk
complex(KIND=8), dimension(neigenvalues,neigenvalues)::Hamiltonian
complex(kind=8)::Expon_KRalpha, Expon_KRbeta,Expon_KRgemma 

!======================================================================================

xRalpha   = -0.5*avalue          ; yRalpha   =    0.5*sqrt(3.0_DP)*avalue
xRbeta    =   avalue                 ; yRbeta    =    0.0
xRgemma= -0.5*avalue          ; yRgemma=   -0.5*sqrt(3.0_DP)*avalue

CS_KRalpha      =   COS(xk*xRalpha+yk*yRalpha)
CS_KRbeta       =   COS(xk*xRbeta+yk*yRbeta)
CS_KRgemma  =   COS(xk*xRgemma+yk*yRgemma)

Hamiltonian(1,1)=2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP+v23
Hamiltonian(1,2)=2.0_DP*( tpp(1,2,6)*CS_KRalpha+tpp(1,2,4)*CS_KRbeta+tpp(1,2,2)*CS_KRgemma )
Hamiltonian(1,3)=2.0_DP*( tpp(1,3,6)*CS_KRalpha+tpp(1,3,4)*CS_KRbeta+tpp(1,3,2)*CS_KRgemma )

Hamiltonian(2,2)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP+v23
Hamiltonian(2,3)=2.0_DP*( tpp(2,3,6)*CS_KRalpha+tpp(2,3,4)*CS_KRbeta+tpp(2,3,2)*CS_KRgemma )
Hamiltonian(3,3)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ+v23


Hamiltonian(9,9)  =2.0_DP*( tpp(1,1,6)*CS_KRalpha+tpp(1,1,4)*CS_KRbeta+tpp(1,1,2)*CS_KRgemma )+DeltaP-v12
Hamiltonian(9,10)=Hamiltonian(1,2)
Hamiltonian(9,11)=Hamiltonian(1,3)

Hamiltonian(10,10)=2.0_DP*( tpp(2,2,6)*CS_KRalpha+tpp(2,2,4)*CS_KRbeta+tpp(2,2,2)*CS_KRgemma )+DeltaP-v12
Hamiltonian(10,11)=Hamiltonian(2,3)
Hamiltonian(11,11)=2.0_DP*( tpp(3,3,6)*CS_KRalpha+tpp(3,3,4)*CS_KRbeta+tpp(3,3,2)*CS_KRgemma )+DeltaZ-v12

Hamiltonian(4,4)=2.0_DP*( tdd(1,1,6)*CS_KRalpha+tdd(1,1,4)*CS_KRbeta+tdd(1,1,2)*CS_KRgemma )+Delta2  !  v2=0 eV
Hamiltonian(4,5)=2.0_DP*( tdd(1,2,6)*CS_KRalpha+tdd(1,2,4)*CS_KRbeta+tdd(1,2,2)*CS_KRgemma )
Hamiltonian(4,6)=2.0_DP*( tdd(1,3,6)*CS_KRalpha+tdd(1,3,4)*CS_KRbeta+tdd(1,3,2)*CS_KRgemma )
Hamiltonian(4,7)=2.0_DP*( tdd(1,4,6)*CS_KRalpha+tdd(1,4,4)*CS_KRbeta+tdd(1,4,2)*CS_KRgemma )
Hamiltonian(4,8)=2.0_DP*( tdd(1,5,6)*CS_KRalpha+tdd(1,5,4)*CS_KRbeta+tdd(1,5,2)*CS_KRgemma )

Hamiltonian(5,5)=2.0_DP*( tdd(2,2,6)*CS_KRalpha+tdd(2,2,4)*CS_KRbeta+tdd(2,2,2)*CS_KRgemma )+Delta1  !  v2=0 eV
Hamiltonian(5,6)=2.0_DP*( tdd(2,3,6)*CS_KRalpha+tdd(2,3,4)*CS_KRbeta+tdd(2,3,2)*CS_KRgemma )
Hamiltonian(5,7)=2.0_DP*( tdd(2,4,6)*CS_KRalpha+tdd(2,4,4)*CS_KRbeta+tdd(2,4,2)*CS_KRgemma )
Hamiltonian(5,8)=2.0_DP*( tdd(2,5,6)*CS_KRalpha+tdd(2,5,4)*CS_KRbeta+tdd(2,5,2)*CS_KRgemma )

Hamiltonian(6,6)=2.0_DP*( tdd(3,3,6)*CS_KRalpha+tdd(3,3,4)*CS_KRbeta+tdd(3,3,2)*CS_KRgemma )+Delta1    !  v2=0 eV
Hamiltonian(6,7)=2.0_DP*( tdd(3,4,6)*CS_KRalpha+tdd(3,4,4)*CS_KRbeta+tdd(3,4,2)*CS_KRgemma )
Hamiltonian(6,8)=2.0_DP*( tdd(3,5,6)*CS_KRalpha+tdd(3,5,4)*CS_KRbeta+tdd(3,5,2)*CS_KRgemma )

Hamiltonian(7,7)=2.0_DP*( tdd(4,4,6)*CS_KRalpha+tdd(4,4,4)*CS_KRbeta+tdd(4,4,2)*CS_KRgemma )+Delta2    !  v2=0 eV
Hamiltonian(7,8)=2.0_DP*( tdd(4,5,6)*CS_KRalpha+tdd(4,5,4)*CS_KRbeta+tdd(4,5,2)*CS_KRgemma )
Hamiltonian(8,8)=2.0_DP*( tdd(5,5,6)*CS_KRalpha+tdd(5,5,4)*CS_KRbeta+tdd(5,5,2)*CS_KRgemma )+Delta0    !  v2=0 eV


xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
xRbeta    =   0                                                         ; yRbeta    =    carbondistance          
xRgemma=  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma=    -0.5*carbondistance 

Expon_KRalpha      =exp(-AI*(xk*xRalpha+yk*yRalpha ))
Expon_KRbeta       =exp(-AI*(xk*xRbeta+yk*yRbeta ))
Expon_KRgemma   =exp(-AI*(xk*xRgemma+yk*yRgemma))

Hamiltonian(1,4)=tdp(1,1,3)*Expon_KRalpha+tdp(1,1,2)*Expon_KRgemma+tdp(1,1,1)*Expon_KRbeta
Hamiltonian(1,5)=tdp(1,2,3)*Expon_KRalpha+tdp(1,2,2)*Expon_KRgemma+tdp(1,2,1)*Expon_KRbeta
Hamiltonian(1,6)=tdp(1,3,3)*Expon_KRalpha+tdp(1,3,2)*Expon_KRgemma+tdp(1,3,1)*Expon_KRbeta 
Hamiltonian(1,7)=tdp(1,4,3)*Expon_KRalpha+tdp(1,4,2)*Expon_KRgemma+tdp(1,4,1)*Expon_KRbeta  
Hamiltonian(1,8)=tdp(1,5,3)*Expon_KRalpha+tdp(1,5,2)*Expon_KRgemma+tdp(1,5,1)*Expon_KRbeta 

Hamiltonian(2,4)=tdp(2,1,3)*Expon_KRalpha+tdp(2,1,2)*Expon_KRgemma+tdp(2,1,1)*Expon_KRbeta
Hamiltonian(2,5)=tdp(2,2,3)*Expon_KRalpha+tdp(2,2,2)*Expon_KRgemma+tdp(2,2,1)*Expon_KRbeta  
Hamiltonian(2,6)=tdp(2,3,3)*Expon_KRalpha+tdp(2,3,2)*Expon_KRgemma+tdp(2,3,1)*Expon_KRbeta
Hamiltonian(2,7)=tdp(2,4,3)*Expon_KRalpha+tdp(2,4,2)*Expon_KRgemma+tdp(2,4,1)*Expon_KRbeta 
Hamiltonian(2,8)=tdp(2,5,3)*Expon_KRalpha+tdp(2,5,2)*Expon_KRgemma+tdp(2,5,1)*Expon_KRbeta 

Hamiltonian(3,4)=tdp(3,1,3)*Expon_KRalpha+tdp(3,1,2)*Expon_KRgemma+tdp(3,1,1)*Expon_KRbeta  
Hamiltonian(3,5)=tdp(3,2,3)*Expon_KRalpha+tdp(3,2,2)*Expon_KRgemma+tdp(3,2,1)*Expon_KRbeta 
Hamiltonian(3,6)=tdp(3,3,3)*Expon_KRalpha+tdp(3,3,2)*Expon_KRgemma+tdp(3,3,1)*Expon_KRbeta 
Hamiltonian(3,7)=tdp(3,4,3)*Expon_KRalpha+tdp(3,4,2)*Expon_KRgemma+tdp(3,4,1)*Expon_KRbeta  
Hamiltonian(3,8)=tdp(3,5,3)*Expon_KRalpha+tdp(3,5,2)*Expon_KRgemma+tdp(3,5,1)*Expon_KRbeta


xRalpha   = -0.5*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   -0.5*carbondistance 
xRbeta    =   0                                                         ; yRbeta    =    carbondistance          
xRgemma=  0.5*sqrt(3.0_DP)*carbondistance          ; yRgemma=    -0.5*carbondistance 

Expon_KRalpha      =exp(AI*(xk*xRalpha+yk*yRalpha ))
Expon_KRbeta       =exp(AI*(xk*xRbeta+yk*yRbeta))
Expon_KRgemma   =exp(AI*(xk*xRgemma+yk*yRgemma))

Hamiltonian(4,9)  =tdp(1,1,6)*Expon_KRalpha+tdp(1,1,5)*Expon_KRgemma+tdp(1,1,4)*Expon_KRbeta
Hamiltonian(5,9)  =tdp(1,2,6)*Expon_KRalpha+tdp(1,2,5)*Expon_KRgemma+tdp(1,2,4)*Expon_KRbeta 
Hamiltonian(6,9)  =tdp(1,3,6)*Expon_KRalpha+tdp(1,3,5)*Expon_KRgemma+tdp(1,3,4)*Expon_KRbeta 
Hamiltonian(7,9)  =tdp(1,4,6)*Expon_KRalpha+tdp(1,4,5)*Expon_KRgemma+tdp(1,4,4)*Expon_KRbeta 
Hamiltonian(8,9)  =tdp(1,5,6)*Expon_KRalpha+tdp(1,5,5)*Expon_KRgemma+tdp(1,5,4)*Expon_KRbeta 

Hamiltonian(4,10)  =tdp(2,1,6)*Expon_KRalpha+tdp(2,1,5)*Expon_KRgemma+tdp(2,1,4)*Expon_KRbeta 
Hamiltonian(5,10)  =tdp(2,2,6)*Expon_KRalpha+tdp(2,2,5)*Expon_KRgemma+tdp(2,2,4)*Expon_KRbeta 
Hamiltonian(6,10)  =tdp(2,3,6)*Expon_KRalpha+tdp(2,3,5)*Expon_KRgemma+tdp(2,3,4)*Expon_KRbeta 
Hamiltonian(7,10)  =tdp(2,4,6)*Expon_KRalpha+tdp(2,4,5)*Expon_KRgemma+tdp(2,4,4)*Expon_KRbeta 
Hamiltonian(8,10)  =tdp(2,5,6)*Expon_KRalpha+tdp(2,5,5)*Expon_KRgemma+tdp(2,5,4)*Expon_KRbeta 

Hamiltonian(4,11)  =tdp(3,1,6)*Expon_KRalpha+tdp(3,1,5)*Expon_KRgemma+tdp(3,1,4)*Expon_KRbeta 
Hamiltonian(5,11)  =tdp(3,2,6)*Expon_KRalpha+tdp(3,2,5)*Expon_KRgemma+tdp(3,2,4)*Expon_KRbeta 
Hamiltonian(6,11)  =tdp(3,3,6)*Expon_KRalpha+tdp(3,3,5)*Expon_KRgemma+tdp(3,3,4)*Expon_KRbeta 
Hamiltonian(7,11)  =tdp(3,4,6)*Expon_KRalpha+tdp(3,4,5)*Expon_KRgemma+tdp(3,4,4)*Expon_KRbeta 
Hamiltonian(8,11)  =tdp(3,5,6)*Expon_KRalpha+tdp(3,5,5)*Expon_KRgemma+tdp(3,5,4)*Expon_KRbeta 

Hamiltonian(1,9)=tpp(1,1,10)
Hamiltonian(1,10)=tpp(1,2,10)
Hamiltonian(1,11)=tpp(1,3,10)

 Hamiltonian(2,9)=tpp(2,1,10)
Hamiltonian(2,10)=tpp(2,2,10)
Hamiltonian(2,11)=tpp(2,3,10)

 Hamiltonian(3,9)=tpp(3,1,10)
Hamiltonian(3,10)=tpp(3,2,10)
Hamiltonian(3,11)=tpp(3,3,10)


!====================================================================================================

!write(*, '(2E16.8, 2E16.8, /, 2E16.8, 2E16.8)')   ((Hamiltonian(M,N), M=1,2),N=1,2)
!write(*, *)  Hamiltonian


!open (unit =UNIT4, file = FILENAME4,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)

K=0
do i=1,11
   do j=1,i
   H(K)=Hamiltonian(j,i)
   K=K+1
   end do
end do   

  ! write(UNIT4, '(2E16.8)')  H(K)

        call ZHPEVX('V','A','U',neigenvalues,H,Vlower,Vupper,1,neigenvalues,abstol,Mfound,eig,&
          &vec,neigenvalues,cwork,aux,iwork,ifail,info)
        if(info.ne.0) then
        write(6,*) 'info: ',info
        write(6,*) 'Mfound ',Mfound
        do i=0,3
        write(6,*) 'i,ifail: ',i,ifail(i)
        enddo
        stop
        endif
	
!        write(6,*) eig
!	write(6,*)

end subroutine Hamiltonianofxs2




!~ ==========================================================

                                                       !~ Hamiltonian of Graphene

!~ ==========================================================

subroutine Hamiltonianof1layergraphene(XK, Yk, eig, vec) 
use moduleone
Implicit real (kind=8) (a-h,o-z) 
!----------------------------------------------------------------------------------------------------------------------
        integer,parameter:: neigenvalues=2
        real(kind=8),dimension(0:neigenvalues-1):: eig
        complex(kind=8):: vec(0:neigenvalues-1,0:neigenvalues-1)
        complex(kind=8):: H(0:(neigenvalues*(neigenvalues+1))/2-1)
        complex(kind=8):: cwork(0:2*neigenvalues-1)
        real(kind=8):: aux(0:7*neigenvalues-1)
        integer:: iwork(0:5*neigenvalues-1),ifail(0:neigenvalues-1)
        real(kind=8),parameter:: abstol=2*tiny(abstol)  ! see ZHPEVX
	
complex(kind=8):: tdp(1:3,1:5,1:6),tdd(1:5,1:5,1:6),tpp(1:3,1:3,1:10)
real(KIND=8)::Delta0,Delta1,Delta2,DeltaP,DeltaZ
real(KIND=8):: v12, v23
!----------------------------------------------------------------------------------------------------------------------
real(KIND=8):: XK, Yk, tp
complex(KIND=8), dimension(neigenvalues,neigenvalues)::Hamiltonian
complex(kind=8)::Expon_KRalpha, Expon_KRbeta,Expon_KRgemma 

tp=2.7d0

!~ avalue=sqrt(3)*carbondistance

!======================================================================================

xRalpha   = 0.5d0*sqrt(3.0_DP)*carbondistance          ; yRalpha   =   0.5d0*carbondistance       
xRbeta    =   0.d0                                                         ; yRbeta    =    -carbondistance             
xRgemma=  -0.5d0*sqrt(3.0_DP)*carbondistance         ; yRgemma=  0.5d0*carbondistance     

Expon_KRalpha      =exp(-AI*(xk*xRalpha+yk*yRalpha ))
Expon_KRbeta       =exp(-AI*(xk*xRbeta+yk*yRbeta))
Expon_KRgemma   =exp(-AI*(xk*xRgemma+yk*yRgemma))

Hamiltonian(1,1)  =cmplx(0.0d0, 0.0d0)
Hamiltonian(1,2)  =-tp*Expon_KRalpha-tp*Expon_KRgemma-tp*Expon_KRbeta 
Hamiltonian(2,2)  =cmplx(0.0d0, 0.0d0)

!====================================================================================================

!write(*, '(2E16.8, 2E16.8, /, 2E16.8, 2E16.8)')   ((Hamiltonian(M,N), M=1,2),N=1,2)
!write(*, *)  Hamiltonian
!open (unit =UNIT4, file = FILENAME4,  status="UNKNOWN" , ACTION="WRITE", iostat=ierror)

K=0
do i=1,neigenvalues
   do j=1,i
   H(K)=Hamiltonian(j,i)
   K=K+1
   end do
end do   

  ! write(UNIT4, '(2E16.8)')  H(K)

        call ZHPEVX('V','A','U',neigenvalues,H,Vlower,Vupper,1,neigenvalues,abstol,Mfound,&
          &eig,vec,neigenvalues,cwork,aux,iwork,ifail,info)
        if(info.ne.0) then
        write(6,*) 'info: ',info
        write(6,*) 'Mfound ',Mfound
        !do i=0,3
        !write(6,*) 'i,ifail: ',i,ifail(i)
        !enddo
        stop
        endif
	
    !~ write(6,*) eig
	!~ write(6,*)

end
!end subroutine 
