MODULE moduleone
implicit none
integer, parameter :: DP=SELECTED_REAL_KIND(15,300)
INTEGER (kind=4) :: IERROR

complex (kind=DP):: AI=(0._DP,1.0_DP)
complex (kind=DP):: Re=(1.0_DP,0._DP)
character(len=1024):: line,word,word2, lines

real(kind=DP), parameter :: pi=3.141592653589793238460_DP,  sq3=dsqrt(3.0_DP)
 
! Boltzmann constant =8.6173303(50)×10-5  ev/K
real(kind=DP), parameter :: Boltzmannconstant=8.617330350*10.0**(-5),   temperatureKelvin=300.0d0 

real(kind=DP) :: BoltzmannBETA=1.d0/(Boltzmannconstant* temperatureKelvin)

real(kind=DP) :: hbar=6.58211951440*10.0**(-16)   !  unit  eV.s

!~ BoltzmannBETA=1.d0/(Boltzmannconstant* temperatureKelvin)

INTEGER (kind=4) :: UNIT1=21
INTEGER (kind=4) :: UNIT2=22
INTEGER (kind=4) :: UNIT3=23
INTEGER (kind=4) :: UNIT4=24
INTEGER (kind=4) :: UNIT5=25
INTEGER (kind=4) :: UNIT6=26
INTEGER (kind=4) :: UNIT7=27
INTEGER (kind=4) :: UNIT8=28
INTEGER (kind=4) :: UNIT9=29
CHARACTER (LEN=56) :: FILENAME1='.dat'
CHARACTER (LEN=56) :: FILENAME2='eigenvalues.dat'
CHARACTER (LEN=56) :: FILENAME3='Hoppingts.dat'
CHARACTER (LEN=56) :: FILENAME4='Hamitonilon.dat'
CHARACTER (LEN=56) :: FILENAME5='dos.dat'
CHARACTER (LEN=56) :: FILENAME6='BandKs.dat'
CHARACTER (LEN=56) :: FILENAME7='Bandstructure0.dat'
CHARACTER (LEN=56) :: FILENAME8='gap1.dat'
CHARACTER (LEN=56) :: FILENAME9='information.out'
INTEGER (kind=4) :: UNIT10=30
INTEGER (kind=4) :: UNIT11=31
INTEGER (kind=4) :: UNIT12=32
INTEGER (kind=4) :: UNIT13=33
INTEGER (kind=4) :: UNIT14=34
INTEGER (kind=4) :: UNIT15=35
INTEGER (kind=4) :: UNIT16=36
INTEGER (kind=4) :: UNIT17=37
INTEGER (kind=4) :: UNIT18=38
INTEGER (kind=4) :: UNIT19=39
INTEGER (kind=4) :: UNIT20=40
CHARACTER (LEN=56) :: FILENAME10='input.txt'
CHARACTER (LEN=56) :: FILENAME11='wavenumber.dat'
CHARACTER (LEN=56) :: FILENAME12='qvectors.dat '
CHARACTER (LEN=56) :: FILENAME13='showgap1.plt'
CHARACTER (LEN=56) :: FILENAME14='showBandstructure0.plt'
CHARACTER (LEN=56) :: FILENAME15=' '
CHARACTER (LEN=56) :: FILENAME16=' '
CHARACTER (LEN=56) :: FILENAME17=' '
CHARACTER (LEN=56) :: FILENAME18='DynamicPolyPI6VIII.dat'
CHARACTER (LEN=56) :: FILENAME19=' '

INTEGER (kind=4) :: UNIT200=200
INTEGER (kind=4) :: UNIT201=201
INTEGER (kind=4) :: UNIT202=202
INTEGER (kind=4) :: UNIT203=203
INTEGER (kind=4) :: UNIT204=204
INTEGER (kind=4) :: UNIT205=205
INTEGER (kind=4) :: UNIT206=206

CHARACTER (LEN=56) :: FILENAME200='dos0.dat'
CHARACTER (LEN=56) :: FILENAME201='K_ferimilevel.dat'
CHARACTER (LEN=56) :: FILENAME202='K_ferimilevelplot.plt'
CHARACTER (LEN=56) :: FILENAME203='K_ferimilevel_1.dat'
CHARACTER (LEN=56) :: FILENAME204='K_ferimilevel_2.dat'

CHARACTER (LEN=56) :: FILENAME205='dosU_D.dat'

!================================================================================

real(kind=DP), parameter :: carbondistance=1.0_DP                            ! 1.420_DP*(1.0_DP/10**10)
real(kind=DP), parameter :: avalue= carbondistance*dsqrt(3.0_DP)     ! 2.46_DP 10^{-10}

integer,parameter:: Nlayer=3
integer,parameter:: Hdimension=2

integer(KIND=4), parameter :: firstiseed=12

integer(KIND=8), parameter::  length=2**21
integer(KIND=4), parameter::  Blengths=840
integer(KIND=4), parameter::  Ndeltax=800   ! 20000*4

integer(KIND=4), parameter::  neNumberofWs2=12*Nlayer
integer(KIND=4), parameter::  neNumberofMos2=13*Nlayer
integer(KIND=4), parameter::  neNumberofGraphene=2*Nlayer

integer(KIND=4), parameter::  NE0=Ndeltax/2

!~ real(kind=DP), parameter ::     AU= 10.91866729512_DP
!~ real(kind=DP), parameter ::     AL=-10.91866729512_DP

real(kind=DP), parameter ::     AU= 20.0_DP
real(kind=DP), parameter ::     AL=-20.0_DP

real(kind=DP), parameter ::    deltaX=(AU-AL)/( real(Ndeltax)*1.0d0)


real(kind=DP), parameter ::     Fermibound=1.5_DP*deltaX          !0.01_DP

real(kind=DP), parameter ::  Sunitcellmos2=8.64778*10.0**(-20)            !!   S=SQRT(3)*a^2/2   a=3.16A
real(kind=DP), parameter ::  SunitcellWs2=4.304756*2.0*10.0**(-20)   !!  a=3.153A
real(kind=DP), parameter ::  Sunitcellgraphene=2.619*2.0*10.0**(-20)

real(kind=DP), parameter :: distanceoflayermos2=1.586*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayerWs2=1.571*10.0**(-10)
real(kind=DP), parameter :: distanceoflayerGraphene=3.350*10.0**(-10)

real(kind=DP), parameter :: distanceoflayermos2W=2.968*10.0**(-10) 
real(kind=DP), parameter :: distanceoflayerws2W=3.0195*10.0**(-10) 

END MODULE moduleone
!==============================================================================


