DEST      = .
EXTHDRS         =
HDRS            =
MAKEFILE  = Makefile


FC        = gfortran
FFLAGS    = -O2  -lblas -llapack
DFLAGS    =  -fcheck=
CFLAGS    =  -c
PFLAGS    =  -fopenmp 
MODS      = -mod .


OUT       = -o $(basename $(@F)).exe


.SUFFIXES: .exe .mod .obj .f90
.f90.obj:
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $<

OBJ0 = Hamiltonian.f90 mainprogram.f90 time.f90 wavenumbers.f90
OBJ1 = Tshopping.f90 DOS5plus.f90 plotparameters.f90

all:  graphene

graphene: module1.f90 module2.f90 mainprogram.f90  $(OBJ0) $(OBJ1)  
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS)  $(LDFLAGS) $(OUT) $^

module1: module1.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $^
module2: module2.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(PFLAGS) $(CFLAGS) $^

	  
cleanup:
	@rm -f fort.*
	@rm -f *.o
	@rm -f *.out
	@rm -f *.eps
	@rm -f *.obj
	@rm -f *.mod
	@rm -f *.plt
	@rm -f *.dat    
	@rm -f *.pdf 

cleanall:
	@rm -f fort.*
	@rm -f *.o
	@rm -f *.out
	@rm -f *.eps
	@rm -f *.obj
	@rm -f *.mod
	@rm -f *.exe
	@rm -f *.plt
	@rm -f *.dat    
	@rm -f *.pdf 
