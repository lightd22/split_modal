PROCESSOR := $(shell uname -m)


  F90=gfortran
  FFLAGS=-g -C -O2 -ffree-form -I/opt/local/include #-fbounds-check -Wtabs -fcheck=all 
  FFLAGS2=$(FFLAGS)
  LDFLAGS=-L/opt/local/lib -lnetcdf -lnetcdff -framework vecLib

  .PHONY= 2d_test clean

SOURCES= mDGmod.f90 mDGsweep.f90 tfcn.f90
OBJECTS=$(SOURCES:.f90=.o)

all: $(SOURCES)  test_modsplit_2d

2d_test: test_modsplit_2d
	./test_modsplit_2d #> tmpOut.txt

test_modsplit_2d: $(OBJECTS) split_2d_modal.f90
	$(F90) $(FFLAGS) $(OBJECTS) split_2d_modal.f90 \
	         -o $@ $(LDFLAGS) 

clean:
	rm -f *.o *.mod test_modsplit_2d

%.o : %.f90
	$(F90) -c $(FFLAGS) $<
