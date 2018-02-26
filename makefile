
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none

all: compute_avgForce.f08 stringlib.f90
	gfortran compute_avgForce.f08 stringlib.f90 $(flags) -o compute_avgForce.x

test: test.f08
	gfortran test.f08 $(flags) -o test.x

clean : 
	rm -f *.o *.mod *.x
