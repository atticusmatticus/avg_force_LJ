
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none

hist_2D: compute_avgForce.2Dhist.f08 stringlib.f90
	gfortran compute_avgForce.2Dhist.f08 stringlib.f90 $(flags) -o compute_avgForce.2Dhist.x

sphere: compute_avgForce.1D.sphere.f08 stringlib.f90
	gfortran compute_avgForce.1D.sphere.f08 stringlib.f90 $(flags) -o compute_avgForce.1D.sphere.x

ellipse: compute_avgForce.1D.ellipse.f08 stringlib.f90
	gfortran compute_avgForce.1D.ellipse.f08 stringlib.f90 $(flags) -o compute_avgForce.1D.ellipse.x

test: test.f08
	gfortran test.f08 $(flags) -o test.x

clean : 
	rm -f *.o *.mod *.x
